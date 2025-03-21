import hyperscan, os, sys, argparse, subprocess, shutil, zipfile
import traceback, resource, time, json, hashlib, tarfile, re
import random, signal, tempfile, platform, urllib.request
from multiprocessing import Pool, cpu_count, Manager
from Bio.SeqUtils.MeltingTemp import Tm_NN, DNA_NN4
from Bio import SeqIO
from pathlib import Path
from collections import defaultdict, Counter
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from functools import partial
from contextlib import contextmanager, nullcontext


class ConfigManager:
	"""Manages configuration and paths for AmpliconHunter."""
	def __init__(self, user_base_dir=None):
		# Determine base directory
		if user_base_dir:
			self.base_dir = Path(user_base_dir)
		else:
			self.base_dir = Path(os.environ.get("AMPLICONHUNTER_BASE_DIR", 
											  str(Path.home() / ".ampliconhunter")))
		
		# Create necessary directories
		self.cache_dir = self.base_dir / "cache"
		self.cache_dir.mkdir(parents=True, exist_ok=True)
		
		self.hmm_cache_dir = self.cache_dir / "hmm"
		self.hmm_cache_dir.mkdir(parents=True, exist_ok=True)
		
		self.genera_cache_dir = self.cache_dir / "genera"
		self.genera_cache_dir.mkdir(parents=True, exist_ok=True)
		
		self.db_dir = self.base_dir / "db"
		self.db_dir.mkdir(parents=True, exist_ok=True)
		
		self.results_dir = self.base_dir / "results"
		self.results_dir.mkdir(parents=True, exist_ok=True)
		
		self.bin_dir = self.base_dir / "bin"
		self.bin_dir.mkdir(parents=True, exist_ok=True)
		
		# Add bin directory to PATH for external tools
		os.environ["PATH"] = f"{self.bin_dir}:{os.environ.get('PATH', '')}"

	def get_hmm_cache_path(self, cache_key):
		"""Get path to cached HMM file."""
		return self.hmm_cache_dir / f"{cache_key}.hmm"
	
	def get_genera_cache_path(self, job_id):
		"""Get path to cached genera results."""
		return self.genera_cache_dir / f"{job_id}_genera.json"
	
	def get_db_path(self, db_type):
		"""Get path for database files."""
		return self.db_dir / db_type
	
	def get_taxonomy_path(self, db_type):
		"""Get path for taxonomy files."""
		return self.db_dir / f"{db_type}_taxonomy.tsv"
	
	def get_genome_list_path(self, db_type):
		"""Get path for genome list files."""
		return self.db_dir / f"{db_type}_genomes_files.txt"
	
	def get_job_dir(self, job_id):
		"""Get path for job results."""
		return self.results_dir / job_id

class DependencyManager:
	"""Manages external tool dependencies for AmpliconHunter."""
	def __init__(self, config):
		self.config = config
		self.system = platform.system().lower()
		
	def check_and_install_datasets(self):
		"""Check if NCBI datasets tool is installed and install if needed."""
		if shutil.which("datasets"):
			print("NCBI datasets tool found in PATH")
			return True
		
		datasets_path = self.config.bin_dir / "datasets"
		if datasets_path.exists():
			print(f"NCBI datasets tool found at {datasets_path}")
			# Add execution permission
			os.chmod(datasets_path, os.stat(datasets_path).st_mode | 0o111)
			return True
		
		print("NCBI datasets tool not found, attempting to download...")
		
		# URLs for different platforms
		urls = {
			"linux": "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets",
			"darwin": "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/mac/datasets",
			"windows": "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/win64/datasets.exe"
		}
		
		if self.system not in urls:
			print(f"Unsupported system: {self.system}")
			print("Please install NCBI datasets tool manually: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/")
			return False
		
		try:
			url = urls[self.system]
			print(f"Downloading datasets from {url}")
			urllib.request.urlretrieve(url, datasets_path)
			
			# Make executable
			os.chmod(datasets_path, os.stat(datasets_path).st_mode | 0o111)
			print(f"Successfully installed NCBI datasets tool to {datasets_path}")
			return True
		except Exception as e:
			print(f"Error downloading NCBI datasets tool: {str(e)}")
			print("Please install the tool manually: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/")
			return False
	
	def check_dependencies(self):
		"""Check all required external dependencies."""
		required_cmds = ['nhmmer', 'mafft', 'hmmbuild']
		missing = []
		
		for cmd in required_cmds:
			if not shutil.which(cmd):
				missing.append(cmd)
		
		if missing:
			print(f"Warning: The following required commands are not found in PATH: {', '.join(missing)}")
			print("Some functionality may be limited. Please install these tools for full functionality.")
			return False
		
		# Check and install datasets tool
		datasets_installed = self.check_and_install_datasets()
		
		return not missing and datasets_installed
	def check_hmm_dependencies(self):
		"""Check if HMM-related dependencies are available."""
		required_cmds = ['nhmmer', 'mafft', 'hmmbuild']
		missing = []

		for cmd in required_cmds:
			if not shutil.which(cmd):
				missing.append(cmd)

		if missing:
			print(f"Warning: The following required commands for HMM processing are not found: {', '.join(missing)}")
			return False

		return True


def download_refseq_database(config, db_type="complete", timeout_hours=2):
	"""
	Download and prepare RefSeq database.
	
	Args:
		config: ConfigManager instance
		db_type: Either "complete" or "all"
		timeout_hours: Maximum time to wait for download (hours)
	
	Returns:
		Boolean indicating success
	"""
	if db_type not in ["complete", "all"]:
		print("Error: db_type must be either 'complete' or 'all'")
		return False
	
	# Define paths
	db_dir = config.get_db_path(f"refseq-{db_type}")
	genome_list = config.get_genome_list_path(f"refseq-{db_type}")
	taxonomy_file = config.get_taxonomy_path(f"refseq-{db_type}")
	
	# Create temporary directory for download
	with tempfile.TemporaryDirectory() as temp_dir:
		temp_dir = Path(temp_dir)
		zip_file = temp_dir / f"refseq-{db_type}.zip"
		
		# Build download command
		datasets_cmd = shutil.which("datasets")
		if not datasets_cmd:
			datasets_cmd = str(config.bin_dir / "datasets")
		
		max_attempts = 3
		for attempt in range(1, max_attempts + 1):
			print(f"Download attempt {attempt} of {max_attempts}...")
			
			# Build command with or without assembly-level filter
			if db_type == "complete":
				cmd = [
					datasets_cmd, "download", "genome", "taxon", "2",
					"--reference", "--assembly-level", "complete",
					"--include", "genome", "--filename", str(zip_file)
				]
			else:
				cmd = [
					datasets_cmd, "download", "genome", "taxon", "2",
					"--reference", "--include", "genome", "--filename", str(zip_file)
				]
			
			# Execute download with timeout
			try:
				process = subprocess.run(
					cmd,
					stdout=subprocess.PIPE,
					stderr=subprocess.PIPE,
					text=True,
					timeout=timeout_hours * 3600
				)
				
				if process.returncode == 0:
					print("Download successful")
					break
				else:
					print(f"Download failed with error: {process.stderr}")
					if attempt < max_attempts:
						wait_time = 30 * (2 ** (attempt - 1))  # Exponential backoff
						print(f"Waiting {wait_time} seconds before retry...")
						time.sleep(wait_time)
			except subprocess.TimeoutExpired:
				print(f"Download timed out after {timeout_hours} hours")
				return False
		
		if attempt == max_attempts and process.returncode != 0:
			print("All download attempts failed")
			return False
		
		# Create the database directory
		db_dir.mkdir(parents=True, exist_ok=True)
		
		# Extract the zip file
		print(f"Extracting genome files to {db_dir}")
		with zipfile.ZipFile(zip_file, 'r') as zip_ref:
			zip_ref.extractall(db_dir)
		
		# Process files
		print("Processing genome files...")
		
		# Move .fna files to the database directory and rename to .fa
		for fna_file in Path(db_dir / "ncbi_dataset" / "data").glob("**/*.fna"):
			dest_file = db_dir / f"{fna_file.stem}.fa"
			shutil.copy2(fna_file, dest_file)
		
		# Rename each fasta file to just be GCF#
		for file in db_dir.glob("*.fa"):
			# Extract GCF or GCA accession
			accession_match = re.search(r'(GC[AF]_\d+\.\d+)', file.name)
			if accession_match:
				accession = accession_match.group(1)
				new_filename = db_dir / f"{accession}.fa"
				file.rename(new_filename)
		
		# Find all .fa files and save their absolute paths
		with open(genome_list, 'w') as f:
			for fa_file in db_dir.glob("*.fa"):
				f.write(f"{fa_file.absolute()}\n")
		
		# Generate taxonomy file using the metadata
		metadata_file = db_dir / "ncbi_dataset" / "data" / "assembly_data_report.jsonl"
		if metadata_file.exists():
			print("Generating taxonomy file...")
			_generate_taxonomy_file(metadata_file, taxonomy_file)
		else:
			print(f"Warning: Metadata file not found at {metadata_file}")
		
		# Clean up temporary files
		shutil.rmtree(db_dir / "ncbi_dataset", ignore_errors=True)
		
		# Count downloaded genomes
		genome_count = len(list(db_dir.glob("*.fa")))
		print(f"RefSeq {db_type} database update completed with {genome_count} genomes")
		
		return True

def _generate_taxonomy_file(metadata_file, output_file):
	"""Generate taxonomy file from metadata JSON."""
	# This is a simplified version of Look_Up_RefSeq_tax.py
	try:
		with open(output_file, 'w') as out:
			# Write header
			out.write("accession\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")
			
			with open(metadata_file, 'r') as f:
				for line in f:
					try:
						data = json.loads(line.strip())
						accession = data.get("accession", "")
						if not accession:
							continue
						
						# Get taxonomy information
						taxonomy = data.get("organism", {}).get("taxonomic_lineage", [])
						tax_dict = {"superkingdom": "Bacteria"}  # Default for bacteria
						
						for tax in taxonomy:
							rank = tax.get("rank", "")
							if rank in ["phylum", "class", "order", "family", "genus"]:
								tax_dict[rank] = tax.get("name", "")
						
						# Get species from scientific name
						species = data.get("organism", {}).get("sci_name", "")
						tax_dict["species"] = species
						
						# Write to file
						out.write(f"{accession}")
						for rank in ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]:
							out.write(f"\t{tax_dict.get(rank, '')}")
						out.write("\n")
					except json.JSONDecodeError:
						continue
					except Exception as e:
						print(f"Error processing entry: {str(e)}")
						continue
		print(f"Taxonomy file generated: {output_file}")
		return True
	except Exception as e:
		print(f"Error generating taxonomy file: {str(e)}")
		return False


class TimeoutError(Exception):
	pass

@contextmanager
def timeout(seconds):
	def signal_handler(signum, frame):
		raise TimeoutError(f"Execution timed out after {seconds} seconds")
	
	# Set the signal handler and alarm
	signal.signal(signal.SIGALRM, signal_handler)
	signal.alarm(seconds)
	try:
		yield
	finally:
		# Disable the alarm
		signal.alarm(0)


class TaxonomyManager:
	"""Manages taxonomy mapping and analysis for amplicon sharing patterns."""
	def __init__(self, taxonomy_file=None):
		self.taxonomy_map = {}  # genome_id -> taxonomy info
		if taxonomy_file:
			self._load_taxonomy(taxonomy_file)
	
	def _load_taxonomy(self, taxonomy_file):
		"""Load taxonomy mapping from file."""
		try:
			with open(taxonomy_file) as f:
				header = f.readline()  # Skip header
				for line in f:
					fields = line.strip().split('\t')
					if len(fields) >= 8:  # Ensure we have all taxonomy levels:	 #biosample_id	superkingdom	phylum	class	order	family	genus	species
						accession = fields[0]  # GCF with version (e.g., GCF_000022305.1)
						genus = fields[6].strip()
						species = fields[7].strip()
						
						# Skip if genus is empty or "[" (indicating placeholder)
						if not genus or genus.startswith('['):
							continue
						
						tax_info = {
							'genus': genus,
							'species': species
						}
						
						# Store both exact match and version-less key
						self.taxonomy_map[accession] = tax_info
						base_accession = accession.split('.')[0]
						if base_accession not in self.taxonomy_map:
							self.taxonomy_map[base_accession] = tax_info   #ugh this is very error prone ....
				
				print(f"Loaded {len(self.taxonomy_map)} taxonomy mappings")
				print(f"Found {len(set(info['genus'] for info in self.taxonomy_map.values()))} unique genera")
				
		except Exception as e:
			print(f"Warning: Failed to load taxonomy file: {str(e)}")
			print(f"Traceback: {traceback.format_exc()}")
	
	@staticmethod
	def parse_fasta_header(header):
		"""Parse FASTA header to extract source and other metadata.
		
		Args:
			header (str): FASTA header line (with or without leading '>')
			
		Returns:
			tuple: (source_id, orientation, other_metadata)
		"""
		try:
			# Remove '>' if present
			header = header[1:] if header.startswith('>') else header
			
			# Split on exact field markers
			try:
				#print("# First split on .source= to get everything before it")  #f"{sequence_id}.source={o}.coordinates={m}-{e}.orientation={o}.Tm={t}"
				sequence_id, remainder = header.split(".source=", 1)
				#print("# Then split remainder on .coordinates= to get source")
				source_id, remainder = remainder.split(".coordinates=", 1)
				#print("# Finally split on .orientation= to get coordinates and orientation")
				coordinates, remainder = remainder.split(".orientation=", 1)
				#print("# The remainder might have additional fields (like .Tm=)")
				orientation, temperature = remainder.split(".Tm=",1)
				#print(source_id, sequence_id, coordinates, orientation)
				return source_id, orientation, {
					'sequence_id': sequence_id,
					'coordinates': coordinates,
					'temperature': temperature
				}
			except ValueError as e:
				print(f"Error parsing header format: {str(e)}")
				print(f"Header: {header}")
				return None, None, None

		except Exception as e:
			print(f"Error parsing header: {str(e)}")
			return None, None, None

	def get_species_info(self, source_id):
		"""Get species info for a given genome ID.
		
		Args:
			source_id (str): Genome identifier (GCF_XXXXXX.X or GCF_XXXXXX)
			
		Returns:
			tuple: (species, genus) or ("Unknown", "Unknown") if not found
		"""
		try:
			# Try exact match first
			if source_id in self.taxonomy_map:
				tax_info = self.taxonomy_map[source_id]
				return tax_info['species'], tax_info['genus']
			
			# If source ends in .fa, remove it and try again
			if source_id.endswith('.fa'):
				clean_id = source_id[:-3]
				if clean_id in self.taxonomy_map:
					tax_info = self.taxonomy_map[clean_id]
					return tax_info['species'], tax_info['genus']
				
			
			print(f"No taxonomy found for genome ID: {source_id} or {clean_id}")
			return "Unknown", "Unknown"

		except Exception as e:
			print(f"Error getting species info for {source_id}: {str(e)}")
			print(f"Traceback: {traceback.format_exc()}")
			return "Unknown", "Unknown"

	def process_header(self, header):
		"""Process a FASTA header and extract taxonomy information.

		Args:
			header (str): FASTA header line

		Returns:
			tuple: (source_id, species, genus, orientation) or (None, None, None, None) if parsing fails
		"""
		source_id, orientation, metadata = self.parse_fasta_header(header)
		if not source_id:
			return None, None, None, None

		species, genus = self.get_species_info(source_id)
		return source_id, species, genus, orientation
	

def write_parameters(args):
	# Read actual primer sequences from the primer file
	with open(args.primer_file, 'r') as f:
		primer_line = f.read().strip()
		forward_primer, reverse_primer = primer_line.split()[:2]

	params = {
		"forward_primer": forward_primer,  # Actually parse the primer sequences
		"reverse_primer": reverse_primer,
		"database": os.path.basename(args.input_file).split('_')[0].lower(),
		"mismatches": args.mismatches,
		"melting_temp": args.Tm,
		"clamp": args.clamp,
		"min_length": args.Lmin,
		"max_length": args.Lmax,
		"use_decoy": args.decoy,
		"use_hmm": bool(args.hmm),
		"hmm_file": args.hmm if args.hmm else None,
		"dnac1": args.dnac1,
		"dnac2": args.dnac2,
		"Na": args.Na,
		"Tris": args.Tris,
		"Mg": args.Mg,
		"dNTPs": args.dNTPs,
		"saltcorr": args.saltcorr,
		"job_id": os.path.basename(args.output_directory),
		"timestamp": datetime.utcnow().isoformat()
	}

	params_file = os.path.join(args.output_directory, "parameters.json")
	with open(params_file, 'w') as f:
		json.dump(params, f, indent=2)
	
	print(f"Parameters written to {params_file}")


def parse_hmm_scores(hmm_file):
	"""Parse nhmmer tblout format file and extract scores.
	
	Format:
	# target name	accession   query name	accession   hmmfrom   hmm to	alifrom	ali to	envfrom	env to	 sq len	strand	E-value	 score	bias	description
	
	Args:
		hmm_file (str): Path to the nhmmer tblout output file
		
	Returns:
		dict: Dictionary mapping sequence IDs to their HMM scores and E-values
		Format: {
			'seq_id': {
				'score': float,
				'evalue': float,
				'orientation': str  # Extracted from description
			}
		}
	"""
	scores = {}
	try:
		with open(hmm_file, 'r') as f:
			for line in f:
				# Skip comments, empty lines, and lines starting with #
				if not line.strip() or line.startswith('#'):
					continue
					
				# Split the line on whitespace
				fields = line.strip().split(maxsplit=15)  # Split into 16 fields max
				if len(fields) < 15:  # Need at least 15 fields
					continue
				
				# Get target name (sequence ID), E-value, and score
				target_name = fields[0]
				evalue = float(fields[12])  # E-value is 13th field
				score = float(fields[13])   # Score is 14th field
				
				# Get orientation from description (last field)
				description = fields[15] if len(fields) > 15 else ""
				orientation = None
				if "orientation=" in description:
					orientation = description.split("orientation=")[1].split(".")[0]
				
				scores[target_name] = {
					'score': score,
					'evalue': evalue,
					'orientation': orientation
				}
				
	except Exception as e:
		print(f"Error parsing HMM scores from {hmm_file}: {str(e)}")
		traceback.print_exc()
		return None
		
	return scores if scores else None

def plot_hmm_distribution(hmm_data, amplicon_info, output_file, score_type='score', decoy_type=""):
	"""Create violin plots of HMM scores or E-values by orientation.
	
	Args:
		hmm_data (dict): Dictionary of HMM scores from parse_hmm_scores()
		amplicon_info (list): List of dictionaries containing amplicon orientations
		output_file (str): Path to save the output plot
		score_type (str): Either 'score' or 'evalue'
		decoy_type (str): Optional string to add to plot title (e.g., "Decoy")
	"""
	import matplotlib.pyplot as plt
	import seaborn as sns
	import pandas as pd
	import numpy as np
	
	if not hmm_data:
		print(f"No HMM data available for plotting")
		return
		
	# Set publication style
	set_publication_style()
	
	# Prepare data for plotting
	plot_data = []
	for seq_id, data in hmm_data.items():
		if data['orientation']:  # Only include if we have orientation info
			value = data[score_type]
			if score_type == 'evalue':
				value = -np.log10(value + 1e-300)  # Add small constant to avoid log(0)
			plot_data.append({
				'orientation': data['orientation'],
				'value': value
			})
	
	if not plot_data:
		print(f"No valid data points for HMM distribution plot")
		return
		
	# Convert to DataFrame
	df = pd.DataFrame(plot_data)
	
	# Create the plot
	plt.figure(figsize=(12, 6))
	palette = {'FR':'#ffbb78', 'FF':'#ff7f0e', 'RF':'#aec7e8', 'RR':'#1f77b4'}
	
	sns.violinplot(data=df, x='orientation', y='value', hue='orientation', 
				  order=['FR', 'RF', 'FF', 'RR'],
				  palette=palette)
	
	# Set labels and title
	plt.xlabel('Orientation')
	ylabel = 'HMM Bit Score' if score_type == 'score' else '-log10(E-value)'
	plt.ylabel(ylabel)
	
	title = f'HMM {score_type.title()} Distribution by Orientation'
	if decoy_type:
		title += f" ({decoy_type})"
	plt.title(title)
	
	plt.xticks(rotation=45)
	plt.tight_layout()
	
	# Save plot
	plt.savefig(output_file, dpi=300)
	plt.close()


def set_publication_style():
	"""Set matplotlib style for publication-quality figures."""
	plt.style.use('default')
	plt.rcParams['font.size'] = 14
	plt.rcParams['font.weight'] = 'bold'
	plt.rcParams['axes.labelweight'] = 'bold'
	plt.rcParams['axes.titleweight'] = 'bold'
	plt.rcParams['axes.labelsize'] = 14
	plt.rcParams['axes.titlesize'] = 16
	plt.rcParams['xtick.labelsize'] = 12
	plt.rcParams['ytick.labelsize'] = 12
	plt.rcParams['legend.fontsize'] = 12
	plt.rcParams['figure.titlesize'] = 16
	
class RunningStats:
	"""Class to maintain running statistics"""
	def __init__(self):
		self.n = 0
		self.sum = 0
		self.sum_sq = 0
		self._min = float('inf')
		self._max = float('-inf')
		self.median_values = []  # Store samples for approximate median

	def update(self, value):
		if value is not None and not np.isnan(value):
			self.n += 1
			self.sum += value
			self.sum_sq += value * value
			self._min = min(self._min, value)
			self._max = max(self._max, value)
			# Store subset of values for approximate median
			if len(self.median_values) < 10000:
				self.median_values.append(value)
			elif random.random() < 10000/self.n:  # Reservoir sampling
				idx = random.randint(0, len(self.median_values)-1)
				self.median_values[idx] = value

	@property
	def mean(self):
		return self.sum / self.n if self.n > 0 else 0

	@property
	def std(self):
		if self.n < 2:
			return 0
		return np.sqrt((self.sum_sq - (self.sum * self.sum) / self.n) / (self.n - 1))

	@property
	def median(self):
		if not self.median_values:
			return 0
		return np.median(self.median_values)


	
class StreamingVisualizer:
	def __init__(self, output_dir, args):
		self.output_dir = output_dir
		self.args = args
		self.patterns_file = os.path.join(output_dir, "ribotype_patterns.tsv")
		self.genome_patterns = defaultdict(lambda: defaultdict(int))
		
		# Initialize statistics trackers
		self.length_stats = defaultdict(RunningStats)  # Per orientation
		self.temp_stats = defaultdict(RunningStats)	# Per orientation
		self.orientation_counts = defaultdict(int)
		self.total_count = 0
		
		# Initialize plot data structures
		self.length_data = defaultdict(list)
		self.temp_data = defaultdict(list)
		
		# Sample size for distribution plots
		self.max_samples = 10000  # Maximum samples to keep per orientation
		
		# Tracking for taxonomy
		self.taxonomy_manager = TaxonomyManager(args.taxonomy) if args.taxonomy else None
		self.sequence_sources = defaultdict(lambda: defaultdict(int))
		self.available_genera = set()
		
	
	def write_ribotype_patterns(self, genome_patterns, output_dir):
		"""Write ribotype patterns to a TSV file during processing."""
		patterns_file = os.path.join(output_dir, "ribotype_patterns.tsv")

		with open(patterns_file, 'w') as f:
			f.write("genome_id\tpattern\n")  # Header
			for genome_id, pattern in genome_patterns.items():
				# Sort counts in descending order and join with hyphens
				pattern_str = "-".join(str(count) for count in sorted(pattern, reverse=True))
				f.write(f"{genome_id}\t{pattern_str}\n")
	
	def calculate_jaccard_similarity(self, set1, set2):
		"""Calculate Jaccard similarity between two sets."""
		if not set1 or not set2:
			return 0.0
		intersection = len(set1.intersection(set2))
		union = len(set1.union(set2))
		return intersection / union if union > 0 else 0.0


	def generate_species_heatmap(self, genus):
		"""Generate heatmap of Jaccard similarities between species in a genus with enhanced styling."""
		if not self.taxonomy_manager:
			print("No taxonomy manager available")
			return None

		print(f"Generating heatmap for genus: {genus}")

		# Get unique genomes and their amplicons for each species
		species_data = defaultdict(lambda: {'genomes': set(), 'amplicons': set()})

		# First pass: collect all unique genomes and their amplicons by species
		for sequence, sources in self.sequence_sources.items():
			for source in sources:
				species, gen = self.taxonomy_manager.get_species_info(source)
				if gen == genus:
					species_data[species]['genomes'].add(source)
					species_data[species]['amplicons'].add(sequence)

		# Filter and sort species by number of unique genomes
		top_species_info = []
		for species, data in species_data.items():
			genome_count = len(data['genomes'])
			amplicon_count = len(data['amplicons'])
			if genome_count > 0:  # Only include species with at least one genome
				top_species_info.append((species, genome_count, amplicon_count))

		top_species_info.sort(key=lambda x: x[1], reverse=True)
		top_species_info = top_species_info[:10]

		if not top_species_info:
			print(f"No species found for genus {genus}")
			return None

		# Calculate Jaccard similarity matrix
		species_labels = [s[0] for s in top_species_info]
		n = len(species_labels)
		matrix = np.zeros((n, n))

		for i in range(n):
			for j in range(n):
				matrix[i, j] = self.calculate_jaccard_similarity(
					species_data[species_labels[i]]['amplicons'],
					species_data[species_labels[j]]['amplicons']
				)

		# Set up the plot style for publication quality
		plt.style.use('default')
		plt.rcParams.update({
			'font.size': 14,
			'font.weight': 'bold',
			'axes.labelweight': 'bold',
			'axes.titleweight': 'bold',
			'axes.labelsize': 16,
			'axes.titlesize': 18,
			'xtick.labelsize': 14,
			'ytick.labelsize': 14,
			'figure.titlesize': 18
		})

		# Create figure with adjusted size
		plt.figure(figsize=(16, 14))
		plt.clf()

		# Create custom diverging colormap from white to dark red
		colors = plt.cm.YlOrRd(np.linspace(0.1, 1, 100))
		custom_cmap = plt.matplotlib.colors.LinearSegmentedColormap.from_list('custom', colors)

		# Create heatmap with enhanced styling
		ax = sns.heatmap(matrix,
						xticklabels=[f"{s}\n(n={c}; m={a})" for s, c, a in top_species_info],
						yticklabels=[f"{s}\n(n={c}; m={a})" for s, c, a in top_species_info],
						cmap=custom_cmap,
						square=True,
						fmt='.2f',
						annot=True,
						annot_kws={'size': 12, 'weight': 'bold'},
						cbar_kws={'label': 'Jaccard Similarity', 'shrink': 0.8})

		# Enhance title and labels
		plt.title(f'Amplicon Similarity Between Species in {genus}\n' +
				 '(Jaccard similarity of amplicon sets)\n' +
				 'n = number of genomes; m = number of distinct amplicons',
				 pad=20, fontweight='bold')

		# Rotate x-axis labels for better readability
		plt.xticks(rotation=45, ha='right')
		plt.yticks(rotation=0)

		# Adjust subplot parameters for better label visibility
		plt.tight_layout()

		# Add extra padding at the bottom for rotated labels
		plt.subplots_adjust(bottom=0.2)

		# Save plot with high DPI
		plots_dir = os.path.join(self.output_dir, "species_heatmaps")
		os.makedirs(plots_dir, exist_ok=True)
		out_path = os.path.join(plots_dir, f'species_similarity_{genus}.png')
		plt.savefig(out_path, dpi=300, bbox_inches='tight')
		plt.close()

		print(f"Saved enhanced heatmap to {out_path}")
		return [(s[0], s[1]) for s in top_species_info], out_path
	
	def process_header(self, header, sequence=""):
		"""Process a single FASTA header and its sequence"""
		source_id, orientation, metadata = TaxonomyManager.parse_fasta_header(header)
		if not source_id:
			return

		temperature = metadata.get('temperature')
		if temperature is not None:
			try:
				temperature = float(temperature)
			except (ValueError, TypeError):
				temperature = None
		if source_id.endswith('.fa'):
			source_id = source_id[:-3]

		# Track sequence for ribotype patterns
		if source_id and sequence:
			self.genome_patterns[source_id][sequence] += 1

			# Process taxonomy and track sequences by genome using a dict
			if sequence not in self.sequence_sources:
				self.sequence_sources[sequence] = {}
			self.sequence_sources[sequence][source_id] = 1

			# Add genus tracking
			if self.taxonomy_manager:
				species, genus = self.taxonomy_manager.get_species_info(source_id)
				if genus != "Unknown":
					self.available_genera.add(genus)

		# Update orientation statistics
		if orientation:
			self.orientation_counts[orientation] += 1
			self.total_count += 1
			if sequence:
				length = len(sequence)
				self.length_stats[orientation].update(length)
				if len(self.length_data[orientation]) < self.max_samples:
					self.length_data[orientation].append(length)
			if temperature is not None:
				self.temp_stats[orientation].update(temperature)
				if len(self.temp_data[orientation]) < self.max_samples:
					self.temp_data[orientation].append(temperature)






	def get_available_genera(self):
		"""Return sorted list of available genera."""
		return sorted(list(self.available_genera))


			
	def generate_plots(self, decoy_type=""):
		"""Generate all plots using collected statistics and samples"""
		set_publication_style()
		
		
		patterns = {}
		for genome, sequences in self.genome_patterns.items():
			counts = sorted(sequences.values(), reverse=True)
			if counts:
				patterns[genome] = counts

		patterns_file = os.path.join(self.output_dir, "ribotype_patterns.tsv") 
		with open(patterns_file, 'w') as f:
			f.write("genome_id\tpattern\n")  # Header
			for genome_id, pattern in patterns.items():
				pattern_str = "-".join(str(count) for count in sorted(pattern, reverse=True))
				f.write(f"{genome_id}\t{pattern_str}\n")
		
		palette = {'FR':'#ffbb78', 'FF':'#ff7f0e', 'RF':'#aec7e8', 'RR':'#1f77b4'}

		# Convert data to pandas DataFrames for seaborn
		length_data = []
		temp_data = []
		for orientation in ['FR', 'RF', 'FF', 'RR']:
			if self.length_data[orientation]:
				length_data.extend([(l, orientation) for l in self.length_data[orientation]])
			if self.temp_data[orientation]:
				temp_data.extend([(t, orientation) for t in self.temp_data[orientation]])

		length_df = pd.DataFrame(length_data, columns=['length', 'orientation'])
		temp_df = pd.DataFrame(temp_data, columns=['maxtemp', 'orientation'])

		# Length distribution plot
		if not length_df.empty:
			plt.figure(figsize=(12, 6))
			sns.violinplot(data=length_df, x='orientation', y='length', hue='orientation', order=['FR', 'RF', 'FF', 'RR'], palette=palette)
			title_suffix = f" ({decoy_type})" if decoy_type else ""
			plt.title(f'Amplicon Length Distribution by Orientation{title_suffix}')
			plt.xlabel('Orientation')
			plt.ylabel('Length (bp)')
			plt.ylim(self.args.Lmin, self.args.Lmax)
			plt.xticks(rotation=45)
			plt.tight_layout()

			out_path = os.path.join(self.output_dir, 
								   'length_distribution_decoy.png' if decoy_type else 'length_distribution.png')
			plt.savefig(out_path, dpi=300)
			plt.close()
		else:
			print("No lengths to plot")

		# Temperature distribution plot
		if not temp_df.empty:
			plt.figure(figsize=(12, 6))
			sns.violinplot(data=temp_df, x='orientation', y='maxtemp', hue='orientation', order=['FR', 'RF', 'FF', 'RR'], palette=palette)
			title_suffix = f" ({decoy_type})" if decoy_type else ""
			plt.title(f'Melting Temperature Distribution by Orientation{title_suffix}')
			plt.xlabel('Orientation')
			plt.ylabel('Primer Tm (°C)')
			plt.xticks(rotation=45)
			plt.tight_layout()

			out_path = os.path.join(self.output_dir, 
								   'temp_distribution_decoy.png' if decoy_type else 'temp_distribution.png')
			plt.savefig(out_path, dpi=300)
			plt.close()
		else:
			print("No temps to plot")

		# Orientation distribution pie chart
		if sum(self.orientation_counts.values()) > 0:
			plt.figure(figsize=(10, 8))
			ordered_orientations = ['FR', 'FF', 'RF', 'RR']
			counts = [self.orientation_counts[o] for o in ordered_orientations]
			non_zero = [i for i, c in enumerate(counts) if c > 0]

			if non_zero:
				values = [counts[i] for i in non_zero]
				labels = [f'{ordered_orientations[i]}\n({counts[i]:,} amplicons)' for i in non_zero]
				colors = ['#ffbb78', '#ff7f0e', '#aec7e8', '#1f77b4']  # Use consistent colors
				colors = [colors[i] for i in non_zero]

				plt.pie(values, 
					   labels=labels, 
					   autopct='%1.1f%%', 
					   colors=colors)

				title_suffix = f" ({decoy_type})" if decoy_type else ""
				plt.title(f'Amplicon Orientation Distribution{title_suffix}')
				plt.tight_layout()

				out_path = os.path.join(self.output_dir, 
									   'orientation_distribution_decoy.png' if decoy_type else 'orientation_distribution.png')
				plt.savefig(out_path, dpi=300)
			plt.close()
		else:
			print("No orientations to plot")
		
		

	def write_statistics(self):
		"""Write summary statistics to file"""
		summary = {
			'Total Amplicons': self.total_count,
			'Orientation Distribution': dict(self.orientation_counts),
			'Length Statistics': {
				orient: {
					'Mean': stats.mean,
					'Median': stats.median,
					'Std Dev': stats.std,
				}
				for orient, stats in self.length_stats.items()
			},
			'Temperature Statistics': {
				orient: {
					'Mean': stats.mean,
					'Median': stats.median,
					'Std Dev': stats.std,
				}
				for orient, stats in self.temp_stats.items()
			}
		}
		
		with open(os.path.join(self.output_dir, 'summary_statistics.txt'), 'w') as f:
			for key, value in summary.items():
				f.write(f'{key}:\n')
				if isinstance(value, dict):
					for subkey, subvalue in value.items():
						if isinstance(subvalue, dict):
							f.write(f'  {subkey}:\n')
							for k, v in subvalue.items():
								f.write(f'	{k}: {v:.2f}\n')
						else:
							f.write(f'  {subkey}: {subvalue}\n')
				else:
					f.write(f'  {value}\n')
				f.write('\n')

def create_reverse_visualizations(input_dir, output_dir, args):
	"""Create visualizations for decoy amplicons using streaming processing"""
	try:
		fasta_file = os.path.join(input_dir, 'decoy_amplicons.fa')
		
		# Check if file exists and has content
		if not os.path.exists(fasta_file) or os.path.getsize(fasta_file) == 0:
			print(f"No decoy amplicons found in {fasta_file} or file is empty; skipping decoy plots.")
			return
		
		visualizer = StreamingVisualizer(output_dir, args)
		
		# Process headers AND sequences
		print(f"Processing decoy amplicons from {fasta_file}")
		count = 0
		current_header = ""
		current_sequence = ""
		
		with open(fasta_file, 'r') as f:
			for line in f:
				if line.startswith('>'):
					if current_header and current_sequence:  # Store the previous sequence
						visualizer.process_header(current_header, current_sequence)
						count += 1
					current_header = line.strip()
					current_sequence = ""
				else:
					current_sequence += line.strip()
			
			# Don't forget to process the last sequence
			if current_header and current_sequence:
				visualizer.process_header(current_header, current_sequence)
				count += 1

		if count == 0:
			print("No valid decoy amplicon headers found; skipping decoy plots.")
			return
			
		print(f"Processed {count} decoy amplicon headers")
		
		# Generate plots and statistics with explicit decoy suffix
		visualizer.generate_plots(decoy_type="Decoy")
		visualizer.write_statistics()
		
		# Process HMM scores if they exist
		hmm_file = os.path.join(input_dir, 'decoy_amplicons_hmm_scores.txt')
		if os.path.exists(hmm_file) and os.path.getsize(hmm_file) > 0:
			print(f"Found HMM scores file: {hmm_file}")
			hmm_data = parse_hmm_scores(hmm_file)
			if hmm_data:
				print("Generating HMM distribution plots for decoy amplicons")
				plot_hmm_distribution(
					hmm_data, 
					[{'orientation': k} for k, v in visualizer.orientation_counts.items() for _ in range(v)],
					os.path.join(output_dir, 'hmm_score_distribution_decoy.png'),
					score_type='score', 
					decoy_type="Decoy"
				)
				plot_hmm_distribution(
					hmm_data, 
					[{'orientation': k} for k, v in visualizer.orientation_counts.items() for _ in range(v)],
					os.path.join(output_dir, 'hmm_evalue_distribution_decoy.png'),
					score_type='evalue', 
					decoy_type="Decoy"
				)
				print("HMM distribution plots generated for decoy amplicons")
		
		print(f"Decoy visualizations generated in {output_dir}")
				
	except Exception as e:
		print(f"Error generating decoy visualizations: {str(e)}")
		print(f"Traceback: {traceback.format_exc()}")

def create_visualizations(input_dir, output_dir, args):
	"""Create visualizations using streaming processing"""
	try:
		fasta_file = os.path.join(input_dir, 'amplicons.fa')
		if not os.path.exists(fasta_file) or os.path.getsize(fasta_file) == 0:
			print(f"No amplicons found in {fasta_file} or file is empty; skipping plots.")
			return
		
		visualizer = StreamingVisualizer(output_dir, args)
		
		# Read and process sequences
		current_header = ""
		current_sequence = ""
		with open(fasta_file, 'r') as f:
			for line in f:
				if line.startswith('>'):
					if current_header and current_sequence:  # Store the previous sequence
						visualizer.process_header(current_header, current_sequence)
					current_header = line.strip()
					current_sequence = ""
				else:
					current_sequence += line.strip()
			
			# Don't forget to process the last sequence
			if current_header and current_sequence:
				visualizer.process_header(current_header, current_sequence)
		
		# Generate plots and statistics
		visualizer.generate_plots()
		visualizer.write_statistics()
		
		# Process HMM scores if they exist
		hmm_file = os.path.join(input_dir, 'amplicons_hmm_scores.txt')
		if os.path.exists(hmm_file) and os.path.getsize(hmm_file) > 0:
			hmm_data = parse_hmm_scores(hmm_file)
			if hmm_data:
				plot_hmm_distribution(
					hmm_data, 
					[{'orientation': k} for k, v in visualizer.orientation_counts.items() for _ in range(v)],
					os.path.join(output_dir, 'hmm_score_distribution.png'),
					score_type='score'
				)
				plot_hmm_distribution(
					hmm_data, 
					[{'orientation': k} for k, v in visualizer.orientation_counts.items() for _ in range(v)],
					os.path.join(output_dir, 'hmm_evalue_distribution.png'),
					score_type='evalue'
				)
					
	except Exception as e:
		print(f"Error generating visualizations: {str(e)}")
		print(f"Traceback: {traceback.format_exc()}")

degen_comp = {'A':'T', 'T':'A', 'U':'A', 'G':'C', 'C':'G', 'Y':'R', 'R':'Y', 'S':'S', 'W':'W', 'K':'M', 'M':'K', 'B':'V', 'D':'H', 'H':'D', 'V':'B', 'N':'N'}
rev_comp = lambda y,delta=degen_comp: "".join([delta[x] for x in y[::-1]])
degen_equiv = {'R': '[AG]', 'Y': '[CT]', 'S': '[CG]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'}
create_degenerate_string = lambda x:''.join([degen_equiv.get(c,c) for c in x])
degen_equiv_set = {'A':{'A'}, 'T':{'T'}, 'U':{'T'}, 'G':{'G'}, 'C':{'C'}, 'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'C', 'G'}, 'W': {'A', 'T'}, 'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'}, 'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'}, 'N': {'A', 'C', 'G', 'T'}}

# Global variables for worker processes
_db = None
_patterns = None
_mismatches = None
_primers = None

def init_worker(patterns, mismatches, primers):
	global _db, _patterns, _mismatches, _primers
	_patterns = patterns
	_mismatches = mismatches
	_primers = primers
	
	# Initialize database in each worker
	ext = hyperscan.ExpressionExt(flags=hyperscan.HS_EXT_FLAG_HAMMING_DISTANCE, hamming_distance=mismatches)
	_db = hyperscan.Database()
	_db.compile(expressions=patterns, ids=list(range(4)), elements=4, 
				flags=[hyperscan.HS_FLAG_SOM_LEFTMOST]*4, ext=[ext]*4)

def list_nondegen_equivalents(seq):
	poss = ['']
	for c in seq:
		temp = []
		for q in degen_equiv_set[c]:
			for p in poss:
				temp.append(p+q)
		poss = temp.copy()
	return set(poss)

def read_fasta_sequence(file):
	"""Optimized generator function to yield sequences from a FASTA file."""
	# Use binary mode for faster IO
	with open(file, 'rb') as fasta_file:
		sequence_id = None
		sequence_parts = []
		
		# Read in larger chunks for better performance
		chunk_size = 1024 * 1024  # 1MB chunks
		remainder = ""
		
		while True:
			chunk = fasta_file.read(chunk_size)
			if not chunk:
				break
				
			# Convert to string and combine with remainder
			text = remainder + chunk.decode('ascii')
			lines = text.split('\n')
			
			# Save the last partial line
			remainder = lines[-1]
			lines = lines[:-1]
			
			for line in lines:
				if line.startswith('>'):
					if sequence_id:
						yield sequence_id, ''.join(sequence_parts)
					sequence_id = line[1:].strip()
					sequence_parts = []
				else:
					sequence_parts.append(line.strip())
		
		# Don't forget the last sequence
		if sequence_id:
			yield sequence_id, ''.join(sequence_parts)


def run_hmmsearch(hmm_file, fasta_file, output_file, threads=1):
	"""Run nhmmer on sequences and parse results."""
	cmd = ["nhmmer", "--cpu", str(threads), 
		   "-E", "10", "--tblout", output_file,
		   "--dna", hmm_file, fasta_file]
	subprocess.run(cmd, stdout=subprocess.DEVNULL)

def process_genome(args):
	global _db, _primers
	forward_primer, reverse_primer, forward_rc, reverse_rc, forward_primer_nondegen, reverse_primer_nondegen = _primers
	f, Tm, Lmin, Lmax, clamp, include_offtarget, decoy, dnac1, dnac2, Na, Tris, Mg, dNTPs, saltcorr = args
	
	def on_match(id, from_, to, flags, context):
		id_map = {0: 'F', 1: 'F_rc', 2: 'R', 3: 'R_rc'}
		matches.append((id_map[id], from_, to))
		return 0
	
	hyperscan_time, reading_time = 0, 0
	
	amplicons = {}
	if decoy:
		decoy_amplicons = {}
	
	read_start = time.time()
	for sequence_id, sequence in read_fasta_sequence(f):
		reading_time += time.time()-read_start
		# Process forward sequence
		matches = []
		scan_start = time.time()
		_db.scan(sequence.encode(), match_event_handler=on_match)
		hyperscan_time += time.time() - scan_start
		
		matches.sort(key = lambda x:x[1])
		
		# Clamp and temp filtering
		filtered_matches = []
		for m in matches:
			template = sequence[m[1]:m[2]]
			clamp_mismatches = 0
			if m[0] == 'F':
				for i in range(clamp):
					if template[-(i+1)] not in degen_equiv_set[forward_primer[-(i+1)]]:
						clamp_mismatches += 1
			elif m[0] == 'R':
				for i in range(clamp):
					if template[-(i+1)] not in degen_equiv_set[reverse_primer[-(i+1)]]:
						clamp_mismatches += 1
			elif m[0] == 'F_rc':
				for i in range(clamp):
					if template[i] not in degen_equiv_set[forward_rc[i]]:
						clamp_mismatches += 1
			elif m[0] == 'R_rc':
				for i in range(clamp):
					if template[i] not in degen_equiv_set[reverse_rc[i]]:
						clamp_mismatches += 1
				
			if clamp_mismatches == 0:
				

				if m[0] in {'R_rc', 'F_rc'}:
					template = rev_comp(template)
				if m[0][0]=='F':
					primer_possibilities=forward_primer_nondegen
				else:
					primer_possibilities=reverse_primer_nondegen
				all_template_poss = list_nondegen_equivalents("".join([degen_comp[c] for c in template]))
				
				all_template_temps = [] #take min to get this
				for template_possibility in all_template_poss:
					all_primer_temps = []
					for primer_possibility in primer_possibilities:
						try:
							all_primer_temps.append(Tm_NN(primer_possibility, c_seq=template_possibility, dnac1=dnac1, dnac2=dnac2, Na=Na, Tris=Tris, Mg=Mg, dNTPs=dNTPs, saltcorr=saltcorr, nn_table=DNA_NN4))
							
						except:
							all_primer_temps.append(-1)


					all_template_temps.append(max(all_primer_temps)) # max(all_primer_temps) is the max temp at which that template possibility would anneal - interpretation is *all* primers are present, so only one need be above the threshold
				max_melting_temp = min(all_template_temps) #we need *all* possible degeneracies of the template to meet this threshold, because the interpretation is that we don't know *which* one we actually have present.


				if Tm:
					if max_melting_temp > Tm:
						filtered_matches.append((m[0], m[1], m[2], max_melting_temp))
				else:
					filtered_matches.append((m[0], m[1], m[2], max_melting_temp))
					
		matches = filtered_matches
		
		# Find amplicons
		for i in range(len(matches)-1):
			if matches[i][0] in {'F', 'R'}:
				for j in range(i+1,len(matches)):
					length = matches[j][2] - matches[i][1]
					if length > Lmax:
						break
					if matches[j][0] in {'F', 'R'}:
						break
					if length < Lmin:
						continue
					if matches[j][0] in {'F_rc', 'R_rc'}:
						orientation = matches[i][0]+matches[j][0][0]
						if include_offtarget or orientation in {'FR','RF'}:
							amplicon = sequence[matches[i][1]:matches[j][2]]
							if orientation[0]=='R':
								amplicon = rev_comp(amplicon)
							temp = min(matches[i][3], matches[j][3])
							newID = f"{sequence_id}.source={os.path.basename(f)}.coordinates={matches[i][1]}-{matches[j][2]}.orientation={orientation}.Tm={temp}"
							
							amplicons[newID] = amplicon

		if decoy:
			# Process reverse sequence for decoys
			sequence = sequence[::-1]
			
			matches = []
			_db.scan(sequence.encode(), match_event_handler=on_match)

			matches.sort(key = lambda x:x[1])

			# Clamp and temp filtering
			filtered_matches = []
			for m in matches:
				template = sequence[m[1]:m[2]]
				clamp_mismatches = 0
				if m[0] == 'F':
					for i in range(clamp):
						if template[-(i+1)] not in degen_equiv_set[forward_primer[-(i+1)]]:
							clamp_mismatches += 1
				elif m[0] == 'R':
					for i in range(clamp):
						if template[-(i+1)] not in degen_equiv_set[reverse_primer[-(i+1)]]:
							clamp_mismatches += 1
				elif m[0] == 'F_rc':
					for i in range(clamp):
						if template[i] not in degen_equiv_set[forward_rc[i]]:
							clamp_mismatches += 1
				elif m[0] == 'R_rc':
					for i in range(clamp):
						if template[i] not in degen_equiv_set[reverse_rc[i]]:
							clamp_mismatches += 1

				if clamp_mismatches == 0:
					

					if m[0] in {'R_rc', 'F_rc'}:
						template = rev_comp(template)
					if m[0][0]=='F':
						primer_possibilities=forward_primer_nondegen
					else:
						primer_possibilities=reverse_primer_nondegen
					all_template_poss = list_nondegen_equivalents("".join([degen_comp[c] for c in template]))
					

					all_template_temps = [] #take min to get this
					for template_possibility in all_template_poss:
						all_primer_temps = []
						for primer_possibility in primer_possibilities:

							try:
								all_primer_temps.append(Tm_NN(primer_possibility, c_seq=template_possibility, dnac1=dnac1, dnac2=dnac2, Na=Na, Tris=Tris, Mg=Mg, dNTPs=dNTPs, saltcorr=saltcorr, nn_table=DNA_NN4))
								
							except:
								all_primer_temps.append(-1)


						all_template_temps.append(max(all_primer_temps)) # max(all_primer_temps) is the max temp at which that template possibility would anneal - interpretation is *all* primers are present, so only one need be above the threshold
					max_melting_temp = min(all_template_temps) #we need *all* possible degeneracies of the template to meet this threshold, because the interpretation is that we don't know *which* one we actually have present.


					if Tm:
						if max_melting_temp > Tm:
							filtered_matches.append((m[0], m[1], m[2], max_melting_temp))
					else:
						filtered_matches.append((m[0], m[1], m[2], max_melting_temp))

			matches = filtered_matches

			# Find amplicons
			for i in range(len(matches)-1):
				if matches[i][0] in {'F', 'R'}:
					for j in range(i+1,len(matches)):
						length = matches[j][2] - matches[i][1]
						if length > Lmax:
							break
						if matches[j][0] in {'F', 'R'}:
							break
						if length < Lmin:
							continue
						if matches[j][0] in {'F_rc', 'R_rc'}:
							orientation = matches[i][0]+matches[j][0][0]
							if include_offtarget or orientation in {'FR','RF'}:
								amplicon = sequence[matches[i][1]:matches[j][2]]
								if orientation[0]=='R':
									amplicon = rev_comp(amplicon)
								#newID = f"{sequence_id}.coordinates={matches[i][1]}-{matches[j][2]}.orientation={orientation}"
								newID = f"{sequence_id}.source={os.path.basename(f)}.coordinates={matches[i][1]}-{matches[j][2]}.orientation={orientation}"
								temp = min(matches[i][3], matches[j][3])
								newID += f".Tm={temp}"
								decoy_amplicons[newID] = amplicon
		read_start = time.time()
	
	if decoy:
		return amplicons, decoy_amplicons, reading_time, hyperscan_time
	else:
		return amplicons, reading_time, hyperscan_time

def read_primers(primer_file):
	# Read primers
	try:
		with open(primer_file) as f:
			return f.read().strip().split()
	except Exception as e:
		print(f"Error reading primer file: {e}")
		exit(1)
	
def AmpliconHunter(input_file, primer_file, output_directory, threads=1, mismatches=0, 
				  Lmin=50, Lmax=5000, Tm=None, hmm_file=None, include_offtarget=False, 
				  decoy=False, clamp=5, dnac1=1000, dnac2=25, Na=50, Tris=0, Mg=4, 
				  dNTPs=1.6, saltcorr=5, timeout_hours=None, config=None):
	timeout_seconds = int(timeout_hours * 3600) if timeout_hours else None
	num_amplicons = 0
	try:
		# If timeout is specified, use the context manager
		with timeout(timeout_seconds) if timeout_seconds else nullcontext():
			forward_primer, reverse_primer = read_primers(primer_file)

			forward_primer_nondegen = list_nondegen_equivalents(forward_primer)
			reverse_primer_nondegen = list_nondegen_equivalents(reverse_primer)

			forward_rc, reverse_rc = rev_comp(forward_primer), rev_comp(reverse_primer)
			patterns = [create_degenerate_string(P).encode() for P in [forward_primer, forward_rc, reverse_primer, reverse_rc]]

			# Read genome files
			try:
				with open(input_file) as f:
					All_genome_files = f.read().strip().split()
			except Exception as e:
				print(f"Error reading genomes file: {e}")
				exit(1)

			primers = (forward_primer, reverse_primer, forward_rc, reverse_rc, 
					  forward_primer_nondegen, reverse_primer_nondegen)

			dataset = [[f, Tm, Lmin, Lmax, clamp, include_offtarget, decoy, dnac1, dnac2, Na, Tris, Mg, dNTPs, saltcorr] 
					   for f in All_genome_files]

			with Pool(processes=threads, initializer=init_worker, 
					  initargs=(patterns, mismatches, primers)) as pool:
				result = pool.map(process_genome, dataset, 
								 chunksize=max(1, len(All_genome_files)//(threads*2)))

			print('Average hyperscan time per thread:', sum([a[-1] for a in result])/threads)
			print('Average file read time per thread:', sum([a[-2] for a in result])/threads)
			write_start = time.time()

			if not decoy:
				with open(f"{output_directory}/amplicons.fa", 'w+') as Out:
					for amplicons in result:
						if amplicons[0]:
							Out.write("\n".join([f">{ID}\n{amplicons[0][ID]}" 
											   for ID in amplicons[0]])+"\n")
							num_amplicons += len(amplicons[0])
				write_time = time.time() - write_start
				if hmm_file and num_amplicons != 0:
					run_hmmsearch(hmm_file, f"{output_directory}/amplicons.fa",
										f"{output_directory}/amplicons_hmm_scores.txt", threads=threads)
			else:
				num_decoy_amplicons = 0
				with open(f"{output_directory}/amplicons.fa", 'w+') as Out1, \
					 open(f"{output_directory}/decoy_amplicons.fa", 'w+') as Out2:
					for result_tuple in result:
						regular_amplicons, decoy_amplicons, read_time, scan_time = result_tuple
						
						if regular_amplicons:
							Out1.write("\n".join([f">{ID}\n{regular_amplicons[ID]}" 
												for ID in regular_amplicons])+"\n")
						if decoy_amplicons:
							Out2.write("\n".join([f">{ID}\n{decoy_amplicons[ID]}" 
												for ID in decoy_amplicons])+"\n")

						num_amplicons += len(regular_amplicons)
						num_decoy_amplicons += len(decoy_amplicons)
				write_time = time.time() - write_start
				if hmm_file and num_amplicons != 0:
					run_hmmsearch(hmm_file, f"{output_directory}/amplicons.fa",
								 f"{output_directory}/amplicons_hmm_scores.txt", threads=threads)
					run_hmmsearch(hmm_file, f"{output_directory}/decoy_amplicons.fa",
								 f"{output_directory}/decoy_amplicons_hmm_scores.txt", threads=threads)
			print("Write time:", write_time)

	except TimeoutError:
		print(f"\nExecution timed out after {timeout_hours} hours")
		# Clean up any temporary files or intermediate results
		cleanup_partial_results(output_directory)
		return None
		
	return num_amplicons

def cleanup_partial_results(output_directory):
	"""Clean up partial results when timeout occurs by removing the entire output directory"""
	if os.path.exists(output_directory):
		shutil.rmtree(output_directory)


def check_dependencies():
	"""Check if required external programs are available."""
	required_cmds = ['nhmmer', 'mafft', 'hmmbuild']
	for cmd in required_cmds:
		if subprocess.run(['which', cmd], capture_output=True).returncode != 0:
			print(f"Error: Required command '{cmd}' not found in PATH")
			exit(1)

			

def calculate_summary_statistics(input_file, output_directory, start_time, parameters=None):
	"""Calculate comprehensive summary statistics for the run using source file information."""
	end_time = time.time()
	runtime = end_time - start_time
	
	# Read the input file to get the list of genome files
	with open(input_file) as f:
		genome_files = f.read().strip().split()
	total_genomes = len(genome_files)
	
	# Initialize counters
	amplified_genomes = set()
	amplicons_per_genome = defaultdict(int)
	off_target_count = 0
	total_amplicons = 0
	
	# Track amplicons by genome
	amplicons_by_genome = defaultdict(lambda: defaultdict(int))
	
	# Stream through the amplicons file
	amplicons_file = os.path.join(output_directory, 'amplicons.fa')
	current_sequence = ""
	current_header = ""
	
	with open(amplicons_file) as f:
		for line in f:
			if line.startswith('>'):
				# Process previous sequence if it exists
				if current_sequence and current_header:
					total_amplicons += 1
					# Extract source file from header #newID = f"{sequence_id}.source={os.path.basename(f)}.coordinates={matches[i][1]}-{matches[j][2]}.orientation={orientation}.Tm={temp}"
					source_file = current_header.split('.coordinates=')[0].split('.source=')[1]
					orientation = current_header.split('.Tm=')[0].split('.orientation=')[1]
					
					if source_file:
						amplified_genomes.add(source_file)
						amplicons_per_genome[source_file] += 1
						# Track sequence and its count for this genome
						amplicons_by_genome[source_file][current_sequence] += 1
					
					if orientation and orientation not in ['FR', 'RF']:
						off_target_count += 1
				
				current_header = line[1:].strip()
				current_sequence = ""
			else:
				current_sequence += line.strip()
		
		# Process the last sequence
		if current_sequence and current_header:
			total_amplicons += 1
			source_file = current_header.split('.coordinates=')[0].split('.source=')[1]
			orientation = current_header.split('.Tm=')[0].split('.orientation=')[1]
			
			if source_file:
				amplified_genomes.add(source_file)
				amplicons_per_genome[source_file] += 1
				amplicons_by_genome[source_file][current_sequence] += 1
			
			if orientation and orientation not in ['FR', 'RF']:
				off_target_count += 1
	
	# Calculate basic statistics
	num_amplified_genomes = len(amplified_genomes)
	percent_amplified = (num_amplified_genomes / total_genomes * 100) if total_genomes > 0 else 0
	avg_amplicons_per_genome = (total_amplicons / num_amplified_genomes) if num_amplified_genomes > 0 else 0
	off_target_rate = (off_target_count / total_amplicons * 100) if total_amplicons > 0 else 0
	
	# Calculate ribotype patterns
	ribotypes = {}  # genome -> set of sequences
	multi_ribotypes = {}  # genome -> tuple of (seq, count) pairs
	
	for genome, sequences in amplicons_by_genome.items():
		# Basic ribotype is just the set of sequences
		ribotypes[genome] = frozenset(sequences.keys())
		# Multi-ribotype includes counts, must be sorted for consistency
		multi_ribotypes[genome] = tuple(sorted((seq, count) 
											 for seq, count in sequences.items()))
	
	# Count unique patterns
	unique_ribotype_patterns = set(ribotypes.values())
	unique_multi_ribotype_patterns = set(multi_ribotypes.values())
	
	# Calculate diversity (number of unique patterns total / number of genomes)
	ribotype_diversity = len(unique_ribotype_patterns) / num_amplified_genomes if num_amplified_genomes > 0 else 0
	multi_ribotype_diversity = len(unique_multi_ribotype_patterns) / num_amplified_genomes if num_amplified_genomes > 0 else 0
	
	# Calculate unique rates (count genomes with unique patterns)
	ribotype_counts = Counter(ribotypes.values())
	unique_ribotype_count = sum(1 for count in ribotype_counts.values() if count == 1)
	unique_ribotype_rate = (unique_ribotype_count / num_amplified_genomes * 100) if num_amplified_genomes > 0 else 0
	
	multi_ribotype_counts = Counter(multi_ribotypes.values())
	unique_multi_ribotype_count = sum(1 for count in multi_ribotype_counts.values() if count == 1)
	unique_multi_ribotype_rate = (unique_multi_ribotype_count / num_amplified_genomes * 100) if num_amplified_genomes > 0 else 0
	
	# Compile statistics
	stats = {
		"runtime_seconds": runtime,
		"runtime_formatted": f"{runtime//3600:.0f}h {(runtime%3600)//60:.0f}m {runtime%60:.0f}s",
		"total_genomes": total_genomes,
		"genomes_amplified": num_amplified_genomes,
		"percent_genomes_amplified": round(percent_amplified, 2),
		"total_amplicons": total_amplicons,
		"avg_amplicons_per_genome": round(avg_amplicons_per_genome, 2),
		"max_amplicons_per_genome": max(amplicons_per_genome.values()) if amplicons_per_genome else 0,
		"min_amplicons_per_genome": min(amplicons_per_genome.values()) if amplicons_per_genome else 0,
		"off_target_amplicons": off_target_count,
		"off_target_rate": round(off_target_rate, 2),
		"ribotype_diversity": round(ribotype_diversity * 100, 2),	   # Convert to percentage
		"multi_ribotype_diversity": round(multi_ribotype_diversity * 100, 2),  # Convert to percentage
		"unique_ribotype_rate": round(unique_ribotype_rate, 2),		 # Already percentage, just change precision
		"unique_multi_ribotype_rate": round(unique_multi_ribotype_rate, 2),	# Already percentage, just change precision
		"parameters": parameters if parameters else {}
	}
	
	# Save to JSON file
	stats_file = os.path.join(output_directory, 'run_statistics.json')
	with open(stats_file, 'w') as f:
		json.dump(stats, f, indent=2)
	
	return stats

def check_cached_hmm(args, config):
	"""Check if we have a cached HMM for these parameters."""
	# Read primers
	forward_primer, reverse_primer = read_primers(args.primer_file)
	# Generate cache key from parameters
	params = {
		'forward_primer': forward_primer,
		'reverse_primer': reverse_primer,
		'min_length': args.Lmin,
		'max_length': args.Lmax,
		'clamp': args.clamp
	}
	param_str = json.dumps(params, sort_keys=True)
	cache_key = hashlib.sha256(param_str.encode()).hexdigest()
	
	cached_hmm = config.get_hmm_cache_path(cache_key)
	return str(cached_hmm) if cached_hmm.exists() else None

def cache_hmm(params, hmm_file, config):
	"""Cache an HMM file with its parameters."""
	param_str = json.dumps(params, sort_keys=True)
	cache_key = hashlib.sha256(param_str.encode()).hexdigest()
	cached_path = config.get_hmm_cache_path(cache_key)
	shutil.copy2(hmm_file, cached_path)

	


# This is a simplified version without complex locking

def process_fasta_file(file_path):
    """Process a single file to convert sequences to uppercase.
    
    Args:
        file_path: Path to the FASTA file
        
    Returns:
        Tuple: (file_path, is_empty)
    """
    # Check if file is empty
    if os.path.getsize(file_path) == 0:
        os.remove(file_path)
        return file_path, True
    
    # Try using sed for speed
    try:
        subprocess.run(
            ["sed", "-i", "/^>/!s/[a-z]/\\U&/g", file_path], 
            check=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE
        )
        return file_path, False
    except:
        # Python fallback if sed fails
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
            tmp_path = tmp.name
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        tmp.write(line)
                    else:
                        tmp.write(line.upper())
        shutil.move(tmp_path, file_path)
        return file_path, False

def convert_fasta_uppercase(args):
    """
    Convert FASTA sequences to uppercase, leaving headers untouched.
    Simplified version with cleaner multiprocessing.
    """
    directory = args.directory
    recursive = not args.no_recursive
    threads = args.threads if args.threads else min(20, cpu_count())
    
    # Find all FASTA files
    extensions = ['.fa', '.fasta', '.fna']
    files_to_process = []
    
    print(f"Scanning directory: {directory}")
    
    if recursive:
        for root, _, files in os.walk(directory):
            for file in files:
                if any(file.endswith(ext) for ext in extensions):
                    files_to_process.append(os.path.join(root, file))
    else:
        for file in os.listdir(directory):
            file_path = os.path.join(directory, file)
            if os.path.isfile(file_path) and any(file.endswith(ext) for ext in extensions):
                files_to_process.append(file_path)
    
    if not files_to_process:
        print("No FASTA files found.")
        return
    
    print(f"Found {len(files_to_process)} FASTA files")
    print(f"Using {threads} threads for processing")
    
    # Try to import tqdm
    try:
        from tqdm import tqdm
        has_tqdm = True
    except ImportError:
        has_tqdm = False
    
    # Process files in parallel
    results = []
    with Pool(processes=threads) as pool:
        if has_tqdm:
            # Use tqdm for nice progress bars
            results = list(tqdm(
                pool.imap_unordered(process_fasta_file, files_to_process),
                total=len(files_to_process),
                desc="Converting"
            ))
        else:
            # Simple chunked processing with manual progress updates
            chunk_size = max(1, len(files_to_process) // (threads * 10))
            completed = 0
            
            for i, result in enumerate(pool.imap_unordered(
                process_fasta_file, 
                files_to_process, 
                chunksize=chunk_size
            )):
                results.append(result)
                completed += 1
                
                # Show progress periodically
                if completed % 10 == 0 or completed == len(files_to_process):
                    percent = completed * 100 // len(files_to_process)
                    print(f"Progress: {completed}/{len(files_to_process)} files ({percent}%)")
    
    # Count empty files
    empty_files = sum(1 for _, is_empty in results if is_empty)
    
    print(f"Completed! Processed {len(files_to_process)} files")
    print(f"Removed {empty_files} empty files")



	
def main():
	start_time = time.time()
	
	parser = argparse.ArgumentParser(description="AmpliconHunter: In-Silico PCR sequence extractor")
	
	# Add base directory argument
	parser.add_argument("--base-dir", help="Base directory for AmpliconHunter files (default: ~/.ampliconhunter)")
	
	# Add download-refseq command
	subparsers = parser.add_subparsers(dest="command", help="Command to execute", required=False)
	
	# Download refseq command
	download_parser = subparsers.add_parser("download-refseq", help="Download RefSeq database")
	download_parser.add_argument("--type", choices=["complete", "all"], default="complete",
							  help="Type of RefSeq database to download (default: complete)")
	download_parser.add_argument("--timeout", type=float, default=2,
							  help="Maximum download time in hours (default: 2)")
			  
	# Convert FASTA to uppercase command
	convert_parser = subparsers.add_parser("convert", help="Convert FASTA sequence lines to uppercase")
	convert_parser.add_argument("directory", help="Directory containing FASTA files to process")
	convert_parser.add_argument("--no-recursive", action="store_true", 
							   help="Do not process subdirectories recursively")
	convert_parser.add_argument("--threads", type=int, 
							   help="Number of parallel threads to use (default: auto-detect)")

	# Run amplicon hunt command
	run_parser = subparsers.add_parser("run", help="Run in-silico PCR")
	run_parser.add_argument("input_file", help="Text file containing absolute paths to input FASTA files")
	run_parser.add_argument("primer_file", help="TSV file containing primer information")
	run_parser.add_argument("output_directory", help="Output directory for results")
	run_parser.add_argument("--threads", type=int, help="Number of threads to use (default: all processors)")
	run_parser.add_argument("--Tm", type=float, help="Minimum melting temperature threshold")
	run_parser.add_argument("--mismatches", type=int, default=0, help="Number of allowed mismatches")
	run_parser.add_argument("--clamp", type=int, default=5, help="3'-most CLAMP bases are not allowed mismatches")
	run_parser.add_argument("--Lmin", type=int, default=50, help="Minimum length filter")
	run_parser.add_argument("--Lmax", type=int, default=5000, help="Maximum length filter")
	run_parser.add_argument("--decoy", action="store_true", help="Use reversed genomes as decoy sequences (default: false)")
	run_parser.add_argument("--hmm", nargs='?', const='', help="HMM processing mode (without FILE: build and use new HMM; with FILE: use existing HMM file)")
	run_parser.add_argument("--clobber", action="store_true", help="Overwrite pre-existing run in output directory")
	run_parser.add_argument("--dnac1", type=float, default=1000, help="Concentration of primer strand [nM]")
	run_parser.add_argument("--dnac2", type=float, default=25, help="Concentration of template strand [nM]")
	run_parser.add_argument("--Na", type=float, default=50, help="Sodium ion concentration [mM]")
	run_parser.add_argument("--Tris", type=float, default=0, help="Tris buffer concentration [mM]")
	run_parser.add_argument("--Mg", type=float, default=4, help="Magnesium ion concentration [mM]")
	run_parser.add_argument("--dNTPs", type=float, default=1.6, help="Total deoxynucleotide concentration [mM]")
	run_parser.add_argument("--saltcorr", type=float, default=5, help="Salt correction method (0-5)")
	run_parser.add_argument("--taxonomy", help="Path to taxonomy mapping file")
	run_parser.add_argument("--no-plots", action="store_true", help="Skip generating visualization plots")
	run_parser.add_argument("--timeout", type=float, help="Maximum execution time in hours")
	
	args = parser.parse_args()
	
	# Handle case when no command is specified
	if args.command is None:
		parser.print_help()
		return 0
	
	# Initialize config manager
	config = ConfigManager(args.base_dir)
	
	# Initialize dependency manager and check dependencies
	dep_manager = DependencyManager(config)
	dep_manager.check_dependencies()
	
	if args.command == "download-refseq":
		# Download RefSeq database
		success = download_refseq_database(config, args.type, args.timeout)
		return 0 if success else 1
	elif args.command == "convert":
		convert_fasta_uppercase(args)
		return 0
	elif args.command == "run":
		# Run AmpliconHunter
		print('Welcome to AmpliconHunter!')
		
		if args.threads is None:
			args.threads = max(1, cpu_count())
		
		if args.clobber and os.path.exists(args.output_directory):
			shutil.rmtree(args.output_directory)
			
		os.makedirs(args.output_directory, exist_ok=args.clobber)
		if not args.clobber and len(os.listdir(args.output_directory)) != 0:
			print(f"Error - {args.output_directory} not empty. Please specify a new directory, or run with --clobber to overwrite existing files.")
			exit(1)
		
		write_parameters(args)
		
		if args.hmm is not None:
			if not dep_manager.check_hmm_dependencies():
				print("Warning: Missing dependencies for HMM processing. Continuing without HMM.")
				args.hmm = None
			else:
				if not args.hmm:
					cached_hmm = check_cached_hmm(args, config)
					if cached_hmm:
						print(f"Using cached HMM from previous run")
						args.hmm = cached_hmm
						os.makedirs(os.path.join(args.output_directory, "refseq"), exist_ok=True)
						output_hmm = os.path.join(args.output_directory, "refseq", "hmm.txt")
						shutil.copy2(cached_hmm, output_hmm)
						args.hmm = output_hmm
					else:
						refseq_dir = config.get_db_path("refseq-complete")
						if not os.path.exists(refseq_dir):
							print(f"Error: Required directory {refseq_dir} not found. Run 'ampliconhunter download-refseq' first.")
							exit(1)
						if not os.listdir(refseq_dir):
							print(f"Error: Directory {refseq_dir} is empty")
							exit(1)

						os.makedirs(os.path.join(args.output_directory, "refseq"), exist_ok=args.clobber)
						refseq_files = os.path.join(args.output_directory, "refseq", "all_files.txt")

						with open(refseq_files, 'w') as f:
							for file in os.listdir(refseq_dir):
								if file.endswith(('.fa', '.fasta', '.fna')):
									f.write(os.path.join(refseq_dir, file) + '\n')

						print("Running AmpliconHunter on refseq with no mismatches:")
						num_amps = AmpliconHunter(refseq_files, args.primer_file, 
												os.path.join(args.output_directory, "refseq"),
												threads=args.threads, mismatches=0, 
												Lmin=args.Lmin, Lmax=args.Lmax, 
												include_offtarget=False, clamp=args.clamp,
												config=config)

						if num_amps != 0:
							align_time = time.time()
							print(f"Running mafft on {num_amps} refseq amplicons:")
							os.system(f"mafft --auto --quiet --thread {args.threads} {args.output_directory}/refseq/amplicons.fa > {args.output_directory}/refseq/amplicons-aligned.afa")
							print(f"mafft took: {time.time()-align_time:.3f}s. Running hmmbuild:")
							build_time = time.time()
							os.system(f"hmmbuild -o {args.output_directory}/refseq/hmm-build.log --cpu {args.threads} {args.output_directory}/refseq/hmm.txt {args.output_directory}/refseq/amplicons-aligned.afa")
							print(f"hmmbuild took: {time.time()-build_time:.3f}s")
							args.hmm = f"{args.output_directory}/refseq/hmm.txt"
							
							forward_primer, reverse_primer = read_primers(args.primer_file)
							hmm_params = {
								'forward_primer': forward_primer,
								'reverse_primer': reverse_primer,
								'min_length': args.Lmin,
								'max_length': args.Lmax,
								'clamp': args.clamp
							}
							cache_hmm(hmm_params, args.hmm, config)
						else:
							print("No amplicons generated in refseq with these primers. Proceeding without HMM.")
							args.hmm = None
		
		# Run the main AmpliconHunter analysis
		print("Running AmpliconHunter on target genomes")
		try:
			num_amplicons = AmpliconHunter(args.input_file, args.primer_file, args.output_directory,
						   threads=args.threads, mismatches=args.mismatches, 
						   Lmin=args.Lmin, Lmax=args.Lmax, Tm=args.Tm, 
						   hmm_file=args.hmm, decoy=args.decoy, 
						   include_offtarget=True, clamp=args.clamp, 
						   dnac1=args.dnac1, dnac2=args.dnac2, 
						   Na=args.Na, Tris=args.Tris, Mg=args.Mg, 
						   dNTPs=args.dNTPs, saltcorr=args.saltcorr,
						   timeout_hours=args.timeout,
						   config=config)
			if num_amplicons is None:
				print("Run terminated due to timeout")
				sys.exit(1)
				
		except KeyboardInterrupt:
			print("\nRun interrupted by user")
			cleanup_partial_results(args.output_directory)
			sys.exit(1)

		# Generate plots if requested and if amplicons were found
		if not args.no_plots and num_amplicons > 0:
			print("\nGenerating visualizations...")
			try:
				# Create a plots subdirectory
				plots_dir = os.path.join(args.output_directory, "plots")
				os.makedirs(plots_dir, exist_ok=True)
				
				# Call the visualization functions
				create_visualizations(args.output_directory, plots_dir, args)
				
				if args.decoy:
					decoy_file = os.path.join(args.output_directory, 'decoy_amplicons.fa')
					if os.path.exists(decoy_file):
						print(f"Found decoy amplicons file with size: {os.path.getsize(decoy_file)} bytes")
						try:
							create_reverse_visualizations(args.output_directory, plots_dir, args)
							print("Decoy visualizations generated successfully")
						except Exception as e:
							print(f"Error generating decoy visualizations: {str(e)}")
							print(f"Traceback: {traceback.format_exc()}")
					else:
						print("No decoy_amplicons.fa file found")
					
				print(f"Plots have been generated in {plots_dir}")
			except Exception as e:
				print(f"Warning: Failed to generate plots: {str(e)}")
				print("Analysis completed successfully, but visualization generation failed.")
		
		with open(args.primer_file) as f:
			forward_primer, reverse_primer = f.read().strip().split()
		
		# Generate run statistics
		run_parameters = {
			"forward_primer": forward_primer,
			"reverse_primer": reverse_primer,
			"database": os.path.basename(args.input_file).split('_')[0],
			"mismatches": args.mismatches,
			"melting_temp": args.Tm,
			"clamp": args.clamp,
			"min_length": args.Lmin,
			"max_length": args.Lmax,
			"use_decoy": args.decoy,
			"use_hmm": bool(args.hmm),
			"hmm_file": args.hmm if args.hmm else None,
			"dnac1": args.dnac1,
			"dnac2": args.dnac2,
			"Na": args.Na,
			"Tris": args.Tris,
			"Mg": args.Mg,
			"dNTPs": args.dNTPs,
			"saltcorr": args.saltcorr,
			"output_directory": args.output_directory
		}
		
		stats = calculate_summary_statistics(args.input_file, args.output_directory, start_time, run_parameters)
		print("\nRun Summary:")
		print(f"Runtime: {stats['runtime_formatted']}")
		print(f"Genomes amplified: {stats['genomes_amplified']}/{stats['total_genomes']} ({stats['percent_genomes_amplified']}%)")
		print(f"Total amplicons: {stats['total_amplicons']}")
		print(f"Average amplicons per genome: {stats['avg_amplicons_per_genome']}")
		print(f"Off-target rate: {stats['off_target_rate']}%")
	
	else:
		# No command specified or invalid command
		parser.print_help()

if __name__ == "__main__":
	start_time = time.time()
	start_resources = resource.getrusage(resource.RUSAGE_SELF)
	
	main()
	
	end_resources = resource.getrusage(resource.RUSAGE_SELF)
	end_time = time.time()

	# Wall clock time
	real = end_time - start_time
	# User mode time
	user = end_resources.ru_utime - start_resources.ru_utime
	# System mode time
	sys = end_resources.ru_stime - start_resources.ru_stime

	print(f"real: {real:.3f}s")
	print(f"user: {user:.3f}s")
	print(f"sys: {sys:.3f}s")
		
