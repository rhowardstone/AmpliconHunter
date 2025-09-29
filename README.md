# AmpliconHunter

A scalable tool for accurate PCR amplicon prediction from microbiome samples using degenerate primers. ([Preprint](https://arxiv.org/abs/2509.13300))

([Web Interface](https://ah1.engr.uconn.edu))




## Overview

AmpliconHunter is a high-performance, in-silico PCR tool designed to identify potential amplicons from large collections of microbial genomes. It supports degenerate primers and provides accurate prediction of primer binding and amplification based on melting temperature and other PCR parameters. The tool now includes support for FASTQ input format, and can extract barcodes (or otherwise flanking regions) from sequence data while normalizing output for strand and trimming primers (à la cutadapt).

If you do not need the advanced features (HMM construction, FASTQ support, visualization, decoy analysis, etc.), we recommend you check out our optimized [AmpliconHunter2](https://github.com/rhowardstone/AmpliconHunter2) tool (~2-3x faster) to see if it meets your needs.

## Features

- **High Performance**: Optimized with Hyperscan pattern matching for processing large genome collections
- **Degenerate Primers**: Full support for IUPAC degenerate nucleotides
- **Temperature-based Filtering**: Uses BioPython's nearest-neighbor model for accurate melting temperature calculations
- **Off-target Detection**: Multiple methods to assess potential off-target amplification
- **FASTQ Support**: Process sequencing data directly with quality score preservation
- **Barcode Extraction**: Extract forward and reverse barcodes adjacent to primer sequences
- **Primer Trimming**: Option to remove primer sequences from amplicon output
- **Visualization**: Comprehensive plots for amplicon length, melting temperature, and orientation distribution
- **Taxonomic Analysis**: When taxonomy information is provided, generates species-level heatmaps of amplicon similarity
- **Ribotype Patterns**: Analysis of amplicon copy number patterns across genomes
- **HMM Profiling**: Optional hidden Markov model analysis for amplicon sequence homology

## Installation

### Requirements

- Python 3.7+
- Hyperscan
- BioPython
- NumPy, Pandas, Matplotlib, Seaborn
- Optional: HMMER suite (nhmmer, hmmbuild), MAFFT
- Optional: NCBI datasets tool (will be auto-downloaded if missing)

### Installation

```bash
git clone https://github.com/rhowardstone/AmpliconHunter.git
cd AmpliconHunter
pip install -e .
```

## Quick Start

### 1. Download RefSeq Database (Required for first use)

```bash
ampliconhunter download-refseq --type complete
```

### 2. Prepare Your Input Files

**Genome list file** (e.g., `genome_list.txt`):
```
/path/to/genome1.fa
/path/to/genome2.fa
/path/to/genome3.fa
```

**Primer file** (e.g., `primers.txt`):
```
AGRGTTYGATYMTGGCTCAG	RGYTACCTTGTTACGACTT
```

### 3. Run Analysis

```bash
ampliconhunter run genome_list.txt primers.txt results_dir --plots
```

## Commands

AmpliconHunter provides three main commands:

### 1. download-refseq - Download Reference Genomes

Downloads and prepares the NCBI RefSeq bacterial genome database for use with AmpliconHunter.

```
usage: ampliconhunter download-refseq [-h] [--type {complete,all}] [--timeout TIMEOUT]

optional arguments:
  -h, --help            show this help message and exit
  --type {complete,all} Type of RefSeq database to download (default: complete)
                        - complete: Only complete reference genome assemblies
                        - all: All reference genomes including incomplete assemblies
  --timeout TIMEOUT     Maximum download time in hours (default: 999)
```

Example:
```bash
# Download only complete genomes (recommended)
ampliconhunter download-refseq --type complete

# Download all reference genomes (much larger)
ampliconhunter download-refseq --type all --timeout 24
```

### 2. convert - Convert FASTA to Uppercase

Utility command to prepare FASTA files by converting all sequence lines to uppercase while preserving headers. This is useful when working with databases where sequences might use lowercase letters for masking or low-complexity regions.

```
usage: ampliconhunter convert [-h] [--no-recursive] [--threads THREADS] directory

positional arguments:
  directory            Directory containing FASTA files to process

optional arguments:
  -h, --help          show this help message and exit
  --no-recursive      Do not process subdirectories recursively
  --threads THREADS   Number of parallel threads to use (default: auto-detect)
```

Example:
```bash
# Convert all FASTA files in a directory tree
ampliconhunter convert /path/to/genomes/

# Convert only files in the specified directory (not subdirectories)
ampliconhunter convert --no-recursive /path/to/genomes/

# Use specific number of threads
ampliconhunter convert --threads 8 /path/to/genomes/
```

### 3. run - Perform In-Silico PCR

Main command to run amplicon prediction analysis.

```
usage: ampliconhunter run [-h] [--threads THREADS] [--Tm TM]
                         [--mismatches MISMATCHES] [--clamp CLAMP] 
                         [--Lmin LMIN] [--Lmax LMAX] [--decoy]
                         [--hmm [HMM]] [--clobber] [--dnac1 DNAC1]
                         [--dnac2 DNAC2] [--Na NA] [--Tris TRIS] [--Mg MG]
                         [--dNTPs DNTPS] [--saltcorr SALTCORR]
                         [--taxonomy TAXONOMY] [--plots] [--timeout TIMEOUT]
                         [--fb-len FB_LEN] [--rb-len RB_LEN]
                         [--include-offtarget] [--trim-primers]
                         [--input-fq]
                         input_file primer_file output_directory

positional arguments:
  input_file            Text file containing absolute paths to input genome files
  primer_file           TSV file containing forward and reverse primer sequences
  output_directory      Output directory for results

optional arguments:
  -h, --help            show this help message and exit
  --threads THREADS     Number of threads to use (default: all processors)
  --Tm TM               Minimum melting temperature threshold (°C) (default: do not filter)
  --mismatches MISMATCHES
                        Number of allowed mismatches in primer binding (default: 0)
  --clamp CLAMP         Number of 3'-most bases that must match perfectly (default: 5)
  --Lmin LMIN           Minimum amplicon length in bp (default: 50)
  --Lmax LMAX           Maximum amplicon length in bp (default: 5000)
  --decoy               Also search reversed genome sequences as negative control
  --hmm [HMM]           Enable HMM analysis. Without argument: build HMM from RefSeq.
                        With FILE argument: use existing HMM file
  --clobber             Overwrite existing output directory
  --dnac1 DNAC1         Primer concentration in nM (default: 1000)
  --dnac2 DNAC2         Template concentration in nM (default: 25)
  --Na NA               Sodium ion concentration in mM (default: 50)
  --Tris TRIS           Tris buffer concentration in mM (default: 0)
  --Mg MG               Magnesium ion concentration in mM (default: 4)
  --dNTPs DNTPS         Total dNTP concentration in mM (default: 1.6)
  --saltcorr SALTCORR   Salt correction method 0-5 (default: 5)
  --taxonomy TAXONOMY   Path to taxonomy TSV file for taxonomic analysis
  --plots               Generate visualization plots
  --timeout TIMEOUT     Maximum execution time in hours
  --fb-len FB_LEN       Forward barcode length to extract (default: 0)
  --rb-len RB_LEN       Reverse barcode length to extract (default: 0)
  --include-offtarget   Include off-target (FF/RR) amplicons in output
  --trim-primers        Remove primer sequences from amplicons
  --input-fq            Input files are in FASTQ format (default: FASTA)
```

## Input File Formats

### Genome List File
A text file with one genome file path per line:
```
/data/genomes/GCF_000005845.2.fa
/data/genomes/GCF_000009045.1.fa
/data/genomes/GCF_000013425.1.fa
```

### Primer File
Tab-separated file with forward and reverse primer sequences:
```
AGRGTTYGATYMTGGCTCAG	RGYTACCTTGTTACGACTT
```

### Taxonomy File (Optional)
Tab-separated file with taxonomic information:
```
accession	superkingdom	phylum	class	order	family	genus	species
GCF_000005845.2	Bacteria	Pseudomonadota	Gammaproteobacteria	Enterobacterales	Enterobacteriaceae	Escherichia	Escherichia coli
```

## Output Files

AmpliconHunter generates comprehensive output:

### Core Output Files
- `amplicons.fa` or `amplicons.fq`: Predicted amplicons in FASTA/FASTQ format
- `parameters.json`: Complete record of run parameters
- `run_statistics.json`: Detailed statistics including:
  - Total genomes processed
  - Number of genomes with amplicons
  - Average amplicons per genome
  - Ribotype diversity metrics
  - Off-target rates
- `ribotype_patterns.tsv`: Amplicon copy number patterns per genome

### When Using --decoy Option
- `decoy_amplicons.fa` or `decoy_amplicons.fq`: Amplicons from reversed genomes
- Additional decoy-specific visualization plots

### When Using --hmm Option
- `refseq/`: Directory containing HMM build process files
- `amplicons_hmm_scores.txt`: HMM scores for each amplicon
- HMM score distribution plots

### When Using --plots Option
- `plots/`: Directory containing all visualizations:
  - `length_distribution.png`: Amplicon length distributions by orientation
  - `temp_distribution.png`: Melting temperature distributions
  - `orientation_distribution.png`: Pie chart of FR/RF/FF/RR orientations
  - `hmm_score_distribution.png`: HMM bit scores (if --hmm used)
  - `hmm_evalue_distribution.png`: HMM E-values (if --hmm used)
  - `species_heatmaps/`: Jaccard similarity heatmaps for each genus

## FASTA Header Format

AmpliconHunter uses a structured header format for output sequences:

```
>sequence_id.source=genome_file.coordinates=start-end.orientation=XX.Tm=temp[.fb=barcode][.rb=barcode]
```

Components:
- `sequence_id`: Original sequence identifier from input
- `source`: Source genome filename
- `coordinates`: Start-end positions in source sequence
- `orientation`: Primer orientation (FR, RF, FF, or RR)
- `Tm`: Minimum melting temperature of primer pair
- `fb`: Forward barcode sequence (if --fb-len specified)
- `rb`: Reverse barcode sequence (if --rb-len specified)

## Advanced Features

### FASTQ Support
Process sequencing data directly while preserving quality scores:
```bash
ampliconhunter run --input-fq --trim-primers \
  fastq_files.txt primers.txt results_fastq
```

### Barcode Extraction
Extract barcodes adjacent to primer sequences:
```bash
# Extract 8bp forward and 6bp reverse barcodes
ampliconhunter run --fb-len 8 --rb-len 6 \
  genomes.txt primers.txt results_barcoded
```

### Off-Target Analysis
Include potential off-target amplicons (FF and RR orientations):
```bash
ampliconhunter run --include-offtarget --decoy \
  genomes.txt primers.txt results_offtarget
```

### HMM Profile Analysis
Build and use HMM profiles for homology assessment:
```bash
# Build HMM from RefSeq and apply to analysis
ampliconhunter run --hmm --plots \
  genomes.txt primers.txt results_hmm

# Use existing HMM file
ampliconhunter run --hmm my_profile.hmm \
  genomes.txt primers.txt results_hmm2
```

### Complex Example
Full analysis with all features:
```bash
ampliconhunter run \
  --threads 16 \
  --Tm 58 \
  --mismatches 2 \
  --clamp 3 \
  --Lmin 400 \
  --Lmax 600 \
  --trim-primers \
  --fb-len 10 \
  --rb-len 10 \
  --include-offtarget \
  --decoy \
  --hmm \
  --taxonomy taxonomy.tsv \
  --plots \
  --timeout 12 \
  genome_list.txt \
  primers.txt \
  comprehensive_results
```


### IUPAC Degenerate Base Codes

AmpliconHunter supports all IUPAC degenerate nucleotide codes:

| Code | Meaning      | Base(s)  |
|------|--------------|----------|
| A    | Adenine      | A        |
| C    | Cytosine     | C        |
| G    | Guanine      | G        |
| T    | Thymine      | T        |
| R    | puRine       | A, G     |
| Y    | pYrimidine   | C, T     |
| S    | Strong bond  | G, C     |
| W    | Weak bond    | A, T     |
| K    | Keto         | G, T     |
| M    | aMino        | A, C     |
| B    | not A        | C, G, T  |
| D    | not C        | A, G, T  |
| H    | not G        | A, C, T  |
| V    | not T/U      | A, C, G  |
| N    | aNy base     | A, C, G, T|

## Performance Tips

1. **Multi-threading**: Use `--threads` to parallelize across genomes
2. **Memory**: Each thread requires ~2GB RAM for typical bacterial genomes
3. **Storage**: RefSeq complete requires ~15GB, all requires ~150GB
4. **HMM Caching**: HMM profiles are automatically cached for reuse
5. **Timeout**: Use `--timeout` for very large datasets to prevent indefinite runs

## Troubleshooting

### Missing Dependencies
If you see warnings about missing tools (nhmmer, mafft, hmmbuild):
- These are only required for `--hmm` functionality
- Install HMMER: `conda install -c bioconda hmmer`
- Install MAFFT: `conda install -c bioconda mafft`

### NCBI Datasets Tool
The tool will attempt to auto-download the NCBI datasets tool if missing. If this fails:
- Download manually from: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/
- Place in your PATH or in `~/.ampliconhunter/bin/`

### Memory Issues
For large datasets, consider:
- Reducing `--threads` to lower memory usage
- Processing genomes in batches
- Using `--timeout` to limit runtime

## Citation

If you use AmpliconHunter in your research, please cite:

Howard-Stone, R. & Mandoiu, I.I. (2025). AmpliconHunter: A Scalable Tool for PCR Amplicon Prediction from Microbiome Samples. In Proceedings of the 15th IEEE International Conference on Computational Advances in Bio and Medical Sciences (ICCABS 2025).

## License

AmpliconHunter is licensed under the MIT License. See LICENSE for details.

## Contact

For questions and support, please open an issue on the [GitHub repository](https://github.com/rhowardstone/AmpliconHunter).

