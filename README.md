# AmpliconHunter

A scalable tool for accurate PCR amplicon prediction from microbiome samples using degenerate primers.

## Overview

AmpliconHunter is a high-performance, in-silico PCR tool designed to identify potential amplicons from large collections of microbial genomes. It supports degenerate primers and provides accurate prediction of primer binding and amplification based on melting temperature and other PCR parameters.

## Features

- **High Performance**: Optimized with Hyperscan pattern matching for processing large genome collections
- **Degenerate Primers**: Full support for IUPAC degenerate nucleotides
- **Temperature-based Filtering**: Uses BioPython's nearest-neighbor model for accurate melting temperature calculations
- **Off-target Detection**: Multiple methods to assess potential off-target amplification
- **Visualization**: Comprehensive plots for amplicon length, melting temperature, and orientation distribution
- **Taxonomic Analysis**: When taxonomy information is provided, generates species-level heatmaps of amplicon similarity
- **Amplitype Patterns**: Analysis of amplicon copy number patterns across genomes ("amplitypes")
- **HMM Profiling**: Optional hidden Markov model analysis for amplicon sequence homology

## Installation

### Requirements

- Python 3.7+
- Hyperscan
- BioPython
- NumPy, Pandas, Matplotlib, Seaborn
- Optional: HMMER suite (nhmmer, hmmbuild), MAFFT

### Installation via pip

```bash
pip install ampliconhunter
```

### Manual Installation

```bash
git clone https://github.com/rhowardstone/AmpliconHunter.git
cd AmpliconHunter
pip install -e .
```

## Quick Start

### Download RefSeq Database

```bash
ampliconhunter download-refseq --type complete
```

### Run Analysis

```bash
ampliconhunter run genome_list.txt primers.txt results_dir
```

Example of a primer file:
```
AGRGTTYGATYMTGGCTCAG	RGYTACCTTGTTACGACTT
```

## Usage

```
usage: ampliconhunter [-h] [--base-dir BASE_DIR] {download-refseq,run} ...

AmpliconHunter: In-Silico PCR sequence extractor

optional arguments:
  -h, --help            show this help message and exit
  --base-dir BASE_DIR   Base directory for AmpliconHunter files (default: ~/.ampliconhunter)

Commands:
  {download-refseq,run}  Command to execute
    download-refseq     Download RefSeq database
    run                 Run in-silico PCR
```

### Run Command Options

```
usage: ampliconhunter run [-h] [--threads THREADS] [--Tm TM]
                         [--mismatches MISMATCHES] [--clamp CLAMP] 
                         [--Lmin LMIN] [--Lmax LMAX] [--decoy]
                         [--hmm [HMM]] [--clobber] [--dnac1 DNAC1]
                         [--dnac2 DNAC2] [--Na NA] [--Tris TRIS] [--Mg MG]
                         [--dNTPs DNTPS] [--saltcorr SALTCORR]
                         [--taxonomy TAXONOMY] [--no-plots] 
                         [--timeout TIMEOUT]
                         input_file primer_file output_directory

positional arguments:
  input_file            Text file containing absolute paths to input FASTA files
  primer_file           TSV file containing primer information
  output_directory      Output directory for results

optional arguments:
  -h, --help            show this help message and exit
  --threads THREADS     Number of threads to use (default: all processors)
  --Tm TM               Minimum melting temperature threshold
  --mismatches MISMATCHES
                        Number of allowed mismatches (default: 0)
  --clamp CLAMP         3'-most CLAMP bases are not allowed mismatches (default: 5)
  --Lmin LMIN           Minimum length filter (default: 50)
  --Lmax LMAX           Maximum length filter (default: 5000)
  --decoy               Use reversed genomes as decoy sequences (default: false)
  --hmm [HMM]           HMM processing mode (without FILE: build and use new HMM; with FILE: use existing HMM file)
  --clobber             Overwrite pre-existing run in output directory
  --dnac1 DNAC1         Concentration of primer strand [nM] (default: 1000)
  --dnac2 DNAC2         Concentration of template strand [nM] (default: 25)
  --Na NA               Sodium ion concentration [mM] (default: 50)
  --Tris TRIS           Tris buffer concentration [mM] (default: 0)
  --Mg MG               Magnesium ion concentration [mM] (default: 4)
  --dNTPs DNTPS         Total deoxynucleotide concentration [mM] (default: 1.6)
  --saltcorr SALTCORR   Salt correction method 0-5 (default: 5)
  --taxonomy TAXONOMY   Path to taxonomy mapping file
  --no-plots            Skip generating visualization plots
  --timeout TIMEOUT     Maximum execution time in hours
```

## Output Files

AmpliconHunter generates several output files:

- `amplicons.fa`: FASTA file containing all predicted amplicons
- `parameters.json`: JSON file with all run parameters
- `run_statistics.json`: JSON file with detailed statistics about the run
- `plots/`: Directory containing visualization plots:
  - `length_distribution.png`: Distribution of amplicon lengths by orientation
  - `temp_distribution.png`: Distribution of melting temperatures by orientation
  - `orientation_distribution.png`: Pie chart of amplicon orientations
  - `ribotype_patterns.tsv`: TSV file with amplicon patterns for each genome
- When using `--decoy` option:
  - `decoy_amplicons.fa`: FASTA file containing predicted amplicons from reversed genomes
  - Additional decoy-specific plots

## Examples

### Basic Analysis with Default Parameters

```bash
ampliconhunter run refseq_complete_genomes_files.txt primers.txt results
```

### Analysis with Advanced Parameters

```bash
ampliconhunter run --threads 16 --Tm 55 --mismatches 2 --clamp 3 \
  --Lmin 100 --Lmax 2000 --taxonomy taxonomy.tsv \
  refseq_all_genomes_files.txt primers.txt results
```

### Building and Using an HMM Profile

```bash
ampliconhunter run --hmm --threads 16 --mismatches 2 \
  --Lmin 1500 --Lmax 2000 --decoy \
  GTDB_genomes_files.txt primers.txt results
```

## Taxonomy Support

AmpliconHunter can use taxonomy information to analyze amplicon patterns across taxonomic groups. The expected format for the taxonomy file is:

```
accession	superkingdom	phylum	class	order	family	genus	species
GCF_000022305.1	Bacteria	Pseudomonadota	Betaproteobacteria	Burkholderiales	Comamonadaceae	Diaphorobacter	[Acidovorax] ebreus TPSY
```

## Advanced Usage

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

### Primer Orientations

AmpliconHunter reports four possible primer orientations in amplicons:

- **FR**: Forward primer followed by reverse complement of reverse primer (standard)
- **RF**: Reverse primer followed by reverse complement of forward primer (reversed)
- **FF**: Forward primer followed by reverse complement of forward primer
- **RR**: Reverse primer followed by reverse complement of reverse primer

FR and RF represent typical amplicons, while FF and RR represent potential off-target amplification.

## Citation

If you use AmpliconHunter in your research, please cite:

```
Howard-Stone, R. & Mandoiu, I.I. (2025). AmpliconHunter: A Scalable Tool for PCR Amplicon Prediction from Microbiome Samples. Journal of Example, 0(0), 0-0.
```

## License

AmpliconHunter is licensed under the MIT License. See LICENSE for details.

## Contact

For questions and support, please open an issue on the [GitHub repository](https://github.com/rhowardstone/AmpliconHunter).