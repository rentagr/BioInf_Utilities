# BioInf_Utilities

BioInf_Utilities is a set of Python utilities for working with biological sequences.

## Project Structure
```python

BioInf_Utilities/

├── main.py              # Main script (entry point)

├── README.md            # Documentation

└── modules/             # Functionality modules

    ├── dna_rna_tools.py    # DNA/RNA utilities

    ├── bio_files_processor.py # Bioinformatic file format utilities

    └── fastq_utils.py      # FASTQ processing
``` 

## Installation

git clone (https://github.com/rentagr/BioInf_Utilities.git)

cd BioInf_Utilities

## Usage

## DNA/RNA Sequence Tools

#from modules.dna_rna_tools import run_dna_rna_tools

### Check if sequence is valid nucleic acid
result = run_dna_rna_tools("ATCG", "is_nucleic_acid")

### Transcribe DNA to RNA  
result = run_dna_rna_tools("ATCG", "transcribe")

### Get reverse complement
result = run_dna_rna_tools("ATCG", "reverse_complement")

### Multiple sequences at once
results = run_dna_rna_tools("ATCG", "AUCG", "reverse")

## FASTQ Processing
## For file-based processing:

from modules.fastq_utils import filter_fastq_file

### Filter FASTQ files by GC content, length and quality
```python
filter_fastq_file(
    "input.fastq",
    "filtered/output.fastq", 
    gc_bounds=(40, 60),
    length_bounds=(50, 150),
    quality_threshold=20
)
```
## For dictionary-based processing:
```python
from modules.fastq_utils import filter_fastq_dict

sequences = {
    'read1': ('ATCG', 'IIII'),
    'read2': ('GGCC', '!!!!')
}

filtered = filter_fastq_dict(
    sequences,
    gc_bounds=(20, 80),
    length_bounds=(50, 150), 
    quality_threshold=20
)
```
## Bioinformatics File Processing

from bio_files_processor import convert_multiline_fasta_to_oneline, parse_blast_output

### Convert multi-line FASTA to single-line format
convert_multiline_fasta_to_oneline("input.fasta", "output.fasta")

### Parse BLAST results and extract best matches
parse_blast_output("blast_results.txt", "best_matches.txt")

## Available Functions
## DNA/RNA Tools Module

- is_nucleic_acid() - Validate DNA/RNA sequences
- transcribe() - Transcribe between DNA and RNA  
- reverse() - Reverse sequences
- complement() - Get complementary sequences
- reverse_complement() - Get reverse complement
- run_dna_rna_tools() - Main interface for DNA/RNA operations

## FASTQ Utilities Module

- validate_fastq_record() - Validate FASTQ format and characters
- calculate_gc_content() - Calculate GC percentage
- calculate_average_quality() - Calculate average Phred quality
- filter_fastq_dict() - Filter sequence dictionaries by criteria
- filter_fastq_file() - Filter FASTQ files (reads and writes files)
- read_fastq_file() - Read FASTQ files into dictionary format
- write_fastq_file() - Write sequences to FASTQ files

## Bio Files Processor Module
- convert_multiline_fasta_to_oneline() - Convert multi-line FASTA to single-line
- parse_blast_output() - Extract best matches from BLAST results

## Documentation

Full function documentation is available through built-in docstrings:

```python
from main import filter_fastq
help(filter_fastq)
```
All functions have type annotations and docstrings.

## Dependencies

- Python 3.6+
- No external dependencies required

## Contacts

For questions and suggestions:

Email: renata.tag.isl@gmail.com

BioInf_Utilities was developed as part of an educational Bioinformatics Institute project