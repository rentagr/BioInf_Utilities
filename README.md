# BioInf_Utilities

BioInf_Utilities is a set of Python utilities for working with biological sequences.

## Project Structure
```python

BioInf_Utilities/

├── main.py              # Main script (entry point)

├── README.md            # Documentation

└── modules/             # Functionality modules

    ├── dna_rna_tools.py    # DNA/RNA utilities

    └── fastq_utils.py      # FASTQ processing
``` 

## Installation

git clone (https://github.com/rentagr/BioInf_Utilities.git)

cd BioInf_Utilities

## Usage

### DNA/RNA Tools

from main import run_dna_rna_tools

### Check if sequence is valid nucleic acid
result = run_dna_rna_tools("ATCG", "is_nucleic_acid")

### Transcribe DNA to RNA  
result = run_dna_rna_tools("ATCG", "transcribe")

### Get reverse complement
result = run_dna_rna_tools("ATCG", "reverse_complement")

### FASTQ Processing

from main import filter_fastq

sequences = {
    'read1': ('ATCG', 'IIII'),
    'read2': ('GGCC', '!!!!')
}

### Filter by GC content and quality
filtered = filter_fastq(
    sequences,
    gc_bounds=(20, 80),
    length_bounds=(50, 150),
    quality_threshold=20
)

## Available Functions

## DNA/RNA Tools Module

- is_nucleic_acid() - Validate DNA/RNA sequences
- transcribe() - Transcribe between DNA and RNA  
- reverse() - Reverse sequences
- complement() - Get complementary sequences
- reverse_complement() - Get reverse complement
- run_dna_rna_tools() - Main interface for DNA/RNA operations

## FASTQ Utilities Module

- validate_fastq_record() - Validate FASTQ records
- calculate_gc_content() - Calculate GC percentage
- calculate_average_quality() - Calculate average Phred quality
- filter_fastq() - Filter sequences by multiple criteria

## Documentation

Full function documentation is available through built-in docstrings:

from main import filter_fastq
help(filter_fastq)

All functions have type annotations and detailed docstrings with usage examples.

## Dependencies

- Python 3.6+
- No external dependencies required

## Contacts

For questions and suggestions:

Email: renata.tag.isl@gmail.com

BioInf_Utilities was developed as part of an educational Bioinformatics Institute project