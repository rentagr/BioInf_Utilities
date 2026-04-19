
# BioInf_Utilities

[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Linting: flake8](https://img.shields.io/badge/linting-flake8-yellowgreen)](https://flake8.pycqa.org/)

BioInf_Utilities is a Python package for working with biological sequences (DNA, RNA, proteins) and filtering FASTQ files.  
The project is built using **Object-Oriented Programming (OOP)** principles and leverages **Biopython** for efficient FASTQ processing.




# BioInf_Utilities

BioInf_Utilities is a set of Python utilities for working with biological sequences.

## Project Structure
```python
BioInf_Utilities/
├── bio_utilities_oop.py # Main module with OOP classes and FASTQ filter
├── tests/ # Unit tests (pytest)
│ └── test_fastq.py # 8 tests for filter_fastq
├── README.md # Documentation
├── requirements.txt # Dependencies 
``` 

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/rentagr/BioInf_Utilities.git
   cd BioInf_Utilities
   ``` 
2. Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```
## Usage
### Working with Biological Sequences

All sequence classes supports:
- `len()` - get sequence length
- Indexing and slicing (e.g., `seq[2]`, `seq[1:5]`) - returns an object of the same class
- `__str__` and `__repr__` for pretty printing
- Alphabet validation at initialization

#### DNA Sequences
```python
from bio_utilities_oop import DNASequence

dna = DNASequence("ATGC")
print(dna)                       # ATGC
print(dna.complement())           # TACG
print(dna.reverse())              # CGTA
print(dna.reverse_complement())   # GCAT
print(dna.transcribe())           # AUGC (returns RNASequence object)
```
## RNA Sequences
```python
from bio_utilities_oop import RNASequence

rna = RNASequence("AUCG")
print(rna.complement())           # UAGC
print(rna.reverse_complement())   # CGAU
```
## Protein Sequences
```python
from bio_utilities_oop import AminoAcidSequence

protein = AminoAcidSequence("MVR")
print(protein.molecular_weight()) # Approx. 362.4 g/mol
```
## FASTQ Filtering with Biopython – Programmatic Use

The `filter_fastq` function reads a FASTQ file, filters reads by length, average quality, and GC content, and writes the passing reads to a new file.

```python
from bio_utilities_oop import filter_fastq

filter_fastq(
    input_fastq="reads.fastq",
    output_fastq="filtered.fastq",
    length_bounds=(50, 150),       # keep reads between 50 and 150 bp
    quality_threshold=20,           # minimum average Phred score
    gc_bounds=(40, 60)              # keep reads with 40–60% GC content
)
```

## FASTQ Filtering – Command‑Line Interface (CLI)
You can run the FASTQ filter directly from the terminal:

```bash
python bio_utilities_oop.py input.fastq -o filtered.fastq -l 50 -L 150 -q 20 -g 40 -G 60 --log mylog.txt
```
Arguments:

input – positional, input FASTQ file.

-o, --output – output FASTQ file (default: filtered.fastq).

-l, --len-min – minimum read length (default: 0).

-L, --len-max – maximum read length (default: inf).

-q, --qual – minimum average Phred quality (default: 0).

-g, --gc-min – minimum GC% (default: 0).

-G, --gc-max – maximum GC% (default: 100).

--log – log file (default: filter.log).

#### Logging
If `log_file` is provided, the function writes:

Start of filtering with parameters

Number of records processed and saved

Any errors (e.g., missing file, parse error) – both to the log file and the console.

## Testing
The package includes 8 unit tests for filter_fastq using pytest.

To run all tests:

```bash
pytest tests/ -v
```
Test coverage:

- File I/O (read/write)

- Error handling (missing input file)

- Length filter

- GC% filter

- Quality filter

- Combined filters

- Empty output when nothing passes

- Log file creation and content

All tests pass under Python 3.6+ with Biopython and pytest installed.

### Code Quality

This project follows PEP 8 guidelines and uses modern Python tooling:

- **Formatting:** [Black](https://github.com/psf/black) - ensures consistent code style
- **Linting:** [flake8](https://flake8.pycqa.org/) - checks for PEP 8 compliance and potential errors
- **Testing:** pytest - unit tests for the FASTQ filter.

### Dependencies
- **Python 3.6+**
- **Biopython** >= 1.79 (for FASTQ processing)
- **pytest** (for running tests – optional, only needed for development)


All dependencies are listed in requirements.txt.

### Documentation
Full docstrings are provided for all classes and functions. 

- Use Python's built-in [help()]

```python
from bio_utilities_oop import DNASequence, filter_fastq
help(DNASequence)
help(filter_fastq)
```
##

BioInf_Utilities was developed as part of an educational Bioinformatics Institute project