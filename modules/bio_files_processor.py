"""
Bio Files Processor - utilities for working with bioinformatics file formats.
Contains functions for FASTA conversion and BLAST results parsing.
"""

import os
from pathlib import Path
import re


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str = None
) -> None:
    """
    Convert FASTA file with multi-line sequences to single-line format.

    Many FASTA files break sequences across multiple lines for readability.
    This function converts them to single-line sequences for easier processing.

    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path for output file (optional - auto-generated if not provided)
    """
    # Generate output filename if not provided
    if output_fasta is None:
        input_path = Path(input_fasta)
        output_fasta = str(
            input_path.parent / f"{input_path.stem}_oneline{input_path.suffix}"
        )

    print(f"Converting {input_fasta} to single-line format...")

    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        current_header = None
        current_sequence = []

        for line in infile:
            line = line.strip()

            if line.startswith(">"):  # Header line
                # Save previous sequence if exists
                if current_header is not None:
                    outfile.write(f"{current_header}\n")
                    outfile.write(f"{''.join(current_sequence)}\n\n")

                # Start new sequence
                current_header = line
                current_sequence = []
            else:  # Sequence line
                current_sequence.append(line)

        # Don't forget the last sequence
        if current_header is not None:
            outfile.write(f"{current_header}\n")
            outfile.write(f"{''.join(current_sequence)}\n")

    print(f"Conversion complete! Output: {output_fasta}")


def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Parse BLAST results and extract best matches for each query.

    BLAST output files contain alignment results. This function extracts
    the best (first) description for each query sequence.

    Args:
        input_file: Path to BLAST results file
        output_file: Path for output file with best matches
    """
    print(f"Parsing BLAST results from: {input_file}")

    descriptions = set()

    with open(input_file, "r") as file:
        content = file.read()

    # Find all blocks starting with "Sequences producing significant alignments"
    pattern = r"Sequences producing significant alignments:(.*?)(?=\n\n|\n\s*\n|$)"
    blocks = re.findall(pattern, content, re.DOTALL)

    print(f"Found {len(blocks)} query blocks")

    for block in blocks:
        lines = block.strip().split("\n")
        if lines:
            # First line contains the best match
            first_line = lines[0].strip()
            # Split by multiple spaces to get columns
            parts = re.split(r"\s{2,}", first_line)
            if len(parts) >= 2:
                description = parts[0].strip()
                descriptions.add(description)

    # Sort alphabetically and save
    sorted_descriptions = sorted(descriptions)

    with open(output_file, "w") as outfile:
        for desc in sorted_descriptions:
            outfile.write(f"{desc}\n")

    print(f"Extracted {len(sorted_descriptions)} unique descriptions")
    print(f"Results saved to: {output_file}")


def demonstrate_functions():
    """
    Demonstrate the functionality of bio files processor.
    Creates test files and shows how to use each function.
    """
    print("Bio Files Processor Demonstration")
    print("=" * 40)

    # Create test multiline FASTA file
    test_fasta = "test_multiline.fasta"
    if not os.path.exists(test_fasta):
        print("Creating test FASTA file...")
        with open(test_fasta, "w") as f:
            f.write(">sequence1\n")
            f.write("ATCG\n")
            f.write("ATCG\n")
            f.write(">sequence2\n")
            f.write("GGCC\n")
            f.write("GGCC\n")
        print(f"Created: {test_fasta}")

    # Demonstrate FASTA conversion
    print("\n1. Demonstrating FASTA conversion:")
    convert_multiline_fasta_to_oneline(test_fasta)

    # Create test BLAST results file
    test_blast = "test_blast_results.txt"
    if not os.path.exists(test_blast):
        print("\nCreating test BLAST results file...")
        with open(test_blast, "w") as f:
            f.write("Query: sequence1\n\n")
            f.write("Sequences producing significant alignments:\n")
            f.write("Protein_A  95% identity  E-value: 1e-50\n")
            f.write("Protein_B  90% identity  E-value: 1e-45\n\n")
            f.write("Query: sequence2\n\n")
            f.write("Sequences producing significant alignments:\n")
            f.write("Protein_C  98% identity  E-value: 1e-60\n")
            f.write("Protein_D  85% identity  E-value: 1e-40\n")
        print(f"Created: {test_blast}")

    # Demonstrate BLAST parsing
    print("\n2. Demonstrating BLAST results parsing:")
    parse_blast_output(test_blast, "best_matches.txt")

    print("\nDemonstration completed!")


if __name__ == "__main__":
    demonstrate_functions()
