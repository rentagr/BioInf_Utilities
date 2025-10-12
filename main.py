"""
Main script for BioInf_Utilities bioinformatics tools
"""

# Import main functions from modules
from modules.dna_rna_tools import run_dna_rna_tools
from modules.fastq_utils import filter_fastq_file, filter_fastq_dict
from bio_files_processor import convert_multiline_fasta_to_oneline, parse_blast_output


def main():
    """
    Main function that calls all required bioinformatics functions.
    """
    # Call DNA/RNA tools function
    run_dna_rna_tools("ATCG", "reverse_complement")

    # Call FASTQ filtering functions
    filter_fastq_file("input.fastq", "output.fastq")

    example_sequences = {"seq1": ("ATCG", "IIII")}
    filter_fastq_dict(example_sequences)

    # Call file processing functions
    convert_multiline_fasta_to_oneline("input.fasta", "output.fasta")
    parse_blast_output("blast_results.txt", "best_matches.txt")


if __name__ == "__main__":
    main()
