DNA_NUCLEOTIDES = set("ATCGatcg")
RNA_NUCLEOTIDES = set("AUCGaucg")
ALL_NUCLEOTIDES = DNA_NUCLEOTIDES.union(RNA_NUCLEOTIDES)


def is_nucleic_acid(sequence):
    """
    Check if a sequence is a valid DNA or RNA sequence.
    Returns True if all characters are valid nucleic acid nucleotides.
    """
    set_seq = set(sequence)
    is_dna = set_seq.issubset(DNA_NUCLEOTIDES)
    is_rna = set_seq.issubset(RNA_NUCLEOTIDES)
    return is_dna or is_rna


def transcribe(sequence):
    """
    Transcribe between DNA and RNA sequences.
    DNA to RNA: T → U, RNA to DNA: U → T
    """
    trans_map = {"T": "U", "U": "T", "t": "u", "u": "t"}
    return "".join(trans_map.get(nucl, nucl) for nucl in sequence)


def reverse(sequence):
    """Reverse the input sequence."""
    return sequence[::-1]


def complement(sequence):
    """
    Return the complementary sequence for DNA or RNA.
    DNA: A↔T, C↔G, RNA: A↔U, C↔G
    """
    set_seq = set(sequence)

    # Determine if sequence is DNA or RNA
    if set_seq.issubset(DNA_NUCLEOTIDES):
        comp_map = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
            "a": "t",
            "t": "a",
            "g": "c",
            "c": "g",
        }
    else:
        comp_map = {
            "A": "U",
            "U": "A",
            "G": "C",
            "C": "G",
            "a": "u",
            "u": "a",
            "g": "c",
            "c": "g",
        }

    return "".join(comp_map.get(char, char) for char in sequence)


def reverse_complement(sequence):
    """Return the reverse complement of a DNA or RNA sequence."""
    comp = complement(sequence)
    return reverse(comp)


def run_dna_rna_tools(*args):
    """
    Main function to run DNA/RNA tools operations.

    Usage examples:
    run_dna_rna_tools("ATCG", "reverse")
    run_dna_rna_tools("seq1", "seq2", "complement")
    """
    if not args:
        raise ValueError("No arguments provided")

    *sequences, operation = args

    function_map = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    if operation not in function_map:
        raise ValueError(f"Unknown operation: {operation}")

    selected_function = function_map[operation]

    # Single loop for validation and execution
    results = []
    for seq in sequences:
        if operation != "is_nucleic_acid" and not is_nucleic_acid(seq):
            raise ValueError(f"Sequence '{seq}' is not valid DNA/RNA")

        result = selected_function(seq)
        results.append(result)

    return results[0] if len(results) == 1 else results
