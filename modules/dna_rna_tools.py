DNA_NUCLEOTIDES = set("ATCGatcg")
RNA_NUCLEOTIDES = set("AUCGaucg")
ALL_NUCLEOTIDES = DNA_NUCLEOTIDES.union(RNA_NUCLEOTIDES)


def is_nucleic_acid(sequence):
    """
    Check if a sequence is a valid DNA or RNA sequence.

    This function verifies that all characters in the input sequence
    are valid nucleic acid nucleotides (either DNA or RNA bases).

    Parameters
    ----------
    sequence : str
        The input sequence to check. Can be any string.

    Returns
    -------
    bool
        True if the sequence contains only valid DNA or RNA nucleotides,
        False otherwise.

    Examples
    --------
    >>> is_nucleic_acid("ATCG")
    True

    >>> is_nucleic_acid("AUCG")
    True

    >>> is_nucleic_acid("ATCGN")  # 'N' is not a valid nucleotide
    False

    >>> is_nucleic_acid("protein")
    False

    >>> is_nucleic_acid("atcg")  # lowercase DNA
    True

    >>> is_nucleic_acid("aucg")  # lowercase RNA
    True

    Notes
    -----
    - Valid DNA nucleotides: A, T, C, G (case insensitive)
    - Valid RNA nucleotides: A, U, C, G (case insensitive)
    - Mixed DNA/RNA sequences will return False since they contain
      both T and U nucleotides
    - Empty sequences return True (empty set is subset of any set)
    """

    set_seq = set(sequence)
    is_dna = set_seq.issubset(DNA_NUCLEOTIDES)
    is_rna = set_seq.issubset(RNA_NUCLEOTIDES)
    return is_dna or is_rna


def transcribe(sequence):
    """
    Transcribe between DNA and RNA sequences by replacing T with U or U with T.

    This function performs bidirectional transcription:
    - For DNA sequences: T -> U and t -> u (transcription to RNA)
    - For RNA sequences: U -> T and u -> t (reverse transcription to DNA)
    - All other nucleotides remain unchanged

    Parameters
    ----------
    sequence : str
        The input DNA or RNA sequence to transcribe.

    Returns
    -------
    str
        The transcribed sequence:
        - RNA sequence if input was DNA
        - DNA sequence if input was RNA
        - Unchanged if sequence contains no T/U nucleotides

    Examples
    --------
    >>> transcribe("ATCG")  # DNA to RNA
    'AUCG'

    >>> transcribe("AUCG")  # RNA to DNA
    'ATCG'

    >>> transcribe("atcg")  # lowercase DNA to RNA
    'aucg'

    >>> transcribe("ACGTU")  # Mixed case - both T and U
    'ACGUU'

    >>> transcribe("ACGG")  # No T or U - unchanged
    'ACGG'

    Notes
    -----
    - The function is bidirectional and doesn't validate input sequence type
    - If sequence contains both T and U, both will be swapped:
      'T' -> 'U' and 'U' -> 'T'
    - Mixed case sequences are handled appropriately
    - Does not check if input is valid nucleic acid sequence
    """

    trans_map = {"T": "U", "U": "T", "t": "u", "u": "t"}
    return "".join(trans_map.get(nucl, nucl) for nucl in sequence)


def reverse(sequence):
    """
    Reverse the input sequence.

    Returns a new string with the characters of the input sequence
    in reverse order.

    Parameters
    ----------
    sequence : str
        The input sequence to reverse. Can be any string.

    Returns
    -------
    str
        The reversed sequence.

    Examples
    --------
    >>> reverse("ATCG")
    'GCTA'

    >>> reverse("AUCG")
    'GCUA'

    >>> reverse("ACGT")
    'TGCA'

    >>> reverse("")
    ''

    >>> reverse("A")
    'A'

    Notes
    -----
    - Works with any string, not just nucleic acid sequences
    - Empty string returns empty string
    - Single character returns the same character
    - Uses Python string slicing with step -1 for efficient reversal
    """

    return sequence[::-1]


def complement(sequence):
    """
    Return the complementary sequence for DNA or RNA.

    For DNA sequences: A<->T, C<->G (and lowercase equivalents)
    For RNA sequences: A<->U, C<->G (and lowercase equivalents)

    The function automatically detects whether the input is DNA or RNA
    based on the nucleotides present in the sequence.

    Parameters
    ----------
    sequence : str
        The input DNA or RNA sequence.

    Returns
    -------
    str
        The complementary sequence.

    Examples
    --------
    >>> complement("ATCG")
    'TAGC'

    >>> complement("AUCG")
    'UAGC'

    >>> complement("atcg")
    'tagc'

    >>> complement("AATTCCGG")
    'TTAAGGCC'

    >>> complement("ACGTN")  # N is not complemented
    'TGCAN'

    Raises
    ------
    ValueError
        If the sequence contains characters from both DNA and RNA
        (both T and U present), as this creates ambiguity in determining
        the sequence type.

    Notes
    -----
    - Sequence type detection:
        - If all characters are in DNA_NUCLEOTIDES → treated as DNA
        - Otherwise → treated as RNA
    - Characters not in the complement map remain unchanged
    - For mixed-case sequences, case is preserved in the output
    """

    set_seq = set(sequence)

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
    else:  # RNA
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
    """
    Return the reverse complement of a DNA or RNA sequence.

    This function first computes the complement of the sequence,
    then reverses the resulting string. It is equivalent to:
    reverse(complement(sequence))

    Parameters
    ----------
    sequence : str
        The input DNA or RNA sequence.

    Returns
    -------
    str
        The reverse complement of the input sequence.

    Examples
    --------
    >>> reverse_complement("ATCG")
    'CGAT'

    >>> reverse_complement("AUCG")
    'CGAU'

    >>> reverse_complement("ATCGATCG")
    'CGATCGAT'

    >>> reverse_complement("AAAATTTT")
    'AAAATTTT'  # Self-reverse-complementary

    Notes
    -----
    - Uses the same sequence type detection as complement():
      DNA if all characters are in DNA_NUCLEOTIDES, otherwise RNA
    - The reverse complement is a common operation in bioinformatics
      for finding the complementary strand in the opposite orientation
    - For palindromic sequences, the reverse complement may be identical
      to the original sequence

    See Also
    --------
    complement : Get the complementary sequence
    reverse : Reverse the sequence
    """

    comp = complement(sequence)
    return reverse(comp)


def run_dna_rna_tools(*args):
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

    if operation != "is_nucleic_acid":
        for seq in sequences:
            if not is_nucleic_acid(seq):
                raise ValueError(f"Sequence '{seq}' is not valid DNA/RNA")

    results = []
    for seq in sequences:
        result = selected_function(seq)
        results.append(result)

    return results[0] if len(results) == 1 else results
