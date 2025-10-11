# Validating constants
VALID_NUCLEOTIDES = set("ATCGNnatcgn")
PHRED33_MIN = 33  # Minimum ASCII code for Phred+33 ('!')
PHRED33_MAX = 126  # Maximum ASCII code for Phred+126 ('~')


def validate_fastq_record(sequence: str, quality: str) -> bool:
    """
    Validates a FASTQ record.

    Parameters
    ----------
    sequence : str
        Nucleotide sequence
    quality : str
        Quality string in Phred+33 format

    Returns
    -------
    bool
        True if the record is valid, False otherwise

    Examples
    --------
    >>> validate_fastq_record("ATCG", "IIII")
    True
    >>> validate_fastq_record("ATCG", "III")  # different length
    False
    >>> validate_fastq_record("ATCX", "IIII")  # invalid nucleotide
    False
    """
    # Check if lengths are equal
    if len(sequence) != len(quality):
        return False

    # Check nucleotide validity
    if not set(sequence.upper()).issubset(VALID_NUCLEOTIDES):
        return False

    # Check quality validity (Phred+33 from '!' to '~')
    for q_char in quality:
        if not (PHRED33_MIN <= ord(q_char) <= PHRED33_MAX):
            return False

    return True


def calculate_gc_content(sequence: str) -> float:
    """
    Calculates GC content of a sequence in percentage.

    Parameters
    ----------
    sequence : str
        Nucleotide sequence

    Returns
    -------
    float
        GC content in percentage (0-100)

    Examples
    --------
    >>> calculate_gc_content("ATCG")
    50.0
    >>> calculate_gc_content("GGCC")
    100.0
    >>> calculate_gc_content("ATAT")
    0.0
    >>> calculate_gc_content("")
    0.0
    """
    if not sequence:
        return 0.0

    sequence_upper = sequence.upper()
    g_count = sequence_upper.count("G")
    c_count = sequence_upper.count("C")
    total_bases = len(sequence_upper)

    return (g_count + c_count) / total_bases * 100


def calculate_average_quality(quality_string: str) -> float:
    """
    Calculates average quality of a sequence in Phred+33 format.

    Parameters
    ----------
    quality_string : str
        Quality string in Phred+33 encoding

    Returns
    -------
    float
        Average quality value

    Examples
    --------
    >>> calculate_average_quality("IIII")  # Phred+33: I=40
    40.0
    >>> calculate_average_quality("!!!!")  # Phred+33: !=0
    0.0
    >>> calculate_average_quality("FFFF")  # Phred+33: F=37
    37.0
    """
    if not quality_string:
        return 0.0

    quality_scores = [ord(char) - 33 for char in quality_string]
    return sum(quality_scores) / len(quality_scores)


def parse_bounds(bounds, default_low: float, default_high: float):
    """
    Parses filtering bounds into unified format.

    Parameters
    ----------
    bounds : any
        Filtering bounds. Can be:
        - None: use default values
        - Number: set as upper bound
        - Tuple of 2 numbers: (lower_bound, upper_bound)
    default_low : float
        Default lower bound value
    default_high : float
        Default upper bound value

    Returns
    -------
    tuple
        Tuple (lower_bound, upper_bound)

    Examples
    --------
    >>> parse_bounds((20, 80), 0, 100)
    (20.0, 80.0)
    >>> parse_bounds(44.4, 0, 100)
    (0.0, 44.4)
    >>> parse_bounds(None, 0, 100)
    (0.0, 100.0)

    Raises
    ------
    ValueError
        If bounds format is incorrect
    """
    if bounds is None:
        return (default_low, default_high)
    elif isinstance(bounds, (int, float)):
        return (default_low, float(bounds))
    elif isinstance(bounds, (tuple, list)) and len(bounds) == 2:
        return (float(bounds[0]), float(bounds[1]))
    else:
        raise ValueError(f"Incorrect bounds format: {bounds}")


def check_length_filter(sequence: str, length_bounds) -> bool:
    """
    Checks if sequence satisfies length bounds.

    Parameters
    ----------
    sequence : str
        Nucleotide sequence
    length_bounds : tuple
        Length bounds (min, max)

    Returns
    -------
    bool
        True if sequence length is within bounds, False otherwise

    Examples
    --------
    >>> check_length_filter("ATCG", (1, 10))
    True
    >>> check_length_filter("ATCG", (5, 10))
    False
    """
    length = len(sequence)
    return length_bounds[0] <= length <= length_bounds[1]


def check_gc_filter(sequence: str, gc_bounds) -> bool:
    """
    Checks if sequence satisfies GC content bounds.

    Parameters
    ----------
    sequence : str
        Nucleotide sequence
    gc_bounds : tuple
        GC content bounds (min, max) in percentage

    Returns
    -------
    bool
        True if GC content is within bounds, False otherwise

    Examples
    --------
    >>> check_gc_filter("ATCG", (40, 60))
    True
    >>> check_gc_filter("GGCC", (0, 50))
    False
    """
    gc_content_val = calculate_gc_content(sequence)
    return gc_bounds[0] <= gc_content_val <= gc_bounds[1]


def check_quality_filter(quality_string: str, quality_threshold: float) -> bool:
    """
    Checks if sequence quality satisfies threshold value.

    Parameters
    ----------
    quality_string : str
        Quality string in Phred+33 encoding
    quality_threshold : float
        Average quality threshold value

    Returns
    -------
    bool
        True if average quality >= threshold, False otherwise

    Examples
    --------
    >>> check_quality_filter("IIII", 20)  # quality 40 >= 20
    True
    >>> check_quality_filter("!!!!", 1)   # quality 0 < 1
    False
    """
    avg_quality = calculate_average_quality(quality_string)
    return avg_quality >= quality_threshold


def filter_fastq(
    seqs: dict,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold: float = 0,
) -> dict:
    """
    Filters FASTQ sequences by GC content, length and quality.

    Arguments:
        seqs: dict - dictionary of FASTQ sequences {name: (sequence, quality)}
        gc_bounds: any - GC content bounds (number or tuple). Default (0, 100)
        length_bounds: any - length bounds (number or tuple). Default (0, 2**32)
        quality_threshold: float - minimum average quality. Default 0

    Returns:
        dict - filtered dictionary of sequences

    Examples
    --------
    >>> sequences = {
    ...     'read1': ('ATCG', 'IIII'),
    ...     'read2': ('GGCC', '!!!!')
    ... }
    >>> filter_fastq(sequences, gc_bounds=(0, 50))
    {'read1': ('ATCG', 'IIII')}
    """

    # Parse bounds into unified format
    gc_min, gc_max = parse_bounds(gc_bounds, 0, 100)
    length_min, length_max = parse_bounds(length_bounds, 0, 2**32)

    filtered_sequences = {}

    for seq_name, (sequence, quality) in seqs.items():
        # Check FASTQ record validity
        if not validate_fastq_record(sequence, quality):
            continue

        # Check all filtering conditions
        passes_length = check_length_filter(sequence, (length_min, length_max))
        passes_gc = check_gc_filter(sequence, (gc_min, gc_max))
        passes_quality = check_quality_filter(quality, quality_threshold)

        # If all conditions are met, add to result
        if passes_length and passes_gc and passes_quality:
            filtered_sequences[seq_name] = (sequence, quality)

    return filtered_sequences
