# Constants for validation
VALID_NUCLEOTIDES = set("ATCGNnatcgn")
PHRED33_MIN = 33  # Minimum ASCII for Phred+33 ('!')
PHRED33_MAX = 126  # Maximum ASCII for Phred+126 ('~')


def validate_fastq_record(sequence: str, quality: str) -> bool:
    """
    Validate a FASTQ record for correct format and characters.
    """
    # Check sequence and quality have same length
    if len(sequence) != len(quality):
        return False

    # Check all nucleotides are valid
    if not set(sequence.upper()).issubset(VALID_NUCLEOTIDES):
        return False

    # Check quality scores are in valid range
    return all(PHRED33_MIN <= ord(q_char) <= PHRED33_MAX for q_char in quality)


def calculate_gc_content(sequence: str) -> float:
    """
    Calculate GC content percentage of a sequence.
    """
    if not sequence:
        return 0.0
    seq_upper = sequence.upper()
    return (seq_upper.count("G") + seq_upper.count("C")) / len(seq_upper) * 100


def calculate_average_quality(quality_string: str) -> float:
    """
    Calculate average quality score from Phred+33 encoded string.
    """
    if not quality_string:
        return 0.0
    quality_scores = [ord(char) - 33 for char in quality_string]
    return sum(quality_scores) / len(quality_scores)


def parse_bounds(bounds, default_low: float, default_high: float):
    """
    Parse filtering bounds into standardized format.
    Supports: None, single number, or tuple of two numbers.
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
    """Check if sequence length is within bounds."""
    min_len, max_len = length_bounds
    length = len(sequence)
    return min_len <= length <= max_len


def check_gc_filter(sequence: str, gc_bounds) -> bool:
    """Check if GC content is within bounds."""
    min_gc, max_gc = gc_bounds
    gc_content_val = calculate_gc_content(sequence)
    return min_gc <= gc_content_val <= max_gc


def check_quality_filter(quality_string: str, quality_threshold: float) -> bool:
    """Check if average quality meets threshold."""
    avg_quality = calculate_average_quality(quality_string)
    return avg_quality >= quality_threshold


def filter_fastq_dict(
    seqs: dict,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold: float = 0,
) -> dict:
    """
    Filter FASTQ sequences stored in dictionary format.

    This is the original filter function that works with dictionaries.
    Used internally by the file-based filtering function.
    """
    # Parse bounds
    gc_min, gc_max = parse_bounds(gc_bounds, 0, 100)
    length_min, length_max = parse_bounds(length_bounds, 0, 2**32)

    filtered_sequences = {}

    for seq_name, (sequence, quality) in seqs.items():
        # Validate FASTQ record
        if not validate_fastq_record(sequence, quality):
            continue

        # Apply all filters
        passes_length = check_length_filter(sequence, (length_min, length_max))
        passes_gc = check_gc_filter(sequence, (gc_min, gc_max))
        passes_quality = check_quality_filter(quality, quality_threshold)

        # Add to results if all filters pass
        if passes_length and passes_gc and passes_quality:
            filtered_sequences[seq_name] = (sequence, quality)

    return filtered_sequences
