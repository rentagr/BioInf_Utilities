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


import os
from pathlib import Path


def read_fastq_file(input_path: str) -> dict:
    """
    Read FASTQ file and convert to dictionary format.

    Args:
        input_path: Path to input FASTQ file

    Returns:
        Dictionary with format {header: (sequence, quality)}

    Raises:
        FileNotFoundError: If input file doesn't exist
        ValueError: If file has invalid FASTQ format
    """
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"File {input_path} not found")

    sequences = {}

    with open(input_path, "r") as file:
        lines = [line.strip() for line in file]

    # FASTQ files should have lines divisible by 4
    if len(lines) % 4 != 0:
        raise ValueError("Invalid FASTQ format: number of lines not divisible by 4")

    # Process each sequence record (4 lines each)
    for i in range(0, len(lines), 4):
        header = lines[i][1:]  # Remove '@' from header
        sequence = lines[i + 1]
        quality = lines[i + 3]

        sequences[header] = (sequence, quality)

    return sequences


def write_fastq_file(sequences: dict, output_path: str) -> None:
    """
    Write sequences from dictionary to FASTQ file.

    Args:
        sequences: Dictionary {header: (sequence, quality)}
        output_path: Path for output FASTQ file
    """
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(output_path, "w") as file:
        for header, (sequence, quality) in sequences.items():
            file.write(f"@{header}\n")
            file.write(f"{sequence}\n")
            file.write("+\n")  # separator line
            file.write(f"{quality}\n")


def filter_fastq_file(
    input_fastq: str,
    output_fastq: str,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """
    Main function to filter FASTQ files.

    This function combines file reading, filtering, and writing.

    Args:
        input_fastq: Path to input FASTQ file
        output_fastq: Path for filtered output file
        gc_bounds: GC content bounds (tuple or single number)
        length_bounds: Sequence length bounds (tuple or single number)
        quality_threshold: Minimum average quality score
    """
    print(f"Reading FASTQ file: {input_fastq}")

    # Step 1: Read FASTQ file into dictionary
    sequences = read_fastq_file(input_fastq)
    print(f"Found {len(sequences)} sequences")

    # Step 2: Filter sequences using existing function
    print("Filtering sequences...")
    filtered_sequences = filter_fastq_dict(
        sequences, gc_bounds, length_bounds, quality_threshold
    )
    print(f"After filtering: {len(filtered_sequences)} sequences remain")

    # Step 3: Write filtered sequences to new file
    print(f"Writing results to: {output_fastq}")
    write_fastq_file(filtered_sequences, output_fastq)

    print("FASTQ filtering completed successfully")
