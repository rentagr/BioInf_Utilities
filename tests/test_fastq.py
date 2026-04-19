import pytest
from bio_utilities_oop import filter_fastq
from Bio import SeqIO
import tempfile
import os

"""
Tests for the filter_fastq function.

The tests use artificial FASTQ data created directly in the code.
The quality symbols comply with the Phred+ standard33:
'!' (ASCII 33) -> quality 0
    'I' (ASCII 73) -> quality 40

Each test creates a temporary FASTQ file in the tmp_path fixture,
writing lines like:
    @read_id
    SEQUENCE
    +
    QUALITY_STRING

For example, the string "@r1\\nAGCT\\n+\\nIIII\\n" represents one reed
with the AGCT sequence and a quality of 40 for each nucleotide.

This allows to control the filtration by length, GC composition.,
quality and their combinations without using external files.
"""

# Read/write file test
def test_filter_fastq_writes_output(tmp_path): # Defining the test function
    input_file = tmp_path / "in.fastq" # Creating a temporary file path object
    # Reads: r1 (AGCT sequence, quality IIII — 40), r2 (AAAA, quality !!!! — 0)
    input_file.write_text("@r1\nAGCT\n+\nIIII\n@r2\nAAAA\n+\n!!!!\n") # Writes a string (the contents of the FASTQ file) to this temporary file / Creating a temporary FASTQ input with two reads
    output_file = tmp_path / "out.fastq" # The path to the output file that our filter_fastq function should create
    
    filter_fastq(str(input_file), str(output_file), length_bounds=(2,10), quality_threshold=30) # Calling the function under test. str(input_file) converts the path to a string
    
    assert output_file.exists()
    records = list(SeqIO.parse(output_file, "fastq")) # Reading the FASTQ output file
    assert len(records) == 1
    assert records[0].seq == "AGCT"

# The file does not exist error test
def test_filter_fastq_missing_input():
    with pytest.raises(FileNotFoundError):
        filter_fastq("nonexistent.fastq", "out.fastq")

# Filter by length 
def test_length_filter(tmp_path):
    input_file = tmp_path / "in.fastq"
    # Two reads: r1 is 1 long, r2 is 4 long. Both have high quality (I = 40)
    input_file.write_text("@r1\nA\n+\nI\n@r2\nAAAA\n+\nIIII\n")
    output_file = tmp_path / "out.fastq"
    filter_fastq(str(input_file), str(output_file), length_bounds=(2,10)) # The minimum length is 2, the maximum is 10. The reed r1 (length 1) must be discarded, r2 (length 4) must pass.
    records = list(SeqIO.parse(output_file, "fastq"))
    assert len(records) == 1
    assert records[0].seq == "AAAA"

# GC filter
def test_gc_filter(tmp_path):
    input_file = tmp_path / "in.fastq"
    # r1: GCGC - GC% = 100%, r2: ATAT - GC% = 0%. The quality is excellent (think so)
    input_file.write_text("@r1\nGCGC\n+\nIIII\n@r2\nATAT\n+\nIIII\n")
    output_file = tmp_path / "out.fastq"
    # Leave the reads with GC from 0% to 30%. Only r2 (0%) is suitable.
    filter_fastq(str(input_file), str(output_file), gc_bounds=(0,30))
    records = list(SeqIO.parse(output_file, "fastq"))
    assert len(records) == 1
    assert records[0].seq == "ATAT"

# Filter by quality
def test_quality_filter(tmp_path):
    input_file = tmp_path / "in.fastq"
    # Both reads have the same AGCT sequence, but different qualities: IIII (average quality 40), !!!! (average quality is 0).
    input_file.write_text("@r1\nAGCT\n+\nIIII\n@r2\nAGCT\n+\n!!!!\n")
    output_file = tmp_path / "out.fastq"
    # Leave the reads with an average quality of ≥ 30. Only the first one passes
    filter_fastq(str(input_file), str(output_file), quality_threshold=30)
    records = list(SeqIO.parse(output_file, "fastq"))
    assert len(records) == 1
    assert records[0].seq == "AGCT"

# Filter Combination
def test_combined_filters(tmp_path):
    input_file = tmp_path / "in.fastq"
    # Three reeds: r1: AGCT (length 4, GC 50%, quality 40); r2: GCGC (length 4, GC 100%, quality 40); r3: A (length 1, GC 0%, quality 40)
    input_file.write_text("@r1\nAGCT\n+\nIIII\n@r2\nGCGC\n+\nIIII\n@r3\nA\n+\nI\n")
    output_file = tmp_path / "out.fastq"
    #Conditions: length 3-5, GC 40-60%, quality ≥30.
    #r1: length 4 (yes), GC 50% (yes), quality (yes) - passes
    #r2: GC 100% - does not pass
    #r3: length 1 - does not pass
    filter_fastq(str(input_file), str(output_file), length_bounds=(3,5), gc_bounds=(40,60), quality_threshold=30)
    records = list(SeqIO.parse(output_file, "fastq"))
    assert len(records) == 1
    assert records[0].seq == "AGCT"

# Empty output file if nothing fits.
def test_empty_output(tmp_path):
    input_file = tmp_path / "in.fastq"
    # One read with poor quality (0)
    input_file.write_text("@r1\nAAAA\n+\n!!!!\n")
    output_file = tmp_path / "out.fastq"
    # Require a quality of ≥30, the reed is not suitable.
    filter_fastq(str(input_file), str(output_file), quality_threshold=30)
    records = list(SeqIO.parse(output_file, "fastq"))
    assert len(records) == 0

# Checking that the log file is being created (if transmitted)
def test_log_file_created(tmp_path):
    input_file = tmp_path / "in.fastq"
    # One good reed
    input_file.write_text("@r1\nAGCT\n+\nIIII\n")
    output_file = tmp_path / "out.fastq"
    log_file = tmp_path / "log.txt"
    # Passing the log_file parameter with the path to the temporary log file
    filter_fastq(str(input_file), str(output_file), log_file=str(log_file))
    assert log_file.exists()
    log_content = log_file.read_text()
    assert "Starting filtering" in log_content
    assert "Records processed" in log_content