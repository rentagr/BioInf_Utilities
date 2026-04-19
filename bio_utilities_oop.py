"""
Object-oriented implementation of BioUtilities for working with
biological sequences and FASTQ filtering.

Implemented classes:

BiologicalSequence
NucleicAcidSequence
DNASequence
RNASequence
AminoAcidSequence
"""

from abc import ABC, abstractmethod
from typing import Union, Tuple, Any, Optional
from Bio import SeqIO
from Bio.SeqUtils import GC
import logging

"""
#from abc import ABC, abstractmethod
from the module for working with abstract classes (abc) we import the ABC class to inherit from,
and also the decorator to mark methods as abstract.
--------------------------------------------------------------------------------------------------------------------------------------
#from typing import Union, Tuple, Any
import the module for type annotations (improves code readability), helps to understand what data types are expected.
Union - one or several types (integer or slice); Tuple - tuple; Any - any type.
--------------------------------------------------------------------------------------------------------------------------------------
#from Bio import SeqIO (Biopython)
From the Biopython package we import SeqIO. It is needed for reading and writing files with biological sequences
(FASTA, FASTQ, GenBank, etc.). It provides a convenient interface: the parse() function reads records from a file,
and the write() function writes them.
--------------------------------------------------------------------------------------------------------------------------------------
#from Bio.SeqUtils import GC
From the Biopython submodule Bio.SeqUtils we import the module that calculates the percentage of GC (guanine and cytosine)
in a sequence. Returns a number from 0 to 100.
"""


class BiologicalSequence(ABC):
    """
    Abstract class for any sequence.
    Requires implementation of the check_alphabet method.
    """

    def __init__(self, sequence: str):
        self._sequence = sequence.upper()
        if not self.check_alphabet():
            raise ValueError(f"Invalid sequence: {sequence}")

    """
    The constructor takes a sequence string.
    Converts it to uppercase (to unify representation) and stores it in the protected attribute _sequence.
    Then calls the check_alphabet() method in the child class to validate characters.
    If validation fails, raises a ValueError with a message indicating the user's input is invalid. 
    """

    def __len__(self) -> int:
        return len(self._sequence)

    """
    The __len__ method is called by the len() function. Returns the length of the sequence.
    """

    def __getitem__(self, index: Union[int, slice]) -> Any:
        if isinstance(index, slice):
            # Returning an object of the same class for slices
            return self.__class__(self._sequence[index])
        return self._sequence[index]

    """
    __getitem__ allows indexing: obj[i] or obj[i:j].
    Type annotation: index can be int or slice.
    Return type Any, because it may return either a string (if int) or a new instance of the same class (if slice).
    If index is a slice, we take a slice of the string self._sequence[index] and create a new object of the same class
    (self.__class__) with that substring. Thus, a slice of DNA remains DNA.
    If index is an integer, we return a single character (string).
    """

    def __str__(self) -> str:
        return self._sequence

    """
    __str__ is called when converting to a string.
    Returns the sequence itself (string). Thus we see the string,
    not something like <main.BiologicalSequence object at 0x...>.
    """

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self._sequence})"

    """
    __repr__ provides a formal representation of the object, useful for debugging.
    Returns a string like: DNASequence(ATGC)
    """

    @abstractmethod
    def check_alphabet(self) -> bool:
        pass

    """
    Checks the correctness of the sequence alphabet.
    This is an abstract method with no implementation in this class.
    Any non-abstract child class must implement it.
    The method should return bool: True if all characters are valid for this sequence type,
    validating the alphabet. Abstractness ensures that each specific class (DNA, RNA, protein)
    performs its own validation.
    """


# Nucleic acids (DNA/RNA)
class NucleicAcidSequence(BiologicalSequence):
    """
    Common class for nucleic acids. Contains methods:
    complement, reverse, reverse_complement.
    """

    _complement_map: dict = {}
    """
    class attribute (dictionary) that will store complementarity rules for a specific nucleic acid type
    (overridden in descendants). We will use the map defined in the object's class.
    """

    def complement(self):
        """Returns the complementary sequence (object of the same class)."""
        complemented = "".join(
            self._complement_map.get(base, base) for base in self._sequence
        )
        return self.__class__(complemented)

    def reverse(self):
        """Returns the reversed sequence"""
        return self.__class__(self._sequence[::-1])

    def reverse_complement(self):
        """Returns the reverse complementary sequence."""
        return self.reverse().complement()


# DNA
class DNASequence(NucleicAcidSequence):
    """Class for working with DNA"""

    _complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}
    _alphabet = set("ATGC")

    def check_alphabet(self) -> bool:
        return all(base in self._alphabet for base in self._sequence)

    def transcribe(self):
        """Transcription of DNA to RNA (replace T with U)."""
        rna_seq = self._sequence.replace("T", "U")
        return RNASequence(rna_seq)


# RNA
class RNASequence(NucleicAcidSequence):
    """Class for working with RNA."""

    _complement_map = {"A": "U", "U": "A", "G": "C", "C": "G"}
    _alphabet = set("AUGC")

    def check_alphabet(self) -> bool:
        return all(base in self._alphabet for base in self._sequence)


# Squirrels from Pushchino
class AminoAcidSequence(BiologicalSequence):
    """Class for working with amino acid sequences."""

    _alphabet = set("ACDEFGHIKLMNPQRSTVWY")

    def check_alphabet(self) -> bool:
        return all(ama in self._alphabet for ama in self._sequence)

    def molecular_weight(self) -> float:
        """Example method: approximate molecular weight."""
        # Table of amino acid residue masses (g/mol, approximate)
        weights = {
            "A": 90.1,
            "R": 174.1,
            "N": 132.1,
            "D": 133.1,
            "C": 121.2,
            "Q": 146.2,
            "E": 147.1,
            "G": 75.1,
            "H": 155.2,
            "I": 131.2,
            "L": 131.2,
            "K": 146.2,
            "M": 149.2,
            "F": 165.2,
            "P": 115.0,
            "S": 105.1,
            "T": 119.1,
            "W": 204.0,
            "Y": 181.2,
            "V": 117.1,
        }
        return sum(weights.get(ama, 0) for ama in self._sequence)


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    length_bounds: Tuple[float, float] = (0, float("inf")),
    quality_threshold: float = 0,
    gc_bounds: Tuple[float, float] = (0, 100),
    log_file: Optional[str] = None,
) -> None:
    """
    Filters a FASTQ file by length, average quality, and GC content.

    Parameters:
        input_fastq (str): path to input FASTQ file
        output_fastq (str): path to save filtered records
        length_bounds (tuple): (min, max) allowed length
        quality_threshold (float): minimum average quality (Phred)
        gc_bounds (tuple): (min, max) allowed GC content in percent
    """
    logger = None
    if log_file:
        logger = logging.getLogger('FilterFastQ')
        logger.setLevel(logging.INFO)
        # Cleaning up old handlers so as not to duplicate
        if logger.handlers:
            logger.handlers.clear()
        fh = logging.FileHandler(log_file) # Creating a handler in a file. He will write messages to the file.
        fh.setLevel(logging.INFO)
        ch = logging.StreamHandler() # Creating a handler in the console (standard output stream). It will display the same messages on the screen.
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s') # Setting the message format
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        logger.addHandler(fh)
        logger.addHandler(ch)
        logger.info(f"Starting filtering: {input_fastq} -> {output_fastq}")
        logger.info(f"Parameters: length={length_bounds}, quality={quality_threshold}, GC={gc_bounds}")
    # Counters: total- how many reads have been read, passed- how many reads have passed
    total = 0
    passed = 0
    try:
        with open(input_fastq, "r") as in_handle, open(output_fastq, "w") as out_handle:
            for record in SeqIO.parse(in_handle, "fastq"):
                total += 1
                seq = record.seq
                # Length
                if not (length_bounds[0] <= len(seq) <= length_bounds[1]):
                    continue
                # GC content
                gc_content = GC(seq)
                if not (gc_bounds[0] <= gc_content <= gc_bounds[1]):
                    continue
                # Quality
                qualities = record.letter_annotations["phred_quality"]
                avg_quality = sum(qualities) / len(qualities)
                if avg_quality < quality_threshold:
                    continue
                # Passed all filters
                SeqIO.write(record, out_handle, "fastq")
                passed += 1

        if logger:
            logger.info(f"Records processed: {total}, saved: {passed}")
    # Exception handling
    except FileNotFoundError as e: # If the input file does not exist
        if logger:
            logger.error(f"File is not find: {e}")
        raise
    except Exception as e: # Catches any other possible error.
        if logger:
            logger.error(f"Error during processing: {e}", exc_info=True)
        raise
            
if __name__ == "__main__":
    import argparse
    """
    The command line interface for terminal operations.       
    
    For optional arguments:    
    type=float - converts the entered string to a floating-point number. By default, type=str.
    default=... - the default value if the user has not specified an argument. Float("inf") is used for --len-max, meaning there is no upper bound.
    for --log, the type is a string.
    """

    parser = argparse.ArgumentParser(description="FASTQ filter by length, quality and GC%") # Creating a parser
    parser.add_argument("input", help="FASTQ input file") # Adding the required positional argument
    parser.add_argument("-o", "--output", default="filtered.fastq", help="FASTQ output file") # Here and below -adding optional arguments.
    parser.add_argument("-l", "--len-min", type=float, default=0, help="Minimum reed length")
    parser.add_argument("-L", "--len-max", type=float, default=float("inf"), help="Maximum reed length")
    parser.add_argument("-q", "--qual", type=float, default=0, help="Minimum Average Quality (Phred)")
    parser.add_argument("-g", "--gc-min", type=float, default=0, help="Minimum GC%")
    parser.add_argument("-G", "--gc-max", type=float, default=100, help="Maximum GC%")
    parser.add_argument("--log", default="filter.log", help="File for logs")

    args = parser.parse_args()

    filter_fastq(
        input_fastq=args.input,
        output_fastq=args.output,
        length_bounds=(args.len_min, args.len_max),
        quality_threshold=args.qual,
        gc_bounds=(args.gc_min, args.gc_max),
        log_file=args.log   
    )
