import argparse
from typing import List, Optional, Dict
from itertools import product


class FASTA_Parser:

    pathname: str = ''
    file_string: str = ''
    readingframe: int = 0
    records: Optional[int] = None
    sequences: Optional[Dict] = None

    def __init__(self, pathname, readingframe) -> None:
        """
        Initializes a new instance of the class with the given record.

        Args:
            record: A genomic data record.

        Returns:
            None.
        """
        self.pathname = pathname
        self.readingframe = int(readingframe)
        self._process_file()

    def _process_file(self) -> int:
        """
        For internal use on instantiation. 
        Counts the number of records in the FASTA file and stores the sequences in a dictionary.

        Args:
            None.

        Returns:
            The number of records in the FASTA file.
        """
        record_count = 0
        identifier, seq_length, seq = None, 0, ''
        d = {}
        with open(self.pathname, 'r') as f:
            for line in f:

                if line.startswith('>'):
                    if identifier is None:
                        # skip for first iteration
                        pass
                    else:
                        d[identifier] = {'seq_length': seq_length, 'seq': seq}
                        seq_length = 0
                        seq = ''

                    record_count += 1
                    identifier = line.split(' ')[0].replace('>', '')
                else:
                    seq += line.strip()
                    seq_length += len(line.strip())
                self.file_string += line

        self.records = record_count
        self.sequences = d

        self._generate_orfs()

    def get_longest_sequence(self) -> str:
        """
        Returns the longest sequence in the dictionary of sequences.

        Returns:
            A string representing the longest sequence in the dictionary of sequences.
            The string is formatted as "<sequence_id>: <sequence>".
        """
        highest_key = max(
            self.sequences, key=lambda k: self.sequences[k]['seq_length'])
        return f'{highest_key}: {self.sequences[highest_key]["seq_length"]}'

    def get_shortest_sequence(self) -> str:
        """
        Returns the shortest sequence in the dictionary of sequences, along with its key.

        Returns:
            A string in the format "key: sequence", where key is the key of the shortest sequence
            and sequence is the shortest sequence itself.
        """
        lowest_key = min(
            self.sequences, key=lambda k: self.sequences[k]['seq_length'])
        return f'{lowest_key}: {self.sequences[lowest_key]["seq_length"]}'

    def _generate_orfs(self) -> None:
        """
        Internal method. Generates all open reading frames (ORFs) in all sequences.
        """
        for identifier in self.sequences:
            print(f'Generating ORFs for {identifier}')
            orf_parser = ORF_Parser(
                self.sequences[identifier]['seq'], self.readingframe)
            self.sequences[identifier]['orfs'] = orf_parser.orfs

    def get_longest_orf(self) -> str:
        """
        Returns the longest open reading frame (ORF) found in the sequences of this object.

        An ORF is a sequence of nucleotides that starts with a start codon (ATG) and ends with a stop codon
        (TAA, TAG, or TGA), without any other stop codons in between. This method searches for all ORFs in
        all sequences of this object, and returns the longest one found, along with its identifier and
        nucleotide sequence.

        Returns:
            A tuple with three elements:
            - The identifier of the sequence containing the longest ORF.
            - The length (in nucleotides) of the longest ORF.
            - The nucleotide sequence of the longest ORF.
        """
        max_length = 0
        max_identifier = ''
        max_seq = ''
        for identifier in self.sequences:
            for orf in self.sequences[identifier]['orfs']:
                if orf['len'] > max_length:
                    max_length = orf['len']
                    max_identifier = identifier
                    max_seq = orf['orf']
        return (max_identifier, max_length, max_seq)

    def get_longest_orf_by_identifier(self, identifier) -> str:
        """
        Returns the longest open reading frame (ORF) sequence and its length for a given sequence identifier.

        Args:
            identifier (str): The identifier of the sequence to search for ORFs.

        Returns:
            Tuple[int, str]: A tuple containing the length of the longest ORF found and its sequence.
        """
        max_length = 0
        max_seq = ''
        for orf in self.sequences[identifier]['orfs']:
            if orf['len'] > max_length:
                max_length = orf['len']
                max_seq = orf['orf']
        return (max_length, max_seq)

    def find_all_substrings(self, n):
        """
        Finds all substrings of length n in the file string and returns a dictionary
        with the count of occurrences of each substring, sorted by increasing count.

        Args:
            n (int): the length of the substrings to search for.

        Returns:
            dict: a dictionary with the count of occurrences of each substring of length n
                found in the file string, sorted by increasing count.
        """
        combinations = [''.join(i) for i in product('ATCG', repeat=n)]
        counts = dict()

        for c in combinations:
            counts[c] = len(list(
                find_all_overlapping_occurrences(self.file_string, c)))
        sorted_counts = dict(sorted(counts.items(), key=lambda item: item[1]))
        return sorted_counts


def find_all_overlapping_occurrences(string, substring):
    """
    Find all overlapping occurrences of a substring in a string.

    Args:
        string (str): The string to search in.
        substring (str): The substring to search for.

    Yields:
        int: The index of each overlapping occurrence of the substring in the string.
    """
    start = 0
    while start < len(string):
        start = string.find(substring, start)
        if start == -1:
            break
        yield start
        start += 1


class ORF_Parser:
    """
    A class to parse open reading frames (ORFs) from a given DNA sequence.

    Attributes:
    -----------
    seq : str
        The DNA sequence to parse ORFs from.
    readingframe : int
        The reading frame to start parsing ORFs from.
    orfs : list
        A list of dictionaries containing information about each ORF found.
    """

    seq: str = ''
    readingframe: int = 0
    orfs = []

    def __init__(self, seq, readingframe) -> None:
        """
        Initializes the ORF_Parser object.

        Parameters:
        -----------
        seq : str
            The DNA sequence to parse ORFs from.
        readingframe : int
            The reading frame to start parsing ORFs from.
        """
        self.seq = seq
        self.readingframe = int(readingframe)
        self._generate_orfs()

    def _generate_orfs(self) -> None:
        """
        Generates a list of ORFs from the given DNA sequence and reading frame.
        """
        pos = 1
        codon = ''
        orf = ''
        in_orf = False
        for c in self.seq[self.readingframe - 1:]:
            codon += c
            pos += 1

            if len(codon) == 3:
                if codon == 'ATG':
                    in_orf = True
                    orf += codon
                elif codon in ['TAA', 'TAG', 'TGA']:
                    if not in_orf:
                        continue
                    else:
                        orf += codon
                        self.orfs.append({'start_pos': pos - len(orf),
                                          'orf': orf,
                                          'len': len(orf)})
                        orf = ''
                        in_orf = False

                if in_orf:
                    orf += codon
                codon = ''
        print(self.orfs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='FASTA Parser')
    parser.add_argument(
        '-c', '--count', action='store_true', help='count number of records', required=False, default=True)
    parser.add_argument(
        '-r', '--readingframe', help='reading frame to use for ORF detection', required=False, default='1')
    parser.add_argument(
        '-n', help='used for calculating substrings', required=True)
    args = parser.parse_args()

    if args.readingframe not in ['1', '2', '3']:
        print('Invalid reading frame. Must be 1, 2, or 3.')
        exit(1)
    try:
        int(args.n)
    except ValueError:
        print('Invalid value for n. Must be an integer.')
        exit(1)

    if args.count:
        print(f'Counting records with options {args.readingframe}...')
        parser = FASTA_Parser('examples/dna2.fasta', args.readingframe)
        print('====================================================')
        print(f'Longest sequence: {parser.get_longest_sequence()}')
        print('====================================================')
        print(f'Shortest sequence: {parser.get_shortest_sequence()}')
        print('====================================================')
        print(f'Longest ORF: {parser.get_longest_orf()}')
        print('====================================================')
        print(f'Combinations: {parser.find_all_substrings(int(args.n))}')
        print('====================================================')
        print(
            f'Longest ORF by identifier: {parser.get_longest_orf_by_identifier("gi|142022655|gb|EQ086233.1|16")}')
