from Bio.SeqIO.FastaIO import FastaIterator

class SlidingParser:
    
    def __init__(self, window_size):
        self.window_size = window_size
    
    def parse_sequence(self, sequence):
        records_letters = {}
        records_letters = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'Y': 0, 'M': 0, 'S': 0, 'R': 0, 'W': 0,
                                      'K': 0, 'N': 0, 'D': 0, 'B': 0, 'H': 0, 'V': 0, 'all': 0, 'a': 0, 'c': 0, 'g': 0, 't': 0, 'y': 0, 'm': 0, 's': 0, 'r': 0, 'w': 0,
                                      'k': 0, 'n': 0, 'd': 0, 'b': 0, 'h': 0, 'v': 0, 'all_small': 0, 'all_big': 0}
        for letter in sequence:
            if letter.islower():
                records_letters['all_small'] += 1
            else:
                records_letters['all_big'] += 1
            records_letters[letter] += 1
            records_letters['all'] += 1
        return records_letters  

    def parse_file(self, file_path):
        data = {}
        print("Analysing: " + file_path)
        with open(file_path) as file:
            for record in FastaIterator(file):
                data[record.id] = {}
                start_index = 0
                end_index = len(record.seq) - 1

                while start_index + self.window_size < end_index:
                    data[record.id][start_index] = self.parse_sequence(record.seq[start_index : (start_index + self.window_size)])
                    start_index += self.window_size

                data[record.id][start_index] = self.parse_sequence(record.seq[start_index : end_index])
                return data