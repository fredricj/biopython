from . import _twoBitIO


class TwoBitIterator:

    def __init__(self, handle):
        self.index = 0
        isByteSwapped, sequences, twobit_sequences = _twoBitIO.TwoBitIterator(handle)
        self.isByteSwapped = isByteSwapped
        self.sequences = sequences
        self.twobit_sequences = twobit_sequences

    def __iter__(self):
        return self

    def __next__(self):
        index = self.index
        if index == len(self.sequences):
            raise StopIteration
        sequence = self.sequences[index]
        twobit_sequence = self.twobit_sequences[index]
        self.index += 1
        return sequence, twobit_sequence

    def __len__(self):
        return len(self.sequences)
