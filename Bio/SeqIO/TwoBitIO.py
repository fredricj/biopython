from . import _twoBitIO


class TwoBitIterator:

    def __init__(self, handle):
        self.index = 0
        isByteSwapped, sequences = _twoBitIO.TwoBitIterator(handle)
        self.isByteSwapped = isByteSwapped
        self.sequences = sequences

    def __iter__(self):
        return self

    def __next__(self):
        index = self.index
        if index == len(self.sequences):
            raise StopIteration
        sequence = self.sequences[index]
        self.index += 1
        return sequence

    def __len__(self):
        return len(self.sequences)
