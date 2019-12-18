from . import _twoBitIO


class Seq(_twoBitIO.Seq):

    def __new__(cls, data):
        self = super().__new__(cls, data)
        self.data = data
        return self

    def __hash__(self):
        return hash(self.data[:])


class TwoBitIterator:

    def __init__(self, handle):
        self.index = 0
        isByteSwapped, names, sequences = _twoBitIO.TwoBitIterator(handle)
        self.isByteSwapped = isByteSwapped
        self.names = names
        self.sequences = []
        for name, sequence in zip(names, sequences):
            sequence = Seq(sequence)
            sequence.name = name
            self.sequences.append(sequence)

    def __iter__(self):
        return self

    def __next__(self):
        index = self.index
        if index == len(self.sequences):
            raise StopIteration
        sequence = self.sequences[index]
        self.index += 1
        return sequence

    def __getitem__(self, key):
        try:
            index = self.names.index(key)
        except ValueError:
            raise KeyError(key) from None
        sequence = self.sequences[index]
        return sequence

    def keys(self):
        return self.names

    def __len__(self):
        return len(self.sequences)
