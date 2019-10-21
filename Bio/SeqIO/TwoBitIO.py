from . import _twoBitIO


class Seq:

    def __init__(self, data):
        self.data = data

    def __getitem__(self, key):
        data = self.data[key]
        return Seq(data)

    def __eq__(self, other):
        data = self.data
        if not isinstance(data, bytes):
            # then ask the data provider to provide bytes
            data = data[:]
            if not isinstance(data, bytes):
                return NotImplemented
        if isinstance(other, bytes):
            pass
        elif isinstance(other, str):
            other = other.encode('ascii')
        else:
            return NotImplemented
        return (data == other)


class TwoBitIterator:

    def __init__(self, handle):
        self.index = 0
        isByteSwapped, names, sequences = _twoBitIO.TwoBitIterator(handle)
        self.isByteSwapped = isByteSwapped
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

    def __len__(self):
        return len(self.sequences)
