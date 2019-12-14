from . import _twoBitIO


class Seq(_twoBitIO.Seq):

    def __new__(cls, data):
        self = super().__new__(cls, data)
        self.data = data
        return self

    def __getitem__(self, key):
        data = self.data[key]
        if isinstance(key, int):
            try:
                letter = chr(data)
            except TypeError:
                letter = data.decode("ascii")
            return letter
        else:
            return Seq(data)

    def __hash__(self):
        return hash(self.data[:])

    def _convert_for_comparison(self, other):
        data = self.data
        if isinstance(data, bytes):
            left = data
        else:
            # then ask the data provider to provide bytes
            left = data[:]
            if not isinstance(left, bytes):
                raise NotImplementedError
        if isinstance(other, Seq):
            right = other.data[:]
        elif isinstance(other, bytes):
            right = other
        elif isinstance(other, str):
            right = other.encode('ascii')
        else:
            raise NotImplementedError
        return left, right

    def __eq__(self, other):
        left, right = self._convert_for_comparison(other)
        return left == right

    def __lt__(self, other):
        left, right = self._convert_for_comparison(other)
        return left < right

    def __le__(self, other):
        left, right = self._convert_for_comparison(other)
        return left <= right

    def __gt__(self, other):
        left, right = self._convert_for_comparison(other)
        return left > right

    def __ge__(self, other):
        left, right = self._convert_for_comparison(other)
        return left >= right

    def __repr__(self):
        length = len(self.data)
        if length > 60:
            data = self.data[:54] + b'...' + self.data[-3:]
        else:
            data = self.data[:]
        sequence = data.decode('ascii')
        return 'Seq("%s")' % sequence


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
