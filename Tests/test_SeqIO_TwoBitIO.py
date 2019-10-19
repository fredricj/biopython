"""Tests for SeqIO TwoBitIO module."""
  

import os
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqIO import TwoBitIO

path = "TwoBit/sequence.littleendian.2bit"
handle = open(path)
f = TwoBitIO.TwoBitIterator(handle, 0, 1)
assert len(f) == 5
assert f.isByteSwapped is False
handle.close()

path = "TwoBit/sequence.bigendian.2bit"
handle = open(path)
f = TwoBitIO.TwoBitIterator(handle, 0, 1)
assert len(f) == 5
assert f.isByteSwapped is True
handle.close()

def perform_test(length=50, start=0, end=None, n=10):
    if end is None:
        end = length
    records = []
    for i in range(10):
        nucleotides = ['ACGTNacgtn'[random.randint(0,9)] for i in range(length)]
        sequence = ''.join(nucleotides)
        seq = Seq(sequence)
        record = SeqRecord(seq, id='name_%d' % i)
        records.append(record)
    handle = open("test.fa", 'w')
    SeqIO.write(records, handle, 'fasta')
    handle.close()
    os.system("faToTwoBit test.fa test.2bit")
    handle = open("test.2bit")
    sequences = TwoBitIO.TwoBitIterator(handle, start, end)
    handle.close()
    for sequence, record in zip(sequences, records):
        assert sequence == record.seq[start:end]

for length in range(1,21):
    for start in range(length):
        for end in range(start+1, length+1):
            print("Testing sequence length %d start %d end %d" % (length, start, end))
            for i in range(10):
                perform_test(length, start, end)
