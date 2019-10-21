"""Tests for SeqIO TwoBitIO module."""
  

import os
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqIO import TwoBitIO

path = "TwoBit/sequence.littleendian.2bit"
handle = open(path)
f = TwoBitIO.TwoBitIterator(handle)
assert len(f) == 5
assert f.isByteSwapped is False
handle.close()

path = "TwoBit/sequence.bigendian.2bit"
handle = open(path)
f = TwoBitIO.TwoBitIterator(handle)
assert len(f) == 5
assert f.isByteSwapped is True
handle.close()


for length in range(1,21):
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
    sequences = TwoBitIO.TwoBitIterator(handle)
    for sequence, record in zip(sequences, records):
        for start in range(length):
            for end in range(start+1, length+1):
                print("Testing sequence length %d start %d end %d" % (length, start, end))
                for i in range(10):
                    seq1 = sequence[start:end]
                    seq2 = record.seq[start:end]
                    assert seq1 == seq2
                    for step in range(1, end-start+1):
                        seq1 = sequence[start:end:step]
                        seq2 = record.seq[start:end:step]
                        assert seq1 == seq2
                        seq1 = sequence[end:start:-step]
                        seq2 = record.seq[end:start:-step]
                        assert seq1 == seq2
