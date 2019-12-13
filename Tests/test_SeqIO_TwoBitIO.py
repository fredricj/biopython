"""Tests for SeqIO TwoBitIO module."""
  

import os
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqIO import TwoBitIO

path = "TwoBit/sequence.fa"
handle = open(path)
records = SeqIO.parse(handle, 'fasta')
records = list(records)
handle.close()

path = "TwoBit/sequence.littleendian.2bit"
handle = open(path)
sequences = TwoBitIO.TwoBitIterator(handle)
assert len(sequences) == 5
assert sequences.isByteSwapped is False
for sequence, record in zip(sequences, records):
    assert sequence == str(record.seq)
    assert sequence.name == record.id
    assert len(sequence) == len(record.seq)
handle.close()

path = "TwoBit/sequence.bigendian.2bit"
handle = open(path)
sequences = TwoBitIO.TwoBitIterator(handle)
assert len(sequences) == 5
assert sequences.isByteSwapped is True
for sequence, record in zip(sequences, records):
    assert sequence == str(record.seq)
    assert sequence.name == record.id
    assert len(sequence) == len(record.seq)
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
        seq1 = sequence
        seq2 = str(record.seq)
        assert seq1 == seq2
        for start in range(length):
            for end in range(start+1, length+1):
                print("Testing sequence length %d start %d end %d" % (length, start, end))
                for i in range(10):
                    assert seq1[start:end] == seq2[start:end]
                    for step in range(1, end-start+1):
                        assert seq1[start:end:step] == seq2[start:end:step]
                        assert seq1[end:start:-step] == seq2[end:start:-step]
