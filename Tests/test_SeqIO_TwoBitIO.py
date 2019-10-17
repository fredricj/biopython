"""Tests for SeqIO TwoBitIO module."""
  

from Bio import SeqIO

path = "TwoBit/sequence.littleendian.2bit"
records = SeqIO.parse(path, '2bit')
record = next(records)
