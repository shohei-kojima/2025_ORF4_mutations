#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Python 3.12.2
pandas 2.2.2
'''

import os,sys,glob,gzip,collections
import pandas as pd

samples_str = '''146.Sendai-H/1121/4_4d
147.Sendai-H/1121/4_7d
149.Sendai-H/1121/4_18d
150.Sendai-H/1121/4_25d
151.Sendai-H/1121/4_32d
152.Sendai-H/1121/4_39d
153.Sendai-H/1121/4_46d
154.Sendai-H/1121/4_53d
155.Sendai-H/1121/4_60d
156.Sendai-H/1121/4_67d
157.Sendai-H/1121/4_74d
158.Sendai-H/1121/4_81d
160.Sendai-H/1948/04_4d
161.Sendai-H/1948/04_7d
163.Sendai-H/1948/04_18d
164.Sendai-H/1948/04_25d
166.Sendai-H/1948/04_39d
167.Sendai-H/1948/04_46d
169.Sendai-H/1948/04_60d
170.Sendai-H/1948/04_67d
171.Sendai-H/1948/04_74d
172.Sendai-H/1948/04_81d
173.Sendai-H/1948/04_88d
174.Sendai-H/1948/04_95d
175.Sendai-H/1948/04_102d'''

samples = {}
for line in samples_str.split('\n'):
    id, sample = line.split('.')
    sample = sample.replace('Sendai-H/', '').replace('/4_', '-').replace('/04_', '-')
    samples[id] = sample


MAX_N = 3
MIN_QUAL = 20

seq1_CAA = 'CATTGCAATTTGT'
seq1_TAA = 'CATTGTAATTTGT'
seq2_ori = 'TGAGAGATAA'
seq2_del = 'TGAGAAAA'


def complement(seq):
    return seq.translate(str.maketrans('ATGCatgc', 'TACGtacg'))[::-1]


seq1_muts = ['CAA', 'TAA', 'other']
seq2_muts = ['ori', 'ATT', 'del', 'other']
combinations = []
for s1 in seq1_muts:
    for s2 in seq2_muts:
        combinations.append('%s:%s' % (s1, s2))



data = []
for id in samples:
    sample = samples[id]
    fq1 = '/path/to/R1.fastq.gz'
    fq2 = '/path/to/R2.fastq.gz'
    
    infile1 = gzip.open(fq1, 'rt')
    infile2 = gzip.open(fq2, 'rt')
    counter = collections.Counter()
    for line1, line2 in zip(infile1, infile2):
        if line1[:7] == '@M02121':
            seq1 = next(infile1).strip()
            seq2 = next(infile2).strip()
            next(infile1)
            next(infile2)
            qual1 = next(infile1).strip()
            qual2 = next(infile2).strip()
            seq1 = seq1[15:-15]
            seq2 = seq2[15:-15]
            qual1 = qual1[15:-15]
            qual2 = qual2[15:-15]
            # check N
            if seq1.count('N') >= MAX_N:
                continue
            if seq2.count('N') >= MAX_N:
                continue
            # check seq qual
            skip = False
            for c in qual1:
                qual = ord(c) - 33
                if qual < MIN_QUAL:
                    skip = True
                    break
            for c in qual2:
                qual = ord(c) - 33
                if qual < MIN_QUAL:
                    skip = True
                    break
            if skip:
                continue
            # analyze
            seq1_mut = 'other'
            if seq1_CAA in seq1:
                seq1_mut = 'CAA'
            elif seq1_TAA in seq1:
                seq1_mut = 'TAA'
            seq2_mut = 'other'
            if seq2_ori in seq2:
                seq2_mut = 'ori'
            elif seq2_del in seq2:
                seq2_mut = 'del'
            res = '%s:%s' % (seq1_mut, seq2_mut)
            counter[res] += 1
    
    tmp = []
    for c in combinations:
        tmp.append(counter[c])
    data.append(tmp)

df = pd.DataFrame(data, index = samples.values(), columns = combinations)
print(df)
df.to_csv('read_counts.tsv', sep = '\t')
