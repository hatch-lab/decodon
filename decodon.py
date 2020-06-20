# coding=utf-8

"""
Given a DNA sequence, give the most divergent sequence that produces the same polypeptide

Usage:
  decodon.py INPUT [--N=1]

Arguments:
  INPUT  The input sequence

Options:
  -h --help  Show this screen.
  --N=<int>  [default: 1] The number of degenerate sequences to print

Output:
  Degenerate DNA sequence
"""

import sys
import os
from pathlib import Path
from docopt import docopt
from schema import Schema, And, Or, Use, SchemaError, Optional

import string, re

## Constants
def CODON_FILTER(x):
  codon_table = {
    ord(c): None for c in string.printable
  }
  codon_table.pop(ord('A'))
  codon_table.pop(ord('C'))
  codon_table.pop(ord('T'))
  codon_table.pop(ord('G'))
  return x.translate(codon_table)

def aa_to_codon(AA_MAP):
  CODON_MAP = {}
  for aa, codons in AA_MAP.items():
    for codon in codons:
      CODON_MAP[codon] = aa

  return CODON_MAP

AA_MAP = {
  'A': [ 'GCT', 'GCC', 'GCA', 'GCG' ],
  'R': [ 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG' ],
  'N': [ 'AAT', 'AAC' ],
  'D': [ 'GAT', 'GAC' ],
  'C': [ 'TGT', 'TGC' ],
  'Q': [ 'CAA', 'CAG' ],
  'E': [ 'GAA', 'GAG' ],
  'G': [ 'GGT', 'GGC', 'GGA', 'GGG' ],
  'H': [ 'CAT', 'CAC' ],
  'I': [ 'ATT', 'ATC', 'ATA' ],
  'L': [ 'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG' ],
  'K': [ 'AAA', 'AAG' ],
  'M': [ 'ATG' ],
  'F': [ 'TTT', 'TTC' ],
  'P': [ 'CCT', 'CCC', 'CCA', 'CCG' ],
  'S': [ 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC' ],
  'T': [ 'ACT', 'ACC', 'ACA', 'ACG' ],
  'W': [ 'TGG' ],
  'Y': [ 'TAT', 'TAC' ],
  'V': [ 'GTT', 'GTC', 'GTA', 'GTG' ],
  '*': [ 'TAA', 'TGA', 'TAG' ]
}

# Build reverse map (COON->AA)
CODON_MAP = aa_to_codon(AA_MAP)

## Process inputs
arguments = docopt(__doc__, version='1.0')

schema = Schema({
  'INPUT': And(Use(CODON_FILTER), lambda x: len(x) % 3 == 0 and len(x) > 0, error='INPUT must be a DNA sequence of codons'),
  '--N': And(Use(int), lambda n: n > 0, error='-N must be greater than 0')
})

try:
  arguments = schema.validate(arguments)
except SchemaError as error:
  print(error)
  exit(1)

sequence = CODON_FILTER(arguments['INPUT'])
output_num = arguments['--N']

# Get the codons
codons = re.findall('...', sequence)

# We want each codon to be as different from the start as possible

# From: https://en.wikipedia.org/wiki/Hamming_distance
def hamming_distance(string1, string2):
  dist_counter = 0
  for n in range(len(string1)):
    if string1[n] != string2[n]:
      dist_counter += 1
  return dist_counter

outputs = [('', 0)] * output_num

for codon in codons:
  # Get the alternative codons
  residue = CODON_MAP[codon]
  alts = AA_MAP[residue]
  
  # Calc the distance between this codon and all possible alternatives
  distances = {}
  for alt in alts:
    distances[alt] = hamming_distance(codon, alt)

  # Keep just the N highest scoring prefixes
  distances = sorted(distances.items(), key=lambda x: x[1], reverse=True)
  distances = distances[:output_num]

  diff = len(outputs) - len(distances)
  if diff > 0:
    distances = distances + [distances[0]]*diff

  for k,val in enumerate(distances):
    outputs[k] = ( outputs[k][0]+val[0], outputs[k][1]+val[1] )

# for codon in codons:
#   # Keep just the N highest scoring prefixes
#   outputs = sorted(outputs, key=lambda x: x[1], reverse=True)
#   outputs = outputs[:output_num]

#   # Get the alternative codons
#   residue = CODON_MAP[codon]
#   alts = AA_MAP[residue]
  
#   # Calc the distance between this codon and all possible alternatives
#   distances = {}
#   for alt in alts:
#     distances[alt] = hamming_distance(codon, alt)

#   # Keep just the N highest scoring prefixes
#   distances = sorted(distances.items(), key=lambda x: x[1], reverse=True)
#   distances = distances[:output_num]
  
#   # Keep a copy of the current prefixes
#   copy = outputs.copy()

#   # Rebuild the prefixes by appending the highest N scoring alt codons to each existing prefix
#   outputs.clear()
#   for distance in distances:
#     outputs = outputs + [ (output[0] + distance[0], output[1] + distance[1]) for output in copy ]

# # Keep just the N highest scoring prefixes
# outputs = sorted(outputs, key=lambda x: x[1], reverse=True)
# outputs = outputs[:output_num]

print("{} ({})".format('Sequence', 'Score'))
for output in outputs:
  print("{} ({})".format(output[0], output[1]))