#!/usr/bin/python

import os
import sys

def get_chain(line):
  return line[21]

def get_amino_acid_number(line):
  pass

def put_atom_number(line, atom_number):
  number_str = str(atom_number)
  number_str = ' ' * (5-len(number_str)) + number_str
  return line[:6] + number_str + line[11:]  

def ubq_cut_tail_filter(line):
  if get_amino_acid_number(line) <= 72:
    return True

def is_header_line(line):
  return line.split()[0] in ['HEADER','TITLE','COMPND','SOURCE','KEYWDS','EXPDATA','AUTHOR','REVDAT','JRNL','REMARK','DBREF','SEQADV', 'SEQRES','FORMUL',
                         'HELIX','SHEET','CRYST1', 'ORIGX1', 'ORIGX2','ORIGX3','SCALE1','SCALE2','SCALE3']

def filter_chains(input_file, chains, addtional_filter_func=None):
  atom_number = 1

  for line in input_file:
    if is_header_line(line):
     print line,
    elif line.startswith('ATOM') or line.startswith('TER'):
      if chains == ['_'] or get_chain(line) in chains:
	if addtional_filter_func:
          if addtional_filter_func(line):
            output_line = put_atom_number(line, atom_number)
            print output_line,
            atom_number += 1
        else:
          output_line = put_atom_number(line, atom_number)
          print output_line,
          atom_number += 1

  print 'END'

if __name__ == '__main__':
  input_file = open(sys.argv[1])
  chains = sys.argv[2:] 
  filter_chains(input_file, chains)
