#!/usr/bin/python

import os

import consts
from exec_silent import exec_silent

def _get_pdb_location(pdb):
  return '%s/%s.pdb' % (consts.PDB_CACHE_DIR, pdb)

def _download_pdb_file(pdb):
  # Donwload PDB file from PDB site using wget
  exec_silent('wget -O %s.gz http://www.rcsb.org/pdb/files/%s.pdb.gz' % (_get_pdb_location(pdb), pdb))
  # Extract file using gzip
  exec_silent('gzip -d %s.gz' % _get_pdb_location(pdb))

def get_pdb_file(pdb):
  if not os.path.exists(_get_pdb_location(pdb)):
    _download_pdb_file(pdb)

  return _get_pdb_location(pdb)
   

def all_pdb_ids():
  return [x.strip() for x in open(consts.ALL_PDB_IDS_FILE).readlines()]
