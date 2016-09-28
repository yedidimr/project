#!/usr/bin/python

import os 
import re
import pdb_cache
import random

ATOM_LINE_RE = "^ATOM *\d* *\w* *[A-Z]{3,4} %s .*$|^TER *\d* *\w* *[A-Z]{3,4} %s .*$"#|^HETATM *\d* *\w* *[A-Z]{3,4} %s .*$|^CONECT .*$"
AMINO_ACID_INDEX = "^ATOM *\d* *\w* *[A-Z]{3,4} %s *(-?\d*)"
AMINO_ACID_TER_INDEX = "^TER *\d* *\w* *[A-Z]{3,4} %s *(-?\d*)"
UB_LEN_LIM = (68,82)
UB_SPECIAL_AMINO_ACID_NAME = "ILE"
UB_SPECIAL_AMINO_ACID_INDEX = 44
AMINO_ACID_NAME = "^ATOM *\d* *\w* *([A-Z]{3,4}) %s"

def is_chain_ubiquitin(lines, chain_id):
  chain_start_index = int(re.findall(AMINO_ACID_INDEX %chain_id, lines[0])[0])
  chain_end_index = int(re.findall(AMINO_ACID_TER_INDEX %chain_id, lines[-1])[0])
  chain_len = chain_end_index - chain_start_index
  print "Chain %s len is %d" % (chain_id, chain_len)

  if chain_len > UB_LEN_LIM[0] and chain_len < UB_LEN_LIM[1]:
      for line in lines:
          current_index = int(re.findall(AMINO_ACID_INDEX % chain_id, line)[0])
          
          if current_index - chain_start_index == UB_SPECIAL_AMINO_ACID_INDEX - 1:
              amino_acid_name = re.findall(AMINO_ACID_NAME % chain_id, line, re.MULTILINE)[0]
              if UB_SPECIAL_AMINO_ACID_NAME == amino_acid_name:
                  return True

  return False



def get_pdb_normalized_name(pdb_file):
  basename = os.path.basename(pdb_file)
  file_name, file_extension = os.path.splitext(basename)

  if file_extension.lower() == '.pdb':
    return file_name
  return basename

def get_surface_name(pdb_file, chains):
  return get_pdb_normalized_name(pdb_file) + '_' + ''.join(chains)

def get_interface_name(pdb_file, ligand_chain, receptor_chains):
  return get_pdb_normalized_name(pdb_file) + '_' + ''.join(receptor_chains) + '_' + ligand_chain

def get_solutions_filename(interface, surface):
  return '%s_%s.res' % (surface.get_name(), interface.get_name())

class Surface:
  def __init__(self, pdb, chains):
    if pdb:
      self.pdb = pdb
      self.pdb_file = pdb_cache.get_pdb_file(pdb)
    self.chains = chains
    self.name = None
  
  def get_name(self):
    if self.name:
      return self.name
    return get_pdb_normalized_name(self.pdb_file) + '_' + ''.join(self.chains)

  def output_file(self):
    return '%s/%s.pdb.gz' % (self.get_name(), self.get_name())

  def surface_file(self):
    return '%s/%s.pdb'

  def surface_pdb(self):
    return '%s/%s.pdb' % (self.get_name(), self.get_name())

  def surface_pdb_with_hydrogens(self):
    return '%s/%s.pdb.HB' % (self.get_name(), self.get_name())

  def surface_file(self):
    return '%s/%s.pdb.ms.gz' % (self.get_name(), self.get_name())

  def set_name(self, name):
    self.name = name

class KnownLigandSurface(Surface):
  def __init__(self, pdb, chains, known_ligand):
    Surface.__init__(self, pdb, chains)
    self.known_ligand = known_ligand

class SurfaceFromFile(Surface):
  def __init__(self, pdb_file, chains):
    Surface.__init__(self, None, chains)
    self.pdb_file = pdb_file
    self.pdb = get_pdb_normalized_name(pdb_file)

class PatchDockLigand(Surface):
  def __init__(self, pdb_file):
    Surface.__init__(self, None, None)
    self.pdb_file = pdb_file
    self.pdb = self.get_name()

  def get_name(self):
    return os.path.basename(self.pdb_file).replace(".pdb","")

class KnownUBQSurface(Surface):
  def __init__(self, pdb, chains, ubq_chain):
    Surface.__init__(self, pdb, chains)
    self.ubq_chain = ubq_chain

  def get_name(self):
    bla =  get_pdb_normalized_name(self.pdb_file) + '_' + ''.join(self.chains) + '_ubq_' + self.ubq_chain
    return bla

class KnownUBQSurfaceFromFile(KnownUBQSurface):
   def __init__(self, pdb_file, chains, ubq_chain):
      KnownUBQSurface.__init__(self, None, chains, ubq_chain)
      self.pdb_file = pdb_file
      self.pdb = get_pdb_normalized_name(pdb_file)


class Interface:
  def __init__(self, pdb, ligand_chain, receptor_chains):
    if pdb:
      self.pdb = pdb
      self.pdb_file = pdb_cache.get_pdb_file(pdb)
      if not self.validate_ub_chain(ligand_chain):
        print("chain %s is not ubiquitin!" % ligand_chain)
        exit(2)
      else: # delete else
        print("chain %s is ubiquitin!" % ligand_chain)
    self.ligand_chain = ligand_chain
    self.chains = [self.ligand_chain]
    self.receptor_chains = receptor_chains
    self.receptor_pdb_file = None
    self.ligand_pdb_file = None

  def validate_ub_chain(self, chain_id):
    """
    filepath - a path of a pdb file
    ids - a list of the chainsid
    open_outputfile - a path to output file (if not specified output file will be filepath_ID_ID...out.pdb)
    """
    # extract atom lines from pdb file
    with open(self.pdb_file) as f:
      pdb_data = f.read()
    f.close()

    # chain_atoms = re.findall(ATOM_LINE_RE %( chain_id), pdb_data, re.MULTILINE)
    chain_atoms = re.findall(ATOM_LINE_RE %(chain_id, chain_id), pdb_data, re.MULTILINE)
    if is_chain_ubiquitin(chain_atoms, chain_id):
        return True
    return False

# class Interface:
#   def __init__(self, pdb, ligand_chain, receptor_chains):
#     if pdb:
#       self.pdb = pdb
#       self.pdb_file = pdb_cache.get_pdb_file(pdb)
#     self.ligand_chain = ligand_chain
#     self.chains = [self.ligand_chain]
#     self.receptor_chains = receptor_chains
#     self.receptor_pdb_file = None
#     self.ligand_pdb_file = None

  def get_name(self):
    return get_pdb_normalized_name(self.pdb_file) + '_' + ''.join(self.receptor_chains) + '_' + self.ligand_chain
  def output_file(self):
    return '%s/%s.pdb.gz' % (self.get_name(), self.get_name())

  def ligand_file(self):
    return '%s/%s_ligand.pdb' % (self.get_name(), self.get_name())

  def ligand_pdb(self):
    return '%s/%s_ligand.pdb' % (self.get_name(), self.get_name())
 
  def receptor_surface_file(self):
    return '%s/%s.pdb.ms.gz' % (self.get_name(), self.get_name())

  def site_file(self):
    return '%s/%s_site.ms.gz' % (self.get_name(), self.get_name())
   
  def ligand_pdb_with_hydrogens(self):
    return '%s/%s_ligand.pdb.HB' % (self.get_name(), self.get_name())

  def receptor_pdb_with_hydrogens(self):
    return '%s/%s.pdb.HB' % (self.get_name(), self.get_name())

class PatchDockReceptor(Surface):
  def __init__(self, receptor_pdb_file):
    Surface.__init__(self, None, None)
    self.pdb_file = receptor_pdb_file
    self.pdb = self.get_name()

  def get_name(self):
    return os.path.basename(self.pdb_file).replace(".pdb","")

  def ligand_pdb(self):
    return self.surface_pdb()


    

# class InterfaceFromFile(Interface):
#   def __init__(self, pdb_file, ligand_chain, receptor_chains):
#     Interface.__init__(self, None, ligand_chain, receptor_chains)
#     self.pdb_file = pdb_file
#     self.pdb = get_pdb_normalized_name(pdb_file)


class InterfaceFromFile(Interface):
  def __init__(self, pdb_file, ligand_chain, receptor_chains):
    Interface.__init__(self, None, ligand_chain, receptor_chains)
    self.pdb_file = pdb_file
    self.pdb = get_pdb_normalized_name(pdb_file)

    if not self.validate_ub_chain(ligand_chain):
      print ("chain %s is not ubiquitin!" % ligand_chain)
      exit(2)
    else: # delete else
      print("chain %s is ubiquitin!" % ligand_chain)

class DualInterface:
  def __init__(self, pdb, chains1, chains2):
    if pdb:
      self.pdb = pdb
      self.pdb_file = pdb_cache.get_pdb_file(pdb)
    self.chains1 = chains1
    self.chains2 = chains2

  def get_name(self):
    return get_pdb_normalized_name(self.pdb_file) + '_' + ''.join(self.chains1) + '_' + ''.join(self.chains2)

  def output_file1(self):
    return '%s/%s_rec1.pdb.gz' % (self.get_name(), self.get_name())

  def output_file2(self):
    return '%s/%s_rec2.pdb.gz' % (self.get_name(), self.get_name())

  def surface_file1(self):
    return '%s/%s_rec1.pdb.ms.gz' % (self.get_name(), self.get_name())

  def surface_file2(self):
    return '%s/%s_rec2.pdb.ms.gz' % (self.get_name(), self.get_name())

  def site_file1(self):
    return '%s/%s_site1.ms.gz' % (self.get_name(), self.get_name())

  def site_file2(self):
    return '%s/%s_site2.ms.gz' % (self.get_name(), self.get_name())
