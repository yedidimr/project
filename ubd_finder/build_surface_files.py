#!/usr/bin/python

import os
import shutil
import sys

from common import get_surface_name
import consts
from exec_silent import exec_silent

def build_surface_files(surface):
  cwd = consts.SITE_ENGINE_UTILS_DIR

  surface_name = surface.get_name()

  # Create a directory with the surface name
  if os.path.exists(surface_name):
    shutil.rmtree(surface_name)
  os.mkdir(surface_name)

  # Output file names. 
  pdb_output_file = '%s/%s.pdb' % (surface_name, surface_name)
  surface_output_file = '%s/%s.pdb.ms' % (surface_name, surface_name)
 

  # Create the ouput PDB file.
  # Output PDB will contain all atoms from a specific protein chain.
  print surface
  if surface.chains:
    exec_silent('python %s/get_chains.py %s %s > %s' % (cwd, surface.pdb_file, ' '.join(surface.chains), pdb_output_file))
  else:
    exec_silent('cp %s %s' % (surface.pdb_file, pdb_output_file))

  
  # Seems like two runMSPoints.pl cannot run on the same time. WTF?
  # Create the output surface with call to runMSPoint.pl
  exec_silent('cd %s ; perl %s %s 10 1.4' % (surface_name, consts.RUN_MS_POINTS_EXE, os.path.basename(pdb_output_file)))

  # Create pdbs with hydrogens
  if consts.CALC_ENERGY or consts.REFINE_ENERGY or consts.ADD_HYDROGENS:
    exec_silent('%s %s %s' % (
        #consts.REDUCE_EXE,
        consts.ADD_HYDROGEN_EXE,
        pdb_output_file,
        surface.surface_pdb_with_hydrogens()))
#    exec_silent('%s -OH -HIS -NOADJust -NOROTMET %s > %s' % (
#        consts.REDUCE_EXE,
#        pdb_output_file,
#        surface.surface_pdb_with_hydrogens()))
 
  # gzip the output files
  exec_silent('cp %s %s.tmp' % (pdb_output_file, pdb_output_file))
  exec_silent('cp %s %s.tmp' % (surface_output_file, surface_output_file))

  exec_silent('gzip %s' % pdb_output_file)
  exec_silent('gzip %s' % surface_output_file)

  exec_silent('mv %s.tmp %s' % (pdb_output_file, pdb_output_file))
  exec_silent('mv %s.tmp %s' % (surface_output_file, surface_output_file))

  # Output PDB will contain all atoms from a specific protein chain.
  if surface.chains:
    exec_silent('python %s/get_chains.py %s %s > %s' % (cwd, surface.pdb_file, ' '.join(surface.chains), pdb_output_file))
  else:
    exec_silent('cp %s %s' % (surface.pdb_file, pdb_output_file))
  
  # Check for success 
  pdb_output_file_gz = '%s/%s.pdb.gz' % (surface_name, surface_name)
  surface_output_file_gz = '%s/%s.pdb.ms.gz' % (surface_name, surface_name)

  if not os.path.exists(pdb_output_file_gz):
    print 'Failure to create %s' % pdb_output_file_gz
    sys.exit(1)

  if not os.path.exists(surface_output_file_gz):
    print 'Failure to create %s' % surface_output_file_gz
    sys.exit(1)
  
  

if __name__ == '__main__':
  main(sys.argv)
