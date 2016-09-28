#!/usr/bin/python

import os
import shutil
import sys

from common import get_interface_name
from exec_silent import exec_silent
import consts

def build_interface_files(interface, receptor_get_chains_cmd = 'get_chains.py', ligand_get_chains_cmd = 'get_chains.py'):
  cwd = consts.SITE_ENGINE_UTILS_DIR

  interface_name = interface.get_name()

  # Create a directory with the surface name
  if os.path.exists(interface_name):
    shutil.rmtree(interface_name)
  os.mkdir(interface_name)

  # Output file names. 
  receptor_pdb_output_file = '%s/%s.pdb' % (interface_name, interface_name)
  ligand_pdb_output_file = '%s/%s_ligand.pdb' % (interface_name, interface_name)
  receptor_surface_output_file = '%s/%s.pdb.ms' % (interface_name, interface_name)
  lignad_surface_output_file = '%s/%s_ligand.pdb.ms' % (interface_name, interface_name)
  site_output_file = '%s/%s_site.ms' % (interface_name, interface_name)
 

  print interface.receptor_chains
  # Create the ouput PDB file.
  if interface.receptor_pdb_file:
    exec_silent('cp %s %s' % (interface.receptor_pdb_file, receptor_pdb_output_file))
  else:
    exec_silent('python %s/%s %s %s > %s' % (cwd, receptor_get_chains_cmd, interface.pdb_file, ' '.join(interface.receptor_chains), receptor_pdb_output_file))

  if interface.ligand_pdb_file:
    exec_silent('cp %s %s' % (interface.ligand_pdb_file, ligand_pdb_output_file))
  else:
    exec_silent('python %s/%s %s %s > %s' % (cwd, ligand_get_chains_cmd,  interface.pdb_file, interface.ligand_chain, ligand_pdb_output_file))

  # Create the output surface with call to runMSPoint.pl
  exec_silent('cd %s && perl %s %s 10 1.4' % (interface_name, consts.RUN_MS_POINTS_EXE, os.path.basename(receptor_pdb_output_file)))
  exec_silent('cd %s && perl %s %s 10 1.4' % (interface_name, consts.RUN_MS_POINTS_EXE, os.path.basename(ligand_pdb_output_file)))

  # Create site file with call to surf_inter
  exec_silent('%s/surf_inter %s %s %s %s > %s' % (cwd,
                                                ligand_pdb_output_file,
                                                lignad_surface_output_file,
                                                receptor_surface_output_file,
                                                consts.SITE_ENGINE_SITE_THRESHOLD,
                                                site_output_file))

  # Create pdbs with hydrogens
  if consts.CALC_ENERGY or consts.REFINE_ENERGY or consts.ADD_HYDROGENS:
    exec_silent('%s %s %s' % (
        consts.ADD_HYDROGEN_EXE,
        receptor_pdb_output_file,
        interface.receptor_pdb_with_hydrogens()))
    exec_silent('%s %s %s' % (
        consts.ADD_HYDROGEN_EXE,
        ligand_pdb_output_file,
        interface.ligand_pdb_with_hydrogens()))
#    exec_silent('%s -OH -HIS -NOADJust -NOROTMET  %s > %s' % (
#        consts.REDUCE_EXE,
#        receptor_pdb_output_file,
#        interface.receptor_pdb_with_hydrogens()))
#    exec_silent('%s -OH -HIS -NOADJust -NOROTMET %s > %s' % (
#        consts.REDUCE_EXE,
#        ligand_pdb_output_file,
#        interface.ligand_pdb_with_hydrogens()))
 
  # gzip the output files
  exec_silent('gzip %s' % receptor_pdb_output_file)
  exec_silent('gzip %s' % ligand_pdb_output_file)
  exec_silent('gzip %s' % receptor_surface_output_file)
  exec_silent('gzip %s' % site_output_file)

  # Create again the ouput PDB file.
  if interface.receptor_pdb_file:
    exec_silent('cp %s %s' % (interface.receptor_pdb_file, receptor_pdb_output_file))
  else:
    exec_silent('python %s/%s %s %s > %s' % (cwd, receptor_get_chains_cmd, interface.pdb_file, ' '.join(interface.receptor_chains), receptor_pdb_output_file))

  if interface.ligand_pdb_file:
    exec_silent('cp %s %s' % (interface.ligand_pdb_file, ligand_pdb_output_file))
  else:
    exec_silent('python %s/%s %s %s > %s' % (cwd, ligand_get_chains_cmd,  interface.pdb_file, interface.ligand_chain, ligand_pdb_output_file))

  # Check for success 
  pdb_output_file_gz = '%s/%s.pdb.gz' % (interface_name, interface_name)
  surface_output_file_gz = '%s/%s.pdb.ms.gz' % (interface_name, interface_name)

  if not os.path.exists(pdb_output_file_gz):
    print 'Failure to create %s' % pdb_output_file_gz
    sys.exit(1)

  if not os.path.exists(surface_output_file_gz):
    print 'Failure to create %s' % surface_output_file_gz
    sys.exit(1)

def build_ubq_interface_files(interface):
  build_interface_files(interface, ligand_get_chains_cmd="get_chains_cut_ubq_tail.py")
  
if __name__ == '__main__':
  main(sys.argv)
