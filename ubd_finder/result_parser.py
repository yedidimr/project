#!/usr/bin/python

import commands
import time

from common import Surface, KnownUBQSurface
import consts
from exec_silent import exec_silent

INF = 999999999999999999999999999999999

class Result:
  score = None
  num_matches = None
  transform = None

  def __init__(self, interface, surface, result_raw_data, rank):
    # Parse result data from a raw data of a single result in the solution file
    self.interface = interface
    self.surface = surface

    score_line = result_raw_data[7].strip()
    matches_lines = result_raw_data[9].strip()
    transform_line = result_raw_data[2].strip()
    aa_matches_lines = result_raw_data[13:] #amino acid

    self.rank = rank
    self.score = float(score_line.split(':')[1].strip())
    self.matches = int(matches_lines.split(':')[1].strip())
    assert len(aa_matches_lines) == self.matches
    self.all_aa = [aa_match.split()[0] for aa_match in aa_matches_lines]
    self.transform = str(transform_line.split(':')[1].strip())
    self.rmsd = INF
    self.energy = INF
    self.refined_energy = INF

  def _to_str(self):
    result_str = "%s %s Rank: %s Score: %s Matches: %s Transform: %s" % (self.surface.pdb, "".join(self.surface.chains), self.rank, self.score, self.matches, self.transform)
    if self.rmsd != INF:
      result_str += ' RMSD: %s' % self.rmsd
    if self.energy != INF:
      result_str += ' Energy: %s' % self.energy
    if self.refined_energy != INF:
      result_str += ' FireDock Refinemet: %s' % self.refined_energy

    return result_str

  def __repr__(self):
    return self._to_str();

  def __str__(self):
    return self._to_str()

class ResultParser:
  def __init__(self):
    self.results = []
    self.job_done = False

  def add_results(self, interface, surface, sol_file):
    self._parse_results(interface, surface, sol_file)

  def _parse_results(self, interface, surface, sol_file):
    solution_data = [x.strip() for x in open(sol_file).readlines()]

    # Read solution by solution from the solutions file
    solution_num = 0
    while True:
      try:
        # Next solution 
        start = solution_data.index('Solution Num : %d' % solution_num)
        end = solution_data.index('End of Match List')
      except ValueError:
        # Done. No more solutions.
        break

      # Add solution to the list of results
      self.parse_result(interface, surface, solution_data[start:end], solution_num)
     
      solution_data = solution_data[end+1:]
      solution_num += 1

  def parse_result(self, interface, surface, result_raw_data, solution_num):
      result = Result(interface, surface, result_raw_data, solution_num)
      self.results.append(result)

  def get_results_by_score(self, num_results=0):
    if num_results:
      return sorted(self.results, key=lambda result: result.score, reverse=True)[:num_results]
    else:
      return sorted(self.results, key=lambda result: result.score, reverse=True)

  def is_job_done(self):
    return self.job_done

  def set_job_done(self):
    self.job_done = True


class UBQResult(Result):
  def __init__(self, interface, ubq_surface, result_raw_data, solution_num):
    Result.__init__(self, interface, ubq_surface, result_raw_data, solution_num)

class UBQResultParser(ResultParser):
  def __init__(self):
    ResultParser.__init__(self)

  def _parse_results(self, interface, surface, sol_file):
    solution_data = [x.strip() for x in open(sol_file).readlines()]

    results = []
    
    # Read solution by solution from the solutions file
    solution_num = 0
    while True:
      try:
        # Next solution 
        start = solution_data.index('Solution Num : %d' % solution_num)
        end = solution_data.index('End of Match List')
      except ValueError:
        # Done. No more solutions.
        break

      # Add solution to the list of results
      result = None
      if isinstance(surface, KnownUBQSurface): 
        result = UBQResult(interface, surface, solution_data[start:end], solution_num)
      else:
        result = Result(interface, surface, solution_data[start:end], solution_num)
      results.append(result)

      solution_data = solution_data[end+1:]
      solution_num += 1
 
    if consts.CALC_RMSD: 
      if isinstance(surface, KnownUBQSurface):
        self._calc_rmsd_all_results(interface, surface, results)

    if consts.CALC_ENERGY:
      self._calc_energy_all_results(interface, surface, results)

    if consts.REFINE_ENERGY:
      for result in results[:50]:
        self._calc_refined_energy_all_results(interface, surface, [result])

    if consts.OUTPUT_COMPLEXES:
      self._build_all_complexes(interface, surface, results[:50])
    
    for result in results:
      self.results.append(result)

  def _calc_rmsd_all_results(self, interface, ubq_surface, results):
    transforms_tmp_file = '/tmp/transforms_%s' % time.time()
    open(transforms_tmp_file,'w').writelines([res.transform + '\n' for res in results])
    rmsds = commands.getoutput('%s %s %s %s %s %s' % (consts.TRANS_RMSD_MANY_EXE,
                                                      interface.pdb_file,
                                                      interface.ligand_chain,
                                                      ubq_surface.pdb_file,
                                                      ubq_surface.ubq_chain,
                                                      transforms_tmp_file))
    rmsds = [float(x.strip()) for x in rmsds.splitlines()]

    for i in range(len(results)):
      results[i].rmsd = rmsds[i]

  def _calc_energy_all_results(self, interface, ubq_surface, results):
    transforms_tmp_file = '/tmp/transforms_%s' % time.time()
    open(transforms_tmp_file,'w').writelines([res.transform + '\n' for res in results])
   
    energies = commands.getoutput('%s %s %s %s' % (consts.CALC_MULTI_ENERGIES_EXE,
                                                   ubq_surface.surface_pdb_with_hydrogens(),
                                                   interface.ligand_pdb_with_hydrogens(),
                                                   transforms_tmp_file))
    energies = [float(x.split('|')[0].strip()) for x in energies.splitlines()[1:]]

    for i in range(len(results)):
      results[i].energy = energies[i]

  def _calc_refined_energy_all_results(self, interface, ubq_surface, results):
    transforms_tmp_file = '/tmp/transforms_%s' % time.time()
    open(transforms_tmp_file,'w').writelines(['%s %s' % (i, results[i].transform) + '\n' for i in range(len(results))])

    params_tmp_file = '/tmp/firedock_params_%s' % time.time()
    result_tmp_file = '/tmp/firedock_res_%s' % time.time()

    exec_silent('%s %s %s U U AA %s %s 0 50 0.8 0 %s' % (
        consts.FIREDOCK_BUILD_PARAMS_EXE,
        ubq_surface.surface_pdb_with_hydrogens(),
        interface.ligand_pdb_with_hydrogens(),
        transforms_tmp_file,
        result_tmp_file,
        params_tmp_file))
   
    exec_silent('%s %s' % (consts.FIREDOCK_RUN_EXE, params_tmp_file)) 
  def _build_all_complexes(self, interface, surface, results):
    for idx, result in enumerate(results):
      output_pdb = "res_%d.pdb" % (idx+1)
      print "%s %s %s %s %s %s %s             " % (consts.CREATE_COMPLEX_EXE,
                                                   output_pdb,
                                                   interface.pdb_file,
                                                   interface.ligand_chain,
                                                   result.transform,
                                                   surface.pdb_file,
                                                   " ".join(surface.chains))
      commands.getoutput("%s %s %s %s %s %s %s" % (consts.CREATE_COMPLEX_EXE,
                                                   output_pdb,
                                                   interface.pdb_file,
                                                   interface.ligand_chain,
                                                   result.transform,
                                                   surface.pdb_file,
                                                   " ".join(surface.chains)))
 

  def get_results_by_rmsd(self, num_results=0):
    if num_results:
      return sorted(self.results, key=lambda result: result.rmsd)[:num_results]
    else:
      return sorted(self.results, key=lambda result: result.rmsd)

  def get_results_by_energy(self, num_results=0):
    if num_results:
      return sorted(self.results, key=lambda result: result.energy)[:num_results]
    else:
      return sorted(self.results, key=lambda result: result.energy)
 
  def get_results_by_refinement_energy(self, num_results=0):
    if num_results:
      return sorted(self.results, key=lambda result: result.refined_energy)[:num_results]
    else:
      return sorted(self.results, key=lambda result: result.refined_energy)

class ReverseResult:
  score = None
  num_matches = None
  transform = None

  def __init__(self, interface, surface, result_raw_data, rank):
    # Parse result data from a raw data of a single result in the solution file
    self.interface = interface
    self.surface = surface

    score_line = result_raw_data[7].strip()
    matches_lines = result_raw_data[9].strip()
    transform_line = result_raw_data[2].strip()

    self.rank = rank
    self.score = float(score_line.split(':')[1].strip())
    self.matches = int(matches_lines.split(':')[1].strip())
    self.transform = str(transform_line.split(':')[1].strip())
    self.rmsd = 10000

  def __repr__(self):
    return "Interface: %s %s_%s Rank: %s Score: %s Matches: %s" % (self.interface.pdb, "".join(self.interface.receptor_chains), self.interface.ligand_chain, self.rank, self.score, self.matches)

  def __str__(self):
    return "Interface: %s %s_%s Rank: %s Score: %s Matches: %s" % (self.interface.pdb, "".join(self.interface.receptor_chains), self.interface.ligand_chain, self.rank, self.score, self.matches)


class ReverseResultParser:
  def __init__(self):
    self.results = []

  def add_results(self, interface, surface, sol_file):
    self._parse_results(interface, surface, sol_file)

  def _parse_results(self, interface, surface, sol_file):
    solution_data = [x.strip() for x in open(sol_file).readlines()]

    # Read solution by solution from the solutions file
    solution_num = 0
    while True:
      try:
        # Next solution 
        start = solution_data.index('Solution Num : %d' % solution_num)
        end = solution_data.index('End of Match List')
      except ValueError:
        # Done. No more solutions.
        break

      # Add solution to the list of results
      self.parse_result(interface, surface, solution_data[start:end], solution_num)
     

      solution_data = solution_data[end+1:]
      solution_num += 1

  def parse_result(self, interface, surface, result_raw_data, solution_num):
      result = ReverseResult(interface, surface, result_raw_data, solution_num)
      self.results.append(result)

  def get_results_by_score(self, num_results=0):
    if num_results:
      return sorted(self.results, key=lambda result: result.score, reverse=True)[:num_results]
    else:
      return sorted(self.results, key=lambda result: result.score, reverse=True)

  def get_results_by_rmsd(self, num_results=0):
    if num_results:
      return sorted(self.results, key=lambda result: result.rmsd)[:num_results]
    else:
      return sorted(self.results, key=lambda result: result.rmsd)


class ReverseUBQResult(Result):
  def __init__(self, interface, ubq_surface, result_raw_data, solution_num):
    Result.__init__(self, interface, ubq_surface, result_raw_data, solution_num)
    self.rmsd = 999999999999
    self.energy = 99999999999

  def __repr__(self):
    return "Interface: %s %s_%s Rank: %s Score: %s Matches: %s RMSD: %s" % (self.interface.pdb, "".join(self.interface.receptor_chains), self.interface.ligand_chain, self.rank, self.score, self.matches, self.rmsd)

  def __str__(self):
    return "Interface: %s %s_%s Rank: %s Score: %s Matches: %s RMSD: %s" % (self.interface.pdb, "".join(self.interface.receptor_chains), self.interface.ligand_chain, self.rank, self.score, self.matches, self.rmsd)

class ReverseUBQResultParser(ResultParser):
  def __init__(self):
    ResultParser.__init__(self)

  def _parse_results(self, interface, surface, sol_file):
    solution_data = [x.strip() for x in open(sol_file).readlines()]

    results = []
    # Read solution by solution from the solutions file
    solution_num = 0
    while True:
      try:
        # Next solution 
        start = solution_data.index('Solution Num : %d' % solution_num)
        end = solution_data.index('End of Match List')
      except ValueError:
        # Done. No more solutions.
        break

      # Add solution to the list of results
      result = None
      if isinstance(surface, KnownUBQSurface): 
        result = ReverseUBQResult(interface, surface, solution_data[start:end], solution_num)
      else:
        result = ReverseResult(interface, surface, solution_data[start:end], solution_num)
      results.append(result) 

      solution_data = solution_data[end+1:]
      solution_num += 1
  
    if consts.CALC_RMSD:
      if isinstance(surface, KnownUBQSurface):
        self._calc_rmsd_all_results(interface, surface, results)

    if consts.CALC_ENERGY:
      self._calc_energy_all_results(interface, surface, results)

    if consts.OUTPUT_COMPLEXES:
      self._build_all_complexes(interface, surface, results[:10])
    
    for result in results:
      self.results.append(result)

  def _calc_rmsd_all_results(self, interface, ubq_surface, results):
    transforms_tmp_file = '/tmp/transforms_%s' % time.time()
    open(transforms_tmp_file,'w').writelines([res.transform + '\n' for res in results])
    rmsds = commands.getoutput('%s %s %s %s %s %s' % (consts.TRANS_RMSD_MANY_EXE,
                                                      interface.pdb_file,
                                                      interface.ligand_chain,
                                                      ubq_surface.pdb_file,
                                                      ubq_surface.ubq_chain,
                                                      transforms_tmp_file))
    rmsds = [float(x.strip()) for x in rmsds.splitlines()]

    for i in range(len(results)):
      results[i].rmsd = rmsds[i]

  def _calc_energy_all_results(self, interface, ubq_surface, results):
    transforms_tmp_file = '/tmp/transforms_%s' % time.time()
    open(transforms_tmp_file,'w').writelines([res.transform + '\n' for res in results])
   
    energies = commands.getoutput('%s %s %s %s' % (consts.CALC_MULTI_ENERGIES_EXE,
                                                   interface.ligand_pdb_with_hydrogens(),
                                                   surface.surface_pdb_with_hydrogens(),
                                                   transforms_tmp_file))
    print 'XXXXXXXXXXXXXXXXXXXXX', energies
  
  def _build_all_complexes(self, interface, surface, results):
    for idx, result in enumerate(results):
      output_pdb = "res_%d.pdb" % (idx+1)
      print "%s %s %s %s %s %s %s             " % (consts.CREATE_COMPLEX_EXE,
                                                   output_pdb,
                                                   interface.pdb_file,
                                                   interface.ligand_chain,
                                                   result.transform,
                                                   surface.pdb_file,
                                                   " ".join(surface.chains))
      commands.getoutput("%s %s %s %s %s %s %s" % (consts.CREATE_COMPLEX_EXE,
                                                   output_pdb,
                                                   interface.pdb_file,
                                                   interface.ligand_chain,
                                                   result.transform,
                                                   surface.pdb_file,
                                                   " ".join(surface.chains)))


  def get_results_by_rmsd(self, num_results=0):
    if num_results:
      return sorted(self.results, key=lambda result: result.rmsd)[:num_results]
    else:
      return sorted(self.results, key=lambda result: result.rmsd)

  def get_results_by_energy(self, num_results=0):
    if num_results:
      return sorted(self.results, key=lambda result: result.energy)[:num_results]
    else:
      return sorted(self.results, key=lambda result: result.energy)
 

def run_firedock_refinement(result):
  transforms_tmp_file = '/tmp/transforms_%s' % time.time()
  open(transforms_tmp_file,'w').writelines(['1 %s\n' % result.transform])

  params_tmp_file = '/tmp/firedock_params_%s' % time.time()
  result_tmp_file = '/tmp/firedock_res_%s' % time.time()

  exec_silent('%s %s %s U U AA %s %s 0 50 0.8 0 %s' % (
      consts.FIREDOCK_BUILD_PARAMS_EXE,
      result.surface.surface_pdb_with_hydrogens(),
      result.interface.ligand_pdb_with_hydrogens(),
      transforms_tmp_file,
      result_tmp_file,
      params_tmp_file))

  exec_silent('%s %s' % (consts.FIREDOCK_RUN_EXE, params_tmp_file)) 

  result_line = open('%s.ref' % result_tmp_file).readlines()[40]
  result.refined_energy = float(result_line.split('|')[5].strip())
