#!/usr/bin/python

import os
import Queue
import threading
import traceback

from build_interface_files import build_interface_files
from build_surface_files import build_surface_files
import common
import consts
from exec_silent import exec_silent

def get_siteengine_sol_filename(surface, interface):
    x1 = surface.pdb[:5]
    if len(x1) < 5:
      x1 += '_' * (5 - len(x1))
    x2 = interface.pdb[:5]
    if len(x2) < 5:
      x2 += '_' * (5 - len(x2))
  
    return '%s_%s_2_sol.res' % (x1,x2)
  

class SimpleRunner:
  def execute_single_candidate(self, interface, surface, result_parser):
    # Run site engine
    intreface_params = '%s %s %s' % (interface.output_file(), interface.receptor_surface_file(), interface.site_file())
    surface_params = '%s %s %s' % (surface.output_file(), surface.surface_file(), surface.surface_file())
    exec_silent('%s %s %s' % (consts.SITE_ENGINE_EXE,
                            surface_params,
                            intreface_params))
    
    # Resolve solutions filename
#    tmp_result_filename = '%s__%s__2_sol.res' % (surface.pdb, interface.pb)
    tmp_result_filename = get_siteengine_sol_filename(surface, interface)
    result_filename = common.get_solutions_filename(interface, surface)
    exec_silent('mv %s %s' % (tmp_result_filename, result_filename))
    
    if result_parser:
      result_parser.add_results(interface, surface, result_filename)

  def execute(self, interface, candidate_surfaces, result_parser=None):
    build_interface_files(interface)
    for surface in candidate_surfaces:
      # Build input files for siteengine
      build_surface_files(surface)
      self.execute_single_candidate(interface, surface, result_parser)



class ParallelExectureThread(threading.Thread):
  def __init__(self, runner, result_parser):
    threading.Thread.__init__(self)
    self.runner = runner
    self.result_parser = result_parser

  def run(self):
    while not self.runner.queue.empty():
      interface, surface = self.runner.queue.get()
      build_surface_files(surface)
      self.runner.execute_single_candidate(interface, surface, self.result_parser)

      
### Runner for parallel executions (with threads) on the same machine
class ParallelRunner(SimpleRunner):
  def __init__(self, num_threads=10):
    self.num_threads = num_threads
    self.queue = Queue.Queue()

  def execute(self, interface, candidate_surfaces, result_parser=None):
    # Build input files for siteengine
    build_interface_files(interface)

    for surface in candidate_surfaces:
      self.queue.put((interface, surface)) 

    # Create executing threads
    threads = [ParallelExectureThread(self, result_parser) for i in range(self.num_threads)]

    # Start threads
    for thread in threads:
      thread.start()

    # Wait for threads done
    for thread in threads:
      thread.join()

class ReverseParallelExectureThread(threading.Thread):
  def __init__(self, runner, interfaces, surface, result_parser):
    threading.Thread.__init__(self)
    self.runner = runner
    self.interfaces = interfaces
    self.surface = surface
    self.result_parser = result_parser

  def run(self):
    for interface in self.interfaces:
      build_interface_files(interface)
      self.runner.execute_single_candidate(interface, self.surface, self.result_parser)

      
### Runner for parallel executions (with threads) on the same machine
class ReverseParallelRunner(SimpleRunner):
  def __init__(self, num_threads=10):
    self.num_threads = num_threads

  def execute(self, interfaces, surface, result_parser=None):
    # Build input files for siteengine
    build_surface_files(surface)

    # Create executing threads
    threads = [ReverseParallelExectureThread(self, interfaces[i::self.num_threads], surface, result_parser) for i in range(self.num_threads)]

    # Start threads
    for thread in threads:
      thread.start()

    # Wait for threads done
    for thread in threads:
      thread.join()
