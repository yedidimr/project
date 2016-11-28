
import os
import sys
import gst
from ubd_finder import consts
from ubd_finder.build_surface_files import build_surface_files
from ubd_finder.build_interface_files import build_interface_files
from ubd_finder.common import Interface, Surface, KnownUBQSurface, SurfaceFromFile
from ubd_finder.result_parser import ResultParser, UBQResultParser
from ubd_finder.runners import SimpleRunner, ParallelRunner 

consts.OUTPUT_COMPLEXES = True


def get_results(interface, surface):

  result_parser = UBQResultParser()
  runner = ParallelRunner(1)
  runner.execute(interface, [surface], result_parser)
  
  output_filnes = [] 
  for x in result_parser.get_results_by_score(10):
    output_filnes.append(repr(x) + '\n')
  open('RESULTS_%s' % surface.get_name(), "w").writelines(output_filnes)

  all_aa=[]
  for res in result_parser.results:
    all_aa.extend(res.all_aa)

  return all_aa, result_parser.filename



if __name__ == '__main__':
  
  if len(sys.argv) != 2 and len(sys.argv) != 3:
    print "Usage: %s <PDB-id/PDB-path> [chains (default: _ for entrie file)]" % sys.argv[0]
    print "PDB-id/PDB-path is the pdb and chaing to search the ubiquitin binding domain"
    print "search is done the reference ubd pdbs which are hard-coded in the code"
    print "\nfor HHARI run:\n\t %s 4kbl A" % sys.argv[0]
    sys.exit(1)

  pdb = sys.argv[1]

  chains = None
  if len(sys.argv) == 2:
    chains = "_"
  else:
    chains = sys.argv[2]

  surf = None
  if len(pdb) == 4:
    surf = Surface(pdb, list(chains)) 
  else:
    if not os.path.exists(pdb):
       print "Error: not such file %s" % pdb
       sys.exit(1)

    surf = SurfaceFromFile(pdb, list(chains)) #pdb_file, chains):
  

  #                         Interface(pdb, ligand_chain, receptor_chains)
  interfaces = [ Interface('3K9P', 'B', ['A']), 
                          Interface('1yd8', 'U', ['H']),
                          Interface('1nbf', 'D', ['A']),
                          Interface('2qho', 'G', ['H']),
                          Interface('3i3t', 'B', ['A']),
                          Interface('3ihp', 'C', ['A']),
                          Interface('2ayo', 'B', ['A']),
                          Interface('2d3g', 'B', ['P']),
                          Interface('3nhe', 'B', ['A']),
                          Interface('2hd5', 'B', ['A']),
                          Interface('3tmp', 'B', ['A']),
                          Interface('2ibi', 'B', ['A']),
                        ]

  all_res = []
  num_of_appearence_per_pdb = dict()
  aa_to_pdb = dict()
  max_val = ""
  max_count = 0
  all_res_files = []
  for interface in interfaces:
    all_aa, filename = get_results(interface, surf)
    all_res_files.append(filename)
    all_res.append(all_aa)
    for aa in all_aa:
      if aa in aa_to_pdb:
        aa_to_pdb[aa].append(i)
      else:
        aa_to_pdb[aa] = [i]
      if aa in num_of_appearence_per_pdb:
        num_of_appearence_per_pdb[aa] +=1
        if num_of_appearence_per_pdb[aa] > max_count:
          max_count = num_of_appearence_per_pdb[aa]
          max_val = aa
      else:
        num_of_appearence_per_pdb[aa] = 1
        if num_of_appearence_per_pdb[aa] > max_count:
          max_count = num_of_appearence_per_pdb[aa]
          max_val = aa

  pdbs_per_num_of_appearence = dict()
  for pdb, count in num_of_appearence_per_pdb.items():
    if count not in pdbs_per_num_of_appearence:
      pdbs_per_num_of_appearence[count] = [pdb]
    else:
      pdbs_per_num_of_appearence[count].append(pdb)

  num_of_appearence = sorted(pdbs_per_num_of_appearence.keys(),reverse = True)
  top_num_of_appearence = num_of_appearence[:len(num_of_appearence)/4]
  most_popular_aa = [] # by descending order
  [most_popular_aa.extend(pdbs_per_num_of_appearence[count]) for count in top_num_of_appearence]



  print "num of results", len(all_res)

  print ""
  tree = gst.STree(all_res)
  lcs_list , lcs_set, depth = tree.lcs() # longest common ancestor
                                  # lcs(self, stringIdxs=-1), param stringIdxs: Optional: List of indexes of strings.
  
  print "!!!!", depth
  print "all results are saved to the following files. next run you can supply them"
  print "\n".join(all_res_files)

  import pdb
  pdb.set_trace()