
import os
import sys

from ubd_finder import consts
from ubd_finder.build_surface_files import build_surface_files
from ubd_finder.build_interface_files import build_interface_files
from ubd_finder.common import Interface, Surface, KnownUBQSurface, SurfaceFromFile
from ubd_finder.result_parser import ResultParser, UBQResultParser
from ubd_finder.runners import SimpleRunner, ParallelRunner 

consts.OUTPUT_COMPLEXES = True

if __name__ == '__main__':
  site_3k9p = Interface('3K9P', 'B', ['A'])

  if len(sys.argv) != 2 and len(sys.argv) != 3:
    print "Usage: %s <PDB-id/PDB-path> [chains (default: _ for entrie file)]" % sys.argv[0]
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

    surf = SurfaceFromFile(pdb, list(chains))
  


  result_parser = UBQResultParser()
  runner = ParallelRunner(1)
  runner.execute(site_3k9p, [surf], result_parser)
  
  output_filnes = [] 
  for x in result_parser.get_results_by_score(10):
    output_filnes.append(repr(x) + '\n')
  open('RESULTS_%s' % surf.get_name(), "w").writelines(output_filnes)

def get_results(interface, surface):

  result_parser = UBQResultParser()
  runner = ParallelRunner(1)
  runner.execute(interface, [surface], result_parser)
  
  output_filnes = [] 
  for x in result_parser.get_results_by_score(10):
    output_filnes.append(repr(x) + '\n')
  open('RESULTS_%s' % surface.get_name(), "w").writelines(output_filnes)

  all_aa=[]
  alldict = dict()
  for res in result_parser.results:
    all_aa.extend(res.all_aa)
    for aa in res.all_aa:
      if aa in alldict:
        alldict[aa] +=1
      else:
        alldict[aa] = 1
  return alldict, all_aa


def fingerprint(string, basis=2**16, r=2**32-3):
    """ used to computes karp-rabin fingerprint of the pattern
    employs Horner method (modulo r) """
    partial_sum = 0
    for x in string:
        partial_sum = (partial_sum*basis + ord(x)) % r
    return partial_sum

def text_fingerprint(string, length, basis=2**16, r=2**32-3):
    """ used to computes karp-rabin fingerprint of the text """
    f = []
    b_power = pow(basis,length-1,r)
    list.append(f, fingerprint(string[0:length], basis, r))
    # f[0] equals first text fingerprint 
    for s in range(1,len(string)-length+1):
        new_fingerprint = ((f[s-1] - ord(string[s-1])*b_power)*basis
                         +ord(string[s+length-1])) % r
            # compute f[s], based on f[s-1]
        list.append(f,new_fingerprint)# append f[s] to existing f       
    return f

def find_matches_KR(pattern, text, basis=2**16, r=2**32-3):
    if len(pattern) > len(text):
        return []
    p = fingerprint(pattern,basis,r)
    f = text_fingerprint(text,len(pattern),basis,r)
    matches = [s for s, f_s in enumerate(f) if f_s == p]
    # list comprehension 
    return matches

def longest_common_substring(s1, s2):
    m = [[0] * (1 + len(s2)) for i in xrange(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in xrange(1, 1 + len(s1)):
        for y in xrange(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return s1[x_longest - longest: x_longest]

def long_substr(data):
    substr = ''
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                if j > len(substr) and all(data[0][i:i+j] in x for x in data):
                    substr = data[0][i:i+j]
    return substr

if __name__ == '__main__':
  site_3k9p = Interface('3K9P', 'B', ['A']) # (pdb, ligand_chain, receptor_chains)
  site_1yd8 = Interface('1yd8', 'U', ['H']) # (pdb, ligand_chain, receptor_chains)
  # site_3k9p = Interface('1nbf', 'D', ['A']) # (pdb, ligand_chain, receptor_chains)
  # name = '1nbf'
  # chains =  ['D', 'A']
  # site_3k9p = InterfaceFromFile('../pdb/1nbf.pdb', 'D', ['A']) # (pdb, ligand_chain (ubiquitin), receptor_chains (ubiquiting binding domain)
  # site_3k9p = InterfaceFromFile('/tmp/1nbf.pdb', 'D', ['A']) # (pdb, ligand_chain, receptor_chains)
  # site_3k9p = InterfaceFromFile('HHARI_4kbl_A.pdb', 'A', 'A') # (pdb, ligand_chain, receptor_chains)
  interfaces = [site_1yd8,site_3k9p]                  
  # interfaces = [site_1yd8]                      
  # interfaces = [site_3k9p]                  

  
  if len(sys.argv) != 2 and len(sys.argv) != 3:
    print "Usage: %s <PDB-id/PDB-path> [chains (default: _ for entrie file)]" % sys.argv[0]
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
  
  all_res = []
  all_dict = dict()
  max_val = ""
  max_count = 0
  for interface in interfaces:
    _, all_aa = get_results(interface, surf)
    all_res.append(all_aa)
    for aa in all_aa:
      if aa in all_dict:
        all_dict[aa] +=1
        if all_dict[aa] > max_count:
          max_count = all_dict[aa]
          max_val = aa
      else:
        all_dict[aa] = 1
        if all_dict[aa] > max_count:
          max_count = all_dict[aa]
          max_val = aa

  print "END"
  import pdb
  pdb.set_trace()
  longest_common_substring(all_res[0], all_res[1])