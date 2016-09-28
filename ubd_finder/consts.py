#!/usr/bin/python

import os

CALC_ENERGY = False
CALC_RMSD = True
REFINE_ENERGY = False
ADD_HYDROGENS = True
OUTPUT_COMPLEXES = False

CWD = os.path.dirname(__file__) 

PDB_CACHE_DIR = '/tmp'

SITE_ENGINE_UTILS_DIR = os.path.join(CWD, 'utils/')

SITE_ENGINE_EXE = os.path.join(CWD, 'utils/siteEngine')
SITE_ENGINE_SITE_THRESHOLD = 5.0

#TRANS_RMSD_EXE = '/home/moldisk/astrinmi/UBD_Search/SiteEngineRunner/utils/transRMSD'
#TRANS_RMSD_MANY_EXE = '/home/moldisk/astrinmi/UBD_Search/SiteEngineRunner/utils/transRMSDMany'
CREATE_COMPLEX_EXE = os.path.join(CWD, 'utils/createComplex')
REDUCE_EXE = os.path.join(CWD, 'utils/reduce.3.14.080821.linuxi386')
ADD_HYDROGEN_EXE = os.path.join(CWD, 'utils/addHydrogens.pl')
RUN_MS_POINTS_EXE = os.path.join(CWD, 'utils/runMSPoints.pl')


