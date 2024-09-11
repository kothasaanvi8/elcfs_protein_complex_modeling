import argparse, datetime
from functools import reduce
import glob
import itertools as it
from multiprocess import cpu_count, Event, Pool, Process
import numpy as np
import os
import pandas as pd
import pickle, random, re, requests, time, shutil, subprocess, sys
from tqdm import tqdm as progressMonitor
print(cpu_count())

#useful to keep track of sys vars so as to better monitor space remaining
sysVars = list(globals().keys())

rootDir = '/Users/wilkinsbusiness/Library/CloudStorage/GoogleDrive-elcfs.proteincomplex.modeling@gmail.com/My Drive/'
workDir = './' #phase2 directory
workDir_elcfs = rootDir + 'proteinPairs_complexMaps/' #modeling library
workDir_ph1 = rootDir + 'Primary Research/JLMwSCBC_notebook/' #phase 1 directory
workDir_other = rootDir + 'otherStudies/'

sys.path.insert(0, rootDir)
for p in workDir, workDir_elcfs, workDir_ph1, workDir_other: sys.path.append(p)
from util import modelEvaluating
from utils import operationsLocal as operations
from utils import reference, alertMe
pushoverKey_user = 'uith8rmy2npjj1oqpjwcanow3un984'
pushoverAPI = 'aw4v3424kaznrw598r6qge9icddwg7'

setupDir = workDir + 'setup/'
refDir = workDir + 'srcData/'
ph1Model_perfDir = workDir_ph1 + 'modelPerformance/'