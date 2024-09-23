import argparse, datetime
from functools import reduce
import glob
import itertools as it
from multiprocess import cpu_count, Event, Pool, Process
import numpy as np
import os
import pandas as pd
import pickle, random, re, requests, time
import shutil
from sklearn.metrics import auc, average_precision_score, precision_recall_curve
import subprocess, sys
from tqdm import tqdm as progressMonitor

rootDir = '/Users/wilkinsbusiness/Documents/GitHub/'
workDir = rootDir + 'elcfs_protein_complex_modeling/'
workDir_elcfs = workDir + 'proteinPairs_complexMaps/' #modeling library
workDir_ph1 = workDir + 'Primary Research/JLMwSCBC_notebook/' #phase 1 directory
workDir_other = workDir + 'otherStudies/'

for p in rootDir, workDir, workDir_elcfs, workDir_ph1, workDir_other: sys.path.append(p)
from util import modelEvaluating
from utils import operationsLocal as operations
from utils import reference, alertMe
print(cpu_count())
sysVars = list(globals().keys())
pushoverKey_user = 'uith8rmy2npjj1oqpjwcanow3un984'
pushoverAPI = 'aw4v3424kaznrw598r6qge9icddwg7'


setupDir = workDir + 'setup/'
refDir = workDir + 'srcData/'
ph1Model_perfDir = workDir_ph1 + 'modelPerformance/'