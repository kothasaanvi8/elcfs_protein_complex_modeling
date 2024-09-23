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
workDir_lib = workDir + 'utils/'

for p in rootDir, workDir, workDir_lib: sys.path.append(p)
from utils import modelEvaluating
from utils import operationsLocal as operations
from utils import reference, alertMe
print(cpu_count())
sysVars = list(globals().keys())
pushoverKey_user = 'uith8rmy2npjj1oqpjwcanow3un984'
pushoverAPI = 'aw4v3424kaznrw598r6qge9icddwg7'


setupDir = workDir + 'setup/'
refDir = workDir + 'srcData/'
ph1Model_perfDir = workDir_ph1 + 'modelPerformance/'