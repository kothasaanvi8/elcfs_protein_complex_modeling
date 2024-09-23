#Collection of functions for simple tasks
from bs4 import BeautifulSoup
import glob
import itertools as it
import json
from multiprocessing import cpu_count, Pool
import numpy as np
import os
import pandas as pd
import pickle, re, requests
from requests.adapters import HTTPAdapter, Retry
from sklearn.metrics import auc, precision_recall_curve
from tqdm import tqdm as progressMonitor


def fmtTbl(x):
    if isinstance(x, (int, np.int64)):
        return '{:,d}'.format(x)

    elif isinstance(x, float):
        return '{:.2f}'.format(x)

def flatten(lst):
  return [item for sublist in lst for item in sublist]

def mapProts(ids, idsFmt, targetFmt, printFormats=False):
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    uniprotURL = 'https://rest.uniprot.org'

    if printFormats:
        xRef_dbsURL = \
            'https://rest.uniprot.org/database/' + \
            'stream?format=json&query=%28*%29+AND+%28category_exact%3A%22Genome+annotation+databases%22%29'
        data = json.loads(requests.get(xRef_dbsURL).text)
        xRef_dbsTable = pd.DataFrame(data['results'])
        xRef_dbsTable.drop(
            columns=['id', 'pubMedId', 'doiId', 'linkType', 'servers', 'dbUrl', 'category', 'statistics'],
            inplace=True)
        xRef_dbsTable.drop_duplicates(subset=['name'], inplace=True, ignore_index=True)

        return xRef_dbsTable

    # filter results by organism (Homo-Sapiens) and return in tab-delimited form
    requestSuccess = False
    while not requestSuccess:
        try:
            condAppend = '?format=tsv&query=%28model_organism%3A9606%29'
            respMap = requests.post(
                '{0}/idmapping/run'.format(uniprotURL),
                data={'from': idsFmt, 'to': targetFmt, 'ids': ','.join(ids)})
            job = respMap.json()['jobId']
            jobURL = '{0}/idmapping/details/{1}'.format(uniprotURL, job)
            respJob = session.get(jobURL)
            respJob_streaming = \
                respJob.json()['redirectURL'].replace('/results/', '/results/stream/')
            protMapping_raw = requests.get(respJob_streaming + condAppend).text
            protMapping = \
                pd.DataFrame(
                    [(map[0], map[1]) for map in [line.split('\t')
                                                  for line in protMapping_raw.split('\n')[1:] if line]],
                columns=[idsFmt, targetFmt])
            requestSuccess = True
        except:
            print('resubmitting job...')

    protMapping = protMapping[[idsFmt, targetFmt]].groupby(idsFmt, as_index=False).agg(','.join)
    protMapping[targetFmt] = [e if ',' not in e else '(' + e + ')' for e in protMapping[targetFmt]]

    return protMapping

#useful for mergers and other database-like operations in which ID order doesn't matter,
# e.g., [id1, id2] vs [id2, id1]
def freezePairs(df, col1='id1', col2='id2', pool=False):
    if pool:
        print('pooling...')
        numWorkers = cpu_count() - 1
        with Pool(numWorkers) as pool:
            pairs = [[row[col1], row[col2]] for _, row in df.iterrows()]
            pairsFrozen = list(progressMonitor(pool.imap(frozenset, pairs), total=len(pairs)))
        df.insert(2, 'pairsFrozen', pairsFrozen)

    else:
        df.insert(2, 'pairsFrozen', df.loc[:, [col1, col2]].apply(frozenset, axis=1))

    return df

#generates n choose 2 pairs from list of input IDs (n choose 2)
# Note: as n grows, function becomes more expensive
def generatePairs(ids):
    return set([frozenset(pair) for pair in list(it.combinations(ids, 2))])

#review this function (below)
def setupSpaces_pairs(pairsLsts, dir='./', names=['drew2021', 'expandedPairs']):

    mainDir = dir + 'otherStudies/'
    suffix = '-set/'

    #generate pickle file for union set of all pairs
    unionSet = set().union(*pairsLsts)
    uDir = mainDir + 'unionPairs' + suffix

    #generate pickle file for intersection set of all pairs
    intersectionSet = set.intersection(*pairsLsts)
    iDir = mainDir + 'intersectionSet_' + '_'.join(names) + suffix

    pickle.dump(unionSet, open(uDir + 'pairs.pkl', 'wb'))
    pickle.dump(intersectionSet, open(iDir + 'pairs.pkl', 'wb'))

def download(workDir, src, url, overwrite=False):
    path = workDir + 'srcData/{0}/'.format(src)
    os.makedirs(path, exist_ok=True)
    filename = path + url.split('/')[-1]
    if not glob.glob(filename) or overwrite:
        res = requests.get(url)
        with open(filename, 'wb') as f:
          f.write(res.content)
        print(filename)

    return filename

def findDownloadables(url, keywords=None):
    resp = requests.get(url)
    respParse = BeautifulSoup(resp.content, 'html.parser')
    htmlOutput = [link.get('href') for link in respParse.find_all('a')]
    filenames = [link for link in htmlOutput if '.' in link.split('/')[-1]]
    if keywords:
        filenames = [link for link in filenames if any([kw in link for kw in keywords])]

    return filenames

def process_pattern(args):
    patt, line = args
    match = patt.search(line)
    if match:
        return [e for e in match.groups() if e and e != "'"]
    else:
        return []

def find_saved_output_files(script_path):

    if not os.path.isfile(script_path):
        print(f"Script not found: {script_path}")

    """
    patts:
    open('filename') or open(var)
    pandas' to_csv('filename') or to_csv(var)
    pandas' to_excel('filename') or to_excel(var)
    pickle.dump(data, 'filename') or pickle.dump(data, var)
    numpy's np.save('filename') or np.save(var)
    numpy's np.savetxt('filename') or np.savetxt(var)
    matplotlib's plt.savefig('filename') or plt.savefig(var)
    imageio's imsave('filename') or imsave(var)
    PIL's Image.save('filename') or Image.save(var)
    """

    # patts for common file-saving methods in .py
    patts = [
        re.compile(r"open\(([\'\"]?)(.*?)\1[,\)]"),
        re.compile(r"to_csv\(([\'\"]?)(.*?)\1[,\)]"),
        re.compile(r"to_excel\(([\'\"]?)(.*?)\1[,\)]"),
        re.compile(r"pickle\.dump\([^,]+,\s*(['\"]?)(.*?)\1[,\)]"),
        re.compile(r"np\.save\(([\'\"]?)(.*?)\1[,\)]"),
        re.compile(r"np\.savetxt\(([\'\"]?)(.*?)\1[,\)]"),
        re.compile(r"plt\.savefig\(([\'\"]?)(.*?)\1[,\)]"),
        re.compile(r"imageio\.imsave\(([\'\"]?)(.*?)\1[,\)]"),
        re.compile(r"Image\.save\(([\'\"]?)(.*?)\1[,\)]")
    ]

    """
    # patts for common file-saving methods in .json (jupyter notebook files)
    patts = [
        re.compile(r".*open\((\"|\')?(.*?)\1[,\)]"),
        re.compile(r".*to_csv\((\"|\')?(.*?)\1[,\)]"),
        re.compile(r".*to_excel\((\"|\')?(.*?)\1[,\)]"),
        re.compile(r".*pickle\.dump\([^,]+,\s*(\"|\')?(.*?)\1[,\)]"),
        re.compile(r".*np\.save\((\"|\')?(.*?)\1[,\)]"),
        re.compile(r".*np\.savetxt\((\"|\')?(.*?)\1[,\)]"),
        re.compile(r".*plt\.savefig\((\"|\')?(.*?)\1[,\)]"),
        re.compile(r".*imageio\.imsave\((\"|\')?(.*?)\1[,\)]"),
        re.compile(r".*Image\.save\((\"|\')?(.*?)\1[,\)]")
    ]
    """

    # Checks for jupyter notebook file
    if script_path.endswith('.ipynb'):
        with open(script_path, 'r') as f:
            notebook_data = json.load(f)

        # Access different parts of the notebook
        cells = notebook_data['cells']
        lines = \
            flatten([[l for l in cell['source']]
                          for cell in cells if cell['cell_type'] == 'code'])

  # Checks for .py file
    else:
        with open(script_path, 'r') as file:
            lines = file.readlines()

    print('pooling...')
    numWorkers = cpu_count() - 1
    with Pool(numWorkers) as pool:
        output_files = list(progressMonitor(
            pool.imap(process_pattern,
                      [(patt, line) for patt in patts for line in lines]),
            total=len([(patt, line) for patt in patts for line in lines])))

    return  flatten(output_files)
