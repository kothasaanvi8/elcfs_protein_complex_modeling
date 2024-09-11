#Constrain data and subsequent predictions to ground truth (CORUM) and relevant proteins
import glob
import numpy as np
import pandas as pd
import pickle
from functools import reduce
import math
from multiprocessing import cpu_count, Pool
from tqdm import tqdm as progressMonitor
from . import operations


#Downloads CORUM dataset (or loads if file already exists),
# generates all labeled pairs and converts complexes into unique searchable sets, and
# counts proteins, complexes, etc. Note: functionality for downloading latest CORUM repository is not possible
# because the manager(s) for the website, http://mips.helmholtz-muenchen.de/corum/, are recovering from a
# cyberattack. Refer to Elisabeth Gasteiger (https://www.biostars.org/p/9563995/)
class CORUM(object):

  def __init__(self, corumCore_excelFilename,
               nameCol='ComplexName', organismCol='Organism',
               subunitsCol='subunits(Entrez IDs)'):

    corumCols_select = [nameCol, organismCol, subunitsCol]
    self.nameCol, self.organismCol, self.subunitsCol = corumCols_select
    self.rawData = pd.read_excel(corumCore_excelFilename, usecols=corumCols_select)

    self.allCplx_names, self.allCplx_frozen, self.allCplx_species, \
      self.allCplx_total, self.allProts, self.organisms = \
        self.formatDF(
          self.rawData.dropna(subset=subunitsCol))

    self.humanComplexes_names, self.humanCplx_frozen, _, self.humanCplx_total, \
      self.humanProts, _ = self.formatDF(
          self.rawData.dropna(subset=subunitsCol).loc[
          self.rawData.dropna(subset=subunitsCol)[self.organismCol]=='Human', :])

    self.allPairs, self.allPairs_human = \
      operations.generatePairs(self.allProts), operations.generatePairs(self.humanProts)

  def dropEmpties_freeze2List(self, line):
    return frozenset([ele for ele in line.replace('None', '').split(';') if ele])

  def formatDF(self, df):
    cplxNames = [row[self.nameCol] for _, row in df.iterrows()]
    cplxFrozen = \
     [self.dropEmpties_freeze2List(row[self.subunitsCol])
     for _, row in df.iterrows()]
    cplxOrganisms = [row[self.organismCol] for _, row in df.iterrows()]
    numComplexes = sum([1 for _, row in df.iterrows()])
    prots = [str(ele) for ele in list(set().union(*
     [row[self.subunitsCol].replace('None', '').split(';')
     for _, row in df.iterrows()])) if ele]
    species = list(set([row[self.organismCol] for _, row in df.iterrows()]))

    return cplxNames, cplxFrozen, cplxOrganisms, numComplexes, prots, species


#Circumscribe space of protein-coding genes from latest data
class geneSpace(object):

  def __init__(self, geneInfo_filename, downloadLatest=False, workDir=None, geneInfo_url=None):

    if downloadLatest:
      print('proceeding will overwrite any existing file...')

    self.geneInfo_filename = geneInfo_filename
    self.downloadLatest = downloadLatest
    self.workDir = workDir
    self.geneInfo_url = geneInfo_url

  def pullReference(self):
    geneIDs_info = pd.read_csv(self.geneInfo_filename, sep='\t')
    geneIDs_info.GeneID = geneIDs_info.GeneID.astype(str)
    geneIDs_info = geneIDs_info
    geneIDs_all = list(geneIDs_info.GeneID.unique())
    print(len(geneIDs_all))

    geneIDs_proteinCoding = set(
      geneIDs_info.loc[geneIDs_info.type_of_gene=='protein-coding',
      'GeneID'].to_list())
    print(len(geneIDs_proteinCoding))

    self.geneIDs_info = geneIDs_info
    self.geneIDs_all = geneIDs_all
    self.geneIDs_proteinCoding = geneIDs_proteinCoding

  def make(self):
    if self.downloadLatest or not self.geneInfo_filename or not glob.glob(self.geneInfo_filename):
      if not self.geneInfo_url:
        self.geneInfo_url = \
          'https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'
      src = '{0}-geneData_global'.format(self.geneInfo_url.split('//')[1].split('/')[0])
      self.geneInfo_filename = operations.download(self.workDir, src, self.geneInfo_url, overwrite=True)
    self.pullReference()


class cmpStudies_props(object):
  def __init__(self, humap1Dir, humap2Dir, bcb19Dir_LM, stringdbDir, geneidDir, workDir='../otherStudies/',
               humap1Patt=None, humap2Patt=None, low_memory=False):
    if not low_memory:
      self.humap1Dir, self.humap2Dir, self.bcb19Dir_LM, self.stringdbDir, self.geneidDir, self.workDir = \
          humap1Dir, humap2Dir, bcb19Dir_LM, stringdbDir, geneidDir, workDir

      self.humap1Patt = humap1Patt if humap1Patt else ['*ppis*', '*.gz', '*prob*', '*clusters*']
      self.humap2Patt = humap2Patt if humap2Patt else ['*ppis*geneid*txt', '*featmat*', '*geneid*prob*', '*_complexes_*']

      self.bcb19Cplxs, self.bcb19Prots_cplxs = dict(), dict()

      entrezGene_data = self.geneidDir + 'Homo_sapiens.gene_info.gz'
      self.geneSpace = geneSpace(entrezGene_data, workDir=self.stringdbDir)
      self.geneSpace.make()

      self.numWorkers = cpu_count() - 1

      (self.humap1Mat_train, self.humap1Mat_feats, self.humap2Mat_train, self.humap2Mat_feats,
        self.humap1Prots_train, self.humap1Prots_feats, self.humap2Prots_train, self.humap2Prots_feats,
        self.humap1Pairs_labeled, self.humap2Pairs_labeled, self.humap1Pairs_feats, self.humap2Pairs_feats,
        self.humap1PPIs, self.humap2PPIs, self.humap1Preds, self.humap2Preds, self.humap1Cplxs, self.humap2Cplxs,
        self.humap1Prots_cplxs, self.humap2Prots_cplxs, self.pairsHumap1p2, self.protsHumap1p2, self.predsHumap1p2,
        self.bcb19Pairs, self.bcb19Pairs_feat, self.stringPreds, self.stringProts,
        self.stringCplxs, self.stringProts_cplxs, self.stringPPIs, self.summaryDesc, self.protsPairs_summary) = \
          None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, \
          None, None, None, None, None, None, None, None, None, None, None, None, None, None
    else:
      self.load()

  def load(self,):
    dir = './drive/MyDrive/elcfs_protein_complex_modeling/google_colab_notebooks/notebookQC/'
    return pd.read_pickle(dir + 'cmpStudies.pkl')

  def getProts(self, mat, cols):
    return set(mat[cols[0]].to_list()).union(set(mat[cols[1]].to_list()))

  def setupTraining(self, dir, patt):
    print('getting training pairs data...')
    pairs = {f: pd.read_csv(f, sep='\t| ', engine='python',
                            dtype='str', names=['id1', 'id2'])
             for f in glob.glob(dir + patt)}

    for f in pairs.keys():
      if 'neg' in f.split('/')[-1]:
        pairs[f].insert(2, 'label', 0)
      else:
        pairs[f].insert(2, 'label', 1)
    trainMat = pd.concat([df for df in pairs.values()], axis=0, ignore_index=True)

    print('generating geneid pairs frozensets...')
    trainMat = operations.freezePairs(trainMat, 'id1', 'id2', pool=True)
    trainProts = self.getProts(trainMat, ['id1', 'id2'])
    trainPairs = trainMat.pairsFrozen.to_list()

    return trainMat, trainProts, trainPairs

  def getFeat(self, dir, patt, idCols=None):
    print('getting feature pairs data...')
    idCols = idCols if idCols else ['id1', 'id2']
    print(idCols)
    print(dir + patt)
    print(glob.glob(dir + patt)[0])
    pairsMat = pd.read_csv(glob.glob(dir + patt)[0], sep=',', usecols=idCols, dtype='str')
    print('generating geneid pairs frozensets...')
    pairsMat = operations.freezePairs(pairsMat, idCols[0], idCols[1], pool=True)
    featProts = self.getProts(pairsMat, idCols)
    featPairs = pairsMat.pairsFrozen.to_list()

    return pairsMat, featProts, featPairs

  def getPPIs(self, trainMat, pairs):
    print('get PPIs...')
    ppis = trainMat.loc[trainMat.pairsFrozen.isin(pairs)].copy()
    print('Total PPIs: {0}'.format(ppis.label.sum()))

    return ppis

  def getPredictions(self, dir, patt):
    print('get predictions...')
    print(glob.glob(dir + patt))
    predsFile = glob.glob(dir + patt)[0]
    preds = pd.read_csv(predsFile, sep='\t',
                        header=None, names=['id1', 'id2', 'prob'], dtype={'id1': 'str', 'id2': 'str'})
    preds = operations.freezePairs(preds, 'id1', 'id2', pool=True)

    return preds

  def getCplx(self, dir, patt, tbl=True):
    print(glob.glob(dir + patt))
    cplxsFilename = glob.glob(dir + patt)[0]
    if not tbl:
      with open(cplxsFilename) as f:
        cplxs = [frozenset([str(ele) for ele in line.rstrip().split(' ')]) for line in f]
      protsCplxs = set().union(*cplxs)

    else:
      cplxsMat = pd.read_csv(cplxsFilename, sep=',')
      cplxsUni = [frozenset(cplx.rstrip().split(' ')) for cplx in cplxsMat['Uniprot_ACCs'].to_list()]
      cplxsGenename = [frozenset(cplx.rstrip().split(' ')) for cplx in cplxsMat['genenames'].to_list()]
      protsCplxs_uni = set().union(*[cplxsUni])
      protsCplxs_genename = set().union(*[cplxsGenename])
      cplxs = {'cplxsUni': cplxsUni, 'cplxsGenename': cplxsGenename}
      protsCplxs = {'protsCplxs_uni': protsCplxs_uni, 'protsCplxs_genename': protsCplxs_genename}

    return cplxs, protsCplxs

  def getBCB19(self):
    #generating pairs for the complete BCB19 and/or ELCFS feature matrix (not just training) could be prohibitive...
    bcb19Pairs = pd.read_csv(self.bcb19Dir_LM + 'featMat_bcb19Train.tsv', sep='\t',
                                  usecols=['idi', 'idii', 'label'], dtype={'idi': 'str', 'idii': 'str'})
    self.bcb19Pairs = operations.freezePairs(bcb19Pairs, 'idi', 'idii', pool=True)

    bcb19Pairs_feat = pd.read_csv(self.bcb19Dir_LM + 'featMat_bcb19.txt', sep=' ', header=None,
                                  names=['id1', 'id2', 'prob'], dtype={'id1': 'str', 'id2': 'str'})
    self.bcb19Pairs_feat = operations.freezePairs(bcb19Pairs_feat, 'id1', 'id2', pool=True)

    bcb19Cplxs_filenames = glob.glob(self.bcb19Dir_LM + '*%.txt')
    for bcb19Cplxs_file in bcb19Cplxs_filenames:
      k = bcb19Cplxs_file.split('/')[-1].split('.txt')[0].split('_')[-1]
      print(k)
      with open(bcb19Cplxs_file) as f:
        self.bcb19Cplxs[k] = [frozenset([str(ele) for ele in line.rstrip().split('\t')]) for line in f]
        self.bcb19Prots_cplxs[k] = set().union(*self.bcb19Cplxs[k])

  def mergePreds_huMAP(self):
    pairsHumap1p2 = \
      set(self.humap1Preds.pairsFrozen.to_list()).union(set(self.humap2Preds.pairsFrozen.to_list()))
    pairsDF_humap1p2 = pd.DataFrame({'pairsFrozen': list(pairsHumap1p2)})
    self.pairsHumap1p2 = pairsHumap1p2
    print('# total unique pairs in hu.MAP 1+2: {0}'.format(len(self.pairsHumap1p2)))

    self.protsHumap1p2 = set().union(*list(pairsHumap1p2))
    print('# total unique proteins in hu.MAP 1+2: {0}'.format(len(self.protsHumap1p2)))

    predsHumap1 = self.humap1Preds.groupby('pairsFrozen', as_index=False).mean(numeric_only=True)
    predsHumap2 = self.humap2Preds.groupby('pairsFrozen', as_index=False).mean(numeric_only=True)

    predsHumap1p2 = pairsDF_humap1p2.merge(predsHumap2, on=['pairsFrozen'], how='left')
    predsHumap1p2.loc[predsHumap1p2.prob.isna(), 'prob'] = \
      predsHumap1p2.loc[predsHumap1p2.prob.isna(), 'pairsFrozen'].map(
        dict(zip(predsHumap1.pairsFrozen, predsHumap1.prob)))
    predsHumap1p2.rename(columns={'prob': 'hu.MAP1+2'}, inplace=True)
    self.predsHumap1p2 = predsHumap1p2.copy()

  def getSTRING(self):
    stringPreds = \
      pd.read_csv(self.stringdbDir + '9606.protein.links.full.v12.0.txt.gz', sep=' ', engine='python')
    stringPreds['prob'] = stringPreds.combined_score.div(1000)
    stringIDs = set(stringPreds.protein1.to_list()).union(set(stringPreds.protein2.to_list()))
    print('STRING pair predictions dims: {0}'.format(stringPreds.shape))
    print('# string IDs: {0}'.format(len(stringIDs)))

    stringAliases = pd.read_csv(self.stringdbDir + '9606.protein.aliases.v12.0.txt.gz', sep='\t')
    stringAliases.dropna(subset=['alias'], inplace=True)
    stringIDs_aliases = set(stringAliases.loc[:, '#string_protein_id'].to_list())
    print('# string ID aliases: {0}'.format(len(stringIDs_aliases)))

    stringAliases_geneids = \
      set(stringAliases.alias.to_list()).intersection(self.geneSpace.geneIDs_proteinCoding)
    stringIDs_geneidPresent = \
      stringIDs.intersection(set(
        stringAliases.loc[stringAliases.alias.isin(stringAliases_geneids),
        '#string_protein_id'].to_list()))
    print('# stringID-mapped geneids: {0}'.format(len(stringAliases_geneids)))
    print('# geneID-mapped stringIDs: {0}'.format(len(stringIDs_geneidPresent)))

    stringPreds_geneidMappings = \
      stringPreds.loc[
        stringPreds.protein1.isin(stringIDs_geneidPresent) &
        stringPreds.protein2.isin(stringIDs_geneidPresent)].copy()
    print('GeneID-mapped STRING pair predictions dims: {0}'.format(stringPreds_geneidMappings.shape))

    stringAliases_aggSTRING = \
      stringAliases.loc[stringAliases.alias.isin(stringAliases_geneids)]
    stringAliases_aggSTRING = \
      stringAliases_aggSTRING.groupby(['#string_protein_id'], as_index=False).agg('; '.join)
    stringAliases_aggSTRING.loc[:, 'alias'] = \
      stringAliases_aggSTRING.loc[:, 'alias'].apply(lambda x: '; '.join(list(set(x.split('; ')))))
    stringID_geneidMapping = \
      dict(zip(stringAliases_aggSTRING['#string_protein_id'], stringAliases_aggSTRING['alias']))

    stringPreds_geneidMappings.insert(
      1, 'id1', stringPreds_geneidMappings.protein1.map(stringID_geneidMapping))
    stringPreds_geneidMappings.insert(
      3, 'id2', stringPreds_geneidMappings.protein2.map(stringID_geneidMapping))
    stringPreds_geneidMappings.loc[:, 'id1'] = \
      stringPreds_geneidMappings.loc[:, 'id1'].apply(lambda x: x.split('; ') if not x.isnumeric() else x)
    stringPreds_geneidMappings.loc[:, 'id2'] = \
      stringPreds_geneidMappings.loc[:, 'id2'].apply(lambda x: x.split('; ') if not x.isnumeric() else x)
    stringPreds_geneidMappings = stringPreds_geneidMappings.explode('id1', ignore_index=True)
    stringPreds_geneidMappings_id1Exp_dims = stringPreds_geneidMappings.shape
    stringPreds_geneidMappings = stringPreds_geneidMappings.explode('id2', ignore_index=True)
    stringPreds_geneidMappings_id2Exp_dims = stringPreds_geneidMappings.shape
    print('GeneID-mapped STRING protein1 expansion dims: {0}'.format(stringPreds_geneidMappings_id1Exp_dims))
    print('GeneID-mapped STRING protein2 expansion dims: {0}'.format(stringPreds_geneidMappings_id2Exp_dims))

    stringPreds_geneidMappings = \
      operations.freezePairs(stringPreds_geneidMappings, 'id1', 'id2', pool=True)
    stringPreds_geneidMappings = stringPreds_geneidMappings

    self.stringProts = \
      set(stringPreds_geneidMappings.id1.to_list()).union(set(stringPreds_geneidMappings.id2.to_list()))
    self.stringPreds = \
      stringPreds_geneidMappings.loc[:,['id1', 'id2', 'pairsFrozen', 'prob']].groupby(
        'pairsFrozen', as_index=False).mean(numeric_only=True)
    print('# geneID-mapped STRING proteins: {0}'.format(len(self.stringProts)))
    print('STRING geneID-mapped pair predictions aggregated dims: {0}'.format(len(self.stringPreds)))

    stringCplxs = pd.read_csv(self.stringdbDir + '9606.clusters.proteins.v12.0.txt.gz', sep='\t')
    stringCplxs = stringCplxs.drop(columns=['#string_taxon_id']).groupby(
      'cluster_id', as_index=False).agg('; '.join)
    stringCplxs.loc[:, 'protein_id'] = \
      stringCplxs.loc[:, 'protein_id'].apply(lambda x: set(x.split('; ')))
    stringCplxs.insert(
      2, 'totalMapping', stringCplxs.protein_id.apply(
        lambda x: True if len(x) == len(x.intersection(set(stringID_geneidMapping.keys()))) else False))
    stringCplxs = stringCplxs.loc[stringCplxs.totalMapping].copy()
    stringCplxs.loc[:, 'protein_id'] = stringCplxs.loc[:, 'protein_id'].apply(lambda x: list(x))
    stringCplxs = stringCplxs.explode('protein_id')
    stringCplxs.insert(2, 'complexes', stringCplxs.loc[:, 'protein_id'].map(stringID_geneidMapping))
    stringCplxs = \
      stringCplxs.drop(columns=['protein_id', 'totalMapping']).groupby(
        'cluster_id', as_index=False).agg('; '.join)
    stringCplxs.loc[:, 'complexes'] = \
      stringCplxs.loc[:, 'complexes'].apply(lambda x: ', '.join(list(set(x.split('; ')))))
    self.stringCplxs = [frozenset(cplx.split(', ')) for cplx in stringCplxs.complexes.to_list()]
    self.stringProts_cplxs = set().union(*self.stringCplxs)
    print('# completely, geneID-mapped STRING protein complexes: {0}'.format(len(self.stringCplxs)))

    with Pool(self.numWorkers) as pool:
      stringCplxs_pairs = \
        list(progressMonitor(pool.imap(operations.generatePairs, self.stringCplxs), total=len(self.stringCplxs)))
    predsPairs = set(self.stringPreds.pairsFrozen.to_list())
    stringPPIs = \
      [predsPairs.intersection(pairs) for pairs in stringCplxs_pairs]
    self.stringPPIs = set().union(*stringPPIs)
    print('# STRING geneID-mapped, predicted-pairs overlapping complex pairs: {0}'.format(self.stringPPIs))

  def go(self):
    #LM 2019/2020
    print('get LM 2019 information...')
    self.getBCB19()

    #hu.MAP 1.0
    print('get hu.MAP 1.0 information...')
    self.humap1Mat_train, self.humap1Prots_train, self.humap1Pairs_labeled = \
      self.setupTraining(self.humap1Dir, self.humap1Patt[0])
    self.humap1Mat_feats, self.humap1Prots_feats, self.humap1Pairs_feats = \
      self.getFeat(self.humap1Dir, self.humap1Patt[1], ['geneid1', 'geneid2'])
    self.humap1PPIs = self.getPPIs(self.humap1Mat_train, self.humap1Pairs_feats)
    self.humap1Preds = self.getPredictions(self.humap1Dir, self.humap1Patt[2])
    self.humap1Cplxs, self.humap1Prots_cplxs = self.getCplx(self.humap1Dir, self.humap1Patt[3], False)

    #hu.MAP 2.0
    print('get hu.MAP 2.0 information...')
    self.humap2Mat_train, self.humap2Prots_train, self.humap2Pairs_labeled = \
      self.setupTraining(self.humap2Dir, self.humap2Patt[0])
    self.humap2Mat_feats, self.humap2Prots_feats, self.humap2Pairs_feats = \
      self.getFeat(self.humap2Dir, self.humap2Patt[1])
    self.humap2PPIs = self.getPPIs(self.humap2Mat_train, self.humap2Pairs_feats)
    self.humap2Preds = self.getPredictions(self.humap2Dir, self.humap2Patt[2])
    self.humap2Cplxs, self.humap2Prots_cplxs = self.getCplx(self.humap2Dir, self.humap2Patt[3])

    #hu.MAP 1+2
    print('merge hu.MAP 1.0 with hu.MAP 2.0 release...')
    self.mergePreds_huMAP()
    #find prots, pairs, preds in the merger of both hu.MAP releases
    #protsHumap1p2, pairsHumap1p2, predsHumap1p2

    #STRING (1/22/24)
    print('get STRING information...')
    self.getSTRING()

    print('done')

  #!!!! check and subsequently count using grand feat mat or concat sets for each source
  def generateSummary(self):
    print('generating predictions\' summary table...')

    sources = ['STRING v.12', 'hu.MAP 1.0', 'hu.MAP 2.0', 'Lugo-Martinez 2019', 'ELCFS']
    #STRING v.12 (07/26/2023)

    prots = [len(self.stringProts), len(self.humap1Prots_feats), len(self.humap2Prots_feats),
             len(self.getProts(self.bcb19Pairs_feat, ['id1', 'id2'])), 19340]
    pairs = [math.comb(n, 2) for n in prots]

    pairsPreds = [len(self.stringPreds),
                  len(self.humap1Mat_feats), len(self.humap2Mat_feats),
                  len(self.bcb19Pairs_feat), 162758569]

    pairsLabeled = \
      ['-',
       math.comb(len(self.humap1Prots_train), 2),
       math.comb(len(self.humap2Prots_train), 2),
       math.comb(1679, 2), math.comb(3375, 2)]
    #BCB19 (LM2019 results) labeled pairs possible is 1,679 only by design;
    # original total from hu.MAP 1.0 CORUM (release 2012 release)-derived labeled data

    pairsLabeled_preds = \
      ['-',
       len(self.humap1Mat_train), len(self.humap2Mat_train),
       1038169, math.comb(3215, 2)]
    #total CORUM (2018 release) labeled pairs available: 5,693,625

    cplxs = [len(self.stringCplxs),
             len(self.humap1Cplxs), 6965, '1,340 (top 0.005%); 7,758 (top 0.02%)', 'Result pending changes']

    cplxsPPIs = [len(self.stringPPIs), len(self.humap1PPIs), len(self.humap2PPIs),
                 '570 (top 0.005%); 1,135 (top 0.02%)', 'Result pending changes']

    cplxProts = [len(self.stringProts_cplxs),
                 len(self.humap1Prots_cplxs), 6957, '412 (top 0.005%); 10,308 (top 0.02%)',
                 'Result pending changes']

    self.summaryDesc = \
      pd.DataFrame({'Source': sources,
                    'Proteins': prots, 'Pairs': pairs, 'Predictions': pairsPreds,
                    'Pairs labeled': pairsLabeled, 'Pair-labeled predictions': pairsLabeled_preds,
                    'Complexes': cplxs, 'Complexes\' PPIs': cplxsPPIs, 'Complexes\' proteins': cplxProts})
    self.summaryDesc.set_index('Source', inplace=True)
    self.summaryDesc = self.summaryDesc.applymap(operations.fmtTbl)

    return self.summaryDesc

  #get rid of later...
  def makeTbl1(self):
    tbl1Dir = './drive/MyDrive/elcfs_protein_complex_modeling/google_colab_notebooks/notebookQC/'
    tbl1 = pd.read_csv(tbl1Dir + 'tbl1.tsv', sep='\t')
    return tbl1.set_index('Unnamed: 0', drop=True).rename_axis(index=None)


class hpaFeats(object):

  def __init__(self, rootDir,
               print=True, save=False, saveDir=False, overwrite=False, numWorkers=False, divy=False, divies=False):
    self.rootDir = rootDir
    self.dataDir = rootDir + 'srcData/HPA/'
    self.names = ['hpa']
    self.filenames = ['subcellular_location.tsv.zip']
    self.paths = {n: self.dataDir + f for n, f in zip(self.names, self.filenames)}
    self.datasets = None
    self.featMat = None
    self.print = print

    self.save = save
    self.saveDir = saveDir
    self.overwrite = overwrite
    self.divy = divy

    if divies:
      self.divies = divies
    elif divy:
      self.divies = 3

    if not self.saveDir:
      targetDir = self.rootDir + 'featureData/proteinLevel/'
    else:
      targetDir = self.saveDir
    self.filename = \
      targetDir + 'hpaFeats.tsv.tar.bz2' if not divy else targetDir + 'hpaFeats_{0}.tsv.tar.bz2'.format(divy)

    if not numWorkers:
      self.numWorkers = cpu_count() - 1
    else:
      self.numWorkers = numWorkers

  def report(self, input=None):
    print('Printing resulting feature matrix dimensions...')
    if (input is None) & (self.print):
      for k, df in self.datasets.items(): print((k, df.shape))
    elif (self.print):
      print(input.shape)
    print('\n')

  def load(self):
    print('loading original data...')
    self.datasets = \
      {k: pd.read_csv(v, sep='\t') for k, v in self.paths.items()}
    self.report()

  def trim(self):
    print('Dropping extraneous columns...')
    for k in self.datasets.keys():
      self.datasets[k].drop(
        columns=['Additional location', 'Cell cycle dependency', 'Enhanced',
                 'Extracellular location', 'GO id', 'Gene name',
                 'Main location', 'Reliability', 'Single-cell variation intensity',
                 'Single-cell variation spatial', 'Uncertain'], inplace=True)
      self.datasets[k].dropna(subset=['Supported', 'Approved'], how='all', inplace=True)
    self.report()

  def multiRec_comb(self):
    print('Combining multiple same-ID observations...')
    for k, v in self.datasets.items():
      v = v.set_index('Gene').stack().groupby(
        level=0, sort=False).agg(';'.join).to_frame().rename(
        columns={0: 'locations'}).reset_index()
      v.loc[:, 'locations'] = v['locations'].str.lower()
      self.datasets[k] = v
    self.report()

  def prep(self):
    self.trim()
    self.multiRec_comb()

  def saveMat(self):
    if self.save:
      print('{0} saving...'.format(self.filename))
      self.featMat.drop(columns=['pairsFrozen']).to_csv(self.filename, sep='\t', index=False)

  def calculate(self, pair):
    id1, id2 = list(pair)
    annots1 = \
      set(self.datasets['hpa'].loc[self.datasets['hpa'].Gene == str(id1), 'locations'].item().split(';'))
    annots2 = \
      set(self.datasets['hpa'].loc[self.datasets['hpa'].Gene == str(id2), 'locations'].item().split(';'))

    overlap = len(annots1.intersection(annots2)) / max([len(x) for x in [annots1, annots2]])
    set_equality = 1 if annots1 == annots2 else 0
    jaccard_similarity = len(annots1.intersection(annots2)) / len(annots1.union(annots2))
    set_inclusion = 1 if any([annots1.issubset(annots2), annots2.issubset(annots1)]) else 0

    return overlap, set_equality, jaccard_similarity, set_inclusion

  def generate(self):

    if glob.glob(self.filename):
      print('loading previously generated features...')
      self.featMat = pd.read_csv(self.filename, sep='\t', dtype={'id1': 'str', 'id2': 'str'})

    else:
      self.load()
      self.prep()
      for k, v in self.datasets.items():
        print('generating list of pairs...')
        pairsList = list(operations.generatePairs(v.Gene.unique()))
        if (self.divy) | (self.divy == 0):
          print('divying up pairs...')
          pairsList = \
            [list(pairsList_sub)
             for i, pairsList_sub in enumerate(list(np.array_split(pairsList, self.divies))) if i==self.divy][0]
        print('# pairs: {0}'.format(len(pairsList)))

        print('generating features...')
        pool = Pool(self.numWorkers)
        hpaFeats = [list(progressMonitor(pool.imap(self.calculate, pairsList), total=len(pairsList)))]
        pool.close()
      self.featMat = hpaFeats

      hpaFeat_cols = ['overlap', 'set_equality', 'jaccard_similarity', 'set_inclusion']
      hpaFeats = pd.concat([pd.DataFrame({'pairsFrozen': pairsList}),
                            pd.DataFrame(hpaFeats[0], columns=hpaFeat_cols)],
                           axis=1)
      self.featMat = hpaFeats


class nci60Feats(object):

  def __init__(self, rootDir, multiRec='mean', print=True, save=False, saveDir=False, overwrite=False):
    self.rootDir = rootDir
    self.dataDir = rootDir + 'srcData/NCI60/'
    self.names = ['rna exp', 'cell-line protein exp', 'tissue protein exp']
    self.filenames = \
      ['nci60_all_expressionPerCellLine.tsv',
       'nci60_all_abundanceLFQPerCellLineProteome.tsv',
       'nci60_all_abundanceLFQPerTissueDeep.tsv']
    self.paths = {n: self.dataDir + f for n, f in zip(self.names, self.filenames)}
    self.datasets = dict()
    self.featMat = None
    self.print = print

    #LM 2019 chose 1st
    self.multiRec = multiRec
    self.save = save
    self.saveDir = saveDir
    self.overwrite = overwrite
    if not self.saveDir:
      targetDir = self.rootDir + 'featureData/proteinLevel/'
    else:
      targetDir = self.saveDir
    self.filename = targetDir + 'nci60Feats.tsv.tar.bz2'

  def report(self, input=None):
    print('Printing resulting feature matrix dimensions...')
    if (input is None) & (self.print):
      for k, df in self.dataMats.items(): print((k, df.shape))
    elif (self.print):
      print(input.shape)
    print('\n')

  def load(self):
    print('loading original data...')
    self.dataMats = \
      {k: pd.read_csv(v, sep='\t', dtype={'GeneID': 'str'})
      for k, v in self.paths.items()}
    self.report()

  def trim(self):
    print('Dropping extraneous non-GeneID desc columns...')
    print('Dropping GeneID absent records...')
    for k in self.dataMats.keys():
      self.dataMats[k].drop(columns=['GeneName', 'UniprotID'], inplace=True)
      self.dataMats[k].dropna(subset=['GeneID'], inplace=True)
    self.report()

  def sepDelimited(self):
    print('Breaking multiple-ID records into individual observations...')
    for k, df in self.dataMats.items():
      df.loc[:, 'GeneID'] = df.GeneID.apply(str.split, args=(','))
      self.dataMats[k] = df.explode('GeneID', ignore_index=True)
    self.report()

  #options for prots with multiple observations: select first, last, or calculate mean
  def multiRec_comb(self):
    print('Combining multiple same-ID observations...')
    if self.multiRec == 'mean':
      self.dataMats = \
        {k: df.groupby('GeneID', as_index=False).mean(numeric_only=True)
         for k, df in self.dataMats.items()}
    elif self.multiRec == 'first':
      self.dataMats = \
        {k: df.groupby('GeneID', as_index=False).first(numeric_only=True)
         for k, df in self.dataMats.items()}
    elif self.multiRec == 'last':
      self.dataMats = \
        {k: df.groupby('GeneID', as_index=False).last(numeric_only=True)
         for k, df in self.dataMats.items()}
    self.report()

  def prep(self):
    self.trim()
    self.sepDelimited()
    self.multiRec_comb()

  def saveMat(self):
    if self.save:
      print('{0} saving...'.format(self.filename))
      self.featMat.drop(columns=['pairsFrozen']).to_csv(self.filename, sep='\t', index=False)

  def generate(self):

    if glob.glob(self.filename):
      print('loading previously generated features...')
      self.featMat = pd.read_csv(self.filename, sep='\t', dtype={'id1': 'str', 'id2': 'str'})

    else:
      self.load()
      self.prep()
      for k, df in self.dataMats.items():
        print(k)
        df.set_index('GeneID', inplace=True, verify_integrity=True)
        corrMat = df.T.corr()
        print('correlation feature matrix dimensions: {0}'.format(corrMat.shape))
        corrMat = corrMat.where(np.triu(np.ones(corrMat.shape), k=1).astype(bool))
        print('correlation feature matrix sample, restructuring: {0}'.format(corrMat.head()))
        corrMat = corrMat.rename_axis('id1').rename_axis('id2', axis=1).stack().reset_index()
        print('correlation feature matrix sample, restructuring: {0}'.format(corrMat.head()))
        corrMat.rename(columns={0: k}, inplace=True)
        print('correlation feature matrix sample, restructuring: {0}'.format(corrMat.head()))
        corrMat.dropna(subset=[k], inplace=True)
        print('correlation feature matrix final dimensions: {0}'.format(corrMat.shape))
        self.datasets[k] = operations.freezePairs(corrMat, 'id1', 'id2', pool=True)
        self.datasets[k].drop(columns=['id1', 'id2'], inplace=True)

      print('merging datasets...')
      self.featMat = \
          reduce(lambda left, right:
                 pd.merge(left, right, on=['pairsFrozen'], how='outer'),
                 [df for df in self.datasets.values()])

      print('adding identifier columns...')
      self.featMat = \
          pd.concat(
            [pd.DataFrame(self.featMat.pairsFrozen.to_list(), columns=['id1', 'id2'], dtype='str'),
             self.featMat],
          axis=1)

      print('averaging protein expression features...')
      self.featMat.insert(
        6, 'protein exp', self.featMat.loc[:, ['cell-line protein exp', 'tissue protein exp']].mean(axis=1))
      self.report(input=self.featMat)


class scbcFeats(object):

  def __init__(self, rootDir, multiRec='mean', print=True, save=False, saveDir=False, overwrite=False):
    self.rootDir = rootDir
    self.dataDir = rootDir + 'srcData/SubCellBarCode/'
    self.path = self.dataDir + 'mmc2.xlsx'
    self.cellLines = ['A431', 'MCF7', 'H322', 'HCC827', 'U251']
    self.dataMats = dict()
    self.multiRec = multiRec
    self.featMat = None
    self.print = print
    self.save = save
    self.saveDir = saveDir
    self.overwrite = overwrite
    if not self.saveDir:
      targetDir = self.rootDir + 'featureData/proteinLevel/'
    else:
      targetDir = self.saveDir
    self.filename = targetDir + 'scbcFeats.tsv.tar.bz2'

  def report(self, input=None):
    print('Printing resulting feature matrix dimensions...')
    if (input is None) & (self.print):
      for k, df in self.dataMats.items(): print((k, df.shape))
    elif (self.print):
      print(input.shape)
    print('\n')

  def load(self):
    print('loading original data...')
    self.dataMats = \
      {k: pd.read_excel(self.path, sheet_name=k, dtype={'Protein': 'str'})
      for k in self.cellLines}
    self.report()

  def trim(self):
    print('Dropping extraneous non-fraction columns...')
    for k in self.dataMats.keys():
      self.dataMats[k].drop(columns=['MIN PSMs for quant'], inplace=True)
    self.report()

  #options for prots with multiple observations: select first, last, or calculate mean
  def multiRec_comb(self):
    print('Combining multiple same-ID observations...')
    if self.multiRec == 'mean':
      self.dataMats = \
        {k: df.groupby('Protein', as_index=False).mean(numeric_only=True)
         for k, df in self.dataMats.items()}
    elif self.multiRec == 'first':
      self.dataMats = \
        {k: df.groupby('Protein', as_index=False).first(numeric_only=True)
         for k, df in self.dataMats.items()}
    elif self.multiRec == 'last':
      self.dataMats = \
        {k: df.groupby('Protein', as_index=False).last(numeric_only=True)
         for k, df in self.dataMats.items()}
    self.report()

  def prep(self):
    self.trim()
    self.multiRec_comb()

  def saveMat(self):
    if self.save:
      self.featMat.drop(columns=['pairsFrozen']).to_csv(self.filename, sep='\t', index=False)

  def generate(self):

    if glob.glob(self.filename):
      print('loading previously generated features...')
      self.featMat = pd.read_csv(self.filename, sep='\t', dtype={'id1': 'str', 'id2': 'str'})

    else:
      self.load()
      self.prep()
      for k, df in self.dataMats.items():
        print(k)
        df.set_index('Protein', inplace=True, verify_integrity=True)
        corrMat = df.T.corr()
        print('correlation feature matrix dimensions: {0}'.format(corrMat.shape))
        corrMat = corrMat.where(np.triu(np.ones(corrMat.shape), k=1).astype(bool))
        print('correlation feature matrix sample, restructuring: {0}'.format(corrMat.head()))
        corrMat = corrMat.rename_axis('id1').rename_axis('id2', axis=1).stack().reset_index()
        print('correlation feature matrix sample, restructuring: {0}'.format(corrMat.head()))
        corrMat.rename(columns={0: k}, inplace=True)
        print('correlation feature matrix sample, restructuring: {0}'.format(corrMat.head()))
        corrMat.dropna(subset=[k], inplace=True)
        print('correlation feature matrix final dimensions: {0}'.format(corrMat.shape))
        self.datasets[k] = operations.freezePairs(corrMat, 'id1', 'id2', pool=True)
        self.datasets[k].drop(columns=['id1', 'id2'], inplace=True)

      print('merging datasets...')
      self.featMat = \
          reduce(lambda left, right:
                 pd.merge(left, right, on=['pairsFrozen'], how='outer'),
                 [df for df in self.datasets.values()])

      print('adding identifier columns...')
      self.featMat = \
          pd.concat(
            [pd.DataFrame(self.featMat.pairsFrozen.to_list(), columns=['id1', 'id2'], dtype='str'),
             self.featMat],
          axis=1)

      self.report(input=self.featMat)