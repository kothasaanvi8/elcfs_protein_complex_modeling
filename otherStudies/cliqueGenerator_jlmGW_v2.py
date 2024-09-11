
# Authors: Jose Lugo-Martinez, Gary Robert Wilkins Jr
# Dates: 01/10/2019
# Advisors: Robert F. Murphy, Ziv Bar-Joseph
# Description: This function computes and reports all the maximal cliques in an undirected graph
# Last Modified: 01/10/2023

import glob
from itertools import combinations as comb
import multiprocessing as mp
import networkx as nx
import pandas as pd
import sys
from tqdm import tqdm as progressMonitor


def main(args):
    print(args)
    if len(args) == 6:
        pairsFile = args[1]
        predictionsFile = args[2]
        thresholdPrime = float(args[3])
        thresholdSecondary = float(args[4])
        outputPath = args[5]
    else:
        print('<<ERROR>>  Invalid number of parameters!')
        return

    buildComplexes(pairsFile, predictionsFile, thresholdPrime, thresholdSecondary, outputPath)


class buildComplexes(object):

    def __init__(self, pairsInfo_filename, predictionsFilename,
                 thresholdPrime, thresholdSecondary, outputDir,
                 pairsRequired=1):

        self.predictionsFilename_root = predictionsFilename.split('.txt')[0].split('/')[-1]
        self.thresholdPrime = thresholdPrime
        self.thresholdSecondary = thresholdSecondary
        self.outputDir = outputDir
        self.pairsRequired = pairsRequired
        self.workers = mp.cpu_count()-1

        # Open input files
        self.pairsInfo = pd.read_csv(pairsInfo_filename, sep='\t', usecols=[1, 2], names=['id1', 'id2'], dtype='str')
        self.proteins = set(self.pairsInfo.loc[:, 'id1'].to_list() + self.pairsInfo.loc[:, 'id2'].to_list())
        self.predictedScores = pd.read_csv(predictionsFilename, names=['score'])

        self.pairsMatrix_root = pd.concat([self.pairsInfo, self.predictedScores], axis=1)
        self.pairsMatrix_root.insert(3, 'pairsFrozen',
                                     self.pairsMatrix_root.loc[:, ['id1', 'id2']].apply(frozenset, axis=1))
        self.pairsRoot = self.pairsMatrix_root.pairsFrozen.to_list()
        self.pairsMatrix = self.pairsMatrix_root.loc[self.pairsMatrix_root.score >= thresholdPrime, :].copy()
        self.pairsMatrix.reset_index(drop=True, inplace=True)

        self.proteinsGraph = None
        self.predictedComplexes = None
        self.predictedComplexes_refinedCounted = None
        self.predictedComplexes_refinedCounted_stats = None
        self.predictedComplexes_refinedCounted_noSubsets_indices = None
        self.predictedComplexes_refinedOrganized = None

    def refineCliques(self, cplxInit):

        # ...expensive line of code
        pairsFrozen_cplxSpecific = list(
            set([frozenset(pair) for pair in list(comb(cplxInit, 2))]).intersection(
                set(self.pairsMatrix.pairsFrozen.to_list())))
        # ...
        scores = \
            self.pairsMatrix.loc[
                self.pairsMatrix.pairsFrozen.isin(pairsFrozen_cplxSpecific), 'score'].to_list()

        remainingPairs = list(set(self.pairsRoot).difference(set(pairsFrozen_cplxSpecific)))
        remainingPairs_thresholdSecondary = list(set(
            self.pairsMatrix_root.loc[(
                                          (self.pairsMatrix_root.pairsFrozen.isin(remainingPairs) &
                                           (self.pairsMatrix_root.score >= self.thresholdSecondary))),
                                      'pairsFrozen'].to_list()))
        remainingPairs_disjoint = \
            [pair for pair in remainingPairs_thresholdSecondary if len(pair.intersection(set(cplxInit))) == 1]
        disjointNodes = set().union(*remainingPairs_disjoint) - set(cplxInit)

        retainedNodes = []
        retainedNodes_pairs = []
        for node in list(disjointNodes):
            nodeSpecific_pairs = [pair for pair in remainingPairs_disjoint if node in list(pair)]
            if len(nodeSpecific_pairs) >= self.pairsRequired:
                retainedNodes.append(node)
                retainedNodes_pairs.append(nodeSpecific_pairs)

        # couple of directions we could take:
        # (1) add every protein of the remaining pairs to the base complex if there exists a pair between that
        # ...protein and any comprising the base complex
        # (2) "   " if there exists some minimum number of pairs between that protein and those comprising
        # ...the base complex
        # (3) add the proteins corresponding to those with the highest number of pairs and above 1, or to those
        # corresponding to the top fraction of pairs

        cplx = cplxInit.copy()
        cplx.extend(list(retainedNodes))
        cplx.sort()
        cplxID = ','.join(cplx)

        scores.extend(
            self.pairsMatrix_root.loc[
                self.pairsMatrix_root.pairsFrozen.isin(retainedNodes_pairs), 'score'].to_list())
        complexScore = min(scores)
        complexScore_mean = sum(scores) / len(scores)

        return frozenset(cplx), cplxID, complexScore, complexScore_mean

    def getCliques(self):
        self.predictedComplexes = {}
        self.predictedComplexes_refinedCounted = []
        self.predictedComplexes_refinedCounted_stats = []
        self.predictedComplexes_refinedCounted_noSubsets_indices = []
        self.predictedComplexes_refinedOrganized = {}
        self.proteinsGraph = nx.from_pandas_edgelist(self.pairsMatrix, 'id1', 'id2', 'score').to_undirected()
        print("# of nodes", len(self.proteinsGraph.nodes))
        print("# of edges", len(self.proteinsGraph.edges))

        cliques = nx.find_cliques(self.proteinsGraph)
        idx = 0
        for cliq in cliques:
            self.predictedComplexes[idx] = cliq
            idx += 1
        print(len(self.predictedComplexes))

        print('beginning refinement of predicted complexes...')
        pool = mp.Pool(self.workers)
        predictedComplexes_refined = \
            [list(progressMonitor(
                pool.imap(self.refineCliques, self.predictedComplexes.values()),
                total=len(self.predictedComplexes.values())))]
        pool.close()
        predictedComplexes_refinedList = \
            [pd.DataFrame([ele for ele in row], columns=['cplxFrozen', 'cplxID', 'complexScore', 'complexScore_mean'])
             for row in predictedComplexes_refined]
        predictedComplexes_refinedDF = \
            pd.concat(predictedComplexes_refinedList).drop_duplicates(subset='cplxFrozen')
        print('total refined predicted complexes...{0}'.format(predictedComplexes_refinedDF.shape[0]))
        print('completed refinement of predicted complexes...')

        self.predictedComplexes_refinedCounted_noSubsets_indices = \
            [cplxIdx for cplxIdx in predictedComplexes_refinedDF.index if
             not any([predictedComplexes_refinedDF.loc[cplxIdx, 'cplxFrozen'].issubset(cplxRemaining)
                      for cplxRemaining in list(set(predictedComplexes_refinedDF.cplxFrozen.to_list()) -
                                                {predictedComplexes_refinedDF.loc[cplxIdx, 'cplxFrozen']})])]
        print('completed identification of non-subset complex predictions...')

        for cplxIdx in self.predictedComplexes_refinedCounted_noSubsets_indices:
            cplxLen = len(predictedComplexes_refinedDF.loc[cplxIdx, 'cplxFrozen'])
            if cplxLen in self.predictedComplexes_refinedOrganized:
                self.predictedComplexes_refinedOrganized[cplxLen].append(
                    predictedComplexes_refinedDF.loc[cplxIdx,
                                                     ['cplxID', 'complexScore', 'complexScore_mean']].to_list())
            else:
                self.predictedComplexes_refinedOrganized[cplxLen] = \
                    [predictedComplexes_refinedDF.loc[cplxIdx,
                                                      ['cplxID', 'complexScore', 'complexScore_mean']].to_list()]

    def saveCliques(self):

        for complexSize, complexesSet in self.predictedComplexes_refinedOrganized.items():
            outfilename = self.outputDir + self.predictionsFilename_root + str(complexSize) + '.txt'
            outfile = open(outfilename, "w")
            for cplxID, complexScore, complexScore_mean in complexesSet:
                outline = cplxID + '\t' + str(complexScore) + '\t' + str(complexScore_mean) + '\n'
                outfile.writelines(outline)

            # Close output file
            outfile.close()

    def consolidateCliques(self):

        assemblyDir = glob.glob(self.outputDir + self.predictionsFilename_root + '*.txt')
        assemblyDir.sort()

        predictions = \
            pd.concat(
                [pd.read_csv(predictionsFile, sep='\t',
                             names=['complex', 'score', 'score_pairsAvg-derived'],
                             dtype={'complex': 'str', 'score': 'float64',
                                    'score_pairsAvg-derived': 'float64'})
                 for predictionsFile in assemblyDir])
        predictions.insert(1, 'complexFrozen',
                           [frozenset([str(prot) for prot in cplx.split(',')])
                            for cplx in predictions.loc[:, 'complex'].to_list()])

        predictions.to_csv(self.outputDir + self.predictionsFilename_root + 'consolidatedComplexes_detailed.txt',
                           sep='\t', index=None)

        outfilename = self.outputDir + self.predictionsFilename_root + 'consolidatedComplexes.txt'
        outfile = open(outfilename, "w")
        for cplx in predictions.loc[:, 'complex'].to_list():
            line = ' '.join(cplx.split(',')) + '\n'
            outfile.writelines(line)

        # Close output file
        outfile.close()


if __name__ == '__main__':
    main(sys.argv)
