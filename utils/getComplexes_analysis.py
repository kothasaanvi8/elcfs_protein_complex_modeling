# Useful miscellaneous functions for analysing complexes

import glob, itertools
import math
import matplotlib.pyplot as plt
import numpy as np
import os, re
import pandas as pd
import pickle
import seaborn as sns
import statistics

from util.getComplexes import buildComplexes

def generateInputs_complexBuilder_testCase(proteinSets, outputPrefix,
                                           fakeConfidence_value=0.9,
                                           complexBuilder_threshold=0.9,
                                           runBuilder=False,
                                           *outputDir):

    allPairs_theoreticalCombined = np.empty((0, 2), dtype=str)
    for proteinSet in proteinSets:
        allPairs_theoreticalCombined = \
            np.vstack(
                (allPairs_theoreticalCombined,
                 np.vstack(list(itertools.combinations(proteinSet, 2))))
            )
        print("The new shape of iteratively generated pairs is: " +
              str(np.vstack(list(itertools.combinations(proteinSet, 2))).shape))

    allPairs_theoreticalCombined = \
        np.insert(
            allPairs_theoreticalCombined, 0,
            np.arange(len(allPairs_theoreticalCombined)), axis=1)

    allPairs_theoreticalCombined_confidenceVals = \
        np.empty((len(allPairs_theoreticalCombined), 1))
    allPairs_theoreticalCombined_confidenceVals[:] = str(fakeConfidence_value)

    os.mkdir('./' + outputPrefix)
    np.savetxt('./' + outputPrefix + '/' + outputPrefix + '_pairsInfo.tsv',
               allPairs_theoreticalCombined, delimiter='\t', fmt='%s')
    np.savetxt('./' + outputPrefix + '/' + outputPrefix + '_pairsPredictions.txt',
               allPairs_theoreticalCombined_confidenceVals, delimiter=',',
               fmt='%.4f')

    # build complexes
    if runBuilder:
        buildComplexes(outputPrefix, complexBuilder_threshold, outputDir)

def generateInputs_complexBuilder_testCase_realData(pairsDF, proteinSets, idNames,
                                                    confidenceCol_name, outputPrefix,
                                                    complexBuilder_threshold=0.9,
                                                    runBuilder=False,
                                                    *outputDir):
    allProteins = set()
    for proteinSet in proteinSets:
        allProteins = allProteins.union(proteinSet)

    idNames_plusConf = list(set(idNames).union({confidenceCol_name}))
    pairsAvailable = pairsDF.loc[
        (
                (pairsDF[idNames[0]].isin(allProteins)) &
                (pairsDF[idNames[1]].isin(allProteins))
        ), idNames_plusConf]
    pairsAvailable[idNames[0]] = pairsAvailable[idNames[0]].astype('str')
    pairsAvailable[idNames[1]] = pairsAvailable[idNames[1]].astype('str')

    pairIDs = pairsAvailable.loc[:, idNames].to_numpy()
    pairIDs = np.insert(pairIDs, 0, np.arange(len(pairIDs)), axis=1)
    pairConfidences = pairsAvailable.loc[:, confidenceCol_name]

    os.mkdir('./' + outputPrefix)
    np.savetxt('./' + outputPrefix + '/' + outputPrefix + '_pairsInfo.tsv',
               pairIDs, delimiter='\t', fmt='%s')
    np.savetxt('./' + outputPrefix + '/' + outputPrefix + '_pairsPredictions.txt',
               pairConfidences, delimiter=',', fmt='%.4f')

    # build complexes
    if runBuilder:
        buildComplexes(outputPrefix, complexBuilder_threshold, outputDir)

def generateInputs_complexBuilder(pairsDF, idNames, confidenceCol_name, outputPrefix,
                                  complexBuilder_threshold=0.9, runBuilder=False, *outputDir):

    pairsDF[idNames[0]] = pairsDF[idNames[0]].astype('str')
    pairsDF[idNames[1]] = pairsDF[idNames[1]].astype('str')
    pairIDs = pairsDF.loc[:, idNames].to_numpy()
    pairIDs = np.insert(pairIDs, 0, np.arange(len(pairIDs)), axis=1)
    pairConfidences = pairsDF.loc[:, confidenceCol_name]

    os.mkdir('./' + outputPrefix)
    np.savetxt('./' + outputPrefix + '/' + outputPrefix + '_pairsInfo.tsv',
               pairIDs, delimiter='\t', fmt='%s')
    np.savetxt('./' + outputPrefix + '/' + outputPrefix + '_pairsPredictions.txt',
               pairConfidences, delimiter=',', fmt='%.4f')

    # build complexes
    if runBuilder:
        buildComplexes(outputPrefix, complexBuilder_threshold, outputDir)

def acquirePredictions(dataDir_patt, modelResults_dirPatt, outputPrefix, suffix='',
                       complexBuilder_threshold=0.9, runBuilder=False, *outputDir):
    numFolds = np.arange(5)
    allLabels = {fold: np.loadtxt(predictionsDir, delimiter='\t', skiprows=1, dtype=np.int) for fold, predictionsDir in
                 zip(numFolds, glob.glob(dataDir_patt + 'testLabels.tsv'))}
    probsPos_combinedWeighted_byFolds = \
        {fold: np.loadtxt(predictionsDir, delimiter=',') for fold, predictionsDir in
         zip(numFolds, glob.glob(modelResults_dirPatt + 'probsPos_weighted' + suffix + '.csv'))}
    labels_acrossFolds = []
    combinedWeighted_probsPos_kCV = []
    for combinedWeighted_probsPos, labels_byFolds in zip(probsPos_combinedWeighted_byFolds.values(), allLabels.values()):
        labels_acrossFolds.append(labels_byFolds)
        combinedWeighted_probsPos_kCV.append(combinedWeighted_probsPos)

    labels_acrossFolds = np.concatenate(labels_acrossFolds)
    combinedWeighted_probsPos_kCV = np.concatenate(combinedWeighted_probsPos_kCV)
    print('The lengths of the labels\'s and predictions\'s datasets are ' +
          str(len(labels_acrossFolds)) + ' and ' +
          str(len(combinedWeighted_probsPos_kCV)) + ' respectively.')

    os.mkdir('./' + outputPrefix)
    np.savetxt('./' + outputPrefix + '/' + outputPrefix + '_pairsInfo.tsv',
               labels_acrossFolds, delimiter='\t', fmt='%d')
    np.savetxt('./' + outputPrefix + '/' + outputPrefix + '_pairsPredictions.txt',
               combinedWeighted_probsPos_kCV, delimiter=',', fmt='%.4f')

    # build complexes
    if runBuilder:
        buildComplexes(outputPrefix, complexBuilder_threshold, outputDir)

    return labels_acrossFolds, combinedWeighted_probsPos_kCV

def transformComplexes_variableUse(complexesVar):
    if isinstance(complexesVar, dict):
        listComplexes = list(complexesVar.values())
    elif isinstance(complexesVar, set):
        listComplexes = list(complexesVar)
    elif isinstance(complexesVar, list):
        listComplexes = complexesVar
    else:
        print('input complexes type not accepted.')
        return

    return listComplexes

def condenseComplexes_singleProteins(complexes, verbose=False):
    proteinSet = set()
    listComplexes = transformComplexes_variableUse(complexes)
    for complexSet in listComplexes:
        proteinSet = proteinSet.union(complexSet)
    if verbose:
        print('The protein set contains ' + str(len(proteinSet)) + ' members.')

    return proteinSet

def confineComplexes_proteinSet(complexes, proteinSet, verbose=False):
    listComplexes = transformComplexes_variableUse(complexes)
    confinedComplexes = set(
        [frozenset(proteinSet.intersection(ele)) for ele in listComplexes if proteinSet.intersection(ele)]
    )
    if verbose:
        print('The confined complexes set contains ' +
              str(len(condenseComplexes_singleProteins(confinedComplexes))) + ' members.')

    return confinedComplexes

def getComplex_predictions(outputDir, outputPrefix, realDir=True):

  if not realDir:
    predictionsPath = \
      '/content/' + outputPrefix + '/' + outputPrefix + '_pairsPredictions*.txt'
  else:
    predictionsPath = \
      outputDir + outputPrefix + '/' + outputPrefix + '_pairsPredictions*.txt'

  predictedCliques_sets = dict()
  for filename in glob.glob(predictionsPath):
    filenameExtension = filename.split('.txt')[0].split('_pairsPredictions')[1]
    if filenameExtension.isnumeric():
      if int(filenameExtension) != 2:
        clique = dict()
        data = pd.read_csv(filename, sep='\t',
                           header=None, names=['complex', 'score'])
        for idx in data.index:
          line = [ele for ele in re.split(',', data['complex'].loc[idx]) if ele]
          clique[idx] = set(list(map(str, line)))

        predictedCliques_sets[int(filenameExtension)] = \
          [frozenset(ele) for ele in list(clique.values())]

  return predictedCliques_sets

def assembleComplex_predictionsTogether(predictedComplexes_filesDirectory):
    predictedComplexes_filenames = \
        sorted(glob.glob(predictedComplexes_filesDirectory+'complexPredictions_partitionModel_thresh=0.*csv'))
    partitionModels = \
        {
            str(round(float(predictedComplexes_filename.split('=')[1].split('.csv')[0]), 4)):
                pd.read_csv(predictedComplexes_filename, header=None) for predictedComplexes_filename in
            predictedComplexes_filenames
        }

    partitionModels = \
        {
            modelKey:
                complexesPredicted.apply(
                    lambda line:
                    [re.sub(r'\s', '', prot).replace('\'', '') for prot in line[0].split(',')], axis=1) for
            (modelKey, complexesPredicted) in partitionModels.items()
        }

    partitionModels = \
        {
            modelKey: set([frozenset(lineDelimited) for lineDelimited in complexesPredicted.tolist()]) for
            (modelKey, complexesPredicted) in partitionModels.items()
        }

    pickle.dump(partitionModels, open(predictedComplexes_filesDirectory+'assembledPredictions.pkl', 'wb'))

    return partitionModels

def getComplexes_setPlus_proteins(complexesDF, colName=0):
    complexesCounter = 0

    proteinsSet = set()
    complexesSet = dict()
    complexesLengths = dict()
    for idx in complexesDF.index:
        placeholder = \
          frozenset([str(int(ele)) for ele in re.split(',\(|\(|\),|\)|\t|\s|;|,',
                                             complexesDF.loc[idx, colName]) if ele])
        if len(placeholder) != 0:
          complexesSet[idx] = placeholder
          complexesLengths[idx] = len(complexesSet[idx])
          proteinsSet = proteinsSet.union(complexesSet[idx])
          complexesCounter+=1

    print('There is a total of ' + str(len(complexesSet)) +
          ' predicted complexes and a total of ' + str(len(proteinsSet)) + ' proteins.')
    print('The minimum and maximum lengths of predicted complexes were ' +
          str(np.amin(list(complexesLengths.values()))) + ' and ' +
          str(np.amax(list(complexesLengths.values()))) + ' respectively.')

    return proteinsSet, complexesSet

def getRef_baselines(proteinSet, complexes, returnResults=True, verbose=False):

    if not isinstance(proteinSet, set):
        print('Input proteins confining complexes elements as a set')

    if not all(isinstance(ele, str) for ele in proteinSet):
        print('Limiting labels passed include non-string entries. For comparison with CORUM, ensure all limiting '
              'are string type.')
        return

    complexes = transformComplexes_variableUse(complexes)
    complexesLengths = [len(ele) for ele in complexes]
    complexesConfined = confineComplexes_proteinSet(complexes, proteinSet)
    complexesConfined_lengths = [len(ele) for ele in complexesConfined]
    complexesConfined_origOverlap = \
        [len(newEle.intersection(origEle))/len(origEle) for origEle, newEle in zip(complexes, complexesConfined)]

    complexesProteins_confiningProteins_overlap = \
        len(condenseComplexes_singleProteins(complexes).intersection(proteinSet))

    results = dict(maxComplex_len=max(complexesLengths),
                   meanComplex_len=statistics.mean(complexesLengths),
                   medianComplex_len=statistics.median(complexesLengths),
                   modeComplex_len=statistics.mode(complexesLengths),
                   minComplex_len=min(complexesLengths),
                   maxComplex_confinedLen=max(complexesConfined_lengths),
                   meanComplex_confinedLen=statistics.mean(complexesConfined_lengths),
                   medianComplex_confinedLen=statistics.median(complexesConfined_lengths),
                   modeComplex_confinedLen=statistics.mode(complexesConfined_lengths),
                   minComplex_confinedLen=min(complexesConfined_lengths),
                   origComplexes_confinedComplexes_meanOverlap=statistics.mean(complexesConfined_origOverlap),
                   origComplexes_confinedComplexes_unchanged=len(complexesConfined_origOverlap[complexesConfined_origOverlap == 1]),
                   complexesProteins_confiningProteins_overlap=complexesProteins_confiningProteins_overlap)

    if verbose:
        print('Mean overlap between complexes before and after protein confinement: {0}.'.format(
            results['origComplexes_confinedComplexes_meanOverlap']))
        print('Total complexes unchanged after confinement: {0}.'.format(
            results['origComplexes_confinedComplexes_unchanged']))
        print('____________________________________________________________')

        print('Original complexes stats')
        print('Max complex length: {0}.'.format(results['maxComplex_len']))
        print('Mean complex length: {0}.'.format(results['meanComplex_len']))
        print('Median complex length: {0}.'.format(results['medianComplex_len']))
        print('Mode complex length: {0}.'.format(results['modeComplex_len']))
        print('Min complex length: {0}.'.format(results['minComplex_len']))
        print('____________________________________________________________')

        print('Confined complexes stats')
        print('Max complex confined length: {0}.'.format(results['maxComplex_confinedLen']))
        print('Mean complex confined length: {0}.'.format(results['meanComplex_confinedLen']))
        print('Median complex confined length: {0}.'.format(results['medianComplex_confinedLen']))
        print('Mode complex confined length: {0}.'.format(results['modeComplex_confinedLen']))
        print('Min complex confined length: {0}.'.format(results['minComplex_confinedLen']))
        print('____________________________________________________________')

    if returnResults:
        return results

def compareContrast_referenceExperimental_complexes(refComplexes, expComplexes, verbose=False):
    refComplexes = transformComplexes_variableUse(refComplexes)
    expComplexes = transformComplexes_variableUse(expComplexes)

    refComplexes_proteinSet = condenseComplexes_singleProteins(refComplexes, verbose=True)
    expComplexes_proteinSet = condenseComplexes_singleProteins(expComplexes, verbose=True)

    refComplexes = list(confineComplexes_proteinSet(refComplexes, expComplexes_proteinSet, verbose=True))
    expComplexes_confinedRef_proteinSet = \
        list(confineComplexes_proteinSet(expComplexes, refComplexes_proteinSet, verbose=True))

    lengths = dict()
    exactMatches = dict()
    subsets = dict()
    supersets = dict()
    maxOverlap = dict()
    meanOverlap = dict()
    minOverlap = dict()
    maxPurity = dict()
    meanPurity = dict()
    minPurity = dict()
    #for each reference complex...
    for refEle, refEle_idx in zip(refComplexes, np.arange(len(refComplexes))):

        if verbose:
            print(str(refEle_idx/len(refComplexes)))

        lengths[refEle_idx] = len(refEle)
        #count exact matches, if any.  Note: not penalizing predicted new partners
        exactMatches[refEle_idx] = [1 if refEle in expComplexes_confinedRef_proteinSet else 0][0]
        #count subsets, if any.
        subsets[refEle_idx] = sum([expEle.issubset(refEle) for expEle in expComplexes_confinedRef_proteinSet])
        #count supersets, if any,
        supersets[refEle_idx] = \
            sum([expEle.issuperset(refEle) for expEle in expComplexes])
        #find max, mean, min overlap across predictions (-> recapitulation).  Note: not penalizing predicted new partners
        overlap = [len(expEle.intersection(refEle))/len(refEle) for expEle in expComplexes]
        maxOverlap[refEle_idx] = max(overlap)
        meanOverlap[refEle_idx] = statistics.mean(overlap)
        minOverlap[refEle_idx] = min(overlap)
        #find max, mean, min purity across predictions.  Note: not penalizing predicted new partners
        purity = [len(expEle.intersection(refEle))/len(expEle) for expEle in expComplexes_confinedRef_proteinSet]
        maxPurity[refEle_idx] = max(purity)
        meanPurity[refEle_idx] = statistics.mean(purity)
        minPurity[refEle_idx] = min(purity)

    results = {'complexLengths': lengths, 'exactMatches': exactMatches,
               'subsets': subsets, 'supersets': supersets,
               'maxOverlap': maxOverlap, 'meanOverlap': meanOverlap, 'minOverlap': minOverlap,
               'maxPurity': maxPurity, 'meanPurity': meanPurity, 'minPurity': minPurity}
    pickle.dump(results, open('./results.pkl', 'wb'))

    plt.clf()
    recapPurity_max = pd.DataFrame(list(zip(results['maxOverlap'].values(), results['maxPurity'].values())),
                                  columns=['Max recapitulation', 'Max purity'])
    recapPurity_max.to_csv('./recapPurity_max.csv', sep='\t', index=False)
    sns.scatterplot(data=recapPurity_max, x="Max recapitulation", y="Max purity")
    plt.title('Max recapitulation of reference complexes vs max composition integrity.')
    plt.savefig('./maxRecap_maxPurity_scatter.svg', bbox_inches='tight', dpi=600)

    plt.clf()
    recapPurity_mean = pd.DataFrame(list(zip(results['meanOverlap'].values(), results['meanPurity'].values())),
                                  columns=['Mean recapitulation', 'Mean purity'])
    recapPurity_mean.to_csv('./recapPurity_mean.csv', sep='\t', index=False)
    sns.scatterplot(data=recapPurity_mean, x="Mean recapitulation", y="Mean purity")
    plt.title('Mean recapitulation of reference complexes vs mean composition integrity.')
    plt.savefig('./meanRecap_meanPurity_scatter.svg', bbox_inches='tight', dpi=600)

    plt.clf()
    recapPurity_min = pd.DataFrame(list(zip(results['minOverlap'].values(), results['minPurity'].values())),
                                  columns=['Min recapitulation', 'Min purity'])
    recapPurity_min.to_csv('./recapPurity_min.csv', sep='\t', index=False)
    sns.scatterplot(data=recapPurity_min, x="Min recapitulation", y="Min purity")
    plt.title('Min recapitulation of reference complexes vs min composition integrity.')
    plt.savefig('./minRecap_minPurity_scatter.svg', bbox_inches='tight', dpi=600)

    cliques = []
    cliquesList = list(set(lengths.values()))
    cliquesList.sort()
    for cliqueSize in cliquesList:
        cliqueSize_indices = list(np.flatnonzero(np.array(list(results['complexLengths'].values()))==cliqueSize))
        cliques.append([
            cliqueSize,
            sum([results['exactMatches'][idx] for idx in cliqueSize_indices]),
            sum([results['subsets'][idx] for idx in cliqueSize_indices]),
            sum([results['supersets'][idx] for idx in cliqueSize_indices]),
            sum([1 for expEle in expComplexes if len(expEle)==cliqueSize]),
            sum([not any(expEle.intersection(refEle) for refEle in refComplexes) for expEle in expComplexes if len(expEle)==cliqueSize])
        ])

    cliqueTable = pd.DataFrame(data=cliques,
                               columns=['Clique length',
                                        'Exact matches',
                                        'Subsets of Reference', 'Supersets of Reference',
                                        'Total Predictions', 'Total Novel Predictions'])
    cliqueTable.to_csv('./cliqueTable.csv', index=False, sep='\t')

    '''
    x1 = cliqueTable.loc[:, ['Clique length', 'Matches']]
    x1.insert(1, 'Count Type', 'Matches', allow_duplicates=True)
    x1.rename(columns={'Matches': 'Count'}, inplace=True)

    x2 = cliqueTable.loc[:, ['Clique length', 'Reference Subsets']]
    x2.insert(1, 'Count Type', 'Reference Subsets', allow_duplicates=True)
    x2.rename(columns={'Reference Subsets': 'Count'}, inplace=True)

    x3 = cliqueTable.loc[:, ['Clique length', 'Reference Supersets']]
    x3.insert(1, 'Count Type', 'Reference Supersets', allow_duplicates=True)
    x3.rename(columns={'Reference Supersets': 'Count'}, inplace=True)

    x4 = cliqueTable.loc[:, ['Clique length', 'Total Predictions']]
    x4.insert(1, 'Count Type', 'Total Predictions', allow_duplicates=True)
    x4.rename(columns={'Total Predictions': 'Count'}, inplace=True)

    x5 = cliqueTable.loc[:, ['Clique length', 'Total Novel Predictions']]
    x5.insert(1, 'Count Type', 'Total Novel Predictions', allow_duplicates=True)
    x5.rename(columns={'Total Novel Predictions': 'Count'}, inplace=True)

    cliqueTable_simplified = pd.concat([x1, x2, x3, x4, x5], ignore_index=True)
    '''

    cliqueTable_simplified = pd.melt(cliqueTable,
                                     id_vars=['Clique length'],
                                     value_vars=['Clique length',
                                                 'Exact matches',
                                                 'Subsets of Reference', 'Supersets of Reference',
                                                 'Total Predictions', 'Total Novel Predictions'],
                                     var_name='Count type', value_name='Count')
    cliqueTable_simplified.to_csv('./cliqueTable_simplified.csv', index=False, sep='\t')

    for countType in list(cliqueTable_simplified.loc[:, 'Count type'].unique()):
        try:
            countType_order = \
                math.floor(math.log(cliqueTable_simplified.loc[cliqueTable_simplified['Count Type']==countType, 'Count'].max(), 10))
            if countType_order != 1:
                cliqueTable_simplified.loc[cliqueTable_simplified['Count Type'] == countType, 'Count'] = \
                    cliqueTable_simplified.loc[cliqueTable_simplified['Count Type'] == countType, 'Count']/(10 ** countType_order)
                newLabel = countType + ' (x' + str((10 ** countType_order)) + ')'
                cliqueTable_simplified.loc[cliqueTable_simplified['Count Type']==countType, 'Count Type'] = newLabel
        except:
            pass

    plt.clf()
    sns.catplot(data=cliqueTable_simplified, kind="bar", x="Clique length", y="Count", hue='Count type', height=8, aspect=15/8)
    plt.title('Comparisons of predictions to CORUM (Sep 2018) by clique size.')
    plt.savefig('./cliqueGrouped_CORUMcomparisons.svg', bbox_inches='tight', dpi=600)

def generateTable_cliqueStats_allThresholds(dirPatt, prefix2ID):

    iterIDs = [directory.split(prefix2ID)[1] for directory in glob.glob(dirPatt)]
    iterPaths = [directory for directory in glob.glob(dirPatt)]
    cliqueStats_byThresholds = dict()
    for iterID, iterPath in zip(iterIDs, iterPaths):
        cliqueStats_byThresholds[iterID] = pd.read_csv(iterPath+'cliqueTable_simplified.csv', sep='\t')
        cliqueStats_byThresholds[iterID].insert(0, 'Threshold', iterID, allow_duplicates=True)

    cliqueStats_allThresholds = pd.concat(list(cliqueStats_byThresholds.values()), ignore_index=True)
    cliqueStats_allThresholds.to_csv('./cliqueStats_allThresholds.csv', sep='\t', index=False)
    cliqueStats_allThresholds_agg = cliqueStats_allThresholds.groupby(['Threshold', 'Count Type']).mean()

    plt.title('Comparisons of predictions to CORUM (Sep 2018) by clique size.')
    plt.savefig('./cliqueGrouped_CORUMcomparisons.svg', bbox_inches='tight', dpi=300)

    plt.clf()
    sns.catplot(data=cliqueStats_allThresholds_agg, kind="bar", x="Clique Length", y="Count", hue='Count Type', height=8, aspect=15/8)
    plt.title('Clique stats by threshold.')
    plt.savefig('./cliqueStats_allThresholds.svg', bbox_inches='tight', dpi=300)

def plotScatter_recapPurity_thresholdRange(dirPatt, prefix2ID):
    dirs = glob.glob(dirPatt)
    dirs.sort()
    dirIDs = [directory.split(prefix2ID)[1] for directory in dirs]

    plt.clf()
    plt.title('Recapitulation of reference complexes vs composition integrity.')
    #sns.color_palette("Paired")
    for ident, path in zip(dirIDs, dirs):
        results = pickle.load(open(path+'/results/results.pkl', 'rb'))
        recapPlus_purity = pd.DataFrame(
            list(
                zip(list(results['maxOverlap'].values()), list(results['maxPurity'].values()))
            )
        )
        recapPlus_purity.insert(0, 'Threshold', 'Threshold', allow_duplicates=True)
        sns.scatterplot(data=recapPlus_purity, x="Recapitulation", y="Purity", label=ident)

    plt.legend(title='Thresholds')
    plt.savefig('referenceRecap_vPurity_scatter_rangeThresholds.svg', bbox_inches='tight', dpi=1200, format='svg')

def plotScatter_meanMaxima_thresholdRange(dirPatt, prefix2ID):
    dirs = glob.glob(dirPatt)
    dirs.sort()
    dirIDs = [directory.split(prefix2ID)[1] for directory in dirs]

    recapMean_maxima = []
    purityMean_maxima = []
    for path in dirs:
        results = pickle.load(open(path + '/results/results.pkl', 'rb'))
        recapMean_maxima.append(statistics.mean(list(results['maxOverlap'].values())))
        purityMean_maxima.append(statistics.mean(list(results['maxPurity'].values())))

    meanMaxima_recap = pd.DataFrame(list(zip(dirIDs, recapMean_maxima)), columns=['Threshold', 'Reference Likeness'])
    meanMaxima_purity = pd.DataFrame(list(zip(dirIDs, purityMean_maxima)), columns=['Threshold', 'Reference Likeness'])

    plt.clf()
    plt.title('Mean recapitulation and purity maxima vs threshold.')
    sns.color_palette("Paired")
    sns.scatterplot(data=meanMaxima_recap, x="Threshold", y="Reference Likeness", label="Recapitulation")
    sns.scatterplot(data=meanMaxima_purity, x="Threshold", y="Reference Likeness", label="Purity")
    plt.legend(title='Measures')
    plt.savefig('./referenceRecap_vPurity_scatter_rangeThresholds.svg', bbox_inches='tight', dpi=300, format='svg')

def writelist_ofComplexes_toFile(listComplexes, filename):
    with open(filename, 'a') as the_file:
        for cplx in listComplexes:
            for prot in list(cplx):
                the_file.write(str(int(prot)) + '\t')
            the_file.write('\n')

def writeRestricted_refCplx_toFile(predsSeries, refComplexes, filedir):
    for upperkey in predsSeries.keys():
        for lowerkey in predsSeries[upperkey].keys():
            tempVar_geneBase = set().union(*[set(cplx) for cplx in list(predsSeries[upperkey][lowerkey])])
            tempVar_complexesRestricted = list(
                set([frozenset(set(cplx).intersection(tempVar_geneBase)) for cplx in list(refComplexes) if
                     set(cplx).intersection(tempVar_geneBase)]))
            writelist_ofComplexes_toFile(
                tempVar_complexesRestricted,
                filedir + 'test_complexes_' + str(upperkey) + '_' + str(lowerkey) + '_restricted.txt')

class evaluatePredictions(object):

    def __init__(self, preds,
                 ref, altPreds=None,
                 experiment_name='partition model', experiment_cell_type=None, experiment_threshold=None,
                 penalize_missing_ref_prots=False):

        self.ref = ref
        self.refCplx_lens = np.array([len(cplx) for cplx in ref])
        self.refLens, self.refLens_dist = np.unique(self.refCplx_lens, return_counts=True)
        #probably won't use...
        self.refProts_base = set().union(*[set(cplx) for cplx in list(ref)])

        self.preds = preds
        if type(preds) is not dict:
            self.predCplx_lens = np.array([len(cplx) for cplx in preds])
            self.lens, self.lensDist = np.unique(self.predCplx_lens, return_counts=True)
        else:
            self.predCplx_lens = dict()
            self.lens = dict()
            self.lensDist = dict()
            for thresh in preds.keys():
                self.predCplx_lens[thresh] = np.array([len(cplx) for cplx in preds[thresh]])
                self.lens[thresh], self.lensDist[thresh] = \
                    np.unique(self.predCplx_lens[thresh], return_counts=True)

        self.altPreds = altPreds
        if self.altPreds:
            if type(altPreds) is not dict:
                print('Only input type dictionary for comparative predictions.')
                return
            else:
                self.altPreds_cplxLens = dict()
                self.altPreds_lens = dict()
                self.altPreds_lensDist = dict()
                for method in altPreds.keys():
                    self.altPreds_cplxLens[method] = np.array([len(cplx) for cplx in altPreds[method]])
                    self.altPreds_lens[method], self.altPreds_lensDist[method] = \
                        np.unique(self.altPreds_cplxLens[method], return_counts=True)

        self.experiment_name = experiment_name
        self.experiment_cell_type = experiment_cell_type
        self.experiment_threshold = experiment_threshold
        self.penalize_missing_ref_prots = penalize_missing_ref_prots
        self.predSummary = None
        self.predSummary_compact = None
        self.fidelityStats = None
        self.fidelityStats_compact = None

    def describe_predictions(self, imgDest=None, corumVersion='(version 2.0)'):
        #describe predictions and compare to reference
        cmp_cats = ['exact matches', 'subsets', 'supersets', 'completely novel--no pairs', 'novel--no common triplets', 'total predicted']
        if type(self.preds) is not dict:
            self.predSummary = pd.DataFrame(columns=cmp_cats, index=self.lens)
            self.predSummary.index.name = 'clique length'
            for cliqueLen in self.predSummary.index:
                cplxSets = set([cplx for cplx in list(self.preds) if len(cplx)==cliqueLen])
                self.predSummary.loc[cliqueLen, 'exact matches'] = self.findMatches(self.preds, cplxSets)
                self.predSummary.loc[cliqueLen, 'subsets'] = self.findSubsets(self.preds, cliqueLen)
                self.predSummary.loc[cliqueLen, 'supersets'] = self.findSupersets(self.preds, cliqueLen)
                self.predSummary.loc[cliqueLen, ['completely novel--no pairs', 'novel--no common triplets']] = \
                    self.findNovel(self.preds, cliqueLen)
                self.predSummary.loc[cliqueLen, 'total predicted'] = self.countAll(self.preds, cliqueLen)
            self.predSummary_compact = pd.melt(self.predSummary.reset_index(),
                                               id_vars=['clique length'],
                                               value_vars=cmp_cats, var_name='count type', value_name='count')

            if imgDest:
                self.predSummary_compact['count'] = self.predSummary_compact['count'].astype('int')
                fig, ax = plt.subplots(figsize=(20, 10))
                sns.barplot(data=self.predSummary_compact.loc[self.predSummary_compact['count type'].isin(
                    ['completely novel--no pairs', 'novel--no common triplets', 'total predicted'])],
                            x='clique length', y='count', hue='count type', palette='pastel', ax=ax)
                plt.legend(bbox_to_anchor=(1.0, 1.0))
                ax_twin = ax.twinx()
                sns.lineplot(data=self.predSummary_compact.loc[self.predSummary_compact['count type'].isin(
                    ['exact matches', 'subsets', 'supersets'])],
                             x='clique length', y='count', hue='count type', palette='crest' ,ax=ax_twin)
                plt.legend(bbox_to_anchor=(1.0, 0.5))
                plt.title('Comparison of predicted complexes against \n gold standard CORUM ' + corumVersion + '.')
                plt.savefig(imgDest+'modelPredictions_noveltySummary.svg',
                            bbox_inches='tight', dpi=1200, format='svg')

        else:
            multiCell_index = []
            for upperKey in self.preds.keys():
                for lowerKey in self.preds[upperKey].keys():
                    cliqueLengths = np.unique([len(ele) for ele in list(self.preds[upperKey][lowerKey])])
                    for lenth in cliqueLengths:
                        multiCell_index.append([upperKey, lowerKey, lenth])
            multiCell_index = np.squeeze(np.array(multiCell_index, ndmin=2))
            multiCell_index = {'cell type': multiCell_index[:, 0], 'threshold': multiCell_index[:, 1],
                                 'clique length': multiCell_index[:, 2]}
            self.predSummary = pd.DataFrame(data=multiCell_index)
            self.predSummary = \
                self.predSummary.reindex(columns=[*self.predSummary.columns.tolist(), *cmp_cats], fill_value=np.nan)
            self.predSummary.set_index(['cell type', 'threshold', 'clique length'],
                                        drop=True, inplace=True, verify_integrity=True)
            for independentIndex, jointIndex in enumerate(self.predSummary.index):
                cplxSets = set([cplx for cplx in list(self.preds[jointIndex[0]][jointIndex[1]]) if len(cplx)==int(jointIndex[2])])
                self.predSummary.iloc[independentIndex, list(self.predSummary.columns).index('exact matches')] = \
                    self.findMatches(self.preds[jointIndex[0]][jointIndex[1]], cplxSets)
                self.predSummary.iloc[independentIndex, list(self.predSummary.columns).index('subsets')] = \
                    self.findSubsets(self.preds[jointIndex[0]][jointIndex[1]], int(jointIndex[2]))
                self.predSummary.iloc[independentIndex, list(self.predSummary.columns).index('supersets')] = \
                    self.findSupersets(self.preds[jointIndex[0]][jointIndex[1]], int(jointIndex[2]))
                self.predSummary.iloc[independentIndex,
                                       [list(self.predSummary.columns).index('completely novel--no pairs'),
                                        list(self.predSummary.columns).index('novel--no common triplets')]] = \
                    self.findNovel(self.preds[jointIndex[0]][jointIndex[1]], int(jointIndex[2]))
                self.predSummary.iloc[independentIndex, list(self.predSummary.columns).index('total predicted')] = \
                    self.countAll(self.preds[jointIndex[0]][jointIndex[1]], int(jointIndex[2]))
            self.predSummary_compact = pd.melt(self.predSummary.reset_index(),
                                               id_vars=['cell type', 'threshold', 'clique length'],
                                               value_vars=cmp_cats,
                                               var_name='count type', value_name='count')

            if imgDest:
                os.makedirs(imgDest+'/comparisonsBy_prtnModel_specVariations/')
                self.predSummary_compact['count'] = self.predSummary_compact['count'].astype('int')
                for cellType in self.predSummary_compact.loc[:, 'cell type'].unique():
                    for threshold in self.predSummary_compact.loc[
                        self.predSummary_compact['cell type'] == cellType, 'threshold'
                    ].unique():
                        predSummary_select = \
                            self.predSummary_compact.loc[
                                (self.predSummary_compact['cell type'] == cellType) &
                                (self.predSummary_compact['threshold'] == threshold)].drop(
                                columns=['cell type', 'threshold'])
                        plt.cla()
                        fig, ax = plt.subplots(figsize=(20, 10))
                        sns.barplot(data=predSummary_select.loc[predSummary_select['count type'].isin(
                            ['completely novel--no pairs', 'novel--no common triplets', 'total predicted'])],
                                    x='clique length', y='count', hue='count type', palette='pastel', ax=ax)
                        plt.legend(bbox_to_anchor=(1.0, 1.0))
                        ax_twin = ax.twinx()
                        sns.lineplot(data=predSummary_select.loc[predSummary_select['count type'].isin(
                            ['exact matches', 'subsets', 'supersets'])],
                                     x='clique length', y='count', hue='count type', palette='crest' ,ax=ax_twin)
                        plt.legend(bbox_to_anchor=(1.0, 0.5))
                        plt.suptitle('Comparison of predicted complexes against \n gold standard CORUM ' +
                                  corumVersion + '.')
                        plt.title('Partition model -- Cell type: ' + cellType + '-- ' + threshold)
                        plt.savefig(imgDest+'/comparisonsBy_prtnModel_specVariations/' +
                                    'modelPredictions_noveltySummary__prtnModel-' + cellType +
                                    '@thresh' + threshold + '.svg', bbox_inches='tight', dpi=1200, format='svg')
                        plt.cla()

    def assess_fidelity(self, imgDest=None):

        measure_columns = ['recapitulation max', 'recapitulation mean', 'purity max', 'purity mean']
        if type(self.preds) is not dict:
            if (self.experiment_cell_type is None) | (self.experiment_threshold is None):
                print('please specify experiment specs, (e.g., model name, cell type, threshold)')
                return

            else:
                fidelity_index_preds = \
                    [[list(cplx), self.experiment_name, self.experiment_cell_type, self.experiment_threshold]
                        for cplx in list(self.ref)]
        else:
            fidelity_index_preds = []
            for cplx in list(self.ref):
                for upperKey in self.preds.keys():
                    for lowerKey in self.preds[upperKey].keys():
                        fidelity_index_preds.append([list(cplx), self.experiment_name, upperKey, lowerKey])

        fidelity_index_preds = np.squeeze(np.array(fidelity_index_preds, ndmin=2, dtype=object))
        fidelity_index_preds = {'reference complex': fidelity_index_preds[:, 0],
                                'source': fidelity_index_preds[:, 1],
                                'cell type': fidelity_index_preds[:, 2],
                                'threshold': fidelity_index_preds[:, 3]}
        fidelityStats_preds = pd.DataFrame(data=fidelity_index_preds)
        fidelityStats_preds = \
            fidelityStats_preds.reindex(columns=[*fidelityStats_preds.columns.tolist(), *measure_columns],
                                       fill_value=np.nan)

        if type(self.preds) is not dict:
            for idx in fidelityStats_preds.index:
                fidelityStats_preds.loc[idx, ['recapitulation max', 'recapitulation mean']] = \
                    self.measureRecapitulation(fidelityStats_preds.loc[idx, 'reference complex'], self.preds)
                fidelityStats_preds.loc[idx, ['purity max', 'purity mean']] = \
                    self.measurePurity(fidelityStats_preds.loc[idx, 'reference complex'], self.preds)
        else:
            for idx in fidelityStats_preds.index:
                fidelityStats_preds.loc[idx, ['recapitulation max', 'recapitulation mean']] = \
                    self.measureRecapitulation(fidelityStats_preds.loc[idx, 'reference complex'],
                                               self.preds[fidelityStats_preds.loc[idx, 'cell type']][fidelityStats_preds.loc[idx, 'threshold']])
                fidelityStats_preds.loc[idx, ['purity max', 'purity mean']] = \
                    self.measurePurity(fidelityStats_preds.loc[idx, 'reference complex'],
                                       self.preds[fidelityStats_preds.loc[idx, 'cell type']][fidelityStats_preds.loc[idx, 'threshold']])

        if self.altPreds:
            fidelity_index_altPreds = []
            for cplx in list(self.ref):
                for method in self.altPreds.keys():
                    fidelity_index_altPreds.append([list(cplx), method, '--', '--'])

            fidelity_index_altPreds = np.squeeze(np.array(fidelity_index_altPreds, ndmin=2, dtype=object))
            fidelity_index_altPreds = {'reference complex': fidelity_index_altPreds[:, 0],
                                       'source': fidelity_index_altPreds[:, 1],
                                       'cell type': fidelity_index_altPreds[:, 2],
                                       'threshold': fidelity_index_altPreds[:, 3]}
            fidelityStats_altPreds = pd.DataFrame(data=fidelity_index_altPreds)
            fidelityStats_altPreds = \
                fidelityStats_altPreds.reindex(columns=[*fidelityStats_altPreds.columns.tolist(), *measure_columns],
                                           fill_value=np.nan)

            for idx in fidelityStats_altPreds.index:
                fidelityStats_altPreds.loc[idx, ['recapitulation max', 'recapitulation mean']] = \
                    self.measureRecapitulation(fidelityStats_altPreds.loc[idx, 'reference complex'],
                                               self.altPreds[fidelityStats_altPreds.loc[idx, 'source']])
                fidelityStats_altPreds.loc[idx, ['purity max', 'purity mean']] = \
                    self.measurePurity(fidelityStats_altPreds.loc[idx, 'reference complex'],
                                       self.altPreds[fidelityStats_altPreds.loc[idx, 'source']])

            self.fidelityStats = pd.concat([fidelityStats_preds, fidelityStats_altPreds])
            self.fidelityStats_compact = pd.melt(self.fidelityStats,
                                                 id_vars=['source', 'cell type', 'threshold', 'reference complex'],
                                                 value_vars=measure_columns,
                                                 var_name='fidelity type', value_name='fidelity measure')

            performancesMeasured_byMethod = dict()
            for modelMethod in self.fidelityStats_compact['source'].unique():

                if modelMethod != 'partition model':

                    fidelities = dict()
                    for fidelityMeasure in self.fidelityStats_compact['fidelity type'].unique():
                        simpleDF = \
                            self.fidelityStats_compact.loc[
                                (self.fidelityStats_compact['source'] == modelMethod) &
                                (self.fidelityStats_compact['fidelity type'] == fidelityMeasure)].copy()
                        simpleDF['k length'] = simpleDF.loc[:, 'reference complex'].apply(lambda x: len(x))
                        kStats = np.array([[k,
                                            simpleDF.loc[simpleDF['k length'] == k, 'fidelity measure'].mean(),
                                            simpleDF.loc[simpleDF['k length'] == k, 'fidelity measure'].max()] for k in
                                           simpleDF['k length'].unique()], ndmin=2)
                        kStats_sorted = kStats[kStats[:, 0].argsort()]
                        fidelities[fidelityMeasure] = kStats_sorted
                    performancesMeasured_byMethod[modelMethod] = fidelities

                elif modelMethod == 'partition model':

                    for cellType in self.fidelityStats_compact.loc[
                        self.fidelityStats_compact['source'] == modelMethod, 'cell type'].unique():

                        for threshold in self.fidelityStats_compact.loc[
                            self.fidelityStats_compact['cell type'] == cellType, 'threshold'].unique():
                            modelMethod_keyAmended = modelMethod + ' ' + cellType + ' @' + threshold

                            fidelities = dict()
                            for fidelityMeasure in self.fidelityStats_compact['fidelity type'].unique():
                                simpleDF = \
                                    self.fidelityStats_compact.loc[
                                        (self.fidelityStats_compact['source'] == modelMethod) &
                                        (self.fidelityStats_compact['threshold'] == threshold) &
                                        (self.fidelityStats_compact['cell type'] == cellType) &
                                        (self.fidelityStats_compact['fidelity type'] == fidelityMeasure)].copy()
                                simpleDF['k length'] = simpleDF.loc[:, 'reference complex'].apply(lambda x: len(x))
                                kStats = np.array([[k,
                                                    simpleDF.loc[simpleDF['k length'] == k, 'fidelity measure'].mean(),
                                                    simpleDF.loc[simpleDF['k length'] == k, 'fidelity measure'].max()]
                                                   for k in
                                                   simpleDF['k length'].unique()], ndmin=2)
                                kStats_sorted = kStats[kStats[:, 0].argsort()]
                                fidelities[fidelityMeasure] = kStats_sorted

                            performancesMeasured_byMethod[modelMethod_keyAmended] = fidelities

        else:
            self.fidelityStats = fidelityStats_preds
            self.fidelityStats_compact = pd.melt(self.fidelityStats,
                                                 id_vars=['source', 'cell type', 'threshold', 'reference complex'],
                                                 value_vars=measure_columns,
                                                 var_name='fidelity type', value_name='fidelity measure')
            performancesMeasured_byMethod = dict()
            for modelMethod in self.fidelityStats_compact['source'].unique():

                if modelMethod != 'partition model':

                    fidelities = dict()
                    for fidelityMeasure in self.fidelityStats_compact['fidelity type'].unique():
                        simpleDF = \
                            self.fidelityStats_compact.loc[
                                (self.fidelityStats_compact['source'] == modelMethod) &
                                (self.fidelityStats_compact['fidelity type'] == fidelityMeasure)].copy()
                        simpleDF['k length'] = simpleDF.loc[:, 'reference complex'].apply(lambda x: len(x))
                        kStats = np.array([[k,
                                            simpleDF.loc[simpleDF['k length'] == k, 'fidelity measure'].mean(),
                                            simpleDF.loc[simpleDF['k length'] == k, 'fidelity measure'].max()] for k in
                                           simpleDF['k length'].unique()], ndmin=2)
                        kStats_sorted = kStats[kStats[:, 0].argsort()]
                        fidelities[fidelityMeasure] = kStats_sorted
                    performancesMeasured_byMethod[modelMethod] = fidelities

                elif modelMethod == 'partition model':

                    for cellType in self.fidelityStats_compact.loc[
                        self.fidelityStats_compact['source'] == modelMethod, 'cell type'].unique():

                        for threshold in self.fidelityStats_compact.loc[
                            self.fidelityStats_compact['cell type'] == cellType, 'threshold'].unique():
                            modelMethod_keyAmended = modelMethod + ' ' + cellType + ' @' + threshold

                            fidelities = dict()
                            for fidelityMeasure in self.fidelityStats_compact['fidelity type'].unique():
                                simpleDF = \
                                    self.fidelityStats_compact.loc[
                                        (self.fidelityStats_compact['source'] == modelMethod) &
                                        (self.fidelityStats_compact['threshold'] == threshold) &
                                        (self.fidelityStats_compact['cell type'] == cellType) &
                                        (self.fidelityStats_compact['fidelity type'] == fidelityMeasure)].copy()
                                simpleDF['k length'] = simpleDF.loc[:, 'reference complex'].apply(lambda x: len(x))
                                kStats = np.array([[k,
                                                    simpleDF.loc[simpleDF['k length'] == k, 'fidelity measure'].mean(),
                                                    simpleDF.loc[simpleDF['k length'] == k, 'fidelity measure'].max()]
                                                   for k in
                                                   simpleDF['k length'].unique()], ndmin=2)
                                kStats_sorted = kStats[kStats[:, 0].argsort()]
                                fidelities[fidelityMeasure] = kStats_sorted

                            performancesMeasured_byMethod[modelMethod_keyAmended] = fidelities

            return

        if imgDest:

            # complete later
            if type(self.preds) is not dict:
                print('error, undeveloped ....')
                return

            else:
                os.makedirs(imgDest + '/comparisonsBy_prtnModel_fidelityMeasures/')
                for fidelityMeasure in ['recapitulation max', 'recapitulation mean', 'purity max', 'purity mean']:
                    plt.cla()
                    fig, (ax_max, ax_mean) = plt.subplots(figsize=(20, 10), nrows=2)
                    plt.suptitle(fidelityMeasure)
                    ax_max.set_title('Max')
                    ax_mean.set_title('Mean')
                    for performanceKey in sorted(performancesMeasured_byMethod):
                        if 'partition model' not in performanceKey:
                            sns.lineplot(x=performancesMeasured_byMethod[performanceKey][fidelityMeasure][:, 0],
                                         y=performancesMeasured_byMethod[performanceKey][fidelityMeasure][:, 2],
                                         linewidth=5, label=performanceKey, ax=ax_max)
                            sns.lineplot(x=performancesMeasured_byMethod[performanceKey][fidelityMeasure][:, 0],
                                         y=performancesMeasured_byMethod[performanceKey][fidelityMeasure][:, 1],
                                         linewidth=5, label=performanceKey, ax=ax_mean)
                        elif 'partition model' in performanceKey:
                            sns.lineplot(x=performancesMeasured_byMethod[performanceKey][fidelityMeasure][:, 0],
                                         y=performancesMeasured_byMethod[performanceKey][fidelityMeasure][:, 2],
                                         label=performanceKey.replace('partition model@', ''), ax=ax_max)
                            sns.lineplot(x=performancesMeasured_byMethod[performanceKey][fidelityMeasure][:, 0],
                                         y=performancesMeasured_byMethod[performanceKey][fidelityMeasure][:, 1],
                                         label=performanceKey.replace('partition model@', ''), ax=ax_mean)
                    ax_max.legend(ncol=2, bbox_to_anchor=(1.0, 1.0))
                    ax_mean.legend(ncol=2, bbox_to_anchor=(1.0, 1.0))
                    plt.savefig(imgDest + '/comparisonsBy_prtnModel_fidelityMeasures/' +
                                'modelPredictions_fidelityMeasures__' + fidelityMeasure + '.svg',
                                bbox_inches='tight', dpi=1200, format='svg')

    def findMatches(self, predictedComplexes, kCliques):
        if not self.penalize_missing_ref_prots:
            restrictingProts = set().union(*[set(cplx) for cplx in list(predictedComplexes)])
            matches = \
                kCliques.intersection(
                    set([frozenset(set(cplx).intersection(restrictingProts)) for cplx in list(self.ref)])
                )
        else:
            matches = kCliques.intersection(predictedComplexes)
        numMatches = len(matches)

        return numMatches

    #find subsets of reference complexes among predictions
    def findSubsets(self, predictedComplexes, cliqueLen):
        if not self.penalize_missing_ref_prots:
            restrictingProts = set().union(*[set(cplx) for cplx in list(predictedComplexes)])
            refCplx = set([frozenset(set(cplx).intersection(restrictingProts)) for cplx in list(self.ref) if set(cplx).intersection(restrictingProts)])
        else:
            refCplx = self.ref

        selectPredictions = [cplx for cplx in list(predictedComplexes) if len(cplx)==cliqueLen]
        numSubsets = 0
        for ref in list(refCplx):
            for cplx in selectPredictions:
                if len(ref) > len(cplx):
                    numSubsets+=set(ref).issuperset(set(cplx))

        return numSubsets

    #find supersets of reference complexes among predictions
    def findSupersets(self, predictedComplexes, cliqueLen):
        if not self.penalize_missing_ref_prots:
            restrictingProts = set().union(*[set(cplx) for cplx in list(predictedComplexes)])
            refCplx = set([frozenset(set(cplx).intersection(restrictingProts)) for cplx in list(self.ref) if set(cplx).intersection(restrictingProts)])
        else:
            refCplx = self.ref

        selectPredictions = [cplx for cplx in list(predictedComplexes) if len(cplx)==cliqueLen]
        numSupersets = 0
        for ref in list(refCplx):
            for cplx in selectPredictions:
                if (len(ref) < len(cplx)) & (len(ref) > 2):
                    numSupersets+=set(cplx).issuperset(set(ref))

        return numSupersets

    #find all completely novel predictions against ref (:= not even a single pair in common)
    def findNovel(self, predictedComplexes, cliqueLen):
        if not self.penalize_missing_ref_prots:
            restrictingProts = set().union(*[set(cplx) for cplx in list(predictedComplexes)])
            refCplx_adj = set([frozenset(set(cplx).intersection(restrictingProts)) for cplx in list(self.ref) if set(cplx).intersection(restrictingProts)])
            refCplx = set([cplx for cplx in list(refCplx_adj) if len(cplx)>2])
        else:
            refCplx = self.ref

        selectPredictions = [cplx for cplx in list(predictedComplexes) if len(cplx)==cliqueLen]
        numNovel_completely = 0
        numNovel = 0
        for cplx in selectPredictions:
            numNovel_completely+=all([1 if len(set(ref).intersection(set(cplx)))<2 else 0 for ref in list(refCplx)])
            numNovel+=all([1 if len(set(ref).intersection(set(cplx)))<3 else 0 for ref in list(refCplx)])

        return numNovel_completely, numNovel

    #count all predictions
    @staticmethod
    def countAll(predictedComplexes, cliqueLen):
        total = len([cplx for cplx in list(predictedComplexes) if len(cplx)==cliqueLen])

        return total

    @staticmethod
    def measureRecapitulation(refCplx, preds):
        #refCplx: single complex, list of proteins (strings)
        #preds: set of multiple complexes, each a frozenset comprised from a list of proteins (strings)
        #compare each reference complex to each prediction
        comparisons = [len(set(refCplx).intersection(set(cplx)))/len(set(refCplx)) for cplx in list(preds)]
        maxRecapitulation = max(comparisons)
        meanRecapitulation = sum(comparisons)/len(comparisons)
        return maxRecapitulation, meanRecapitulation

    @staticmethod
    def measurePurity(refCplx, preds):
        #refCplx: single complex, list of proteins (strings)
        #preds: set of multiple complexes, each a frozenset comprised from a list of proteins (strings)
        #compare each reference complex to each prediction
        comparisons = [len(set(refCplx).intersection(set(cplx)))/len(set(cplx)) for cplx in list(preds)]
        maxPurity = max(comparisons)
        meanPurity = sum(comparisons)/len(comparisons)
        return maxPurity, meanPurity
