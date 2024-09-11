import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
from sklearn.model_selection import KFold
import time


def dataframeSummary(dataframe, dataframeName):
    print('Printing {0} summary...'.format(dataframeName))
    print('Printing columns of data...')
    print(dataframe.columns)
    print('Printing the shape of data...')
    print(dataframe.shape)
    print('Printing head of data...')
    print(dataframe.head())
    print('____________')


def acquire(dataPath, labelsPath, headings=True, chunkSize=None):

    print('Loading data...')
    if headings:
        data = pd.read_csv(dataPath, sep='\t', chunksize=chunkSize)
        try:
            labels = pd.read_csv(labelsPath, sep='\t', names=['label', 'IDi', 'IDii'], header=0, chunksize=chunkSize)
        except:
            labels = pd.read_csv(labelsPath, sep='\t', names=['label'], header=0, chunksize=chunkSize)
    else:
        data = pd.read_csv(dataPath, sep='\t', header=None, chunksize=chunkSize)
        labels = pd.read_csv(labelsPath, sep='\t', header=None, chunksize=chunkSize)

    if hasattr(data, 'shape') and hasattr(labels, 'shape'):
        print('The shape of the data is {0}'.format(data.shape))
        print('The shape of labels is {0}'.format(labels.shape))

    print('Completed loading data')
    print('________________________________________________')

    return data, labels


def createCV_trainTest_splits(workDir,
                              trainingData_path, trainingLabels_path, testData_path, testLabels_path,
                              numFolds, expName, headings=False):

    trainingData, trainingLabels = acquire(trainingData_path, trainingLabels_path, headings=headings)
    testData, testLabels = acquire(testData_path, testLabels_path, headings=headings)

    combinedData = trainingData.append(testData)
    print('combinedData shape: {0}'.format(combinedData.shape))
    combinedLabels = trainingLabels.append(testLabels)
    print('combinedLabels shape: {0}'.format(combinedData.shape))

    kf = KFold(n_splits=numFolds, shuffle=True, random_state=42)
    CVdatasets_path = workDir + expName + '_5CV/'
    os.mkdir(CVdatasets_path)
    fold = 0
    for trainingIndices, testIndices in kf.split(np.arange(len(combinedData))):

        dir = CVdatasets_path + '{0}_5CV_{1}/'.format(expName, fold)
        os.mkdir(dir)

        trainingData_fold = combinedData.iloc[trainingIndices, :]
        print('trainingData_fold shape: {0}'.format(trainingData_fold.shape))
        trainingData_fold.to_csv(dir + 'trainingData.tsv', sep='\t', header=False, index=False)

        trainingLabels_fold = combinedLabels.iloc[trainingIndices, :]
        print('trainingLabels_fold shape: {0}'.format(trainingLabels_fold.shape))
        trainingLabels_fold.to_csv(dir + 'trainingLabels.tsv', sep='\t', header=False, index=False)

        testData_fold = combinedData.iloc[testIndices, :]
        print('testData_fold shape: {0}'.format(testData_fold.shape))
        testData_fold.to_csv(dir + 'testData.tsv', sep='\t', header=False, index=False)

        testLabels_fold = combinedLabels.iloc[testIndices, :]
        print('testLabels_fold shape: {0}'.format(testLabels_fold.shape))
        testLabels_fold.to_csv(dir + 'testLabels.tsv', sep='\t', header=False, index=False)

        fold += 1


def getProt_mapping(gene_2uni_path, uniprot_2geneID_path):

    blockTime_start = time.time()

    print('Establishing mapping...')
    genecard_to_uniprot = pd.read_excel(gene_2uni_path, sheet_name=None)
    genecard_to_uniprot = genecard_to_uniprot['Sheet0']
    genecard_to_uniprot = genecard_to_uniprot[["yourlist:M20200428A94466D2655679D1FD8953E075198DA8845C47J", "Entry"]]
    genecard_to_uniprot = genecard_to_uniprot.rename(
        columns={"yourlist:M20200428A94466D2655679D1FD8953E075198DA8845C47J": "GeneCardID", "Entry": "uniprot"})
    uniprot_to_geneID = pd.read_excel(uniprot_2geneID_path, sheet_name=None)
    uniprot_to_geneID = uniprot_to_geneID['scbc_uniprot->geneIDs']
    uniprot_to_geneID['geneid'] = uniprot_to_geneID['geneid'].astype(str)
    geneIDmapping = uniprot_to_geneID.merge(genecard_to_uniprot, how='inner', on="uniprot")
    geneIDmapping = geneIDmapping[["GeneCardID", "geneid"]]
    geneIDmapping = dict(zip(geneIDmapping["GeneCardID"], geneIDmapping["geneid"]))

    blockTime_end = time.time()
    print('Time elapsed: {elapsed}'.format(elapsed=blockTime_end - blockTime_start))

    print('Completed establishing mapping...')
    print('________________________________________________')

    return geneIDmapping


def makeSCBC_dict(outputDir, scbcData_path, scbc_cellLines, gene_2uni_path, uniprot_2geneID_path):

    blockTime_start = time.time()
    print('Making SCBC data dictionary...')
    print('____________________________________________________________')

    print('     mapping genes in SCBC to IDs in JLM original data matrix...')
    corrMatsDict = {cellLine: [] for cellLine in scbc_cellLines}
    geneIDmapping = getProt_mapping(gene_2uni_path, uniprot_2geneID_path)
    print('     completed mapping genes from SCBC to JLM original data matrix')
    print('     ________________________________________________')

    for cellLine in scbc_cellLines:
        print('Cell line')
        print(cellLine)
        print('____________')
        placeholder = pd.read_excel(scbcData_path, sheet_name=cellLine).iloc[:, :-1]
        dataframeSummary(placeholder, cellLine)
        try:
            placeholder.set_index('Protein', inplace=True, verify_integrity=True)
        except:
            placeholder.drop_duplicates(subset='Protein', keep=False, inplace=True)
            placeholder.set_index('Protein', inplace=True, verify_integrity=True)  # possibly problematic
        placeholder = placeholder.T.corr().rename(
            columns=geneIDmapping,
            index=geneIDmapping).rename_axis(None).rename_axis(None, axis=1).stack().reset_index()
        dataframeSummary(placeholder, '{0} PPI correlation feature'.format(cellLine))
        corrMatsDict[cellLine] = placeholder.rename(columns={'level_0': 'prot1', 'level_1': 'prot2', 0: 'corr'})
        corrMatsDict[cellLine]['prot1'] = corrMatsDict[cellLine]['prot1'].astype('str')
        corrMatsDict[cellLine]['prot2'] = corrMatsDict[cellLine]['prot2'].astype('str')
    print('     Confirming properties of SCBC cell line dictionary, particularly its cell line specific correlation ' +
          'feature...')
    for cellLine in scbc_cellLines:
        placeholder = corrMatsDict[cellLine]
        dataframeSummary(placeholder, cellLine)

    print('     Saving cell line specific correlation data dictionary to disk...')
    filename = outputDir + 'cell-line_specific_correlation_data_dict.pkl'
    pickle.dump(corrMatsDict, open(filename, 'wb'))
    print('     Completed saving cell line specific correlation data dictionary to disk...')

    blockTime_end = time.time()
    print('Time elapsed: {elapsed}'.format(elapsed=blockTime_end - blockTime_start))

    print('Completed making SCBC data dictionary...')
    print('____________________________________________________________')

    return corrMatsDict


def mergeSCBC(data, labels, scbc_cellLines, scbcData_dict):

    blockTime_start = time.time()

    dataframeSummary(data, 'original data')
    dataframeSummary(labels, 'labels')

    print('Merging data with corresponding IDs...')
    ppiIDs = labels.iloc[:, 1:].astype(str)
    data = pd.concat([ppiIDs, data], axis=1)
    dataframeSummary(data, 'ID merged data')

    print('Merging data with SCBC features...')
    for cellLine in scbc_cellLines:
        data = data.merge(scbcData_dict[cellLine].astype(str),
                          how='left',
                          left_on=['IDi', 'IDii'],
                          right_on=['prot1', 'prot2'])
        data = data.rename(columns={'corr': 'corr_' + cellLine})
    print('Completed merging data with SCBC features')
    print('________________________________________________')

    dataframeSummary(data, 'scbc-merged data')

    blockTime_end = time.time()
    print('Time elapsed: {elapsed}'.format(elapsed=blockTime_end - blockTime_start))

    return data


def generateSCBC_data(outputDir, expName,
                      dataPath, labelsPath,
                      scbc_cellLines, scbcData_dictPath=None, scbcData_path=None,
                      gene_2uni_path=None, uniprot_2geneID_path=None):

    print('Generating and saving SCBC merged data to disk...')

    data, labels = acquire(dataPath, labelsPath)
    dataframeSummary(data, 'data')
    dataframeSummary(labels, 'training labels with IDs')

    print('printing scbc_cellLines')
    print(scbc_cellLines)

    if scbcData_dictPath:
        scbcData_dict = pickle.load(open(scbcData_dictPath, 'rb'))
    else:
        scbcData_dict = makeSCBC_dict(outputDir, scbcData_path, scbc_cellLines, gene_2uni_path, uniprot_2geneID_path)

    scbcMerged_data = mergeSCBC(data, labels, scbc_cellLines, scbcData_dict)
    dataframeSummary(scbcMerged_data, 'merged data')
    print('saving SCBC merged data to disk...')

    try:
        filename = outputDir + 'SCBC-merged_{0}_ppiData.pkl'.format(expName)
        pickle.dump(scbcMerged_data, open(filename, 'wb'))
    except:
        filename = outputDir + 'SCBC-merged_{0}_ppiData.tsv'.format(expName)
        scbcMerged_data.to_csv(filename, sep='\t', index=False)

    print('Completed saving SCBC merged data to disk...')


def extractPartitions(inputData, dataDescriptor='input', sortDir=-1, debug=False, debug_numPartitions=256):

    print('Extracting unique feature combinations from {0} data...'.format(dataDescriptor))
    data = inputData.dropna(axis=0, how='all')
    dataSubset_null = data.notnull()
    dataSubset_null_tonumpy = dataSubset_null.to_numpy()
    featureCombinations, mapping_featureCombinations, uniqueCounts = \
        np.unique(dataSubset_null_tonumpy, return_index=True, return_counts=True, axis=0)
    numPartitions = len(featureCombinations)
    uniqueCounts_sorted = np.argsort(sortDir*uniqueCounts)
    featureCombinations_sorted = featureCombinations[uniqueCounts_sorted, :]

    if debug:
        orderedIDs = uniqueCounts[uniqueCounts_sorted]
        plt.scatter(np.arange(len(orderedIDs)), orderedIDs)
        plt.title("Partition vs # Observations")
        plt.xlabel("Partition (Ordered)")
        plt.ylabel("Observation Count")

        print('     Printing counts, a subset, of unique partitions in data...')
        sorted_uniqueCounts = uniqueCounts[uniqueCounts_sorted]
        print(sorted_uniqueCounts[:10])

        print('     The top {0} partitions cover {1} of the interactions pairs...'.format(
            debug_numPartitions, sorted_uniqueCounts[:debug_numPartitions].sum()/sorted_uniqueCounts.sum()))

        indices = np.arange(len(sorted_uniqueCounts))
        sorted_uniqueCounts_cumsumRatio = np.cumsum(sorted_uniqueCounts)/sorted_uniqueCounts.sum()

        fig1, (ax1, ax2) = plt.subplots(2)
        fig1.suptitle('Ordered Feature Combinations vs Counts and Cumulative Ratio')
        ax1.plot(indices, sorted_uniqueCounts, label=dataDescriptor)
        plt.xlabel('Ordered Feature Combinations')
        ax1.set_ylabel('Counts')
        ax1.legend()

        ax2.plot(indices, sorted_uniqueCounts_cumsumRatio, label=dataDescriptor)
        ax2.set_ylabel('Cumulative Ratio')
        ax2.legend()

    print('     The number of unique feature combinations are {0} and the shape of unique counts is {1}...'.format(
        numPartitions, uniqueCounts.shape))
    print('     Completed extracting unique feature combinations from {0} data...'.format(dataDescriptor))
    print('________________________________________________')

    return featureCombinations, mapping_featureCombinations, uniqueCounts, uniqueCounts_sorted, numPartitions, \
        featureCombinations_sorted


def loadResults_precRecall(path, name):

    precRecall_results = pd.read_csv(path, sep='\t', header=None)
    prec = precRecall_results.iloc[:, 0]
    recall = precRecall_results.iloc[:, 1]
    precRecall_pair = [prec, recall, name]

    return precRecall_pair
