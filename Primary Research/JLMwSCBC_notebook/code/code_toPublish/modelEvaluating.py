import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
from sklearn.metrics import auc
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve


def evalPartitions(trainingPartitions_sorted, trainedModels_properties,
                   testData, testLabels, workdir, expName, suffix='', storeModels=False):
    print('Evaluating models by training data partitions on testing data')
    print('____________________________________________________________')
    features = testData.columns.tolist()
    testData_notnull = testData.notnull()
    addProperties = []
    columnNames = ['training partition ID', 'testing accuracy', 'average precision', 'Pr-R AUC', 'est model accuracy',
                   'num features', 'num test-data observations', 'num train-data observations', 'features']
    os.makedirs(workdir + 'modelPerformance/modelsPerformance_{0}'.format(expName), exist_ok=True)
    if storeModels:
        modelID = 0
        models = {}
        for prtn in np.arange(len(trainingPartitions_sorted)):
            print('model {0} of total {1} models loaded...'.format(prtn, len(trainingPartitions_sorted)))
            filename = workdir + 'models/models_{0}/model_{0}_Partition_{1}'.format(expName, prtn)
            models[modelID] = pickle.load(open(filename + '.pkl', 'rb'))
            setattr(models[modelID], 'verbose', 0)
            modelID += 1

    candidatePartitions = np.empty((len(testData), len(trainingPartitions_sorted)), dtype=np.float16)
    candidatePartitions[:] = np.nan
    candidatePartitions_preds = np.empty((len(testData), len(trainingPartitions_sorted)), dtype=np.float16)
    candidatePartitions_preds[:] = np.nan
    candidatePartitions_predsProbs = np.empty((len(testData), len(trainingPartitions_sorted)), dtype=np.float16)
    candidatePartitions_predsProbs[:] = np.nan
    for partition in np.arange(len(trainingPartitions_sorted)):

        progressRatio = partition / len(trainingPartitions_sorted)
        if np.floor(int(progressRatio*100)) % 5 == 0:
            print('progress: {0}'.format(progressRatio*100))
        selectFeatures = list(np.argwhere(trainingPartitions_sorted[partition, :])[:, 0])
        # when num of training features exceeds # of available testing features
        selectFeatures = [i for i in selectFeatures if i < len(features)]
        selectFeatures_names = [features[i] for i in selectFeatures]
        allRows_selectFeatures = testData_notnull.iloc[:, selectFeatures].all(axis=1)
        rowIndices = allRows_selectFeatures.index[allRows_selectFeatures].tolist()

        if len(rowIndices) != 0:
            candidatePartitions[rowIndices, partition] = 1
            interactionsSubset = testData.iloc[rowIndices, selectFeatures]
            labelsSubset = testLabels.iloc[rowIndices, 0]

            if not storeModels:
                partitionString = '{0:4d}'.format(partition)
                filenamePatt = \
                    glob.glob(workdir + 'models/models_{0}/model_{0}_Partition_{1}_*.pkl'.format(expName, partition))[0]
                model = pickle.load(open(filenamePatt, 'rb'))
                setattr(model, 'verbose', 0)
                candidatePartitions_preds[rowIndices, partition] = model.predict(interactionsSubset)
                try:
                    candidatePartitions_predsProbs[rowIndices, partition] = \
                        model.predict_proba(interactionsSubset)[:, 1]
                except:
                    candidatePartitions_predsProbs[rowIndices, partition] = model.predict_proba(interactionsSubset)[0]
                score = model.score(interactionsSubset, labelsSubset)

            else:
                candidatePartitions_preds[rowIndices, partition] = models[partition].predict(interactionsSubset)
                try:
                    candidatePartitions_predsProbs[rowIndices, partition] = \
                        models[partition].predict_proba(interactionsSubset)[:, 1]
                except:
                    candidatePartitions_predsProbs[rowIndices, partition] = models[partition].predict_proba(interactionsSubset)[0]
                score = models[partition].score(interactionsSubset, labelsSubset)

            if len(np.unique(labelsSubset)) > 1:
                precisionScore = average_precision_score(labelsSubset,
                                                         candidatePartitions_predsProbs[rowIndices, partition],
                                                         average='weighted')
                noSkill_probs = [0 for _ in range(len(labelsSubset))]
                noSkill_prec, noSkill_recall, _ = precision_recall_curve(labelsSubset, noSkill_probs)
                modelPartition_prec, modelPartition_recall, _ = \
                    precision_recall_curve(labelsSubset,
                                           candidatePartitions_predsProbs[rowIndices, partition])
                modelPartition_prrauc = auc(modelPartition_recall, modelPartition_prec)

            else:
                precisionScore = np.nan
                modelPartition_prrauc = np.nan

        else:
            score = np.nan
            precisionScore = np.nan
            modelPartition_prrauc = np.nan

        properties = [partition, score, precisionScore, modelPartition_prrauc,
                      trainedModels_properties['est model accuracy'].loc[partition], len(selectFeatures),
                      len(rowIndices), trainedModels_properties['num train-data observations'].loc[partition],
                      selectFeatures_names]
        addProperties.append(pd.DataFrame([properties], columns=columnNames))

    modelsEval_properties = pd.concat(addProperties)
    modelsEval_properties.set_index('training partition ID', inplace=True)
    modelsEval_properties.to_csv(workdir + 'modelPerformance/modelsPerformance_{0}/'.format(expName) +
                                 'Partition-Classifier_Testing_Performance{0}.csv'.format(suffix))

    np.savetxt(workdir + 'modelPerformance/modelsPerformance_{0}/'.format(expName) +
               'candidatePartitions{0}.csv'.format(suffix),
               np.array(candidatePartitions), delimiter=',')
    np.savetxt(workdir + 'modelPerformance/modelsPerformance_{0}/'.format(expName) +
               'candidatePartitions_preds{0}.csv'.format(suffix),
               np.array(candidatePartitions_preds), delimiter=',')
    np.savetxt(workdir + 'modelPerformance/modelsPerformance_{0}/'.format(expName) +
               'candidatePartitions_predsProbs{0}.csv'.format(suffix),
               np.array(candidatePartitions_predsProbs), delimiter=',')

    print('Completed training models per partition of the training data')
    print('____________________________________________________________')

    return candidatePartitions, candidatePartitions_preds, candidatePartitions_predsProbs, modelsEval_properties


def evalTotal(candidateModels, candidateModel_predsProbs, classifiersProps, testLabels, workdir, expName,
              evaluatingModels_limitingThreshold=0.9, evaluatingModels_max=8, evaluatingModels_min=5, suffix=''):
    print('Evaluating total {0}'.format(expName))
    print('____________________________________________________________')
    probsPos_acc = []
    probsPos_weighted = []

    os.makedirs(workdir + 'modelPerformance/modelsPerformance_{0}'.format(expName), exist_ok=True)
    progress = 0
    for interactionsPair in np.arange(len(testLabels)):

        progressRatio = progress / len(testLabels)
        print('progress: {0}'.format(progressRatio))
        memberPartitions = list(np.argwhere(np.isfinite(candidateModels[interactionsPair, :]))[:, 0])
        classifiersProps_relevant = classifiersProps.loc[memberPartitions]
        if len(classifiersProps_relevant) > 0:
            classifierChoice_accuracyInd = np.nanargmax(classifiersProps_relevant.loc[:, 'est model accuracy'])
            classifierChoice_acc = classifiersProps_relevant.index[classifierChoice_accuracyInd]
            classifierProb_accuracy = candidateModel_predsProbs[interactionsPair, int(classifierChoice_acc)]

            allAccs = classifiersProps_relevant[
                          classifiersProps_relevant['est model accuracy'] >
                          evaluatingModels_limitingThreshold].loc[:, 'est model accuracy'].to_numpy()
            if (len(allAccs) > evaluatingModels_min) and (len(allAccs) > evaluatingModels_max):

                classifierAccs = allAccs[
                    allAccs.argsort()[-evaluatingModels_max:]
                ]
                memberPartitions_select = [memberPartitions[i] for i in allAccs.argsort()[-evaluatingModels_max:]]
                classifiersAll_predProbs = candidateModel_predsProbs[interactionsPair, memberPartitions_select]

            else:
                classifierAccs = classifiersProps_relevant.loc[:, 'est model accuracy'].to_numpy()
                classifiersAll_predProbs = candidateModel_predsProbs[interactionsPair, memberPartitions]

            classifiersAll_weightedScore = np.sum(
                np.multiply(
                    classifierAccs,
                    classifiersAll_predProbs)
            ) / np.sum(classifierAccs)
            probsPos_acc.append(classifierProb_accuracy)
            probsPos_weighted.append(classifiersAll_weightedScore)
            progress += 1

    modelsApplicable = np.nansum(candidateModels, axis=1)
    count0_modelsApplicable = len(modelsApplicable) - np.count_nonzero(modelsApplicable)
    print(
        'The # of interaction pairs for which there are 0 applicable partitions is {0}'.format(
            count0_modelsApplicable))
    interactionPairs_noModel = np.argwhere(modelsApplicable < 1)
    print('The number of members to interactions {0}'.format(len(interactionPairs_noModel)))
    testLabels_adj = \
        [testLabels.iloc[i, 0] for i in np.arange(len(testLabels)) if i not in interactionPairs_noModel]
    print('The length of testLabels without the pairs for which there is no model is {0}'.format(
        len(testLabels) - len(testLabels_adj)))

    np.savetxt(workdir + 'modelPerformance/modelsPerformance_{0}/'.format(expName) +
               'probsPos_acc{0}.csv'.format(suffix), np.array(probsPos_acc), delimiter=',')
    np.savetxt(workdir + 'modelPerformance/modelsPerformance_{0}/'.format(expName) +
               'probsPos_weighted{0}.csv'.format(suffix), np.array(probsPos_weighted), delimiter=',')
    np.savetxt(workdir + 'modelPerformance/modelsPerformance_{0}/'.format(expName) +
               'testLabels_adj{0}.csv'.format(suffix), np.array(testLabels_adj), delimiter=',')
    print('Completed evaluating total {0}'.format(expName))
    print('____________________________________________________________')

    return probsPos_acc, probsPos_weighted, testLabels_adj


def evalSummary(probsPos_accuracy, probsPos_weighted, testLabels_adj, workDir, expName, suffix=''):
    print('Summarizing {0} analysis...'.format(expName))
    print('____________________________________________________________')
    precision_acc, recall_acc, _ = precision_recall_curve(testLabels_adj, probsPos_accuracy)
    precision_weighted, recall_weighted, _ = precision_recall_curve(testLabels_adj, probsPos_weighted)
    try:
        modelPartition_prrauc_acc = auc(recall_acc, precision_acc)
        modelPartition_prrauc_weighted = auc(recall_weighted, precision_weighted)
    except:
        modelPartition_prrauc_acc = np.nan
        modelPartition_prrauc_weighted = np.nan
    print('The unit-accuracy model Pr-R AUC is {0}.'.format(modelPartition_prrauc_acc))
    print('The combined-weighted scores model Pr-R AUC is {0}.'.format(modelPartition_prrauc_weighted))

    np.savetxt(workDir + 'modelPerformance/modelsPerformance_{0}/'.format(expName) +
               'precision_acc{0}.csv'.format(suffix), precision_acc, delimiter=',')
    np.savetxt(workDir + 'modelPerformance/modelsPerformance_{0}/'.format(expName) +
               'recall_acc{0}.csv'.format(suffix), recall_acc, delimiter=',')
    np.savetxt(workDir + 'modelPerformance/modelsPerformance_{0}/'.format(expName) +
               'precision_weighted{0}.csv'.format(suffix), precision_weighted, delimiter=',')
    np.savetxt(workDir + 'modelPerformance/modelsPerformance_{0}/'.format(expName) +
               'recall_weighted{0}.csv'.format(suffix), recall_weighted, delimiter=',')

    results_acc = pd.DataFrame(np.vstack((recall_acc, precision_acc)).T)
    results_weighted = pd.DataFrame(np.vstack((recall_weighted, precision_weighted)).T)
    np.savetxt(workDir + 'modelPerformance/modelsPerformance_{0}/'.format(expName) +
               'unit-accuracy_results{0}.csv'.format(suffix), results_acc, delimiter='\t')
    np.savetxt(workDir + 'modelPerformance/modelsPerformance_{0}/'.format(expName) +
               'combined-weighted_results{0}.csv'.format(suffix), results_weighted, delimiter='\t')

    print('Completed summarizing {0} analysis'.format(expName))
    print('____________________________________________________________')

    return precision_acc, recall_acc, precision_weighted, recall_weighted


def disjointAssembly(partitionsSorted, trainedModels_properties,
                     data, labels,
                     workDir, expName,
                     numChunks, evaluatingChunk_incr, evaluatingStart_interval,
                     evaluatingModels_limitingThreshold=0.9, evaluatingModels_max=8,
                     evaluatingModels_min=5, storeModels=False):
    chunkCount = 0
    startChunk = evaluatingChunk_incr * evaluatingStart_interval
    endChunk = evaluatingChunk_incr * (evaluatingStart_interval + 1)
    if endChunk == evaluatingChunk_incr * np.ceil(numChunks / evaluatingStart_interval):
        endChunk += 1
    for dataChunk, labelChunk in zip(data, labels):

        if (chunkCount >= startChunk) and (chunkCount < endChunk):
            print('printing chunkCount {0}'.format(chunkCount))

            candidatePartitions, candidatePartitions_preds, candidatePartitions_predsProbs, \
            modelsEval_properties = evalPartitions(partitionsSorted, trainedModels_properties,
                                                   dataChunk.reset_index(drop=True), labelChunk.reset_index(drop=True),
                                                   workDir, expName, suffix='_{:0>5d}'.format(chunkCount),
                                                   storeModels=storeModels)

            evalTotal(candidatePartitions, candidatePartitions_predsProbs, modelsEval_properties,
                      labelChunk.reset_index(drop=True), workDir, expName, suffix='_{:0>5d}'.format(chunkCount),
                      evaluatingModels_limitingThreshold=evaluatingModels_limitingThreshold,
                      evaluatingModels_max=evaluatingModels_max, evaluatingModels_min=evaluatingModels_min)

        chunkCount += 1


def loadAssembly(workDir, expName):
    modelPerformance_dirPath = workDir + 'modelPerformance/modelsPerformance_{0}/'.format(expName)
    probsPos_accFilenames = glob.glob(modelPerformance_dirPath + 'probsPos_acc.csv')
    probsPos_weightedFilenames = glob.glob(modelPerformance_dirPath + 'probsPos_weighted.csv')
    testLabels_adjFilenames = glob.glob(modelPerformance_dirPath + 'testLabels_adj.csv')
    probsPos_accFilenames = sorted(probsPos_accFilenames)
    probsPos_weightedFilenames = sorted(probsPos_weightedFilenames)
    testLabels_adjFilenames = sorted(testLabels_adjFilenames)
    assembliesLengths = [len(probsPos_accFilenames),
                         len(probsPos_weightedFilenames),
                         len(testLabels_adjFilenames)]
    assembliesLengths_unique = np.unique(assembliesLengths)
    if len(assembliesLengths_unique) != 1:
        print('The # of files for all 3 key assemblies is not equal. Check disjoint assembly completion...')

    else:
        probsPos_acc = np.empty((0,))
        probsPos_weighted = np.empty((0,))
        testLabels_adj = np.empty((0,))
        for file1, file2, file3 in zip(probsPos_accFilenames, probsPos_weightedFilenames, testLabels_adjFilenames):
            probsPos_acc = np.concatenate((probsPos_acc, np.loadtxt(file1)))
            probsPos_weighted = np.concatenate((probsPos_weighted, np.loadtxt(file2)))
            testLabels_adj = np.concatenate((testLabels_adj, np.loadtxt(file3)))

    return probsPos_acc, probsPos_weighted, testLabels_adj


def assembleCV_results(workDir, modelPerformances_dirs, numFolds, expName):
    allLabels = {fold: np.loadtxt(predictionsDir, delimiter=',') for fold, predictionsDir in
                 zip(np.arange(numFolds), glob.glob(modelPerformances_dirs + 'testLabels_adj.csv'))}
    probsPos_unitAccuracy_byFolds = {fold: np.loadtxt(predictionsDir, delimiter=',') for fold, predictionsDir in
                                     zip(np.arange(numFolds), glob.glob(modelPerformances_dirs + 'probsPos_acc.csv'))}
    probsPos_combinedWeighted_byFolds = {fold: np.loadtxt(predictionsDir, delimiter=',') for fold, predictionsDir in
                                         zip(np.arange(numFolds), glob.glob(modelPerformances_dirs + 'probsPos_weighted.csv'))}

    labels_acrossFolds = []
    unitAccuracy_probsPos_kCV = []
    combinedWeighted_probsPos_kCV = []
    for unitAccuracy_probsPos, combinedWeighted_probsPos, labels_byFolds in zip(probsPos_unitAccuracy_byFolds.values(),
                                                                                probsPos_combinedWeighted_byFolds.values(),
                                                                                allLabels.values()):
        labels_acrossFolds.append(labels_byFolds)
        unitAccuracy_probsPos_kCV.append(unitAccuracy_probsPos)
        combinedWeighted_probsPos_kCV.append(combinedWeighted_probsPos)

    labels_acrossFolds = np.concatenate(labels_acrossFolds)
    unitAccuracy_probsPos_kCV = np.concatenate(unitAccuracy_probsPos_kCV)
    combinedWeighted_probsPos_kCV = np.concatenate(combinedWeighted_probsPos_kCV)
    precision_unitAccuracy_kCV, recall_unitAccuracy_kcv, _ = precision_recall_curve(labels_acrossFolds,
                                                                                    unitAccuracy_probsPos_kCV)
    precision_combinedWeighted_kCV, recall_combinedWeighted_kCV, _ = precision_recall_curve(labels_acrossFolds,
                                                                                            combinedWeighted_probsPos_kCV)
    cvPerformance_path = \
        workDir + 'modelPerformance/modelsPerformance_{0}-5CV/'.format(expName)
    os.makedirs(cvPerformance_path, exist_ok=True)
    np.savetxt(cvPerformance_path + 'precision_acc.csv', precision_unitAccuracy_kCV, delimiter=',')
    np.savetxt(cvPerformance_path + 'recall_acc.csv', recall_unitAccuracy_kcv, delimiter=',')
    np.savetxt(cvPerformance_path + 'precision_weighted.csv', precision_combinedWeighted_kCV, delimiter=',')
    np.savetxt(cvPerformance_path + 'recall_weighted.csv', recall_combinedWeighted_kCV, delimiter=',')

    results_acc = pd.DataFrame(np.vstack((recall_unitAccuracy_kcv, precision_unitAccuracy_kCV)).T)
    results_weighted = pd.DataFrame(np.vstack((recall_combinedWeighted_kCV, precision_combinedWeighted_kCV)).T)
    np.savetxt(cvPerformance_path + 'unit-accuracy_results.csv', results_acc, delimiter='\t')
    np.savetxt(cvPerformance_path + 'combined-weighted_results.csv', results_weighted, delimiter='\t')


def prrCalculator(testLabels, preds):
    testLabels = pd.Series(testLabels)
    preds = pd.Series(preds, index=testLabels.index)

    TP = 0
    FP = 0
    FN = 0

    for i in testLabels.index:
        if testLabels[i] == preds[i] == 1:
            TP += 1
        if preds[i] == 1 and testLabels[i] != preds[i]:
            FP += 1
        if preds[i] == 0 and testLabels[i] != preds[i]:
            FN += 1

    try:
        prec = TP / (TP + FP)
    except:
        prec = 1

    try:
        rec = TP / (TP + FN)
    except:
        rec = 1

    return prec, rec


def generateCV_prr(workDir, modelPerformances_dirs, numFolds, expName):
    numFolds = np.arange(numFolds)
    allLabels = {fold: np.loadtxt(predictionsDir, delimiter=',') for fold, predictionsDir in
                 zip(numFolds, glob.glob(modelPerformances_dirs + 'testLabels_adj.csv'))}
    probsPos_unitAccuracy_byFolds = {fold: np.loadtxt(predictionsDir, delimiter=',') for fold, predictionsDir in
                                     zip(numFolds, glob.glob(modelPerformances_dirs + 'probsPos_acc.csv'))}
    probsPos_combinedWeighted_byFolds = {fold: np.loadtxt(predictionsDir, delimiter=',') for fold, predictionsDir in
                                         zip(numFolds, glob.glob(modelPerformances_dirs + 'probsPos_weighted.csv'))}

    labels_acrossFolds = []
    unitAccuracy_probsPos_acrossFolds = []
    combinedWeighted_probsPos_acrossFolds = []
    for unitAccuracy_probsPos, combinedWeighted_probsPos, labels_byFolds in zip(probsPos_unitAccuracy_byFolds.values(),
                                                                                probsPos_combinedWeighted_byFolds.values(),
                                                                                allLabels.values()):
        labels_acrossFolds.append(labels_byFolds)
        unitAccuracy_probsPos_acrossFolds.append(unitAccuracy_probsPos)
        combinedWeighted_probsPos_acrossFolds.append(combinedWeighted_probsPos)

    labels_acrossFolds = np.concatenate(labels_acrossFolds)
    unitAccuracy_probsPos_acrossFolds = np.concatenate(unitAccuracy_probsPos_acrossFolds)
    combinedWeighted_probsPos_acrossFolds = np.concatenate(combinedWeighted_probsPos_acrossFolds)

    precScores_unitAcc = []
    recScores_unitAcc = []
    precScores_combinedWeighted = []
    recScores_combinedWeighted = []
    prob_thresholds = np.linspace(0, 1, num=100)
    for p in prob_thresholds:

        testPreds_unitAcc = []
        testPreds_combinedWeighted = []

        for prob in unitAccuracy_probsPos_acrossFolds:
            if prob > p:
                testPreds_unitAcc.append(1)
            else:
                testPreds_unitAcc.append(0)

        for prob in combinedWeighted_probsPos_acrossFolds:
            if prob > p:
                testPreds_combinedWeighted.append(1)
            else:
                testPreds_combinedWeighted.append(0)

        prec, rec = prrCalculator(labels_acrossFolds, testPreds_unitAcc)
        precScores_unitAcc.append(prec)
        recScores_unitAcc.append(rec)

        prec, rec = prrCalculator(labels_acrossFolds, testPreds_combinedWeighted)
        precScores_combinedWeighted.append(prec)
        recScores_combinedWeighted.append(rec)

    try:
        cvPerformance_path = \
            workDir + 'modelPerformance/modelsPerformance_ppiPrediction_performanceAnalysis-5CV_{0}/'.format(expName)
        os.makedirs(cvPerformance_path)
        np.savetxt(cvPerformance_path + 'precision_acc.csv', precScores_unitAcc, delimiter=',')
        np.savetxt(cvPerformance_path + 'recall_acc.csv', recScores_unitAcc, delimiter=',')
        np.savetxt(cvPerformance_path + 'precision_weighted.csv', precScores_combinedWeighted, delimiter=',')
        np.savetxt(cvPerformance_path + 'recall_weighted.csv', recScores_combinedWeighted, delimiter=',')

        results_acc = pd.DataFrame(np.vstack((recScores_unitAcc, precScores_unitAcc)).T)
        results_weighted = pd.DataFrame(np.vstack((recScores_combinedWeighted, precScores_combinedWeighted)).T)
        np.savetxt(cvPerformance_path + 'unit-accuracy_results.csv', results_acc, delimiter='\t')
        np.savetxt(cvPerformance_path + 'combined-weighted_results.csv', results_weighted, delimiter='\t')

    except:
        print('ERROR: files not saved!!!!')

    return precScores_unitAcc, recScores_unitAcc, precScores_combinedWeighted, recScores_combinedWeighted


def plot_kCV_prr(workDir, modelPerformances_dirs, numFolds, expName):

    #f = plt.figure(figsize=(10, 5))

    folds = np.arange(numFolds)
    allLabels = {fold: np.loadtxt(predictionsDir, delimiter=',') for fold, predictionsDir in
                 zip(folds, glob.glob(modelPerformances_dirs + 'testLabels_adj.csv'))}
    probsPos_unitAccuracy_byFolds = {fold: np.loadtxt(predictionsDir, delimiter=',') for fold, predictionsDir in
                                     zip(folds, glob.glob(modelPerformances_dirs + 'probsPos_acc.csv'))}
    probsPos_combinedWeighted_byFolds = {fold: np.loadtxt(predictionsDir, delimiter=',') for fold, predictionsDir in
                                         zip(folds, glob.glob(modelPerformances_dirs + 'probsPos_weighted.csv'))}

    labels_acrossFolds = []
    unitAccuracy_probsPos_acrossFolds = []
    combinedWeighted_probsPos_acrossFolds = []
    for unitAccuracy_probsPos, combinedWeighted_probsPos, labels_byFolds, fold in zip(
            probsPos_unitAccuracy_byFolds.values(),
            probsPos_combinedWeighted_byFolds.values(),
            allLabels.values(), folds):

        #precAcc, recallAcc, _ = precision_recall_curve(labels_byFolds, unitAccuracy_probsPos)
        #labAcc = 'Unit-Acc Fold {0} AUC={1:.4f}'.format(fold, auc(recallAcc, precAcc))
        #plt.step(recallAcc, precAcc, label=labAcc)

        #precWeighted, recallWeighted, _ = precision_recall_curve(labels_byFolds, combinedWeighted_probsPos)
        #labWeighted = 'Weight Fold {0} AUC={1:.4f}'.format(fold, auc(recallWeighted, precWeighted))
        #plt.step(recallWeighted, precWeighted, label=labWeighted)

        labels_acrossFolds.append(labels_byFolds)
        unitAccuracy_probsPos_acrossFolds.append(unitAccuracy_probsPos)
        combinedWeighted_probsPos_acrossFolds.append(combinedWeighted_probsPos)

    labels_acrossFolds = np.concatenate(labels_acrossFolds)
    unitAccuracy_probsPos_acrossFolds = np.concatenate(unitAccuracy_probsPos_acrossFolds)
    combinedWeighted_probsPos_acrossFolds = np.concatenate(combinedWeighted_probsPos_acrossFolds)

    precAcc, recallAcc, _ = precision_recall_curve(labels_acrossFolds, unitAccuracy_probsPos_acrossFolds)
    #labAcc = '{0}CV AUC={1:.4f}'.format(numFolds, auc(recallAcc, precAcc))
    #plt.step(recallAcc, precAcc, label=labAcc, lw=2)

    precWeighted, recallWeighted, _ = precision_recall_curve(labels_acrossFolds, combinedWeighted_probsPos_acrossFolds)
    #labWeighted = '{0}CV AUC={1:.4f}'.format(numFolds, auc(recallWeighted, precWeighted))
    #plt.step(recallWeighted, precWeighted, label=labWeighted, lw=2)

    cvPerformance_path = \
        workDir + 'modelPerformance/modelsPerformance_{0}-5CV/'.format(expName)
    os.makedirs(cvPerformance_path)

    np.savetxt(cvPerformance_path + 'precision_acc.csv', precAcc, delimiter=',')
    np.savetxt(cvPerformance_path + 'recall_acc.csv', recallAcc, delimiter=',')
    np.savetxt(cvPerformance_path + 'precision_weighted.csv', precWeighted, delimiter=',')
    np.savetxt(cvPerformance_path + 'recall_weighted.csv', recallWeighted, delimiter=',')

    results_acc = pd.DataFrame(np.vstack((recallAcc, precAcc)).T)
    results_weighted = pd.DataFrame(np.vstack((recallWeighted, precWeighted)).T)
    np.savetxt(cvPerformance_path + 'unit-accuracy_results.csv', results_acc, delimiter='\t')
    np.savetxt(cvPerformance_path + 'combined-weighted_results.csv', results_weighted, delimiter='\t')

    #plt.xlabel('Recall')
    #plt.ylabel('Precision')
    #plt.legend(bbox_to_anchor=(1, 1), loc='upper left', ncol=1, fontsize='small')

    #f.tight_layout()
    #f.savefig(workDir + '{0}.png'.format(expName))


def plotAnalysis(workdir, expName, precRecall_pairSets):
    filename = '{0}_Pr-R_Performance'.format(expName)
    plt.clf()
    plt.figure(figsize=(12.8, 9.6))
    plt.title('{0} Pr-R Performance'.format(expName))
    for precRecall_pair in precRecall_pairSets:
        precRecall_auc = auc(precRecall_pair[0], precRecall_pair[1])
        plt.plot(precRecall_pair[0], precRecall_pair[1],
                 label=precRecall_pair[2] + ', AUC: {0: .3f}'.format(precRecall_auc))
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(bbox_to_anchor=(1, 1), loc='upper right', ncol=1)
    filename_wPath = \
      workdir + 'modelPerformance/modelsPerformance_{0}/'.format(expName) + filename + '.svg'
    plt.savefig(filename_wPath, bbox_inches='tight', dpi=600)

def prrPlot_overSTRING_varying(DF, variedTypes, confidenceMetric, thresholdRange, outputDir):

    for varyType in variedTypes:
        plt.clf()
        for t in thresholdRange:
            gtSTRING = DF.copy()
            gtSTRING.loc[gtSTRING[confidenceMetric] >= t, 'labels'] = 1
            gtSTRING.loc[gtSTRING[confidenceMetric] < t, 'labels'] = 0
            precision, recall, _ = precision_recall_curve(
                gtSTRING.labels.to_numpy(),
                gtSTRING.loc[:, 'probsPos_' + varyType].to_numpy())
            plt.plot(recall, precision, label=str(t) + ', AUC: {0: .3f}'.format(auc(recall, precision)))
        plt.title('{0} Pr-R Performance over STRING'.format(varyType))
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.legend(bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
        plt.savefig(outputDir + '{0}_Pr-R_overSTRING'.format(varyType) + '.png', bbox_inches='tight')