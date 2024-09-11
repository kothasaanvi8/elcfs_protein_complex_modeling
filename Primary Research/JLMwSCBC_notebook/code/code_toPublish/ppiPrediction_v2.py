import dataProcessing
import modelTraining
import modelEvaluating
import numpy as np
import os
import pandas as pd
import pickle


def run(workDir,
        trainingData_path, trainingLabels_path, testData_path, testLabels_path,
        expName='test', hasHeadings=False, savedPartitions=False, reset=False,
        evaluatingModels_limitingThreshold=0.9, evaluatingModels_max=8, evaluatingModels_min=5, suffix='',
        trainingData_chunkSize=None, testData_chunkSize=None,
        modelChunk_incr=None, modelStart_interval=None, evaluatingStart_interval=None,
        evaluatingChunk_size=None, evaluatingChunk_incr=None, storeModels=False,
        modelingOnly=False, parallelProcessing=False,
        evaluatingOnly=False, intermediatesGenerating=False
        ):

    trainingData, trainingLabels = dataProcessing.acquire(trainingData_path, trainingLabels_path, headings=hasHeadings,
                                                          chunkSize=trainingData_chunkSize)
    testData, testLabels = dataProcessing.acquire(testData_path, testLabels_path, headings=hasHeadings,
                                                  chunkSize=testData_chunkSize)

    #need to be sure chunk doesn't break this
    if savedPartitions:
        partitionsSorted = pickle.load(open(savedPartitions, 'rb'))
        print('the length of partitions sorted is {0}'.format(len(partitionsSorted)))
    else:
        #check if partitions have been saved
        #savedPartitions must be used if partitions are going to be extracted!!!
        partitionsFile = workDir + 'models/models_{0}/savedPartitions.pkl'.format(expName)
        if os.path.isfile(partitionsFile) and not reset:
            partitionsSorted = pickle.load(open(partitionsFile, 'rb'))
        else:
            _, _, _, _, _, trainingPartitions_sorted = dataProcessing.extractPartitions(trainingData, 'training')
            _, _, _, _, _, testingPartitions_sorted = dataProcessing.extractPartitions(testData, 'testing')
            partitionsConcatenated = np.concatenate((trainingPartitions_sorted, testingPartitions_sorted), axis=0)
            partitionsSorted = np.unique(partitionsConcatenated, axis=0)
            os.makedirs(workDir + 'models/models_{0}'.format(expName), exist_ok=True)
            pickle.dump(partitionsSorted, open(partitionsFile, 'wb'))

    modelsProperties = modelTraining.build(partitionsSorted, trainingData, trainingLabels, workDir, expName,
                                           parallelProcessing=parallelProcessing, modelChunk_incr=modelChunk_incr,
                                           modelStart_interval=modelStart_interval)

    if not modelingOnly:
        candidatePartitions, candidatePartitions_preds, candidatePartitions_predsProbs, modelsEval_properties = \
            modelEvaluating.evalPartitions(partitionsSorted, modelsProperties, testData, testLabels, workDir, expName,
                                           storeModels=storeModels)

        probsPos_acc, probsPos_weighted, testLabels_adj = \
            modelEvaluating.evalTotal(candidatePartitions, candidatePartitions_predsProbs,
                                      modelsEval_properties, testLabels, workDir, expName,
                                      evaluatingModels_limitingThreshold, evaluatingModels_max, evaluatingModels_min,
                                      suffix)

        precAcc, recAcc, precWeighted, recWeighted = modelEvaluating.evalSummary(probsPos_acc, probsPos_weighted,
                                                                                 testLabels_adj, workDir, expName)

        return precAcc, recAcc, precWeighted, recWeighted

    if evaluatingOnly:

        if parallelProcessing:

            modelsProperties = modelTraining.build(partitionsSorted, trainingData, trainingLabels, workDir, expName,
                                                   parallelProcessing=parallelProcessing,
                                                   modelChunk_incr=modelChunk_incr,
                                                   modelStart_interval=modelStart_interval)

            numChunks = np.ceil(
                len(pd.read_csv(trainingData_path, sep='\t', header=None, usecols=[0])) / evaluatingChunk_size)
            modelEvaluating.disjointAssembly(partitionsSorted, modelsProperties, testData, testLabels,
                                             workDir, expName, numChunks,
                                             evaluatingChunk_incr, evaluatingStart_interval,
                                             evaluatingModels_limitingThreshold, evaluatingModels_max,
                                             evaluatingModels_min, storeModels=storeModels)
        elif intermediatesGenerating:
            probsPos_acc, probsPos_weighted, testLabels_adj = modelEvaluating.loadAssembly(workDir, expName)
            precAcc, recAcc, precWeighted, recWeighted = modelEvaluating.evalSummary(probsPos_acc,
                                                                                     probsPos_weighted,
                                                                                     testLabels_adj, workDir,
                                                                                     expName)

        else:
            candidatePartitions, candidatePartitions_preds, candidatePartitions_predsProbs, modelsEval_properties = \
                modelEvaluating.evalPartitions(partitionsSorted, modelsProperties, trainingData, trainingLabels,
                                               workDir, expName, storeModels=storeModels)

            probsPos_acc, probsPos_weighted, testLabels_adj = \
                modelEvaluating.evalTotal(candidatePartitions, candidatePartitions_predsProbs,
                                          modelsEval_properties, testLabels, workDir, expName,
                                          evaluatingModels_max, evaluatingModels_min, suffix)

            precAcc, recAcc, precWeighted, recWeighted = modelEvaluating.evalSummary(probsPos_acc, probsPos_weighted,
                                                                                     testLabels_adj, workDir, expName)

        return precAcc, recAcc, precWeighted, recWeighted
