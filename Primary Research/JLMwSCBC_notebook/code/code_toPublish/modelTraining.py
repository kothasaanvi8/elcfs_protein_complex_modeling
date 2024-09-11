import datetime
import glob
import numpy as np
import os
import pandas as pd
import pickle
from sklearn.ensemble import RandomForestClassifier

global dateSignature
dateSignature = datetime.datetime.now().strftime("%Y-%m-%d")


def build(trainingPartitions, trainingData, trainingLabels,
          workdir, expName,
          otherEstimator=None,
          reset=False,
          parallelProcessing=False, reverse=False,
          modelChunk_incr=30, modelStart_interval=0):

    print('Training models per partition of the training data')
    print('____________________________________________________________')

    columnNames = ['partition ID', 'est model accuracy', 'class imbalance weight', 'num features',
                   'num train-data observations', 'features']
    trainedModels_properties = pd.DataFrame(columns=columnNames)
    os.makedirs(workdir + 'models/models_{0}/'.format(expName), exist_ok=True)
    if parallelProcessing:
        startChunk = modelChunk_incr * modelStart_interval
        endChunk = modelChunk_incr * (modelStart_interval + 1)
        if endChunk >= len(trainingPartitions):
            endChunk = len(trainingPartitions)
    else:
        startChunk = 0
        endChunk = len(trainingPartitions)

    reverse = -1 if reverse else 1
    for partition in np.arange(startChunk, endChunk)[::reverse]:

        print('Progress: {0}/{1} '.format(partition, len(trainingPartitions)))
        selectFeatures = list(np.argwhere(trainingPartitions[partition, :])[:, 0])
        if hasattr(trainingData, 'shape'):
            features = trainingData.columns.tolist()
            selectFeatures_names = [features[i] for i in selectFeatures]
            trainingData_notnull = trainingData.notnull()
            allRows_selectFeatures = trainingData_notnull.iloc[:, selectFeatures].all(axis=1)
            rowIndices = allRows_selectFeatures.index[allRows_selectFeatures].tolist()
            interactionsSubset = trainingData.iloc[rowIndices, selectFeatures]
            labelsSubset = trainingLabels.iloc[rowIndices, 0]
            weight = labelsSubset.sum() / len(labelsSubset)

        else:
            init = 0
            for trainingData_chunk, trainingLabels_chunk in zip(trainingData, trainingLabels):
                print('number of data chunked...{0}'.format(init))
                if init == 0:
                    interactionsSubset = pd.DataFrame(columns=trainingData_chunk.columns.to_list())
                    labelsSubset = pd.DataFrame(columns=trainingLabels_chunk.columns.to_list())
                    features = trainingData_chunk.columns.tolist()
                    selectFeatures_names = [features[i] for i in selectFeatures]
                trainingData_notnull_chunk = trainingData_chunk.notnull()
                allRows_selectFeatures = trainingData_notnull_chunk.loc[:, selectFeatures].all(axis=1)
                rowIndices = allRows_selectFeatures.index[allRows_selectFeatures].tolist()
                interactionsSubset = \
                    interactionsSubset.append(trainingData_notnull_chunk.loc[rowIndices, selectFeatures])
                labelsSubset = labelsSubset.append(trainingLabels_chunk.loc[rowIndices, 0].to_frame())
                init += 1
            weight = labelsSubset.sum() / len(labelsSubset)

        existingFile = \
            glob.glob(workdir + 'models/models_{0}/model_{0}_Partition_{1}_*.pkl'.format(expName, partition))
        if existingFile and not reset:
            print('Found existing model on disk...')
            trainScore = existingFile[0].split('oob_')[-1].split('_')[0]
            weight = existingFile[0].split('weight_')[-1].split('_')[0]
            # model = pickle.load(open(existingFile, 'rb'))
            # print('Completed loading model from disk...')

        else:
            # crude proxy for whether the object is a classifier/model
            if hasattr(otherEstimator, 'random_state'):
                initModel = otherEstimator
            else:
                initModel = RandomForestClassifier(n_estimators=400, criterion='entropy', max_features=None,
                                                   oob_score=True, random_state=42, n_jobs=-1, verbose=2)
            print('Fitting model to data...')
            initModel.class_weight = {1: weight} if weight > 0 else None
            model = initModel.fit(interactionsSubset, labelsSubset)
            filename = workdir + 'models/models_{0}/model_{0}_Partition_{1}_oob_{2}_weight_{3}_numSelectFeatures_{' \
                                 '4}_numRowIndices_{5}.pkl'.format(expName, partition, model.oob_score_, weight,
                                                                   len(selectFeatures), len(rowIndices))
            print('Completed fitting model to data...now saving...')
            pickle.dump(model, open(filename, 'wb'))
            print('Completed saving model to disk...')
            trainScore = model.oob_score_
            with open(workdir +
                      'models/models_{0}/model_{0}_Partition_{1}_selectFeatures_names.txt'.format(expName, partition),
                      'w+') as f:
                for item in selectFeatures_names:
                    f.write("%s\n" % item)
            f.close()

        if not parallelProcessing:
            properties = [partition, trainScore, weight, len(selectFeatures), len(rowIndices), selectFeatures_names]
            addProperties = pd.DataFrame([properties], columns=columnNames)
            trainedModels_properties = trainedModels_properties.append(addProperties)
            print('______________________________')

    if not parallelProcessing:
        trainedModels_properties.set_index('partition ID', inplace=True)
        trainedModels_properties.to_csv(workdir + f'models/models_{expName}/Partition Specific '
                                                  f'{expName} Training Performance_{dateSignature}.csv')
        print('Completed training models per partition of the training data')
        print('____________________________________________________________')

        return trainedModels_properties
