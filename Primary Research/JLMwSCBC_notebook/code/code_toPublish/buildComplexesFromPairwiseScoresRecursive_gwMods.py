#!/usr/bin/env python

#Author: Jose Lugo-Martinez
#File: buildComplexesFromPairwiseScoresRecursive.py
#Date: July 03, 2018
#Advisor Prof. Robert Murphy and Ziv Bar-Joseph

#Sample call: python buildComplexesFromPairwiseScoresRecursive_gwMods.py all_either_model.tsv MM+TSCS+Exp+Abun+Loc+Kaggle_features_all_either_model_predictions_rf_400 0.959012 50.0

import copy
import numpy as np
import sys


def buildComplexes(pairsInfoFilename, pairsPredictionsFilename, cutoff, targetPercentile):

  maxK = 45
  TOP_CUTOFF = cutoff

  #Open input file
  try:
      infilePairsInfo = open(pairsInfoFilename, "r")
  except IOError as e:
      print("<<ERROR>> Unable to open the file", pairsInfoFilename, "\nThis program will be quiting now."), e
      sys.exit()

  #Open input file
  predictionsFilename = pairsPredictionsFilename + '.txt'
  try:
      infileScores = open(predictionsFilename, "r")
  except IOError as e:
      print("<<ERROR>> Unable to open the file", predictionsFilename, "\nThis program will be quiting now."), e
      sys.exit()

  proteins = {}
  restrictedProteins = {}
  currComplexes = {}
  scores = {}
  comparisonMatrix = {}
  outfilename = pairsPredictionsFilename + '_candidates_complexesSize_2.txt'
  outfile = open(outfilename, 'w')
  for line in infileScores:
      line = line.strip()
      tokens = line.split('\t')
      pScore = float(tokens[0])

      infoLine = infilePairsInfo.readline().strip()
      tokens = infoLine.split('\t')
      prot1ID = tokens[1]
      prot2ID = tokens[2]
      pair1 = (prot1ID, prot2ID)
      comparisonMatrix[pair1] = pScore

      if not (prot1ID in proteins):
          proteins[prot1ID] = prot1ID
      if not (prot2ID in proteins):
          proteins[prot2ID] = prot2ID

      temp = []
      if pScore >= TOP_CUTOFF:
          if not (prot1ID in restrictedProteins):
              restrictedProteins[prot1ID] = prot1ID
          if not (prot2ID in restrictedProteins):
              restrictedProteins[prot2ID] = prot2ID

          if not (prot1ID in temp):
              temp.append(prot1ID)
          if not (prot2ID in temp):
              temp.append(prot2ID)
          temp.sort()
          complexID = ','.join(temp)
          if not (complexID in currComplexes):
              outline = complexID + '\t' + str(pScore) + '\n'
              outfile.writelines(outline)
              currComplexes[complexID] = (temp, pScore)
              scores[complexID] = pScore

  #Close files
  infilePairsInfo.close()
  infileScores.close()
  outfile.close()
  outfilename = pairsPredictionsFilename.split('_features_')[0] + '_candidates_complexesSize_2@' + str(cutoff) + '.gv'
  outfile = open(outfilename, 'w')

  #Prepare headers
  outline = 'graph {' + '\n'
  outline += '\t node [shape = circle];' + '\n'
  outline += '\t edge [style = solid];' + '\n'
  outfile.writelines(outline)
  nodeID = 0
  protID2nodeID = {}
  for protID in restrictedProteins.keys():
      if not (protID in protID2nodeID):
          protID2nodeID[protID] = nodeID
          outline = '\t n' + str(nodeID) + ' [label = p' + str(nodeID) + '];' + '\n'
          outfile.writelines(outline)
          nodeID += 1

  outfilenameGraph = pairsPredictionsFilename.split('_features_')[0] + '_candidates_complexesSize_2@' + str(cutoff) + '.graph'
  outfileGraph = open(outfilenameGraph, 'w')
  outlineGraph = 'p edge ' + str(len(restrictedProteins)) + '\n'
  outfileGraph.writelines(outlineGraph)
  for complexID, attributes in currComplexes.items():
      prot1ID, prot2ID = complexID.split(',')
      edgeScore = copy.copy(attributes[1])
      if prot1ID in protID2nodeID:
          node1ID = protID2nodeID[prot1ID]
      if prot2ID in protID2nodeID:
          node2ID = protID2nodeID[prot2ID]
      outline = '\t n' + str(node1ID) + ' -- ' + 'n' + str(node2ID) + ' [label = ' + str(edgeScore) + '];' + '\n'
      outfile.writelines(outline)
      outlineGraph = 'e ' + str(prot1ID) + ' ' + str(prot2ID) + ' ' + str(edgeScore) + '\n'
      outfileGraph.writelines(outlineGraph)

      #Prepare footers
      outline = '}' + '\n'
  outfile.writelines(outline)

  #Close file
  outfile.close()
  outfileGraph.close()

  percentile = np.percentile(list(scores.values()), targetPercentile)
  averageWeight = np.mean(list(scores.values()))
  stdWeight = np.std(list(scores.values()))
  weightThreshold = percentile

  print('# of proteins', len(proteins), len(restrictedProteins))
  print('Candidate pairwise interactions', len(comparisonMatrix))
  print('\t High-scoring pairwise interactions', len(currComplexes))
  print('\t Average pairwise score', averageWeight)
  print('\t Standard deviation of pairwise scores', stdWeight)
  outline = '\t ' + str(targetPercentile) + ' percentile score ' + str(weightThreshold)
  print(outline)

  k = 3
  while len(currComplexes) > 0 and k <= maxK:
      complexesCount = 0
      candidateComplexes = copy.copy(currComplexes)
      print (k - 1), len(restrictedProteins)
      proteinSet = restrictedProteins.keys()
      restrictedProteins = {}
      scores = {}
      currComplexes = {}
      outfilenamePredicted = pairsPredictionsFilename + '_predicted_complexesSize_' + str(k-1) + '.txt'
      outfilePredicted = open(outfilenamePredicted, "w")
      for candidateComplex, weight in candidateComplexes.values():
          augmented = False
          for prot2ID in proteinSet:
              prune = False
              temp = []
              for i in range(len(candidateComplex)):
                  prot1ID = candidateComplex[i]
                  if prot2ID in candidateComplex:
                      break
                  pair1 = (prot1ID, prot2ID)
                  pair2 = (prot2ID, prot1ID)
                  if pair1 in comparisonMatrix:
                      pScore = comparisonMatrix[pair1]
                  elif pair2 in comparisonMatrix:
                      pScore = comparisonMatrix[pair2]
                  else:
                      break
                  if pScore >= TOP_CUTOFF: #PAIRWISE_INTERACTIONS_CUTOFF:
                      if not (prot1ID in temp):
                          temp.append(prot1ID)
                      if not (prot2ID in temp):
                          temp.append(prot2ID)
                  else:
                      prune = True
                      break
              if prune == False and len(temp) == k:
                  augmented = True
                  temp.sort()
                  complexID = ','.join(temp)
                  if not (complexID in currComplexes):
                      complexMembers = complexID.split(',')
                      scoresPerPair = []
                      for x in range((len(complexMembers) - 1)):
                          for y in range(x + 1, len(complexMembers)):
                              pairF = (complexMembers[x], complexMembers[y])
                              pairB = (complexMembers[y], complexMembers[x])
                              pairScore = 0.0
                              if pairF in comparisonMatrix:
                                  pairScore = comparisonMatrix[pairF]
                              elif pairB in comparisonMatrix:
                                  pairScore = comparisonMatrix[pairB]
                              scoresPerPair.append(pairScore)
                      newWeight = min(scoresPerPair)
                      if newWeight < TOP_CUTOFF:
                          continue

                      currComplexes[complexID] = (temp, newWeight)
                      for prot in temp:
                          if not (prot in restrictedProteins):
                              restrictedProteins[prot] = prot
                      scores[complexID] = newWeight
          if not augmented:
              candidateComplex.sort()
              key = ','.join(candidateComplex)
              outline = key + '\t' + str(weight) + '\n'
              outfilePredicted.writelines(outline)
              complexesCount += 1

      #Close outfile
      outfilePredicted.close()

      outline = '\t Total of predicted ' + str(k-1) + '-complexes'
      print(outline, complexesCount, len(candidateComplexes))
      outline = 'Candidate ' + str(k) + '-complexes'
      if len(scores) < 1:
          print(outline, len(currComplexes))
          continue
      averageWeight = np.mean(list(scores.values()))
      stdWeight = np.std(list(scores.values()))
      weightThreshold = averageWeight + (1.0 * stdWeight)
      print(outline, len(currComplexes), averageWeight, stdWeight, weightThreshold)

      k = k + 1

  return


def main(argv):
    if len(argv) == 5:
        pairsInfoFilename = argv[1]
        pairsPredictionsFilename = argv[2]
        cutoff = float(argv[3])
        targetPercentile = float(argv[4])
    else:
        print("<<ERROR>> Invalid number of parameters!")
        return

    buildComplexes(pairsInfoFilename, pairsPredictionsFilename, cutoff, targetPercentile)


if __name__ == '__main__':
    main(sys.args) #Alternate call, if you are using IDLE.
    #main(['buildComplexesFromPairwiseScoresRecursive.py', 'all_either_model.tsv', 'MM+TSCS+Exp+Abun+Loc+Kaggle_features_all_either_model_predictions_rf_400', '0.760267', '50.0']) #top 0.02%
