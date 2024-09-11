#!/usr/bin/env python

#Author: Jose Lugo-Martinez
#File: buildComplexesFromPairwiseScoresRecursive.py
#Date: July 03, 2018
#Advisor Prof. Robert Murphy and Ziv Bar-Joseph

#Sample call: python buildComplexesFromPairwiseScoresRecursive.py all_either_model.tsv MM+TSCS+Exp+Abun+Loc+Kaggle_features_all_either_model_predictions_rf_400 0.959012 50.0 ./results/

import copy
import numpy as np

MAX_K = 45

def buildComplexes(pairsInfoFilename, pairsPredictionsFilename, cutoff, targetPercentile, outputPATH):
	TOP_CUTOFF = cutoff 	    #0.7740 for 0.05% high-scoring pairs; 0.8663 for 0.02% high-scoring pairs
                                #0.9251 for 0.01% high-scoring pairs; 0.9615 for 0.005% high-scoring pairs
                                #0.9251 for 0.02% high-scoring pairs; 0. for 0.005% high-scoring pairs
                                #0.9840 for 0.002% high-scoring pairs; 0.9898 for 0.001% high-scoring pairs

	try:
		#Open input file
		infilePairsInfo = open(pairsInfoFilename, "r")
	except(IOError), e:
		print "<<ERROR>> Unable to open the file", pairsInfoFilename, "\nThis program will be quiting now.", e
		sys.exit()

	predictionsFilename = pairsPredictionsFilename + '.txt'
	try:
		#Open input file
		infileScores = open(predictionsFilename, "r")
	except(IOError), e:
		print "<<ERROR>> Unable to open the file", predictionsFilename, "\nThis program will be quiting now.", e
		sys.exit()

	proteins = {}
	restrictedProteins = {}
	currComplexes = {}
	scores = {}
	comparisonMatrix = {}
	outfilename = outputPATH + pairsPredictionsFilename + '_candidates_complexesSize_2.txt'
	outfile = open(outfilename, 'w')
	for line in infileScores:
		line = line.strip()
		tokens = line.split('\t')
		pScore = float(tokens[0])
		
		infoLine = infilePairsInfo.readline().strip()
		tokens = infoLine.split('\t')
		classLabel = tokens[0]
		prot1ID = tokens[1]
		prot2ID = tokens[2]
		pair1 = (prot1ID, prot2ID)
		pair2 = (prot2ID, prot1ID)
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
	for complexID, attributes in currComplexes.iteritems():
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

#	return

	percentile = np.percentile(scores.values(), targetPercentile)
	averageWeight = np.mean(scores.values())
	stdWeight = np.std(scores.values())
##        weightThreshold = averageWeight + (1.0 * stdWeight)
##        weightThreshold = averageWeight + (2.0 * stdWeight)
	weightThreshold = percentile

	print '# of proteins', len(proteins), len(restrictedProteins)
	print 'Candidate pairwise interactions', len(comparisonMatrix)
	print '\t High-scoring pairwise interactions', len(currComplexes)
	print '\t Average pairwise score', averageWeight
	print '\t Standard deviation of pairwise scores', stdWeight
	outline = '\t ' + str(targetPercentile) + ' percentile score ' + str(weightThreshold)
	print outline

	k = 3
	while (len(currComplexes) > 0 and k <= MAX_K):
		complexesCount = 0
		candidateComplexes = copy.copy(currComplexes)
		print (k - 1), len(restrictedProteins)
##		proteinSet = proteins.keys()
		proteinSet = restrictedProteins.keys()
		restrictedProteins = {}
		scores = {}
		currComplexes = {}
		alreadyPrinted = {}
##		outfilenameCandidates = outputPATH + pairsPredictionsFilename + '_candidates_complexesSize_' + str(k) + '.txt'
##		outfileCandidates = open(outfilenameCandidates, 'w')
		outfilenamePredicted = outputPATH + pairsPredictionsFilename + '_predicted_complexesSize_' + str(k-1) + '.txt'
                outfilePredicted = open(outfilenamePredicted, "w")
		for candidateComplex, weight in candidateComplexes.values():
##			prune = False
			augmented = False
##			temp = []
####			temp = copy.copy(candidateComplex)
####			if weight < TOP_CUTOFF:
####				continue
			for prot2ID in proteinSet:
				prune = False
				temp = []
				for i in xrange(len(candidateComplex)):
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
				if (prune == False and len(temp) == k):
					augmented = True
					temp.sort()
					complexID = ','.join(temp)
					if not (complexID in currComplexes):
						complexMembers = complexID.split(',')
						scoresPerPair = []
						for x in xrange((len(complexMembers) - 1)):
							for y in xrange(x + 1, len(complexMembers)):
								pairF = (complexMembers[x], complexMembers[y])
								pairB = (complexMembers[y], complexMembers[x])
								pairScore = 0.0
								if pairF in comparisonMatrix:
									pairScore = comparisonMatrix[pairF]
								elif pairB in comparisonMatrix:
									pairScore = comparisonMatrix[pairB]
								scoresPerPair.append(pairScore)
						newWeight = min(scoresPerPair)
####						newWeight = np.mean(scoresPerPair)
####						newWeight = np.percentile(scoresPerPair, 50)
						if newWeight < TOP_CUTOFF:
							continue
##						outline = complexID
##						for pairScore in scoresPerPair:
##                                                        outline += '\t' + str(pairScore) 
##						outline += '\n'
##						outfileCandidates.writelines(outline)
						currComplexes[complexID] = (temp, newWeight)
						for prot in temp:
							if not (prot in restrictedProteins):
                                				restrictedProteins[prot] = prot
						scores[complexID] = newWeight
			if augmented == False:
				candidateComplex.sort()
				key = ','.join(candidateComplex)
				outline = key + '\t' + str(weight) + '\n'
				outfilePredicted.writelines(outline)
				complexesCount += 1
##		#Close outfile
##		outfileCandidates.close()
		outfilePredicted.close()

		outline = '\t Total of predicted ' + str(k-1) + '-complexes'
		print outline, complexesCount, len(candidateComplexes)

		outline = 'Candidate ' + str(k) + '-complexes'
		if (len(scores) < 1):
			print outline, len(currComplexes)
			continue
##		percentile = np.percentile(scores.values(), targetPercentile)
		averageWeight = np.mean(scores.values())
		stdWeight = np.std(scores.values())
		weightThreshold = averageWeight + (1.0 * stdWeight)
##		weightThreshold = averageWeight + (2.0 * stdWeight)
##		weightThreshold = percentile

		print outline, len(currComplexes), averageWeight, stdWeight, weightThreshold

		k = k + 1

	return

def main(argv):
	if (len(argv) == 6):
		pairsInfoFilename = argv[1]
		pairsPredictionsFilename = argv[2]
		cutoff = float(argv[3])
		targetPercentile = float(argv[4])
		outputPATH = argv[5]
	else:
		print "<<ERROR>> Invalid number of parameters!"
		return

	buildComplexes(pairsInfoFilename, pairsPredictionsFilename, cutoff, targetPercentile, outputPATH)

if __name__ == '__main__':
#        main(sys.args) #Alternate call, if you are using IDLE.
#	main(['buildComplexesFromPairwiseScoresRecursive.py', 'all_either_model.tsv', 'MM+TSCS+Exp+Abun+Loc+Kaggle_features_all_either_model_predictions_rf_400', '0.0', '50.0', './results/']) #top 0.10% 
#	main(['buildComplexesFromPairwiseScoresRecursive.py', 'all_either_model.tsv', 'MM+TSCS+Exp+Abun+Loc+Kaggle_features_all_either_model_predictions_rf_400', '0.959012', '50.0', './results/top001_min/']) #top 0.001%
#	main(['buildComplexesFromPairwiseScoresRecursive.py', 'all_either_model.tsv', 'MM+TSCS+Exp+Abun+Loc+Kaggle_features_all_either_model_predictions_rf_400', '0.930573', '50.0', './results/top002_min/']) #top 0.002%
#	main(['buildComplexesFromPairwiseScoresRecursive.py', 'all_either_model.tsv', 'MM+TSCS+Exp+Abun+Loc+Kaggle_features_all_either_model_predictions_rf_400', '0.897115', '50.0', './results/top005_min/']) #top 0.005%
#	main(['buildComplexesFromPairwiseScoresRecursive.py', 'all_either_model.tsv', 'MM+TSCS+Exp+Abun+Loc+Kaggle_features_all_either_model_predictions_rf_400', '0.854287', '50.0', './results/top01_min/']) #top 0.01%
	main(['buildComplexesFromPairwiseScoresRecursive.py', 'all_either_model.tsv', 'MM+TSCS+Exp+Abun+Loc+Kaggle_features_all_either_model_predictions_rf_400', '0.760267', '50.0', './results/top02_min/']) #top 0.02%
#	main(['buildComplexesFromPairwiseScoresRecursive.py', 'all_either_model.tsv', 'MM+TSCS+Exp+Abun+Loc+Kaggle_features_all_either_model_predictions_rf_400', '0.655646', '50.0', './results/top05_min/']) #top 0.05%
