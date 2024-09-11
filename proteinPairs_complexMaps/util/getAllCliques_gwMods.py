#!/usr/bin/env python

#Author: Jose Lugo-Martinez
#File: getAllCliques.py
#Date: January 10, 2019
#Advisor Profs. Robert F. Murphy and Ziv Bar-Joseph
#Description: This function computes and reports all the maximal cliques in an undirected graph.

#Last Modified: November 4, 2021

#Example calls: python getAllCliques_gwMods.py integrated_MM+TSCS+Exp+Abun+Loc_candidates_complexesSize_2_top02.graph integrated_MM+TSCS+Exp+Abun+Loc_features_all_pairs_model_either_predictions_rf_predicted_complexesSize_

import copy
import networkx as nx
import sys


def getAllCliques_gwMods(pairsInfoFilename, pairsPredictionsFilename, cutoff, outputPATH):

    TOP_CUTOFF = cutoff

    # Open input file
    try:
        infilePairsInfo = open(pairsInfoFilename, "r")
    except IOError as e:
        print("<<ERROR>> Unable to open the file", pairsInfoFilename, "\n This program will be quiting now."), e
        sys.exit()

    # Open input file
    predictionsFilename = pairsPredictionsFilename + '.txt'
    try:
        infileScores = open(predictionsFilename, "r")
    except IOError as e:
        print("<<ERROR>> Unable to open the file", predictionsFilename, "\n This program will be quiting now."), e
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
        try:
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
        except:
            pass
    # Close files
    infilePairsInfo.close()
    infileScores.close()
    outfile.close()
    outfilename = pairsPredictionsFilename.split('_features_')[0] + '_candidates_complexesSize_2@' + str(cutoff) + '.gv'
    outfile = open(outfilename, 'w')

    # Prepare headers
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

    outfilenameGraph = pairsPredictionsFilename.split('_features_')[0] + '_candidates_complexesSize_2@' + str(
        cutoff) + '.graph'
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

        # Prepare footers
        outline = '}' + '\n'
    outfile.writelines(outline)

    # Close file
    outfile.close()
    outfileGraph.close()

    ##insert code for "get all cliques"
    infile = open(outfilenameGraph, "r")

    headers = infile.readline().strip()

    G = nx.Graph()
    comparisonMatrix = {}
    for line in infile:
        line = line.strip()
        edgeType, prot1ID, prot2ID, pScore = line.split(' ')
        G.add_edge(prot1ID, prot2ID)
        proteinPair = (prot1ID, prot2ID)
        comparisonMatrix[proteinPair] = float(pScore)
    # Close input file
    infile.close()

    print("# of nodes", len(G.nodes))
    print("# of edges", len(G.edges))

    cliques = nx.find_cliques(G)

    predictedComplexes = {}
    for cliqueVerticesList in cliques:
        cliqueVerticesList.sort()
        complexID = ','.join(cliqueVerticesList)
        complexSize = len(cliqueVerticesList)
        scoresPerPair = []
        for i in range((len(cliqueVerticesList) - 1)):
            for j in range(i + 1, len(cliqueVerticesList)):
                pairF = (cliqueVerticesList[i], cliqueVerticesList[j])
                pairB = (cliqueVerticesList[j], cliqueVerticesList[i])
                pairScore = 0.0
                if pairF in comparisonMatrix:
                    pairScore = comparisonMatrix[pairF]
                elif pairB in comparisonMatrix:
                    pairScore = comparisonMatrix[pairB]
                scoresPerPair.append(pairScore)
        complexScore = min(scoresPerPair)
        if complexSize in predictedComplexes:
            predictedComplexes[complexSize].append((complexID, complexScore))
        else:
            predictedComplexes[complexSize] = [(complexID, complexScore)]

    print(len(predictedComplexes.keys()))
    # gw change
    outfilenamePrefix = pairsPredictionsFilename
    for complexSize, complexesSet in predictedComplexes.items():
        outfilename = outfilenamePrefix + str(complexSize) + '.txt'
        outfile = open(outfilename, "w")
        for complexID, complexScore in complexesSet:
            outline = complexID + '\t' + str(complexScore) + '\n'
            outfile.writelines(outline)
        # Close output file
        outfile.close()

    return

def main(args):
    print(args)
    if len(args) == 5:
        pairsInfoFilename = args[1]
        pairsPredictionsFilename = args[2]
        cutoff = float(args[3])
        outputPATH = args[4]
    else:
        print('<<ERROR>>  Invalid number of parameters!')
        return

    getAllCliques_gwMods(pairsInfoFilename, pairsPredictionsFilename, cutoff, outputPATH)


if __name__ == '__main__':
    main(sys.argv)
