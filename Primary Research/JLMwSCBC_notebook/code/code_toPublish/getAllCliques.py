#!/usr/bin/env python

#Author: Jose Lugo-Martinez
#File: getAllCliques.py
#Date: January 10, 2019
#Advisor Profs. Robert F. Murphy and Ziv Bar-Joseph
#Description: This function computes and reports all the maximal cliques in an undirected graph.

#Last Modified: January 11, 2019

# Example calls: python getAllCliques.py integrated_MM+TSCS+Exp+Abun+Loc_candidates_complexesSize_2_top02.graph integrated_MM+TSCS+Exp+Abun+Loc_features_all_pairs_model_either_predictions_rf_predicted_complexesSize_

import networkx as nx
import sys


def main(args):
    print(args)
    if len(args) == 3:
        graphFilename = args[1]
        outfilenamePrefix = args[2]
    else:
        print ("<<ERROR>> Invalid number of parameters!")
        return
   
    infile  = open(graphFilename, "r")
    next(infile)

    G = nx.Graph()
    comparisonMatrix = {}
    for line in infile:
        print(line)
        line = line.strip()
        edgeType, prot1ID, prot2ID, pScore = line.split(' ')
        G.add_edge(prot1ID, prot2ID)  
        proteinPair = (prot1ID, prot2ID)
        comparisonMatrix[proteinPair] = float(pScore)
    #Close input file    
    infile.close()
    
    print ("# of nodes", len(G.nodes))
    print ("# of edges", len(G.edges))
    
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

    print (len(predictedComplexes.keys()))

    for complexSize, complexesSet in predictedComplexes.items():
        outfilename = outfilenamePrefix + str(complexSize) + '.txt'
        outfile  = open(outfilename, "w")
        for complexID, complexScore in complexesSet:
            outline = complexID + '\t' + str(complexScore) + '\n'
            outfile.writelines(outline)
        #Close output file    
        outfile.close()
    
    return


if __name__ == '__main__':
    main(sys.argv)
