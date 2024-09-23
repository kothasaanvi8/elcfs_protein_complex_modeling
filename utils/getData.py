# Load requisite data from original sources

import numpy as np
from os.path import exists
import pandas as pd
import pickle, re, sys

def lowercaseStr(DF):
    dummyDF = pd.DataFrame({'dummyCol': ['a', 'b', 'c']})
    DF.columns = DF.columns.str.lower()
    for col in DF.columns:
        if (DF[col].dtypes == dummyDF.dummyCol.dtypes) and (len(DF[col].unique())!=2):
            DF[col] = DF[col].str.lower()

    return DF

def splitDF_byDelimiters(DF, fieldNames, delimiters):
    #print('original dimensions: ' + str(DF.shape))
    for fieldName in fieldNames:
        for delimiter in delimiters:
            DF = DF.assign(
                temp=DF[fieldName].str.split(delimiter)
            ).explode('temp')
            DF.drop(columns=[fieldName], inplace=True)
            DF.rename(columns={'temp': fieldName}, inplace=True)
            #print('dimensions post changes: ' + str(DF.shape))

    return DF

def dataframeSummary(DF):
    print('Printing summary...')
    print('Printing columns of data...')
    print(DF.columns)
    print('Printing the shape of data...')
    print(DF.shape)
    print('Printing head of data...')
    print(DF.head())
    print('____________')

class LoadDict:

    def source(self, userInput):
        sourceName = userInput[0]
        otherArgs = userInput[1:]
        return getattr(self, sourceName, lambda: 'There is no data for the named source.')(otherArgs)

    def bioplex(self, otherArgs):
        if len(otherArgs) != 1:
            print('An incorrect number of arguments was passed.')
            return

        outputProteins_flag = otherArgs[0]

        bioplexData = {}
        bioplexProtein_interactions = lowercaseStr(
            pd.read_csv('./sourceData/bioplex/BioPlex_interactionList_v2.tsv', sep='\t')
        )
        bioplexProtein_interactions.rename(columns={'gene a': 'geneid1', 'gene b': 'geneid2'}, inplace=True)
        bioplexData['interactions'] = bioplexProtein_interactions

        if outputProteins_flag:
            bioplexInteractions_summaryUnique_geneidsList = set(bioplexProtein_interactions['geneid1'].to_list()).union(
                set(bioplexProtein_interactions['geneid2'].to_list()))

            bioplexInteractions_summaryUnique_geneidsList = set(
                [str(entry) for entry in bioplexInteractions_summaryUnique_geneidsList]
            )

            bioplexData['proteins'] = bioplexInteractions_summaryUnique_geneidsList

        return bioplexData

    def corum(self, otherArgs):
        #published releases updated 2022.01.30
        publishedReleases = {'2018.09.03', '2018.07.01', '2017.07.02',
                             '2017.05.03', '2016.12.15', '2012.02.17'}

        if len(otherArgs) >= 2:
            releaseDate = otherArgs[0]
            if releaseDate in publishedReleases:

                if releaseDate == '2016.12.15':
                    print('The Helmholtz Zentrum München website maps to release 2017.05.03 in error, '
                          'please check back later...')

                else:
                    corumData = dict()
                    sheetName = otherArgs[1]
                    corumDF_all = pd.read_excel('./sourceData/corum/' + releaseDate + '/corum.xls', sheet_name=None)

                    if sheetName in corumDF_all.keys():
                        corumDF = lowercaseStr(corumDF_all[sheetName])
                        corumData['rawData_xls'] = corumDF
                        complexID_col = [col for col in corumDF.columns if \
                                         all(substring for substring in ['complex', 'id'] if substring in col)][0]
                        geneID_col = [col for col in corumDF.columns if 'entrez' in col][0]
                        print('CORUM complexes dataframe dimensions initially: ' + str(corumDF.shape))
                        if len(otherArgs) > 2:

                            if len(otherArgs) != 6:
                                print('Insufficient number of arguments passed. '
                                      'If one flag is used, all flags must be used.')
                                return
                            else:
                                ribosomeOnly_flag = otherArgs[2]
                                missingGeneid_flag = otherArgs[3]
                                explicitNone_flag = otherArgs[4]
                                onlyHumans_flag = otherArgs[5]

                                if ribosomeOnly_flag and sheetName != 'spliceComplexes':
                                    complexName_col = \
                                        [col for col in corumDF.columns if 'name' in col and 'complex' in col][0]
                                    smallRibosomal_subunitLine = corumDF.loc[
                                        corumDF[complexName_col] == '40s ribosomal subunit, cytoplasmic',
                                        geneID_col].iloc[0]
                                    smallRibosomal_subunit = \
                                        [ele for ele in re.split(
                                            ',\(|\(|\),|\)|\t|\s|;|,',
                                            smallRibosomal_subunitLine) if ele]
                                    largeRibosomal_subunitLine = corumDF.loc[
                                        corumDF[complexName_col] == '60s ribosomal subunit, cytoplasmic',
                                        geneID_col].iloc[0]
                                    largeRibosomal_subunit = \
                                        [ele for ele in re.split(
                                            ',\(|\(|\),|\)|\t|\s|;|,',
                                            largeRibosomal_subunitLine) if ele]

                                    return smallRibosomal_subunit, largeRibosomal_subunit

                                if missingGeneid_flag:
                                    if explicitNone_flag:
                                        corumDF = corumDF.loc[corumDF[geneID_col].str.contains('None') == False, :]
                                    else:
                                        corumDF = corumDF.dropna(subset=[geneID_col])
                                    corumData['rawData_nonEmpty_geneIDs'] = corumDF
                                    print('CORUM complexes dataframe dimensions after removal of <None> or NaN '
                                          'entries: ' + str(corumDF.shape))

                                if corumDF.loc[:, geneID_col].isna().sum() > 0:
                                    print('Program exiting. Chosen specs contain NaNs in GeneID column.')

                                    return

                                if onlyHumans_flag:
                                    #!ATTN: Design a util to address formatting variability, re: column name cases
                                    corumDF = corumDF.loc[corumDF.organism=='human', [complexID_col, geneID_col]]
                                    corumData['filteredData_humanSubset'] = corumDF
                                    print('CORUM complexes dataframe dimensions post-organism-filter [Human]: ' +
                                          str(corumDF.shape))

                        complexesCounter = 0
                        corumDF_subsetComplexes = dict()
                        corumDF_subsetComplexes_dataLengths = dict()
                        CORUMproteins_subset = set()
                        for idx in corumDF.index:
                            #print(corumDF_subset.loc[idx, geneID_col])  # for debugging
                            complexLine = \
                                [ele for ele in re.split(',\(|\(|\),|\)|\t|\s|;|,', corumDF.loc[idx, geneID_col]) if ele]
                            complexDelimited = set(list(map(str, complexLine)))
                            corumDF_subsetComplexes[corumDF.loc[idx, complexID_col]] = complexDelimited
                            corumDF_subsetComplexes_dataLengths[corumDF.loc[idx, complexID_col]] = \
                                len(complexDelimited)
                            CORUMproteins_subset = CORUMproteins_subset.union(complexDelimited)

                            complexesCounter += 1

                        corumData['complexes'] = corumDF_subsetComplexes
                        corumData['complexLengths'] = corumDF_subsetComplexes_dataLengths
                        corumData['proteins'] = CORUMproteins_subset

                        print('# curated complexes in set: ' + str(len(corumDF_subsetComplexes)) + '.')
                        print('--min complex length: ' +
                              str(np.amin(np.array(list(set(list(
                                  corumDF_subsetComplexes_dataLengths.values())))))))
                        print('--mean complex length: ' +
                              str(np.mean(np.array(list(set(list(
                                  corumDF_subsetComplexes_dataLengths.values())))))))
                        print('--median complex length: ' +
                              str(np.median(np.array(list(set(list(
                                  corumDF_subsetComplexes_dataLengths.values())))))))
                        print('--max complex length: ' +
                              str(np.amax(np.array(list(set(list(
                                  corumDF_subsetComplexes_dataLengths.values())))))))
                        print('# total unique proteins in set: ' + str(len(CORUMproteins_subset)) + '.')

                        return corumData

                    else:
                        print('The specified worksheet for CORUM data with release date, ' + releaseDate +
                              ', was not found. Available worksheets include <allComplexesCore>, <allComplexes>, ' +
                              'and <spliceComplexes>.')

            else:
                print('The specified release date is unavailable. Confirm the specified release date and '
                      'the expected format (e.g. yyyy.mm.dd).')

        else:
            print('There are too few arguments for loading the CORUM data.')

    # need the real source!!!
    def fantom(self, otherArgs):
        print('...under development.')

    # under development
    def go(self, otherArgs):
        swissProt_plusGO_cellularComponent = lowercaseStr(
            pd.read_csv('./sourceData/go/uniprot-filtered-organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22.tab',
                        sep='\t', low_memory=False)
        )
        swissProt_plusGO_cellularComponent.rename(columns={'cross-reference (geneid)': 'geneid'}, inplace=True)

        return swissProt_plusGO_cellularComponent

    # need the real source!!!
    def gtex(self, otherArgs):
        print('...under development.')

    def hein2015(self, otherArgs):

        if len(otherArgs) != 5:
            print('Incorrect number of arguments passed.')
            return

        rawForm_flag = otherArgs[0]
        ignoreIsoforms_flag = otherArgs[1]
        mergeGeneids_flag = otherArgs[2]
        manMerge_refPlus_intact2Geneids_flag = otherArgs[3]
        outputProteins_flag = otherArgs[4]

        hein2015_data = {}
        hein2015_proteinInteractions = lowercaseStr(
            pd.read_csv('./sourceData/hein2015/query-IM-24272-18012022_0716.txt', sep='\t', low_memory=False)
        )

        if rawForm_flag:
            hein2015_data['rawInteractions'] = hein2015_proteinInteractions
            hein2015_data['rawInteractions']['idsInteractors_frozen'] = \
                hein2015_data['rawInteractions'].loc[:,
                ['#id(s) interactor a', 'id(s) interactor b']].apply(frozenset, axis=1)

        hein2015_proteinInteractions.loc[
            hein2015_proteinInteractions['#id(s) interactor a'].str.contains('uniprotkb:'),
            '#id(s) interactor a'] = \
            hein2015_proteinInteractions.loc[
                hein2015_proteinInteractions['#id(s) interactor a'].str.contains('uniprotkb:'),
                '#id(s) interactor a'].str.replace('uniprotkb:', '')

        hein2015_proteinInteractions.loc[
            hein2015_proteinInteractions['id(s) interactor b'].str.contains('uniprotkb:'),
            'id(s) interactor b'] = \
            hein2015_proteinInteractions.loc[
                hein2015_proteinInteractions['id(s) interactor b'].str.contains('uniprotkb:'),
                'id(s) interactor b'].str.replace('uniprotkb:', '')
        hein2015_data['processedRemoved_entryPrefix_uniprotkb'] = hein2015_proteinInteractions
        hein2015_data['processedRemoved_entryPrefix_uniprotkb']['idsInteractors_frozen'] = \
            hein2015_data['processedRemoved_entryPrefix_uniprotkb'].loc[:,
            ['#id(s) interactor a', 'id(s) interactor b']].apply(frozenset, axis=1)

        if ignoreIsoforms_flag:
            hein2015_proteinInteractions.loc[
                (hein2015_proteinInteractions['#id(s) interactor a'].str.contains('uniprotkb:')) &
                                             (hein2015_proteinInteractions['#id(s) interactor a'].str.contains('-')),
                '#id(s) interactor a'] = \
                hein2015_proteinInteractions.loc[
                    (hein2015_proteinInteractions['#id(s) interactor a'].str.contains('uniprotkb:')) &
                    (hein2015_proteinInteractions['#id(s) interactor a'].str.contains('-')),
                    '#id(s) interactor a'].str.split('-').str[0]

            hein2015_proteinInteractions.loc[
                (hein2015_proteinInteractions['id(s) interactor b'].str.contains('uniprotkb:')) &
                                             (hein2015_proteinInteractions['#id(s) interactor a'].str.contains('-')),
                'id(s) interactor b'] = \
                hein2015_proteinInteractions.loc[
                    (hein2015_proteinInteractions['id(s) interactor b'].str.contains('uniprotkb:')) &
                    (hein2015_proteinInteractions['id(s) interactor b'].str.contains('-')),
                    'id(s) interactor b'].str.split('-').str[0]
            hein2015_data['processedIgnore_isoforms'] = hein2015_proteinInteractions
            hein2015_data['processedIgnore_isoforms']['idsInteractors_frozen'] = \
                hein2015_data['processedIgnore_isoforms'].loc[:,
                ['#id(s) interactor a', 'id(s) interactor b']].apply(frozenset, axis=1)

        if mergeGeneids_flag:
            hein2015_proteinInteraction_uniprotIDs_geneidsConverted_byDAVID = lowercaseStr(
                pd.read_csv('./sourceData/hein2015/hein2015_proteinInteraction_uniprotIDs_geneidsConverted_byDAVID.txt',
                            sep='\t', usecols=[0, 1], dtype='str')
            )

            hein2015_proteinInteractions = hein2015_proteinInteractions.merge(
                hein2015_proteinInteraction_uniprotIDs_geneidsConverted_byDAVID,
                how='left', left_on=['#id(s) interactor a'], right_on=['from'])
            hein2015_proteinInteractions.rename(columns={'to': 'geneid1'}, inplace=True)
            hein2015_proteinInteractions.drop(columns=['from'], inplace=True)

            hein2015_proteinInteractions = hein2015_proteinInteractions.merge(
                hein2015_proteinInteraction_uniprotIDs_geneidsConverted_byDAVID,
                how='left', left_on=['id(s) interactor b'], right_on=['from'])
            hein2015_proteinInteractions.rename(columns={'to': 'geneid2'}, inplace=True)
            hein2015_proteinInteractions.drop(columns=['from'], inplace=True)

            hein2015_data['processedMerged_geneids2Uniprotid'] = hein2015_proteinInteractions
            hein2015_data['processedMerged_geneids2Uniprotid']['geneidsFrozen'] = \
                hein2015_data['processedMerged_geneids2Uniprotid'].loc[:,
                ['geneid1', 'geneid2']].apply(frozenset, axis=1)

            if manMerge_refPlus_intact2Geneids_flag:
                hein2015_proteinInteractions.loc[hein2015_proteinInteractions['#id(s) interactor a'] ==
                                                     'intact:ebi-10964486', 'geneid1'] = '5646'
                hein2015_proteinInteractions.loc[hein2015_proteinInteractions['id(s) interactor b'] ==
                                                     'intact:ebi-10964486', 'geneid2'] = '5646'

                hein2015_proteinInteractions.loc[hein2015_proteinInteractions['#id(s) interactor a'] ==
                                                     'refseq:np_003365.1', 'geneid1'] = '7416'
                hein2015_proteinInteractions.loc[hein2015_proteinInteractions['id(s) interactor b'] ==
                                                     'refseq:np_003365.1', 'geneid2'] = '7416'

                hein2015_proteinInteractions.loc[hein2015_proteinInteractions['#id(s) interactor a'] ==
                                                     'refseq:np_005990', 'geneid1'] = '7257'
                hein2015_proteinInteractions.loc[hein2015_proteinInteractions['id(s) interactor b'] ==
                                                     'refseq:np_005990', 'geneid2'] = '7257'

                hein2015_proteinInteractions.loc[hein2015_proteinInteractions['#id(s) interactor a'] ==
                                                     'a0a286ycx6', 'geneID_A'] = '75964'
                hein2015_proteinInteractions.loc[hein2015_proteinInteractions['id(s) interactor b'] ==
                                                     'a0a286ycx6', 'geneID_B'] = '75964'

                hein2015_proteinInteractions.loc[((hein2015_proteinInteractions['#id(s) interactor a'] == 'p0dpb5') |
                                                  (hein2015_proteinInteractions['#id(s) interactor a'] == 'p0dpb6')),
                                                 'geneid1'] = '51082'
                hein2015_proteinInteractions.loc[((hein2015_proteinInteractions['id(s) interactor b'] == 'p0dpb5') |
                                                  (hein2015_proteinInteractions['id(s) interactor b'] == 'p0dpb6')),
                                                 'geneid2'] = '51082'

                hein2015_data['processedFound_geneids2Nonuniprotids'] = hein2015_proteinInteractions
                hein2015_data['processedFound_geneids2Nonuniprotids']['geneidsFrozen'] = \
                    hein2015_data['processedFound_geneids2Nonuniprotids'].loc[:,
                    ['geneid1', 'geneid2']].apply(frozenset, axis=1)

        if outputProteins_flag:
            if 'geneid1' in hein2015_proteinInteractions.columns and 'geneid2' in hein2015_proteinInteractions.columns:
                hein2015_proteinInteraction_geneidsList = set(hein2015_proteinInteractions.geneid1.to_list()).union(
                    set(hein2015_proteinInteractions.geneid2.to_list()))
            else:
                hein2015_proteinInteraction_geneidsList = \
                        set(hein2015_proteinInteractions['#id(s) interactor a'].to_list()).union(
                            set(hein2015_proteinInteractions['id(s) interactor b'].to_list()))

            hein2015_proteinInteraction_geneidsList = set(
                [str(entry) for entry in hein2015_proteinInteraction_geneidsList]
            )

            hein2015_data['proteins'] = hein2015_proteinInteraction_geneidsList

        return hein2015_data

    def humap1(self, otherArgs):

        if len(otherArgs) != 4:
            print('Insufficient number of arguments passed. Must specify flags for pairs, complexes, and feature '
                  'matrix.')
            return

        returnPairs = otherArgs[0]
        returnComplexes = otherArgs[1]
        returnFeature_matrix = otherArgs[2]
        returnResults = otherArgs[3]
        if not any((returnPairs, returnComplexes, returnFeature_matrix)):
            print('All flags specified False!')
            return

        else:
            humap1 = {}
            if returnPairs:
                # acquire pairs used in original huMAP 1.0 study
                ppiData = {}
                ppiData_path = './sourceData/humap1/pairs/'
                ppiData_names = ['train_ppis', 'train_neg_ppis', 'test_ppis', 'test_neg_ppis']
                for name in ppiData_names:
                    if name != 'test_ppis':
                        ppiData[name] = lowercaseStr(
                            pd.read_csv(ppiData_path + name + '.txt', sep=' ',
                                        header=None, names=['idi', 'idii'], dtype='str')
                        )
                    else:
                        ppiData[name] = lowercaseStr(
                            pd.read_csv(ppiData_path + name + '.txt', sep='\t',
                                        header=None, names=['idi', 'idii'], dtype='str')
                        )
                humap1['pairs'] = ppiData

            if returnComplexes:
                # acquire complexes used in original huMAP 1.0 study
                complexesData = {}
                complexesData_path = './sourceData/humap1/complexes/'
                complexesNames = ['train', 'test']
                for name in complexesNames:
                    complexesData[name] = pd.read_csv(complexesData_path + name + '_complexes.txt', sep='\t', header=None)
                humap1['complexes'] = complexesData

            if returnFeature_matrix:
                # acquire feature matrix
                humap1['featureMatrix'] = lowercaseStr(
                    pd.read_csv('./sourceData/humap1/pairs/'
                                'blake_bioplex_feature_revisitMerge_pairsOnly_preyMerge2_'
                                'heinCollapseMerge_pairsOnly_preyMerge2.txt.bz2',
                                compression='bz2', low_memory=False)
                )

            if returnResults:
                print('under development...check back later')

            return humap1

    def humap2(self, otherArgs):

        if len(otherArgs) != 5:
            print('Insufficient number of arguments passed. Must specify flags for pairs, complexes, and feature '
                  'matrix.')
            return

        returnFmt = otherArgs[0]
        returnPairs = otherArgs[1]
        returnComplexes = otherArgs[2]
        returnFeature_matrix = otherArgs[3]
        returnResults = otherArgs[4]
        if not any((returnFmt, returnPairs, returnComplexes, returnFeature_matrix, returnResults)):
            print('All flags specified False!')
            return

        elif returnFmt not in ['acc', 'geneid', 'genename']:
            print('Available formats include only \'acc\', \'geneid\', and \'genename\' (only for results)')
            return

        else:

            humap2 = {}
            if returnPairs and returnFmt in ['acc', 'geneid']:
                # acquire pairs used in original huMAP 2.0 study
                ppiData = {}
                ppiData_path = './sourceData/humap2/pairs/' + returnFmt + '/'
                prefix = 'humap2_'
                suffix = \
                    [ending for ending in ['_acc_20200818.txt', '_geneid_20200818.txt'] if returnFmt in ending][0]
                ppiData_names = ['neg_test_ppis', 'neg_train_ppis', 'test_ppis', 'train_ppis']
                for name in ppiData_names:
                    ppiData[name] = lowercaseStr(
                        pd.read_csv(ppiData_path + prefix + name + suffix, sep='\t',
                                    header=None, names=['idi', 'idii'], dtype='str')
                    )
                humap2['pairs'] = ppiData

            elif returnPairs and returnFmt not in ['acc', 'geneid']:
                print('huMAP 2.0 study did not provide pairs in the genename format.')
                return

            if returnComplexes and returnFmt in ['acc', 'geneid']:
                # acquire complexes used in original huMAP 2.0 study
                complexesData = {}
                complexesData_path = './sourceData/humap2/complexes/' + returnFmt + '/'
                prefix = 'humap2_'
                suffix = \
                    [ending for ending in ['_acc_20200818.txt', '_geneid_20200818.txt'] if returnFmt in ending][0]
                complexesNames = ['train_complexes', 'test_complexes']
                for name in complexesNames:
                    try:
                        DF = pd.read_csv(complexesData_path + prefix + name + suffix, sep=' ', header=None)

                    except:
                        DF = pd.read_csv(complexesData_path + prefix + name + suffix, sep='\t', header=None)
                    DF = DF.assign(temp0=DF[0].str.split('\t'))

                    for idx in DF.index:
                        DF.loc[idx, 'temp1'] = ' '.join(DF.loc[idx, 'temp0'])
                    DF.drop(columns=[0, 'temp0'], inplace=True)
                    DF.rename(columns={'temp1': 0}, inplace=True)
                    complexesData[name] = DF

                humap2['complexes'] = complexesData

            elif returnComplexes and returnFmt not in ['acc', 'geneid']:
                print('huMAP 2.0 study did not provide complexes in the genename format.')

            if returnFeature_matrix:
                # acquire feature matrix
                humap2['featureMatrix'] = lowercaseStr(
                    pd.read_csv('./sourceData/humap2/pairs/'
                                'orig9k_bioplex2_hygeo_bioid_hygeo_boldt_apms_hygeo_treiber_hygeo_wgt2_youn_hygeo_'
                                'trimCols_groupbyMean.featmat.bz2',
                                compression='bz2', low_memory=False)
                )

            return humap2

    # under development
    def humphreys2021(self, otherArgs):
        print('...under development.')

    # under development
    def kaggle(self, otherArgs):
        print('...under development.')

    def nci60(self, otherArgs):

        if len(otherArgs) != 1:
            print('Wrong number of arguments passed. Consider flags.')
            return

        else:
            maxGeneid_matches = otherArgs[0]

        nci60 = {}

        nci60_cellLine_proteomeData = lowercaseStr(
            pd.read_csv('./sourceData/nci60/nci60_all_abundanceLFQPerCellLineProteome.tsv', sep=',')
        )
        nci60['cellProtein_expression'] = nci60_cellLine_proteomeData

        nci60_tissueData = lowercaseStr(
            pd.read_csv('./sourceData/nci60/nci60_all_abundanceLFQPerTissueDeep.tsv', sep=',')
        )
        nci60['tissueProtein_expression'] = nci60_tissueData

        nci60_cellLine_allData = lowercaseStr(
            pd.read_csv('./sourceData/nci60/nci60_all_expressionPerCellLine.tsv', sep=',')
        )
        nci60['rnaExpression'] = nci60_cellLine_allData

        nci60['cellProt_geneids'] = nci60['cellProtein_expression'].geneid.to_list()
        nci60['cellProt_geneids'].remove(float(np.nan))
        nci60['cellProt_geneids'] = set(nci60['cellProt_geneids'])

        nci60['tissueProt_geneids'] = nci60['tissueProtein_expression'].geneid.to_list()
        nci60['tissueProt_geneids'].remove(float(np.nan))
        nci60['tissueProt_geneids'] = set(nci60['tissueProt_geneids'])

        nci60['rnaGeneids'] = set(nci60['rnaExpression'].geneid.to_list())

        if maxGeneid_matches:

            allUniprot_human = lowercaseStr(
                pd.read_csv('./sourceData/uniprot/uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_.tab', sep='\t')
            )
            allUniprot_human = allUniprot_human[['entry', 'gene names', 'cross-reference (geneid)']]
            allUniprot_human.rename(columns={'gene names': 'genename', 'cross-reference (geneid)': 'geneid'}, inplace=True)

            nci60_cellLine_allData = splitDF_byDelimiters(nci60_cellLine_allData,
                                                          ['genename', 'geneid'],
                                                          [';', ',', ' '])
            nci60['rnaExpressionExploded'] = nci60_cellLine_allData

            nci60_cellLine_proteomeData = nci60_cellLine_proteomeData.loc[nci60_cellLine_proteomeData.geneid.isna(), :]
            nci60_cellLine_proteomeData.drop(columns=['geneid'], inplace=True)
            nci60_cellLine_proteomeData = nci60_cellLine_proteomeData.merge(
                allUniprot_human.loc[allUniprot_human.geneid.notnull(), ['genename', 'geneid']],
                on=['genename'])

            nci60_tissueData = nci60_tissueData.loc[nci60_tissueData.geneid.isna(), :]
            nci60_tissueData.drop(columns=['geneid'], inplace=True)
            nci60_tissueData = nci60_tissueData.merge(
                allUniprot_human.loc[allUniprot_human.geneid.notnull(), ['genename', 'geneid']],
                on=['genename'])

            nci60_cellLine_allData = nci60_cellLine_allData.loc[nci60_cellLine_allData.geneid.isna(), :]
            nci60_cellLine_allData.drop(columns=['geneid'], inplace=True)
            nci60_cellLine_allData = nci60_cellLine_allData.merge(
                allUniprot_human.loc[allUniprot_human.geneid.notnull(), ['genename', 'geneid']],
                on=['genename'])

            nci60['cellProt_geneids'] = nci60['cellProt_geneids'].union(
                set(nci60_cellLine_proteomeData.geneid.to_list()))
            nci60['tissueProt_geneids'] = nci60['tissueProt_geneids'].union(
                set(nci60_tissueData.geneid.to_list()))
            nci60['rnaGeneids'] = nci60['rnaGeneids'].union(
                set(nci60_cellLine_allData.geneid.to_list()))

        nci60['protGeneids'] = nci60['cellProt_geneids'].union(nci60['tissueProt_geneids'])
        nci60['protPlus_rnaGeneids'] = nci60['rnaGeneids'].intersection(nci60['protGeneids'])

        return nci60

    # under development
    def properseq(self, otherArgs):
        print('...under development.')

    def scbc(self, otherArgs):

        if len(otherArgs) != 9:
            print('Incorrect number of arguments passed!')
            return

        reviewSCBC_flag = otherArgs[0]
        reviewMappings_flag = otherArgs[1]
        mappingFlag = otherArgs[2]
        mappingPath = otherArgs[3]
        mergeData_flag = otherArgs[4]
        dataPath = otherArgs[5]
        labelsPath = otherArgs[6]
        dataHeader_flag = otherArgs[7]
        labelsHeader_flag = otherArgs[8]

        output = {}
        #1: get scbc
        print('Making SCBC data dictionary...')
        print('     mapping genes in SCBC to IDs in JLM original data matrix...')
        corrMatsDict = pd.read_excel('./sourceData/scbc_quantitative_ms(S1).xlsx', sheet_name=None)
        if reviewSCBC_flag:
            output['origSCBC_data'] = corrMatsDict
        print('     completed mapping genes from SCBC to JLM original data matrix')
        print('     ________________________________________________')

        #2: make mappings
        genecard_to_uniprot = \
            lowercaseStr(
                pd.read_excel('./sourceData/scbc/all_SCBC_proteins_IDS_genecard-uniprotkbid.xls', sheet_name='Sheet0')
            )
        if reviewMappings_flag:
            output['genecard2uniprot_minProcessed_mapping'] = genecard_to_uniprot

        genecard_to_uniprot = \
            genecard_to_uniprot[["yourlist:m20200428a94466d2655679d1fd8953e075198da8845c47j", "entry"]]
        genecard_to_uniprot = genecard_to_uniprot.rename(
            columns={"yourlist:m20200428a94466d2655679d1fd8953e075198da8845c47j": "genecardid", "entry": "uniprot"})
        if reviewMappings_flag:
            output['genecard2Uniprot_mapping'] = genecard_to_uniprot

        uniprot_to_geneID = lowercaseStr(
            pd.read_excel('./scbc_uniprot_geneIDs.xls', sheet_name='scbc_uniprot->geneIDs')
        )
        uniprot_to_geneID['geneid'] = uniprot_to_geneID['geneid'].astype(str)
        if reviewMappings_flag:
            output['uniprot2Geneid_mapping'] = uniprot_to_geneID

        geneIDmapping = uniprot_to_geneID.merge(genecard_to_uniprot, how='inner', on="uniprot")
        geneIDmapping = geneIDmapping[["genecardid", "geneid"]]
        geneIDmapping = dict(zip(geneIDmapping["genecardid"], geneIDmapping["geneid"]))
        if reviewMappings_flag:
            output['geneidMapping'] = geneIDmapping

        #3: do mappings
        if mappingFlag or mappingPath:

            if mappingFlag and mappingPath:
                print('Alternate path to saved mapping passed. Reconsider duplicative effort as this a very expensive '
                      'process in time and memory.')
                return

            elif not mappingFlag and mappingPath:
                if ~exists(mappingPath):
                    print('Path to mapping is incorrect!')
                    return

                else:
                    corrMatsDict = pickle.load(open(mappingPath, 'rb'))

            elif mappingFlag and not mappingPath:
                if exists('./util/cell-line_specific_correlation_data_dict.pkl'):
                    corrMatsDict = pickle.load(open('./util/cell-line_specific_correlation_data_dict.pkl', 'rb'))
                    print('Loading existing mapping...')

                else:
                    for cellLine in corrMatsDict.keys():
                        print('Cell line: ' + cellLine)
                        placeholder = lowercaseStr(corrMatsDict[cellLine].iloc[:, :-1])
                        dataframeSummary(placeholder)
                        try:
                            placeholder.set_index('protein', inplace=True, verify_integrity=True)
                        except:
                            placeholder.drop_duplicates(subset='protein', keep=False, inplace=True)
                            placeholder.set_index('protein', inplace=True, verify_integrity=True)  # possibly problematic

                        placeholder = \
                            placeholder.T.corr().rename(columns=geneIDmapping,index=geneIDmapping
                                                        ).rename_axis(None).rename_axis(None, axis=1).stack().reset_index()
                        print('DF summary after correlation addition...')
                        dataframeSummary(placeholder)
                        corrMatsDict[cellLine] = placeholder.rename(columns={'level_0': 'geneid1', 'level_1': 'geneid2', 0: 'corr'})
                        corrMatsDict[cellLine]['geneid1'] = corrMatsDict[cellLine]['geneid2'].astype('str')
                        corrMatsDict[cellLine]['geneid1'] = corrMatsDict[cellLine]['geneid2'].astype('str')

                    print('     Confirming properties of SCBC cell line dictionary, particularly its cell line specific '
                          'correlation feature...')
                    for cellLine, data in corrMatsDict.items():
                        print('Cell line: ' + cellLine)
                        dataframeSummary(data)

                    print('     Saving cell line specific correlation data dictionary to disk...')
                    pickle.dump(corrMatsDict, open('./util/cell-line_specific_correlation_data_dict.pkl', 'wb'))
                    print('     Completed saving cell line specific correlation data dictionary to disk...')
                    print('Completed making SCBC data dictionary...')

            #4: map scbc to input data
            if mergeData_flag and dataPath and labelsPath and dataHeader_flag and labelsHeader_flag:
                data = lowercaseStr(pd.read_csv(dataPath, sep='\t', header=dataHeader_flag))
                dataframeSummary(data)
                labels = lowercaseStr(pd.read_csv(labelsPath, sep='\t', header=labelsHeader_flag))
                dataframeSummary(labels)

                print('Merging data with corresponding IDs...')
                ppiIDs = labels.iloc[:, 1:].astype(str)
                data = pd.concat([ppiIDs, data], axis=1)
                dataframeSummary(data)

                print('Merging data with SCBC features...')
                for cellLine in corrMatsDict.keys():
                    print('Cell line: ' + cellLine)
                    data = data.merge(corrMatsDict[cellLine].astype(str),
                                      how='left',
                                      left_on=['idi', 'idii'],
                                      right_on=['geneid1', 'geneid2'])
                    data = data.rename(columns={'corr': 'corr_' + cellLine})
                    dataframeSummary(data)
                print('Completed merging data with SCBC features')
                print('________________________________________________')
                output['dataSCBC_merged'] = data

            else:
                print('In order to merge, the path to the data and labels, as well as the flags specifying whether '
                      'each respectively has headers must also be passed.')
                return

        if mergeData_flag:
            print('Must pass either flag to create and/or load mapping or an alternate path to mapping.')
            return

    # under development
    def string(self, otherArgs):
        stringHuman_ppi = lowercaseStr(
            pd.read_csv('./sourceData/string/results/9606.protein.links.full.v11.5.txt', sep=' ')
        )
        stringHuman_ppi.rename(
            columns={'protein1': 'stringIDi', 'protein2': 'stringIDii', 'combined_score': 'combinedScore'},
            inplace=True)
        stringHuman_ppi['stringIDi'] = stringHuman_ppi['stringIDi'].str.replace('9606.', '')
        stringHuman_ppi['stringIDii'] = stringHuman_ppi['stringIDii'].str.replace('9606.', '')
        stringHuman_ppi = stringHuman_ppi.assign(combinedScore_frac=stringHuman_ppi['combinedScore'] / 1000)
        stringHuman_ppi = stringHuman_ppi[['stringIDi', 'stringIDii', 'combinedScore', 'combinedScore_frac']]

        return stringHuman_ppi

    def uniprot(self, otherArgs):

        if len(otherArgs)==10:
            indEntry_flag = otherArgs[0]
            removeEmpty_geneidFlag = otherArgs[1]
            hpaData_flag = otherArgs[2]
            dropEmpty_locsFlag = otherArgs[3]
            hpaSeparate_geneidSingle_entriesFlag = otherArgs[4]
            hpaSeparate_genenameSingle_entriesFlag = otherArgs[5]
            hpaMerge_geneidUniprot_flag = otherArgs[6]
            hpaMerge_geneidGenename_flag = otherArgs[7]
            hpaMerge_geneidUniprot_plusGenename_flag = otherArgs[8]
            hpaOutput_proteinsFlag = otherArgs[9]
        else:
            print('Incorrect number of arguments passed.')
            return

        results = {}
        allUniprot_human = lowercaseStr(
            pd.read_csv('./sourceData/uniprot/uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_.tab', sep='\t')
        )
        allUniprot_human.rename(columns={'cross-reference (geneid)': 'geneid', 'gene names': 'genename'}, inplace=True)
        allUniprot_human = allUniprot_human[['entry', 'genename', 'subcellular location [cc]', 'geneid']]
        results['uniprotHuman'] = allUniprot_human

        if indEntry_flag:
            allUniprot_human = allUniprot_human.assign(
                genename=allUniprot_human.genename.str.split(';')).explode('genename')
            allUniprot_human = allUniprot_human.assign(
                genename=allUniprot_human.genename.str.split(' ')).explode('genename')
            allUniprot_human = allUniprot_human.assign(
                geneid=allUniprot_human.geneid.str.split(';')).explode('geneid')
            allUniprot_human['geneid'].replace('', np.nan, inplace=True)
            results['uniprotHuman_splitMultiple_entries'] = allUniprot_human

        if removeEmpty_geneidFlag:
            allUniprot_human.dropna(subset=['geneid'], inplace=True)
            results['uniprotHuman_removeEmpty_geneids'] = allUniprot_human

        if hpaData_flag:
            swissProt_data = lowercaseStr(
                pd.read_csv('sourceData/hpa/uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_+AND+review--.tab',
                            sep='\t')
            )
            '''
            swissProt_data = lowercaseStr(
                pd.read_csv('./sourceData/hpa/swiss_prot_subcellular_location.tsv', sep='\t')
            )
            '''

            swissProt_data.rename(
                columns={'cross-reference (geneid)': 'geneid', 'gene names': 'genename'},
                inplace=True)
            results['rawHPA'] = swissProt_data

            if hpaOutput_proteinsFlag:
                results['hpaProteins'] = set(swissProt_data.geneid.to_list())

            if dropEmpty_locsFlag:
                swissProt_data.dropna(subset=['subcellular location [cc]'], inplace=True)
                results['hpaLocs_allPresent'] = swissProt_data
                if hpaOutput_proteinsFlag:
                    results['hpaProteins'] = set(swissProt_data.geneid.to_list())

            if hpaSeparate_geneidSingle_entriesFlag:
                swissProt_data = swissProt_data.loc[swissProt_data.geneid.notnull(), :]
                print('data dimensions (empty GeneIDs removed): ' + str(swissProt_data.shape))
                swissProt_data = swissProt_data.assign(geneid=swissProt_data.geneid.str.split(';')).explode('geneid')

                print('data dimensions (GeneIDs exploded by delimiter): ' + str(swissProt_data.shape))
                swissProt_data['geneid'].replace('', np.nan, inplace=True)
                swissProt_data.dropna(axis=0, subset=['geneid'], inplace=True)
                print('data dimensions (GeneIDs exploded by delimiter, nulls addressed...): ' + str(swissProt_data.shape))
                print('# unique GeneIDs (from nonEmpty GeneID partition): ' + str(len(swissProt_data.geneid.unique())))
                results['hpaProcessed_geneidsSeparated_entries'] = swissProt_data
                if hpaOutput_proteinsFlag:
                    results['hpaProteins'] = set(swissProt_data.geneid.to_list())

            if hpaSeparate_genenameSingle_entriesFlag:
                swissProt_data = swissProt_data.assign(genename=swissProt_data.genename.str.split(';')).explode('genename')
                print(swissProt_data.shape)
                swissProt_data = swissProt_data.assign(genename=swissProt_data.genename.str.split(' ')).explode('genename')
                print(swissProt_data.shape)
                print('data dimensions (only empty GeneIDs, GeneNames exploded by delimiter): ' + str(swissProt_data.shape))
                swissProt_data.drop(columns=['geneid'], inplace=True)
                print(swissProt_data.shape)
                results['hpaProcessed_genenamesSeparated_entries'] = swissProt_data

            if hpaMerge_geneidUniprot_flag:
                results['hpaMerged_geneidsUniprot'] = swissProt_data.merge(
                    allUniprot_human.loc[allUniprot_human.geneid.notnull(), ['entry', 'genename', 'geneid']],
                    on=['entry'])
                results['hpaMerged_geneidsUniprot'].rename(columns={'geneid_y': 'geneid'}, inplace=True)
                if hpaOutput_proteinsFlag:
                    results['hpaProteins'] = results['hpaProteins'].union(
                        set(results['hpaMerged_geneidsUniprot'].geneid.to_list()))

            if hpaMerge_geneidGenename_flag:
                results['hpaMerged_geneidsGenenames'] = swissProt_data.merge(
                    allUniprot_human.loc[allUniprot_human.geneid.notnull(), ['entry', 'genename', 'geneid']],
                    on=['genename'])
                results['hpaMerged_geneidsGenenames'].rename(columns={'geneid_y': 'geneid'}, inplace=True)
                if hpaOutput_proteinsFlag:
                    results['hpaProteins'] = results['hpaProteins'].union(
                        set(results['hpaMerged_geneidsGenenames'].geneid.to_list()))

            if hpaMerge_geneidUniprot_plusGenename_flag:
                results['hpaMerged_geneidsUniprot_plusGenenames'] = swissProt_data.merge(
                    allUniprot_human.loc[allUniprot_human.geneid.notnull(), ['entry', 'genename', 'geneid']],
                    on=['entry', 'genename'])
                results['hpaMerged_geneidsUniprot_plusGenenames'].rename(columns={'geneid_y': 'geneid'}, inplace=True)
                if hpaOutput_proteinsFlag:
                    results['hpaProteins'] = results['hpaProteins'].union(
                        set(results['hpaMerged_geneidsUniprot_plusGenenames'].geneid.to_list()))

        return results

    def wan2015(self, otherArgs):
        wan2015_data = {}
        wan2015_proteinInteractions = lowercaseStr(
            pd.read_csv('./sourceData/wan2015/BIOGRID-PUBLICATION-185267-4.4.204.tab3.txt', sep='\t')
        )
        wan2015_proteinInteractions.rename(columns={'entrez gene interactor a': 'geneid1',
                                                    'entrez gene interactor b': 'geneid2'}, inplace=True)
        wan2015_data['interactions'] = wan2015_proteinInteractions

        if len(otherArgs) == 1:
            wan2015_interactionsSummary_geneidsList = set(wan2015_proteinInteractions['geneid1'].to_list()).union(
                set(wan2015_proteinInteractions['geneid2'].to_list()))

            wan2015_data['proteins'] = set(
                [str(entry) for entry in wan2015_interactionsSummary_geneidsList]
                )

        return wan2015_data

    def dummy(self, otherArgs):
        print('dummy function passed ' + str(len(otherArgs)) + ' arguments.')


if __name__ == '__main__':

    if len(sys.argv) >= 3:
        userInput = sys.argv[1:]
        LoadDict().source(userInput)
    else:
        print('insufficient number of arguments')

def addSCBC_feats(baseDF, scbcPickle_path, outputPath, cellSpecific_flag=False):
     fractions = ['fs1', 'fs2', 'fp1', 'fp2', 'fp3']
     cells = ['A431', 'MCF7', 'H322', 'HCC827', 'U251'] if cellSpecific_flag else ['nonspecific']
     appendedDF = dict()
     for cell in cells:
         print('beginning scbc correlation merger to base data with all fractions for ' + cell + '...')
         appendedDF[cell] = baseDF
         for frac in fractions:
             print('load scbc correlation pickle...')
             print(scbcPickle_path+cell+'-'+frac+'_noDups.pkl')
             scbcPickle = p.load(open(scbcPickle_path+cell+'-'+frac+'_noDups.pkl', 'rb'))
             print('completed loading scbc correlation pickle.')
             print('merging scbc correlation by fraction ' + frac + '...')
             appendedDF[cell] = appendedDF[cell].merge(scbcPickle.drop(columns=['idi', 'idii']), how='left', on=['pairsFrozen'])
             appendedDF[cell].rename(columns={'corr': 'SCBC_'+cell+'-'+frac+'corr'}, inplace=True)
             print('completed merging scbc correlation by fraction ' + frac + '...')
             print('columns length: ' + str(len(appendedDF[cell].columns)))
         appendedDF[cell].drop(columns=['pairsFrozen'], inplace=True)
         appendedDF[cell].to_csv(outputPath+'lmFeats_newFeats_'+cell+'SCBC_corrAppended.tsv', sep='\t', header=False)
         print('completed scbc correlation merger to base data with all fractions for ' + cell + '.')
     return appendedDF