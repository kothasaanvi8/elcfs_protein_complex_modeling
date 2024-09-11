import multiprocessing as mp
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import hmean
import seaborn as sns
import sys
sys.path.append('../../protein_complex_maps_public/')
sys.path.append('./../proteinPairs_complexMaps')

#import protein_complex_maps.complex_comparison as cc
from protein_complex_maps import complex_comparison as cc


def plotF1_vClique_evaluations(referenceComplexes,
                               predictionsFilesnames,
                               predictionsFilesnames_short,
                               outputFilename,
                               procs=1):

    excluded_complexes = []
    compare2goldstandard_input_list = []
    for i, cluster_filename in enumerate(predictionsFilesnames):
        parameter_dict = dict()
        parameter_dict['cluster_filename'] = cluster_filename
        parameter_dict['gs_boot_iteration'] = 1
        parameter_dict['gs_complexes'] = referenceComplexes
        parameter_dict['ex_complexes'] = excluded_complexes
        parameter_dict['predictionSource'] = predictionsFilesnames_short[i]
        parameter_dict['samples'] = 10000
        parameter_dict['exact'] = False
        parameter_dict['max_clique'] = None
        parameter_dict['pseudocount'] = 1
        parameter_dict['normalize_by_combinations'] = False
        compare2goldstandard_input_list.append(parameter_dict)

    p = mp.Pool(procs)
    compare2goldstandard_results = p.map(compare2goldstandard, compare2goldstandard_input_list)

    precision_dict = dict()
    recall_dict = dict()
    f1score_dict = dict()
    grand_f1score_dict = dict()
    mean_f1score_dict = dict()
    mean_f1score_weighted_dict = dict()
    cumulative_precision_dict = dict()
    cumulative_recall_dict = dict()
    numOfClusters_dict = dict()
    for result in compare2goldstandard_results:
        precision_dict[result['id']] = result['precision_list']
        recall_dict[result['id']] = result['recall_list']
        f1score_dict[result['id']] = result['f1score_list']
        grand_f1score_dict[result['id']] = result['grand_f1score']
        mean_f1score_dict[result['id']] = result['mean_f1score']
        mean_f1score_weighted_dict[result['id']] = result['mean_f1score_weighted']
        cumulative_precision_dict[result['id']] = result['cumulative_precision_list']
        cumulative_recall_dict[result['id']] = result['cumulative_recall_list']
        numOfClusters_dict[result['id']] = result['numOfClusters']

    pr_dict = dict()
    pr_dict['precision'] = precision_dict
    pr_dict['recall'] = recall_dict
    pr_dict['f1score'] = f1score_dict
    pr_dict['grand_f1score'] = grand_f1score_dict
    pr_dict['mean_f1score'] = mean_f1score_dict
    pr_dict['mean_f1score_weighted'] = mean_f1score_weighted_dict
    pr_dict['cumulative_precision'] = cumulative_precision_dict
    pr_dict['cumulative_recall'] = cumulative_recall_dict
    pr_dict['numOfClusters'] = numOfClusters_dict

    #generate table of k-clique based F scores for range of cliques' lengths
    fscores_cliqueLength_bySource = dict()
    fgrandScores = []
    f1meanScores = []
    f1weightedMean_scores = []
    for id in pr_dict['f1score'].keys():
        sourceName_dummy = [id for ele in list(pr_dict['f1score'][id].keys())]
        fGrand_dummy = [list(pr_dict['grand_f1score'][id].values())[0] for ele in list(pr_dict['f1score'][id].keys())]
        fgrandScores.append(fGrand_dummy[0])
        f1mean_dummy = [list(pr_dict['mean_f1score'][id].values())[0] for ele in list(pr_dict['f1score'][id].keys())]
        f1meanScores.append(f1mean_dummy[0])
        f1mean_weightedDummy = \
            [list(pr_dict['mean_f1score_weighted'][id].values())[0] for ele in list(pr_dict['f1score'][id].keys())]
        f1weightedMean_scores.append(f1mean_weightedDummy[0])
        fscores_cliqueLength_bySource[id] = pd.DataFrame(
            zip(
                sourceName_dummy,
                list(range(2, len(pr_dict['f1score'][id]) + 2)),
                pr_dict['f1score'][id],
                fGrand_dummy,
                f1mean_dummy,
                f1mean_weightedDummy
            ),
        columns=['Source', 'Clique Length', 'F1', 'F1-Grand', 'F1-Mean', 'Weighted F1-Mean'])

    fscores_cliqueLength_sourceTable = pd.concat(list(fscores_cliqueLength_bySource.values()), ignore_index=True)
    fscores_cliqueLength_sourceTable.to_csv('./f1Scores_cliqueSize_bySource.csv', index=False, sep='\t')
    sources = list(fscores_cliqueLength_sourceTable.Source.unique())
    f1Grand_stats = pd.DataFrame(list(zip(sources, fgrandScores)), columns=['ID', 'score'])
    f1Mean_stats = pd.DataFrame(list(zip(sources, f1meanScores)), columns=['ID', 'score'])
    f1Mean_weightedStats = pd.DataFrame(list(zip(sources, f1weightedMean_scores)), columns=['ID', 'score'])

    plt.clf()
    plt.title('F1-grand, F1-mean, and weighted F1-mean by source.')
    sns.scatterplot(data=f1Grand_stats, x="ID", y="score", label="F1-Grand")
    sns.scatterplot(data=f1Mean_stats, x="ID", y="score", label="Mean F1")
    sns.scatterplot(data=f1Mean_weightedStats, x="ID", y="score", label="Mean F1 Weighted")
    plt.legend(title='F1 Measures')
    plt.savefig('./f1Measures_bySource.svg', bbox_inches='tight', dpi=300)

    plotF1score_vs_cliqueSize(pr_dict, predictionsFilesnames_short, plot_filename=outputFilename+'.svg')

def compare2goldstandard(parameter_dict):
    cluster_filename = parameter_dict['cluster_filename']
    id = parameter_dict['predictionSource']
    gs_complexes = parameter_dict['gs_complexes']
    boot_iteration = parameter_dict['gs_boot_iteration']
    ex_complexes = parameter_dict['ex_complexes']
    samples = parameter_dict['samples']
    exact = parameter_dict['exact']
    max_clique = parameter_dict['max_clique']
    pseudocount = parameter_dict['pseudocount']
    normalize_by_combinations = parameter_dict['normalize_by_combinations']

    predicted_clusters = []
    clpred_f = open(cluster_filename, "rb")
    for line in clpred_f.readlines():
        predicted_clusters.append(line.split())
    clpred_f.close()

    cplx_compare = cc.ComplexComparison(gs_complexes, predicted_clusters, exclusion_complexes=ex_complexes,
                                        samples=samples, exact=exact, max_clique=max_clique,
                                        pseudocount=pseudocount,
                                        normalize_by_combinations=normalize_by_combinations)
    result_dict = cplx_compare.clique_comparison_metric()
    precision_list = [result_dict[x]['precision'] for x in result_dict.keys()]
    recall_list = [result_dict[x]['recall'] for x in result_dict.keys()]
    f1score_list = [result_dict[x]['f1score'] for x in result_dict.keys()]
    cumulative_precision_list = [result_dict[x]['cumulative_precision'] for x in result_dict.keys()]
    cumulative_recall_list = [result_dict[x]['cumulative_recall'] for x in result_dict.keys()]
    numOfClusters_list = [result_dict[x]['numOfClusters'] for x in result_dict.keys()]

    grand_f1score = cplx_compare.clique_comparison_metric_grandf1score()
    ccmm = cplx_compare.clique_comparison_metric_mean()
    wccmm = cplx_compare.clique_comparison_metric_mean(weighted=True)
    mean_f1score = hmean([ccmm['precision_mean'], ccmm['recall_mean']])
    mean_f1score_weighted = hmean([wccmm['precision_mean'], wccmm['recall_mean']])

    return {'id': id,
            'boot_iteration': boot_iteration,
            'precision_list': precision_list, 'recall_list': recall_list, 'f1score_list': f1score_list,
            'grand_f1score': grand_f1score, 'mean_f1score': mean_f1score, 'mean_f1score_weighted': mean_f1score_weighted,
            'cumulative_precision_list': cumulative_precision_list,
            'cumulative_recall_list': cumulative_recall_list,
            'numOfClusters': numOfClusters_list}

def plotF1score_vs_cliqueSize(pr_dict, psweep_indices, plot_filename=None):

    if plot_filename:
        for id in psweep_indices:
            f1score_list = pr_dict['f1score'][id]
            clique_list = range(2, len(f1score_list) + 2)
            plt.plot(clique_list, f1score_list, '.-', linewidth=1, label=id)

        plt.clf()
        plt.legend(loc="upper right")
        plt.xlabel('Clique Size')
        plt.ylabel('F1-score')
        plt.savefig(plot_filename+'.svg', bbox_inches='tight', dpi=300, format='svg')


# for hu.MAP 1.0
# referenceComplexes_filename = '../sourceData/humap1/complexes/test_complexes.txt'

# for hu.MAP 2.0 (geneid)
'''
referenceComplexes_filename = \
   '../sourceData/humap2/complexes/geneid/humap2_test_complexes_geneid_20200818.txt'
'''
# for hu.MAP 2.0 (acc)
'''
referenceComplexes_filename = \
   '../sourceData/humap2/complexes/acc/humap2_test_complexes_ACC_20200818.txt'
'''

# for CORUM (latest 2018)
'''
referenceComplexes_filename = '../allCORUM.csv'

predictionsFilesnames = ['../drew2017_results/clusters.txt',
                         '../lmPredicted_clusters0.02%.csv',
                         '../lmPredicted_clusters0.05%.csv']

predictionsFilesnames_short = ['Drew 2017',
                               'LM 2019 Top 0.02%',
                               'LM 2019 Top 0.05%']

outputFilename = '../F1vClique_clusterPredictions_comparisonEval_Drew2017+LM_bothThreshs_ref=CORUM2018'

plotF1_vClique_evaluations(referenceComplexes_filename,
                           predictionsFilesnames,
                           predictionsFilesnames_short,
                           outputFilename)
'''
