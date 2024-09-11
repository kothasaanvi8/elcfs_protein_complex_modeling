import argparse
import itertools as it
import multiprocessing as mp
import networkx as nx
import os, pickle
from scipy.stats import hmean
import sys
import random
import subprocess as sp
import tempfile as tf

import protein_complex_maps.evaluation.complex_comparison_py3 as cc


def main(scoredPairs_filename, goldStandard_filename, tempDir):
    parser = argparse.ArgumentParser(description="Finds optimal parameters for clustering of ppi network")
    parser.add_argument("--input_network", action="store", dest="input_network", required=True,
                        help="Filename of ppi network with optional edge weights (format: id\tid\tweight)")
    parser.add_argument("--gold_standard", action="store", dest="gold_standard_filename", required=True,
                        help="Filename of gold standard complexes, one cluster per line, ids tab-separated")
    parser.add_argument("--output_fle", action="store", dest="output_file", required=False,
                        default='./drive/My Drive/otherStudies/output_file',
                        help="Filename of output clusters for best set of parameters")
    parser.add_argument("--random_seed", action="store", type=int, dest="random_seed", required=False, default=42,
                        help="Sets random seed (int), default=None")
    parser.add_argument("--bootstrap_iter", action="store", type=int, dest="bootstrap_iter", required=False, default=10,
                        help="Number of bootstrap iterations (int, default=10)")
    parser.add_argument("--bootstrap_fraction", action="store", type=float, dest="bootstrap_fraction", required=False,
                        default=0.5,
                        help="Fraction of edges to sample for bootstrapping, default=0.5")
    parser.add_argument("--clusterone", action="store", dest="clustone_jar", required=False,
                        default="./drive/My Drive/otherStudies/cluster_one-1.0.jar",
                        help="Location of cluster_one jar file, default= ./cluster_one-1.0.jar")
    parser.add_argument("--clusterone_size", action="store", dest="clusterone_size", nargs='+', required=False,
                        default=[2, 3, 4, 5],
                        help="ClusterOne Size parameter sweep, default = 2")
    parser.add_argument("--clusterone_density", action="store", dest="clusterone_density", nargs='+', required=False,
                        default=[.1, .2, .3, .4, .5, .6, .7, .8, .9],
                        help="ClusterOne Density parameter sweep, default = .1 .25 .3 .35")
    parser.add_argument("--clusterone_max_overlap", action="store", dest="clusterone_max_overlap", nargs='+',
                        required=False,
                        default=[.1, .2, .3, .4, .5, .6, .7, .8, .9],
                        help="ClusterOne max-overlap parameter sweep, default = 0.5 0.75 0.9")
    parser.add_argument("--clusterone_seed_method", action="store", dest="clusterone_seed_method", nargs='+',
                        required=False,
                        default=['nodes', 'cliques', 'unused_nodes', 'edges'],
                        help="ClusterOne seed method parameter sweep "
                             "(nodes, cliques, unused_nodes, edges, default = nodes")
    parser.add_argument("--ppi_fraction", action="store", dest="ppi_fraction", nargs='+', required=False,
                        default=[0.005, 0.01, 0.05, .1, .25, .5, .75, 1.0],
                        help="Use top fraction for further clustering, default = 0.005 0.01 0.05 .1 .25 .5 .75 1.0")
    parser.add_argument("--ppi_threshold_score", action="store", dest="ppi_threshold_score", nargs='+', required=False,
                        default=[1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1],
                        help="Use ppis with score higher or equal to threshold score, "
                             "default = 1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1")
    parser.add_argument("--clixo_bin", action="store", dest="clixo_bin", required=False,
                        default='./drive/My Drive/otherStudies/AIM_H_SAPIENS/clixo_0.3/clixo')
    parser.add_argument("--clixo_alpha", action="store", dest="clixo_alpha", nargs='+', required=False,
                        default=[None],
                        help="Clixo Alpha parameter, default = None")
    parser.add_argument("--clixo_beta", action="store", dest="clixo_beta", nargs='+', required=False,
                        default=[None],
                        help="Clixo Beta parameter, default = None")
    parser.add_argument("--mcl", action="store", dest="mcl_bin", required=False,
                        default='mcl',
                        help="Location of mcl binary, default = 'mcl' ")
    parser.add_argument("--mcl_inflation", action="store", dest="mcl_inflation", nargs='+', required=False,
                        default=[None],
                        help="MCL Inflation (-I) parameter, "
                             "default = [None] (no 2-stage clustering), docs suggest = 1.2 - 5.0")
    parser.add_argument("--cfinder", action="store", dest="cfinder_exe", required=False,
                        default="./drive/My Drive/otherStudies/CFinder-2.0.6--1448/CFinder_commandline64",
                        help="Location of CFinder executable, default= ./CFinder-2.0.6--1448/CFinder_commandline64")
    parser.add_argument("--cfinder_license", action="store", dest="cfinder_license", required=False,
                        default="./drive/My Drive/otherStudies/CFinder-2.0.6--1448/licence.txt",
                        help="Location of CFinder license, default= ./CFinder-2.0.6--1448/licence.txt")
    parser.add_argument("--cfinder_cliquesize", action="store", dest="cfinder_cliquesize", nargs='+', required=False,
                        default=[3],
                        help="Cfinder clique size (-k) parameter, "
                             "default = None (use CFinder's default setting, recommended: 3)")
    parser.add_argument("--cfinder_timeout", action="store", dest="cfinder_timeout", nargs='+', required=False,
                        default=[10],
                        help="Cfinder timeout (-t) parameter, "
                             "default = None (use CFinder's default setting, recommended: 10)")
    parser.add_argument("--trim2threshold", action="store_true", dest="trim2threshold", required=False, default=False,
                        help="Trim final clusters of subunits that do not have an edge that passes the threshold_score, "
                             "default=False")
    parser.add_argument("--twostep_combination", action="store", dest="twostep_combination", nargs='+', required=False,
                        default=['clusterone', 'mcl'],
                        help="Combination of two step clustering, "
                             "default = [clusterone,mcl], options=[clusterone,mcl,cfinder,agglomod,clixo]")
    parser.add_argument("--eval_metric", action="store", dest="eval_metric", required=False, default='mmr',
                        help="Evaluation metric used to determine best set of parameters "
                             "(mmr, acc, sensitivity, ppv, clique_precision_mean, clique_recall_mean), default=mmr")
    parser.add_argument("--output_all", action="store_true", dest="output_all", required=False, default=True,
                        help="Output all clusterings, default=False")
    parser.add_argument("--procs", action="store", type=int, dest="procs", required=False, default=9,
                        help="Number processors to use (int), default=1)")
    parser.add_argument("--temp_dir", action="store", dest="temp_dir", required=False, default=None,
                        help="Where to store temporary files generated during processing, "
                             "default=None (defaults to OS level tmp), in memory suggestion = /dev/shm/")
    parser.add_argument("--nodelete", action="store_true", dest="nodelete", required=False, default=True,
                        help="When set, does not delete temporary files")

    args = \
        parser.parse_args(
            ["--input_network", scoredPairs_filename, "--gold_standard", goldStandard_filename,
             "--temp_dir", tempDir]
        )

    if not os.path.isfile(args.clustone_jar):
        raise IOError("File not found: %s" % args.clustone_jar)

    gstd_file = open(args.gold_standard_filename, "rb")
    gold_standard_complexes = []
    for line in gstd_file.readlines():
        gold_standard_complexes.append(line.split())

    with open(args.input_network, "r") as input_network_file:
        input_network_list = input_network_file.readlines()

    ppi_scores = dict()
    for ppi in input_network_list:
        ppi_scores[frozenset([ppi.split()[0], ppi.split()[1]])] = float(ppi.split()[2])

    random.seed(args.random_seed)
    size_sweep = args.clusterone_size
    overlap_sweep = args.clusterone_max_overlap
    seed_method_sweep = args.clusterone_seed_method
    density_sweep = args.clusterone_density
    fraction_sweep = args.ppi_fraction
    threshold_score_sweep = args.ppi_threshold_score
    inflation_sweep = args.mcl_inflation
    clixo_alpha_sweep = args.clixo_alpha
    clixo_beta_sweep = args.clixo_beta
    cliquesize_sweep = args.cfinder_cliquesize
    timeout_sweep = args.cfinder_timeout

    best_size = None
    best_ii = None
    best_density = None
    best_fraction = None
    best_threshold_score = None
    best_overlap = None
    best_inflation = None
    best_eval = None

    p = mp.Pool(args.procs)
    network_input_list = []
    for ii, parameters in enumerate(
            it.product(size_sweep, density_sweep, fraction_sweep, overlap_sweep, seed_method_sweep, inflation_sweep,
                       clixo_alpha_sweep, clixo_beta_sweep, cliquesize_sweep, timeout_sweep, threshold_score_sweep)):
        sys.stdout.flush()
        size, density, fraction, overlap, seed_method, inflation, clixo_alpha, clixo_beta, cliquesize, timeout, \
            threshold_score = parameters
        threshold_score = float(threshold_score)

        if threshold_score is None:
            sizeOfTopNetwork = int(len(input_network_list) * float(fraction))
            network_list = input_network_list[0:sizeOfTopNetwork]
            threshold_score = float(network_list[-1].split()[2])
        else:
            network_list = []
            for x in input_network_list:
                if float(x.split()[2]) >= float(threshold_score):
                    network_list.append(x)
                else:
                    break
            fraction = 1.0 * len(network_list) / len(input_network_list)

        parameter_dict = dict()
        parameter_dict['network_list'] = network_list
        parameter_dict['ppi_scores'] = ppi_scores
        parameter_dict['args'] = args
        parameter_dict['size'] = str(size)
        parameter_dict['density'] = str(density)
        parameter_dict['overlap'] = str(overlap)
        parameter_dict['seed_method'] = str(seed_method)
        parameter_dict['fraction'] = str(fraction)
        parameter_dict['threshold_score'] = threshold_score
        parameter_dict['inflation'] = str(inflation)
        parameter_dict['clixo_alpha'] = str(clixo_alpha)
        parameter_dict['clixo_beta'] = str(clixo_beta)
        parameter_dict['cliquesize'] = str(cliquesize)
        parameter_dict['timeout'] = str(timeout)
        parameter_dict['twostep_combination'] = args.twostep_combination
        parameter_dict['trim2threshold'] = args.trim2threshold
        parameter_dict['i'] = ii
        parameter_dict['nodelete'] = args.nodelete
        network_input_list.append(parameter_dict)

    cluster_predictions = p.map(cluster_helper, network_input_list)
    for cluster_prediction, ii in cluster_predictions:
        network_list = network_input_list[ii]['network_list']
        size = network_input_list[ii]['size']
        density = network_input_list[ii]['density']
        overlap = network_input_list[ii]['overlap']
        seed_method = network_input_list[ii]['seed_method']
        fraction = network_input_list[ii]['fraction']
        threshold_score = network_input_list[ii]['threshold_score']
        inflation = network_input_list[ii]['inflation']
        clixo_alpha = network_input_list[ii]['clixo_alpha']
        clixo_beta = network_input_list[ii]['clixo_beta']
        timeout = network_input_list[ii]['timeout']
        cliquesize = network_input_list[ii]['cliquesize']
        twostep_combination = network_input_list[ii]['twostep_combination']
        trim2threshold = network_input_list[ii]['trim2threshold']
        nodelete = network_input_list[ii]['nodelete']
        cplx_comparison = cc.ComplexComparison(gold_standard_complexes, cluster_prediction)
        cplx_comparison_normalize = cc.ComplexComparison(gold_standard_complexes, cluster_prediction,
                                                         normalize_by_combinations=True, pseudocount=0.00001)

        try:
            metric_dict = dict()
            metric_dict['acc'] = cplx_comparison.acc()
            metric_dict['sensitivity'] = cplx_comparison.sensitivity()
            metric_dict['ppv'] = cplx_comparison.ppv()
            metric_dict['mmr'] = cplx_comparison.mmr()
            metric_dict['precision_recall_product'] = cplx_comparison.precision_recall_product()
            ccmm = cplx_comparison.clique_comparison_metric_mean()
            metric_dict['clique_precision_mean'] = ccmm['precision_mean']
            metric_dict['clique_recall_mean'] = ccmm['recall_mean']
            wccmm = cplx_comparison.clique_comparison_metric_mean(weighted=True)
            metric_dict['clique_weighted_pr_mean'] = wccmm['precision_mean']
            metric_dict['clique_weighted_re_mean'] = wccmm['recall_mean']
            metric_dict['clique_weighted_hmean'] = hmean([wccmm['precision_mean'], wccmm['recall_mean']])
            ccmm_normalize = cplx_comparison_normalize.clique_comparison_metric_mean()
            metric_dict['clique_precision_mean_normalize'] = ccmm_normalize['precision_mean']
            metric_dict['clique_recall_mean_normalize'] = ccmm_normalize['recall_mean']
            ccmm_normalize_weighted = cplx_comparison_normalize.clique_comparison_metric_mean(weighted=True)
            metric_dict['clique_precision_mean_normalize_weighted'] = ccmm_normalize_weighted['precision_mean']
            metric_dict['clique_recall_mean_normalize_weighted'] = ccmm_normalize_weighted['recall_mean']
            metric_dict['clique_hmean_mean_normalize_weighted'] = hmean(
                [ccmm_normalize_weighted['precision_mean'], ccmm_normalize_weighted['recall_mean']])

        except Exception as e:
            print(e)
            continue

        sizeOfBootstrap = int(len(network_list) * args.bootstrap_fraction)
        shuffled_networks = [sorted(network_list, key=lambda k: random.random()) for i in range(args.bootstrap_iter)]
        bootstrapped_networks = [sn[:sizeOfBootstrap] for sn in shuffled_networks]
        bootstrapped_test_networks = [sn[sizeOfBootstrap:] for sn in shuffled_networks]

        multiproc_input = []
        for i, boot_net in enumerate(bootstrapped_networks):
            parameter_dict = dict()
            parameter_dict['ppi_scores'] = ppi_scores
            parameter_dict['network_list'] = boot_net
            parameter_dict['args'] = args
            parameter_dict['size'] = str(size)
            parameter_dict['density'] = str(density)
            parameter_dict['overlap'] = str(overlap)
            parameter_dict['seed_method'] = str(seed_method)
            parameter_dict['fraction'] = str(fraction)
            parameter_dict['threshold_score'] = threshold_score
            parameter_dict['inflation'] = str(inflation)
            parameter_dict['clixo_alpha'] = str(clixo_alpha)
            parameter_dict['clixo_beta'] = str(clixo_beta)
            parameter_dict['timeout'] = str(timeout)
            parameter_dict['cliquesize'] = str(cliquesize)
            parameter_dict['twostep_combination'] = twostep_combination
            parameter_dict['trim2threshold'] = trim2threshold
            parameter_dict['i'] = i
            parameter_dict['nodelete'] = nodelete
            multiproc_input.append(parameter_dict)
        bootstrapped_cluster_predictions = p.map(cluster_helper, multiproc_input)

        multiproc_input = [(cluster_prediction, predicted_clusters, bootstrapped_test_networks[i]) for
                           predicted_clusters, i in bootstrapped_cluster_predictions]
        bootstrap_cplx_cmp_metrics = p.map(comparison_helper, multiproc_input)
        try:
            pickle.dump(bootstrap_cplx_cmp_metrics, open(args.output_file + 'bootstrap_cplx_cmp_metrics.pkl', 'rb'))
        finally:
            print('pickling bootstrap_cplx_cmp_metrics didn\'t work...:\\')

        if best_eval is None or best_eval < metric_dict[args.eval_metric]:
            best_eval = metric_dict[args.eval_metric]
            best_size = size
            best_ii = ii
            best_density = density
            best_overlap = overlap
            best_fraction = fraction
            best_threshold_score = threshold_score
            best_inflation = inflation
            best_cluster_prediction = cluster_prediction

        if args.output_file is not None and args.output_all:
            output_filename = "%s.ii%s.%s" % (
                ".".join(args.output_file.split('.')[:-1]), ii, args.output_file.split('.')[-1])
            with open(output_filename, "w") as output_file:
                for cluster in cluster_prediction:
                    output_file.write(' '.join(cluster))
                    output_file.write("\n")

    if args.output_file is not None:
        with open(args.output_file +
                  '_{0}_{1}_{2}_{3}_{4}_{5}_{6}.txt'.format(
                      best_size, best_ii, best_density, best_overlap, best_fraction,
                      best_threshold_score, best_inflation), "w") as output_file:
            for cluster in best_cluster_prediction:
                output_file.write(' '.join(cluster))
                output_file.write("\n")

    p.close()
    p.join()


def comparison_helper(parameter_tuple):
    gold_std = parameter_tuple[0]
    pred_clst = parameter_tuple[1]
    test_net = parameter_tuple[2]
    ppi_recovered_count = 0
    for ppi in test_net:
        prot1 = ppi.split()[0]
        prot2 = ppi.split()[1]
        for clst in pred_clst:
            if prot1 in clst and prot2 in clst:
                ppi_recovered_count += 1
                break

    cplx_cmp = cc.ComplexComparison(gold_std, pred_clst)
    cplx_cmp_normalize = cc.ComplexComparison(gold_std, pred_clst, normalize_by_combinations=True, pseudocount=0.00001)

    d = dict()
    d['acc'] = cplx_cmp.acc()
    d['sensitivity'] = cplx_cmp.sensitivity()
    d['ppv'] = cplx_cmp.ppv()
    d['mmr'] = cplx_cmp.mmr()
    d['percent_ppi_recovered'] = (1.0 * ppi_recovered_count) / len(test_net)
    d['precision_recall_product'] = cplx_cmp.precision_recall_product()
    ccmm = cplx_cmp.clique_comparison_metric_mean()
    d['clique_precision_mean'] = ccmm['precision_mean']
    d['clique_recall_mean'] = ccmm['recall_mean']
    ccmm_normalize = cplx_cmp_normalize.clique_comparison_metric_mean()
    d['clique_precision_mean_normalize'] = ccmm_normalize['precision_mean']
    d['clique_recall_mean_normalize'] = ccmm_normalize['recall_mean']
    ccmm_normalize_weighted = cplx_cmp_normalize.clique_comparison_metric_mean(weighted=True)
    d['clique_precision_mean_normalize_weighted'] = ccmm_normalize_weighted['precision_mean']
    d['clique_recall_mean_normalize_weighted'] = ccmm_normalize_weighted['recall_mean']

    return d


def cluster_helper(parameter_dict):
    input_network_list = parameter_dict['network_list']
    ppi_scores = parameter_dict['ppi_scores']
    args = parameter_dict['args']
    size = str(parameter_dict['size'])
    density = str(parameter_dict['density'])
    overlap = str(parameter_dict['overlap'])
    seed_method = parameter_dict['seed_method']
    threshold_score = parameter_dict['threshold_score']
    i = parameter_dict['i']
    inflation = parameter_dict.get('inflation', None)
    trim2threshold = parameter_dict.get('trim2threshold', False)
    twostep_combination = parameter_dict['twostep_combination']

    if not os.path.exists(args.temp_dir) and args.temp_dir is not None:
        os.makedirs(args.temp_dir)
    fileTemp = tf.NamedTemporaryFile(delete=False, dir=args.temp_dir, mode="w", encoding='utf-8')
    dirTemp = tf.mkdtemp(dir=args.temp_dir)

    try:
        for bstrap_ppi in input_network_list:
            fileTemp.write(bstrap_ppi)
        fileTemp.close()

        if twostep_combination[0] == 'clusterone':
            proc = sp.Popen(
                ['java', '-jar', args.clustone_jar, fileTemp.name, '-s', str(size), '-d', str(density),
                 '--max-overlap', str(overlap), '--seed-method', seed_method], stdout=sp.PIPE, stderr=sp.PIPE)
            clust_out, err = proc.communicate()

            predicted_clusters = []
            for line in str(clust_out, 'UTF-8').split('\n'):
                if len(line.split()) > 0:
                    predicted_clusters.append(line.split())

            outfilename = args.temp_dir + 'ClusterONE_complexPredictions_' + str(i) + '.txt'
            outfile = open(outfilename, "w")
            for pred in predicted_clusters:
                outline = ' '.join(pred) + '\n'
                outfile.writelines(outline)
            outfile.close()
    finally:
        pass

    print("finished clustering phase1")
    sys.stdout.flush()

    if len(twostep_combination) >= 2:
        if twostep_combination[1] == 'mcl':
            if inflation != str(None):
                mcl_clusters = []
                for clust in predicted_clusters:
                    fileTemp = tf.NamedTemporaryFile(delete=False, dir=args.temp_dir, mode="w", encoding='utf-8')
                    outTemp = tf.NamedTemporaryFile(delete=False, dir=args.temp_dir, mode="w", encoding='utf-8')

                    try:
                        for prot1, prot2 in it.combinations(clust, 2):
                            try:
                                ppi_score = ppi_scores[frozenset([prot1, prot2])]
                                if ppi_score < threshold_score:
                                    ppi_score = 0.0

                            except KeyError:
                                ppi_score = 0.0

                            ppi_str = "%s\t%s\t%s\n" % (prot1, prot2, ppi_score)
                            fileTemp.write(ppi_str)
                        fileTemp.close()

                        proc = sp.Popen([args.mcl_bin, fileTemp.name, '--abc', '-o', outTemp.name,
                                         '-I', str(inflation)], stdout=sp.PIPE, stderr=sp.PIPE)
                        _, _ = proc.communicate()

                        outfile = open(outTemp.name, "rb")
                        for line in outfile.readlines():
                            mcl_clusters.append(line.split())
                        outfile.close()

                        outfilename = args.temp_dir + 'ClusterONE+MCL_complexPredictions_' + str(i) + '.txt'
                        outfile = open(outfilename, "w")
                        for pred in mcl_clusters:
                            outline = ' '.join(pred) + '\n'
                            outfile.writelines(outline)
                        outfile.close()

                    finally:
                        pass

                predicted_clusters = mcl_clusters
    print("finished clustering phase2")

    if trim2threshold:
        predicted_clusters = trim_clusters2threshold(predicted_clusters, threshold_score, ppi_scores)

    return predicted_clusters, i


def trim_clusters2threshold(predicted_clusters, threshold_score, ppi_scores):
    trimed_clusters = []

    for clust in predicted_clusters:
        trimed_clust = []
        for prot1 in clust:
            prot1_max_score = 0.0
            for prot2 in clust:
                if prot1 != prot2:
                    try:
                        prot1_max_score = max(prot1_max_score, ppi_scores[frozenset([prot1, prot2])])
                    except KeyError:
                        continue
            if prot1_max_score < threshold_score:
                print("removing prot1: %s max score: %s" % (prot1, prot1_max_score))
            else:
                trimed_clust.append(prot1)
        if len(trimed_clust) > 1:
            trimed_clusters.append(trimed_clust)

    return trimed_clusters
