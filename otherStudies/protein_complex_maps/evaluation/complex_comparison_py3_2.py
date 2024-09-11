# Compare cluster predictions to gold standard complexes##

import bisect
import itertools as it
import matplotlib as mpl
import numpy as np
import numpy.random as rand
import pandas as pd
import random
import scipy.misc as misc
from scipy.stats import hmean

mpl.use('Agg')
import matplotlib.pyplot as plt


class ComplexComparison(object):

    def __init__(self, gold_standard=[], clusters=[], exclusion_complexes=[],
                 samples=10000, pseudocount=1,
                 exact=False, max_clique=None, remove_non_gold_standard_proteins=False,
                 normalize_by_combinations=False):
        self.gold_standard = [set(x) for x in gold_standard]
        self.gold_standard_proteins = set()
        for x in gold_standard:
            self.gold_standard_proteins = self.gold_standard_proteins.union(x)

        self.clusters = [set(x) for x in clusters]
        if remove_non_gold_standard_proteins:
            self.remove_non_gold_standard_proteins()

        self.normalize_by_combinations = normalize_by_combinations
        self.exclusion_complexes = [set(x) for x in exclusion_complexes]
        self.intersection_table = None
        self.na_table = None
        self.clique_comparison_metric_results = None
        self.samples = samples
        self.pseudocount = pseudocount
        self.exact = exact
        self.max_clique = max_clique

    def remove_non_gold_standard_proteins(self, ):
        filtered_clusters = []
        for c in self.clusters:
            c_filt = c.intersection(self.gold_standard_proteins)
            if len(c_filt) > 1:
                filtered_clusters.append(c_filt)
        self.clusters = filtered_clusters

        return

    def get_gold_standard(self, ):
        return self.gold_standard

    def get_exclusion_complexes(self, ):
        return self.exclusion_complexes

    def get_gold_standard_proteins(self, ):
        return self.gold_standard_proteins

    def get_clusters(self, ):
        return self.clusters

    def get_na_table(self, ):
        if self.na_table is None:
            self.generate_na_table()

        return self.na_table

    def get_intersection_table(self, ):
        if self.intersection_table is None:
            self.generate_intersection_table()

        return self.intersection_table

    def precision_recall_product(self, topN=None):
        pm = self.precision_measure(topN=topN)
        rm = self.recall_measure(topN=topN)

        return hmean([pm, rm])

    def precision_measure(self, topN=None):
        max_df = self.prediction_wise_max_matching_ratio_distribution(topN=topN)
        lengths = [len(x) for x in self.get_clusters()]
        sum_weighted_max = np.sum(lengths * max_df)
        precision_measure = sum_weighted_max / np.sum(lengths)

        return precision_measure

    def recall_measure(self, topN=None):
        max_df = self.max_matching_ratio_distribution(topN=topN)
        lengths = [len(x) for x in self.get_gold_standard()]
        sum_weighted_max = np.sum(lengths * max_df)
        recall_measure = sum_weighted_max / np.sum(lengths)

        return recall_measure

    def mmr_pwmmr_hmean(self, ):
        mmr = self.mmr()
        pwmmr = self.pwmmr()

        return hmean([mmr, pwmmr])

    def mmr(self, topN=None):
        return self.max_matching_ratio(topN)

    def pwmmr(self, topN=None):
        return self.prediction_wise_max_matching_ratio(topN)

    def max_matching_ratio(self, topN=None):
        max_df = self.max_matching_ratio_distribution(topN=topN)
        sum_max_df = max_df.sum()
        num_of_gold_complexes = len(self.get_gold_standard())
        mmr = sum_max_df / num_of_gold_complexes

        return mmr

    def prediction_wise_max_matching_ratio(self, topN=None):
        max_df = self.prediction_wise_max_matching_ratio_distribution(topN=topN)
        sum_max_df = max_df.sum()
        num_of_clusters = len(self.get_clusters())
        pwmmr = sum_max_df / num_of_clusters

        return pwmmr

    def max_matching_ratio_distribution(self, topN=None):
        df = self.get_na_table()
        if topN is not None:
            df = df.ix[:, :topN]
        max_df = df.max(axis=1)

        return max_df

    def prediction_wise_max_matching_ratio_distribution(self, topN=None):
        df = self.get_na_table()
        if topN is not None:
            df = df.ix[:, :topN]
        max_df = df.max(axis=0)

        return max_df

    def sensitivity_distribution(self, topN=None):
        df = self.get_intersection_table()
        if topN is not None:
            df = df.ix[:, :topN]
        sn_df = df.divide(list(map(len, self.get_gold_standard())), axis=0)
        sn_co = sn_df.max(axis=1)

        return sn_co

    def sensitivity(self, topN=None):
        sn_co = self.sensitivity_distribution(topN=topN)
        goldstd_sizes = list(map(len, self.get_gold_standard()))
        Sn = 1.0 * sum(sn_co * goldstd_sizes) / sum(goldstd_sizes)

        return Sn

    def ppv_distribution(self, topN=None):
        df = self.get_intersection_table()
        if topN is not None:
            df = df.ix[:, :topN]
        ppv_df = df.divide(1.0 * df.sum())
        ppv_df = ppv_df.fillna(0.0)
        PPV_cl = ppv_df.max()

        return PPV_cl

    def ppv(self, topN=None):
        df = self.get_intersection_table()
        PPV_cl = self.ppv_distribution(topN=topN)
        try:
            PPV = 1.0 * sum(PPV_cl * df.sum()) / sum(df.sum())
        except ZeroDivisionError:
            PPV = None

        return PPV

    def acc(self, topN=None):
        Sn = self.sensitivity(topN)
        PPV = self.ppv(topN)
        Acc = (Sn * PPV) ** (1.0 / 2)

        return Acc

    def clique_comparison_metric_grandf1score(self, mean_func=hmean):
        result_dict = self.clique_comparison_metric()
        f1score_list = [result_dict[x]['f1score'] for x in result_dict.keys()]
        grand_f1score = 0.0
        if 0.0 not in f1score_list:
            grand_f1score = mean_func(f1score_list)

        return grand_f1score

    def clique_comparison_metric_mean(self, weighted=False):
        result_dict = self.clique_comparison_metric()
        precision_list = [result_dict[x]['precision'] for x in result_dict.keys()]
        recall_list = [result_dict[x]['recall'] for x in result_dict.keys()]

        if weighted:
            sum_complexes = sum([result_dict[x]['numOfClusters'] for x in result_dict.keys()])
            weights_list = [1.0 * result_dict[x]['numOfClusters'] / sum_complexes for x in result_dict.keys()]
            precision_mean = np.average(precision_list, weights=weights_list)
            recall_mean = np.average(recall_list, weights=weights_list)
        else:
            precision_mean = np.mean(precision_list)
            recall_mean = np.mean(recall_list)

        return {'precision_mean': precision_mean, 'recall_mean': recall_mean}

    def clique_comparison_metric(self, force=False):

        if self.clique_comparison_metric_results is not None and not force:
            return self.clique_comparison_metric_results

        clusters = [clust & self.get_gold_standard_proteins() for clust in self.get_clusters()]
        clust_max_len = np.max(list(map(len, clusters)))
        gs_max_len = np.max(list(map(len, self.get_gold_standard())))
        return_dict = dict()

        if self.max_clique is None:
            max_len = clust_max_len
        else:
            max_len = min(self.max_clique, clust_max_len)

        cumulative_gs_tp = 0.0
        cumulative_tp = 0.0
        cumulative_fn = 0.0
        cumulative_fp = 0.0
        for size in range(2, max_len + 1):

            zero_w_pseudocount = 1.0 * self.pseudocount / (self.samples + 2 * self.pseudocount)
            if size > 4 and return_dict[size - 1]['precision'] == zero_w_pseudocount and \
                    return_dict[size - 1]['recall'] == zero_w_pseudocount and \
                    return_dict[size - 2]['precision'] == zero_w_pseudocount and \
                    return_dict[size - 2]['recall'] == zero_w_pseudocount:
                return_dict[size] = {'precision': zero_w_pseudocount, 'f1score': zero_w_pseudocount,
                                     'recall': zero_w_pseudocount, 'cumulative_recall': zero_w_pseudocount,
                                     'cumulative_precision': zero_w_pseudocount, 'numOfClusters': 0}
                continue

            if gs_max_len < size:
                return_dict[size] = {'precision': zero_w_pseudocount, 'f1score': zero_w_pseudocount,
                                     'recall': zero_w_pseudocount, 'cumulative_recall': zero_w_pseudocount,
                                     'cumulative_precision': zero_w_pseudocount, 'numOfClusters': 0}
                continue

            if self.exact:
                result_dict = self.clique_comparison_exact(size)
            else:
                result_dict = self.clique_comparison(size)

            cumulative_gs_tp += result_dict['gs_tp']
            cumulative_tp += result_dict['tp']
            cumulative_fp += result_dict['fp']
            cumulative_fn += result_dict['fn']
            cumulative_recall = 1.0 * cumulative_gs_tp / (cumulative_gs_tp + cumulative_fn)
            cumulative_precision = 1.0 * cumulative_tp / (cumulative_tp + cumulative_fp)
            recall = 1.0 * result_dict['gs_tp'] / (result_dict['gs_tp'] + result_dict['fn'])
            precision = 1.0 * result_dict['tp'] / (result_dict['tp'] + result_dict['fp'])

            f1score = 0.0
            if precision != 0.0 and recall != 0.0:
                f1score = hmean([recall, precision])
            return_dict[size] = {'precision': precision, 'recall': recall, 'f1score': f1score,
                                 'cumulative_recall': cumulative_recall, 'cumulative_precision': cumulative_precision,
                                 'numOfClusters': result_dict['numOfClusters']}

        self.clique_comparison_metric_results = return_dict

        return return_dict

    def clique_comparison_exact(self, clique_size):
        true_positives = self.pseudocount
        gs_true_positives = self.pseudocount
        false_positives = self.pseudocount
        false_negatives = self.pseudocount
        if len(self.get_exclusion_complexes()) > 0:
            raise Exclusion_Complexes_Exception(
                "ERROR: EXCLUSION COMPLEX functionality is not implemented for exact method "
                "(def clique_comparison_exact), use estimated (def clique_comparison)")

        clusters = [clust & self.get_gold_standard_proteins() for clust in self.get_clusters()
                    if len(clust & self.get_gold_standard_proteins()) >= clique_size]

        for clust in clusters:
            is_positive_list = [np.max(list(map(set(group).issubset, self.get_gold_standard())))
                                for group in it.combinations(clust, clique_size)]
            tp_curr = sum(is_positive_list)
            if self.normalize_by_combinations:
                tp_curr = tp_curr * (1.0 * len(clust) / misc.comb(len(clust), clique_size))

            true_positives += tp_curr
            fp_curr = len(is_positive_list) - tp_curr
            if self.normalize_by_combinations:
                fp_curr = fp_curr * (1.0 * len(clust) / misc.comb(len(clust), clique_size))
            false_positives += fp_curr

        for gs_clust in self.get_gold_standard():
            is_positive_list = [np.max(list(map(set(gs_group).issubset, clusters)))
                                for gs_group in it.combinations(gs_clust, clique_size)]
            gs_tp_curr = sum(is_positive_list)
            if self.normalize_by_combinations:
                gs_tp_curr = gs_tp_curr * (1.0 * len(gs_clust) / misc.comb(len(gs_clust), clique_size))
            gs_true_positives += gs_tp_curr

            fn_curr = sum(np.logical_not(is_positive_list))
            if self.normalize_by_combinations:
                fn_curr = fn_curr * (1.0 * len(gs_clust) / misc.comb(len(gs_clust), clique_size))
            false_negatives += fn_curr

        return_dict = dict()
        return_dict['tp'] = true_positives
        return_dict['gs_tp'] = gs_true_positives
        return_dict['fp'] = false_positives
        return_dict['fn'] = false_negatives
        return_dict['numOfClusters'] = len(clusters)

        return return_dict

    def clique_comparison(self, clique_size):
        true_positives = self.pseudocount
        gs_true_positives = self.pseudocount
        false_positives = self.pseudocount
        false_negatives = self.pseudocount
        clusters = [clust & self.get_gold_standard_proteins()
                    for clust in self.get_clusters() if len(clust & self.get_gold_standard_proteins()) >= clique_size]

        weights_array = [misc.comb(len(clust & self.get_gold_standard_proteins()), clique_size) for clust in clusters]
        wrg = WeightedRandomGenerator(weights_array)

        random_cliques = set()
        random_clique_weights = dict()
        for s in range(self.samples):

            clust = clusters[wrg()]
            clust_intersection = clust & self.get_gold_standard_proteins()

            if len(clust_intersection) <= 0:
                continue

            shuffled_l = rand.permutation(list(clust_intersection))
            random_cliques.add(frozenset(shuffled_l[:clique_size]))
            random_clique_weights[frozenset(shuffled_l[:clique_size])] = \
                (1.0 * len(clust) / misc.comb(len(clust), clique_size))

        for random_clique in random_cliques:
            if np.max(list(map(random_clique.issubset, self.get_gold_standard()))):
                if self.normalize_by_combinations:
                    true_positives += 1 * random_clique_weights[random_clique]
                else:
                    true_positives += 1
            else:
                if len(self.get_exclusion_complexes()) > 0 and \
                        np.max(list(map(random_clique.issubset, self.get_exclusion_complexes()))):
                    continue
                else:
                    if self.normalize_by_combinations:
                        false_positives += 1 * random_clique_weights[random_clique]
                    else:
                        false_positives += 1

        gs_clusters = [gs_clust for gs_clust in self.get_gold_standard() if len(gs_clust) >= clique_size]
        gs_wrg = WeightedRandomGenerator([misc.comb(len(gs_clust), clique_size) for gs_clust in gs_clusters])
        random_gs_cliques = set()
        random_gs_clique_weights = dict()
        for s in range(self.samples):
            gs_clust = gs_clusters[gs_wrg()]
            shuffled_l = rand.permutation(list(gs_clust))
            random_gs_cliques.add(frozenset(shuffled_l[:clique_size]))
            random_gs_clique_weights[frozenset(shuffled_l[:clique_size])] = \
                (1.0 * len(gs_clust) / misc.comb(len(gs_clust), clique_size))

        for random_clique in random_gs_cliques:
            if np.max(list(map(random_clique.issubset, clusters))):
                if self.normalize_by_combinations:
                    gs_true_positives += 1 * random_gs_clique_weights[random_clique]
                else:
                    gs_true_positives += 1
            else:
                if self.normalize_by_combinations:
                    false_negatives += 1 * random_gs_clique_weights[random_clique]
                else:
                    false_negatives += 1

        return_dict = dict()
        return_dict['tp'] = true_positives
        return_dict['gs_tp'] = gs_true_positives
        return_dict['fp'] = false_positives
        return_dict['fn'] = false_negatives
        return_dict['numOfClusters'] = len(clusters)

        return return_dict

    def generate_intersection_table(self, ):
        rows_list = []
        for cc in self.get_gold_standard():
            d = dict()
            for i, clst in enumerate(self.get_clusters()):
                d[i] = 1.0 * len(set.intersection(set(clst), set(cc)))
            rows_list.append(d)
        self.intersection_table = pd.DataFrame(rows_list)

    def generate_na_table(self, ):
        rows_list = []
        for cc in self.get_gold_standard():
            d = dict()
            for i, clst in enumerate(self.get_clusters()):
                numerator = (1.0 * len(set.intersection(set(clst), set(cc)))) ** 2
                denominator = 1.0 * (len(set(clst)) * len(set(cc)))
                d[i] = numerator / denominator
            rows_list.append(d)
        self.na_table = pd.DataFrame(rows_list)


def chunks(l, n):
    n = max(1, n)
    return [l[i:i + n] for i in range(0, len(l), n) if i + n < len(l)]


def rand_combinations(l, n, m):
    try:
        i = int(np.ceil(1.0 * m / (len(l) / n)))
        for ii in range(i):
            shuffled_l = rand.permutation(l)
            for ch in chunks(shuffled_l, n):
                yield ch
    except:
        return


class WeightedRandomGenerator(object):
    def __init__(self, weights):
        self.totals = []
        running_total = 0

        for w in weights:
            running_total += w
            self.totals.append(running_total)

    def next(self):
        rnd = random.random() * self.totals[-1]
        return bisect.bisect_right(self.totals, rnd)

    def __call__(self):
        return self.next()


class Exclusion_Complexes_Exception(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class Dict2Class(object):
    def __init__(self, my_dict):
        for key in my_dict:
            setattr(self, key, my_dict[key])
