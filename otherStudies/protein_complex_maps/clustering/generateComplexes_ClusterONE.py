import argparse
import itertools as it
import multiprocessing as mp
import os
import pickle, random
import subprocess as sp
import sys
import tempfile as tf
from tqdm import tqdm as progressMonitor


def cluster_helper(parameter_dict):
    input_network_list = parameter_dict['network_list']
    args = parameter_dict['args']
    size = parameter_dict['size']
    density = parameter_dict['density']
    overlap = parameter_dict['overlap']
    seed_method = parameter_dict['seed_method']
    temp_dir = parameter_dict['temp_dir']
    clustone_jar = parameter_dict['clustone_jar']
    i = parameter_dict['i']

    if not os.path.exists(temp_dir) and temp_dir is not None:
        os.makedirs(temp_dir)

    fileTemp = tf.NamedTemporaryFile(delete=False, dir=temp_dir)
    try:
        for bstrap_ppi in input_network_list:
            fileTemp.write(bstrap_ppi)
        fileTemp.close()

        print("in clusterone")

        proc = sp.Popen(
            ['java', '-jar', clustone_jar, fileTemp.name,
             '-s', str(size), '-d', str(density), '--max-overlap', str(overlap),
             '--seed-method', seed_method],
            stdout=sp.PIPE, stderr=sp.PIPE)

        print('right before execute')
        clust_out, err = proc.communicate()
        print(clust_out)

        predicted_clusters = []
        for line in clust_out.split('\n'):
            if len(line.split()) > 0:
                predicted_clusters.append(line.split())

        print('writing to disk')
        outfilename = open(temp_dir + '/clusterONE/{0}_{1}_{2_{3}_{4}_ClusterONE_complexes.txt'.format(
            i, str(size), str(density), str(overlap), seed_method), "w")
        for pred in predicted_clusters:
            outline = ' '.join(pred) + '\n'
            outfilename.writelines(outline)
        outfilename.close()

    finally:
        print('something still is not working')

    print("finished clustering phase1")
    sys.stdout.flush()


def generateComplexes(predictedPairs_filename, outputDir):
    parser = argparse.ArgumentParser(description="Finds optimal parameters for clustering of ppi network")
    parser.add_argument("--input_network", action="store", dest="input_network", required=True,
                        help="Filename of ppi network with optional edge weights (format: id\tid\tweight)")
    parser.add_argument("--random_seed", action="store", type=int, dest="random_seed",
                        required=False, default=42, help="Sets random seed (int), default=None")
    parser.add_argument("--clusterone", action="store", dest="clustone_jar",
                        required=False, default="./cluster_one-1.0.jar")
    parser.add_argument("--clusterone_size", action="store", dest="clusterone_size",
                        nargs='+', required=False, default=[2, 3, 4, 5])
    parser.add_argument("--clusterone_density", action="store", dest="clusterone_density",
                        nargs='+', required=False, default=[.1, .2, .3, .4, .5, .6, .7, .8, .9])
    parser.add_argument("--clusterone_max_overlap", action="store", dest="clusterone_max_overlap",
                        nargs='+', required=False, default=[.1, .2, .3, .4, .5, .6, .7, .8, .9])
    parser.add_argument("--clusterone_seed_method", action="store", dest="clusterone_seed_method",
                        nargs='+', required=False, default=['nodes', 'cliques', 'unused_nodes', 'edges'],
                        help="ClusterOne seed method parameter sweep "
                             "(nodes, cliques, unused_nodes, edges, default = nodes")
    parser.add_argument("--temp_dir", action="store", dest="temp_dir",
                        required=False, default=None)
    args = parser.parse_args(["--input_network", predictedPairs_filename, "--temp_dir", outputDir])

    predictionsOutput_dir = args.temp_dir + 'clusterONE/'
    if not os.path.exists(args.temp_dir) and args.temp_dir is not None:
        os.makedirs(args.temp_dir)
        os.makedirs(predictionsOutput_dir)

    with open(args.input_network, "r") as input_network_file:
        network_list = input_network_file.readlines()

    random.seed(args.random_seed)
    size_sweep = args.clusterone_size
    overlap_sweep = args.clusterone_max_overlap
    seed_method_sweep = args.clusterone_seed_method
    density_sweep = args.clusterone_density

    network_input_list = []
    for ii, parameters in enumerate(
            it.product(size_sweep, density_sweep, overlap_sweep, seed_method_sweep)):
        sys.stdout.flush()

        size, density, overlap, seed_method = parameters
        parameter_dict = dict()
        parameter_dict['network_list'] = network_list
        parameter_dict['args'] = args
        parameter_dict['size'] = str(size)
        parameter_dict['density'] = str(density)
        parameter_dict['overlap'] = str(overlap)
        parameter_dict['seed_method'] = str(seed_method)
        parameter_dict['temp_dir'] = args.temp_dir
        parameter_dict['clustone_jar'] = args.clustone_jar
        parameter_dict['i'] = ii
        network_input_list.append(parameter_dict)

    p = mp.Pool(mp.cpu_count()-1)
    cluster_predictions = [list(progressMonitor(p.imap(cluster_helper, network_input_list), total=len(network_list)))]
    p.close()
    p.join()

    pickle.dump(cluster_predictions, open(predictionsOutput_dir + 'ClusterONE_complexes.pkl', 'wb'))
