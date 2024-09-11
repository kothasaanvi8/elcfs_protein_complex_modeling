import os
import sys
import subprocess as sp
import tempfile as tf

def cluster_helper(parameter_dict):
    input_network_list = parameter_dict['network_list']
    ppi_scores = parameter_dict['ppi_scores']
    args = parameter_dict['args']
    size = parameter_dict['size']
    density = parameter_dict['density']
    overlap = parameter_dict['overlap']
    seed_method = parameter_dict['seed_method']
    fraction = parameter_dict['fraction']
    threshold_score = parameter_dict['threshold_score']
    i = parameter_dict['i']
    clixo_alpha = parameter_dict.get('clixo_alpha', None)
    clixo_beta = parameter_dict.get('clixo_beta', None)
    inflation = parameter_dict.get('inflation', None)
    cliquesize = parameter_dict.get('cliquesize', None)
    timeout = parameter_dict.get('timeout', None)
    trim2threshold = parameter_dict.get('trim2threshold', False)
    nodelete = parameter_dict.get('nodelete', False)
    twostep_combination = parameter_dict['twostep_combination']

    if not os.path.exists(args.temp_dir) and args.temp_dir is not None:
        os.makedirs(args.temp_dir)
    fileTemp = tf.NamedTemporaryFile(delete=False, dir=args.temp_dir)
    dirTemp = tf.mkdtemp(dir=args.temp_dir)
    try:
        for bstrap_ppi in input_network_list:
            fileTemp.write(bstrap_ppi)
        fileTemp.close()

        print "in clusterone"
        proc = sp.Popen(
            ['java', '-jar', args.clustone_jar, fileTemp.name,
             '-s', str(size), '-d', str(density), '--max-overlap', str(overlap), '--seed-method', seed_method],
            stdout=sp.PIPE, stderr=sp.PIPE)
        clust_out, err = proc.communicate()

        predicted_clusters = []
        for line in clust_out.split('\n'):
            if len(line.split()) > 0:
                predicted_clusters.append(line.split())

    finally:
        if nodelete:
            pass

    outfilename = './FakeTemp/ClusterONE_complexPredictions_norm_' + str(i) + '.txt'
    outfile = open(outfilename, "w")
    for pred in predicted_clusters:
        outline = ' '.join(pred) + '\n'
        outfile.writelines(outline)
    outfile.close()

    print "finished clustering phase1"
    sys.stdout.flush()
