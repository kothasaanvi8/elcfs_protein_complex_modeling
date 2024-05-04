# Helper function to complex builder

import os, shutil, subprocess, sys

def buildComplexes(outputPrefix, complexBuilder_threshold, *otherArgs):

    try:
        os.environ.pop('PYTHONIOENCODING')
    except KeyError:
        pass

    path_to_complex_builder = './util/getAllCliques_gwMods.py '
    path_to_pairsInfo = './' + outputPrefix + '/' + outputPrefix + '_pairsInfo.tsv '
    path_to_pairsPredictions = './' + outputPrefix + '/' + outputPrefix + '_pairsPredictions '
    cmd = 'python ' + path_to_complex_builder + path_to_pairsInfo + \
          path_to_pairsPredictions + str(complexBuilder_threshold) + ' .'

    print('For reference, the command passing to the complex builder is as follows: ')
    print(cmd)
    cmdOut = subprocess.Popen(cmd, cwd='./',
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, shell=True).communicate()[0]
    print(cmdOut)

    if otherArgs[0]:
        shutil.copytree('./' + outputPrefix + '/', otherArgs[0] + outputPrefix + '/')


if __name__ == '__main__':

    if len(sys.argv) >= 3:
        outputPrefix = sys.argv[1]
        complexBuilder_threshold = str(sys.argv[2])
        otherArgs = sys.argv[3:]
        buildComplexes(outputPrefix, complexBuilder_threshold, otherArgs)
    else:
        print('Passed an insufficient number of arguments.')
