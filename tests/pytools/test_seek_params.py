import os
import sys
import time
import glob
import argparse
import utils
from rank_correlation import files_rank_correlation


min_result_correlation = 0.95

# parameters to test
# "-V" ["CV","EQUAL","ORDER_STAT","VAR","USER","AVERAGE_Z","CV_CUSTOM"]
# "-I" ["LOI","LOO","XFOLD"] # for CV mode
# "-z" ["pearson","z_score", "z_score -e"]
# "-m"
# "-M"
# "-S"
# "-C" 1
param_tests = {
    "test1": "",
    "test2": "-m M",  # subtract each gene's avg z-score from the datasets and platforms
    "test3": "-V CV -I LOI",
    "test4": "-V CV -I LOO",
    "test5": "-V CV -I XFOLD",
    "test6": "-V EQUAL",
    "test7": "-V ORDER_STAT",
    # "test8": "-V VAR -U {gene-var-dir}", # don't have a gene-var dir to use
    "test9": "-V USER -J {dweight_filelist}",
    "test10": "-V CV_CUSTOM -K {goldStd}",
    "test11": "-V AVERAGE_Z",
    "test12": "-z pearson",
    "test13": "-z z_score",
    "test14": "-z z_score -e",
    "test15": "-C 1.0",  # fraction of query genes required to be present
    # "test16": "-S",   # generate random ranking score
}


def mktempdir():
    timestamp = str(int(time.time()))
    temp_dir = "/tmp/param_test_" + timestamp
    os.mkdir(temp_dir)
    return temp_dir


if __name__ == "__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--seekdir', '-s', type=str, required=True,
                           help='Seek directory')
    argParser.add_argument('--seekbin', '-b', type=str, required=True,
                           help='Seek binary files')
    argParser.add_argument('--known-good-results', '-g', type=str, required=True,
                           help='Directory with known-good results files')
    argParser.add_argument('--outputdir', '-o', type=str, required=True,
                           help='output directory where results will be written')
    args = argParser.parse_args()

    utils.checkAndMakePath(args.outputdir)

    # check paths exist
    paths = [args.seekdir, args.seekbin, args.known_good_results, args.outputdir]
    for path in paths:
        if not os.path.exists(path):
            print('path {} does not exist'.format(path))
            sys.exit(-1)

    seekMinerBin = os.path.join(args.seekbin, 'SeekMiner')

    goldStdDir = args.known_good_results
    queryfile = os.path.join(goldStdDir, 'param_test.query.txt')
    goldstd_list = os.path.join(goldStdDir, 'param_test.goldStd.txt')
    dweight_filelist = os.path.join(goldStdDir, 'dweight_filelist.txt')
    test1Dir = os.path.join(goldStdDir, 'test1')
    os.system("ls -1 -d {}/*.dweight > {}".format(test1Dir, dweight_filelist))

    correlation_errors = 0
    for test, params in param_tests.items():
        params = params.format(dweight_filelist=dweight_filelist, goldStd=goldstd_list)
        print("{}, {}".format(test, params))
        resultdir = os.path.join(args.outputdir, test)
        os.mkdir(resultdir)
        cmd = "{seekminer} -x {db}/dataset.map -i {db}/gene_map.txt " \
              "-d {db}/db.combined -p {db}/prep.combined -P {db}/platform.combined " \
              "-Q {db}/quant2 -u {db}/sinfo.combined -R {db}/dataset_size " \
              "-n 1000 -b 200 -q {queryfile} -o {resultdir} -O {params}". \
              format(seekminer=seekMinerBin, db=args.seekdir, resultdir=resultdir, 
                     queryfile=queryfile, params=params)
        # print(cmd)
        os.system(cmd)

        # compare output to gold standard results
        filepattern = os.path.join(goldStdDir, test, r'*.results.txt')
        for goldResultFile in glob.iglob(filepattern):
            filename = os.path.basename(goldResultFile)
            testfile = os.path.join(resultdir, filename)
            corrs = files_rank_correlation(goldResultFile, testfile)
            print('Result Correlations: {}'.format(corrs))
            for corr in corrs:
                if corr < min_result_correlation:
                    correlation_errors += 1
                    print('ERROR: Result correlation too low')

    print('Test Complete, {} correlation errors'.format(correlation_errors))
    if correlation_errors > 0:
        sys.exit(-1)
    sys.exit(0)
