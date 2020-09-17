import os
import sys
import time
import glob
import argparse
import utils
from envbash import load_envbash
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
    "test8": "-L",  # dataset size must be > 10
    "test9": "-V USER -J {dweight_filelist}",
    "test10": "-V CV_CUSTOM -K {goldStd}",
    "test11": "-V AVERAGE_Z",
    "test12": "-z pearson",
    "test13": "-z z_score",
    "test14": "-z z_score -e",
    "test15": "-C 1.0",  # fraction of query genes required to be present
    # "test16": "-V VAR -U {gene-var-dir}", # don't have a gene-var dir to use
    # "test17": "-S",   # generate random ranking score
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
    argParser.add_argument('--verbose', '-v', default=False, action='store_true',
                           help='verbose output')
    args = argParser.parse_args()

    utils.checkAndMakePath(args.outputdir)

    # check paths exist
    paths = [args.seekdir, args.seekbin, args.known_good_results, args.outputdir]
    for path in paths:
        if not os.path.exists(path):
            print('path {} does not exist'.format(path))
            sys.exit(-1)

    seekMinerBin = os.path.join(args.seekbin, 'SeekMiner')
    dbDir = args.seekdir

    bashEnvironmentFile = os.path.join(args.seekdir, 'seek_env')
    print('Load bash environment file {}'.format(bashEnvironmentFile))
    load_envbash(bashEnvironmentFile)

    goldStdDir = args.known_good_results
    queryfile = os.path.join(goldStdDir, 'param_test.query.txt')
    goldstd_list = os.path.join(goldStdDir, 'param_test.goldStd.txt')
    dweight_filelist = os.path.join(goldStdDir, 'dweight_filelist.txt')
    test1Dir = os.path.join(goldStdDir, 'test1')
    os.system("ls -1 -d {}/*.dweight > {}".format(test1Dir, dweight_filelist))

    correlation_errors = 0
    for test, params in param_tests.items():
        params = params.format(dweight_filelist=dweight_filelist, goldStd=goldstd_list)
        print("{}, params {}".format(test, params))
        resultdir = os.path.join(args.outputdir, test)
        outfile = os.path.join(resultdir, "seekminer.out")
        utils.checkAndMakePath(resultdir)
        utils.file_truncate(outfile)
        cmd = f"{seekMinerBin} -x {dbDir}/dataset.map -i {dbDir}/gene_map.txt " \
              f"-d {dbDir}/db.combined -p {dbDir}/prep.combined -P {dbDir}/platform.combined " \
              f"-Q {dbDir}/quant2 -u {dbDir}/sinfo.combined " \
              f"-n 1000 -b 200 -q {queryfile} " \
              f"-R {dbDir}/dataset_size " \
              f"-o {resultdir} -O {params} "
        utils.file_appendline(outfile, cmd)
        if args.verbose:
            cmd += " |& tee -a {}".format(outfile)
            print(cmd, flush=True)
        else:
            cmd += "&>> {}".format(outfile)
        os.system(cmd)

        # compare output to gold standard results
        filepattern = os.path.join(goldStdDir, test, r'*.results.txt')
        for goldResultFile in glob.iglob(filepattern):
            filename = os.path.basename(goldResultFile)
            testfile = os.path.join(resultdir, filename)
            corrs = files_rank_correlation(goldResultFile, testfile, remove_substr=".pcl")
            print('Result Correlations: {}'.format(corrs))
            for corr in corrs:
                if corr < min_result_correlation:
                    correlation_errors += 1
                    print('ERROR: Result correlation too low')

    print('Test Complete, {} correlation errors'.format(correlation_errors))
    if correlation_errors > 0:
        sys.exit(-1)
    sys.exit(0)
