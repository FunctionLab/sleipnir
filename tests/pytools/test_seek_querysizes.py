import os
import sys
import time
import glob
import argparse
import utils
from envbash import load_envbash
from rank_correlation import files_rank_correlation


min_result_correlation = 0.95


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

    bashEnvironmentFile = os.path.join(args.seekdir, 'seek_env')
    print('Load bash environment file {}'.format(bashEnvironmentFile))
    load_envbash(bashEnvironmentFile)

    # The query files have the query strings to run (multiple queries per file
    #   one query per line), located in goldStdDir
    filepattern = os.path.join(goldStdDir, r'*.query.txt')
    queryFiles = [queryfile for queryfile in glob.iglob(filepattern)]

    correlation_errors = 0
    # Run SeekMiner for all query files
    for qfile in queryFiles:
        path, filename = os.path.split(qfile)
        # the first part of the query file name will be used for the result directory name
        queryName = filename.split('.')[0]
        resultDir = os.path.join(args.outputdir, queryName)
        utils.checkAndMakePath(resultDir)
        print('SeekMiner run query {}'.format(queryName))
        seekMinerCmd = 'time {seekminer} -x {db}/dataset.map -i {db}/gene_map.txt ' \
                       '-d {db}/db.combined -p {db}/prep.combined ' \
                       '-P {db}/platform.combined -Q {db}/quant2 ' \
                       '-u {db}/sinfo.combined -n 1000 -b 200  ' \
                       '-R {db}/dataset_size -V CV -I LOI -z z_score -m -M -O ' \
                       '-q {queryfile} -o {resultdir} ' \
                       '&> {resultdir}/seekminer.out'.format(
                         seekminer=seekMinerBin, db=args.seekdir,
                         queryfile=qfile, resultdir=resultDir
                       )
        os.system(seekMinerCmd)

        # compare output to gold standard results
        filepattern = os.path.join(goldStdDir, queryName, r'*.results.txt')
        for goldResultFile in glob.iglob(filepattern):
            resultFilename = os.path.basename(goldResultFile)
            testfile = os.path.join(resultDir, resultFilename)
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
