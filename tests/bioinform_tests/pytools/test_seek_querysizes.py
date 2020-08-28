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
    goldStdDir = args.known_good_results
    dbDir = args.seekdir

    bashEnvironmentFile = os.path.join(args.seekdir, 'seek_env')
    print('Load bash environment file {}'.format(bashEnvironmentFile))
    load_envbash(bashEnvironmentFile)

    # The query files have the query strings to run (multiple queries per file
    #   one query per line), located in goldStdDir
    filepattern = os.path.join(goldStdDir, r'*.query.txt')
    queryFiles = [queryfile for queryfile in glob.iglob(filepattern)]
    print('queryFiles: ' + str(queryFiles))

    correlation_errors = 0
    # Run SeekMiner for all query files
    for queryfile in queryFiles:
        path, filename = os.path.split(queryfile)
        # the first part of the query file name will be used for the result directory name
        queryName = filename.split('.')[0]
        resultDir = os.path.join(args.outputdir, queryName)
        utils.checkAndMakePath(resultDir)
        outfile = os.path.join(resultDir, "seekminer.out")
        utils.file_truncate(outfile)
        print('SeekMiner run query {}'.format(queryName))
        seekMinerCmd = f'time {seekMinerBin} -x {dbDir}/dataset.map -i {dbDir}/gene_map.txt ' \
                       f'-d {dbDir}/db.combined -p {dbDir}/prep.combined ' \
                       f'-P {dbDir}/platform.combined -Q {dbDir}/quant2 ' \
                       f'-u {dbDir}/sinfo.combined -n 1000 -b 200  ' \
                       f'-V CV -I LOI -z z_score -m -M -O ' \
                       f'-R {dbDir}/dataset_size ' \
                       f'-q {queryfile} -o {resultDir} '
        utils.file_appendline(outfile, seekMinerCmd)
        if args.verbose:
            seekMinerCmd += " |& tee -a {}".format(outfile)
            print(seekMinerCmd, flush=True)
        else:
            seekMinerCmd += " &>> {}".format(outfile)
        os.system(seekMinerCmd)

        # compare output to gold standard results
        filepattern = os.path.join(goldStdDir, queryName, r'*.results.txt')
        for goldResultFile in glob.iglob(filepattern):
            resultFilename = os.path.basename(goldResultFile)
            testfile = os.path.join(resultDir, resultFilename)
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
