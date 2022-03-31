'''
Run a set of bioinformative queries on a seek database and compare the results with a
known-good run of the same queries.
Config file settings:
seekPath - full path to seek database installation to test
binPath - full path to seek binaries to test
goldStdPath - full path to known-good output files and gold standard files for the test
outputPath - output directory for the test results
queryPath - full path to the query input directory
geneFile - file with list of all genes to include
geneMap - gene map file
recallPct - Recall depth (in percent) for calculating precision
'''
import os
import sys
import re
import glob
import time
import argparse
import utils
from envbash import load_envbash
import create_plot1a_filelists as filelist1a
import create_plot1c_filelists as filelist1c
import plot_seek1a as plot1a
import plot_seek1c as plot1c
import compare_result_files as cmpfiles


if __name__ == "__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--config', '-c', type=str, required=True,
                           help='Test configuration file')
    argParser.add_argument('--seek-dir', '-s', type=str, required=False,
                           help='Seek Directory')
    argParser.add_argument('--seek-bin', '-b', type=str, required=False,
                           help='Seek binary files')
    argParser.add_argument('--outputdir', '-o', type=str, required=False,
                           help='output directory where results will be written')
    argParser.add_argument('--verbose', '-v', default=False, action='store_true',
                           help='verbose output')
    # TODO - put this in config file instead
    argParser.add_argument('--known-good-results', '-g', type=str, required=False,
                           help='Expected results files from a known-good run')
    args = argParser.parse_args()

    # Parse test config file
    cfg = utils.load_config_file(args.config)

    if args.seek_dir is not None:
        cfg.seekPath = args.seek_dir
    if args.seek_bin is not None:
        cfg.binPath = args.seek_bin
    if args.outputdir is not None:
        cfg.outputPath = args.outputdir

    utils.checkAndMakePath(cfg.outputPath)

    # check paths exist
    paths = [cfg.binPath, cfg.seekPath, cfg.queryPath, cfg.goldStdPath]
    for path in paths:
        if not os.path.exists(path):
            print('path {} does not exist'.format(path))
            sys.exit(-1)

    seekMinerBin = os.path.join(cfg.binPath, 'SeekMiner')
    seekEvaluatorBin = os.path.join(cfg.binPath, 'SeekEvaluator')
    dataDir = cfg.seekPath

    bashEnvironmentFile = os.path.join(cfg.seekPath, 'seek_env')
    print('Load bash environment file {}'.format(bashEnvironmentFile))
    load_envbash(bashEnvironmentFile)
    dbDir = os.environ.get('DB')
    platDir = os.environ.get('PLAT')
    prepDir = os.environ.get('PREP')
    sinfoDir = os.environ.get('SINFO')
    if None in (dbDir, platDir, prepDir, sinfoDir):
        print('Please set seek_env with the DB directory names')
        sys.exit(-1)
    dsetFile = os.environ.get('DSET_FILE')
    geneMapFile = os.environ.get('GENE_MAP_FILE')
    dsetSizeFile = os.environ.get('DSET_SIZE_FILE')
    if None in (dsetFile, geneMapFile, dsetSizeFile):
        print('Please set seek_env with the dataset file names')
        sys.exit(-1)

    assert os.path.normpath(cfg.geneMap) == os.path.join(dataDir, geneMapFile)

    # The query files have the query strings to run (multiple queries per file
    #   one query per line), located in query path
    filepattern = os.path.join(cfg.queryPath, r'*.query.txt')
    queryFiles = [queryfile for queryfile in glob.iglob(filepattern)]

    queryNames = []
    # Run SeekMiner for all query files
    for file in queryFiles:
        path, file = os.path.split(file)
        # the first part of the query file name will be used for the result directory name
        queryName = file.split('.')[0]
        queryNames.append(queryName)
        queryFile = os.path.join(cfg.queryPath, file)
        resultDir = os.path.join(cfg.outputPath, queryName)
        utils.checkAndMakePath(resultDir)
        outfile = os.path.join(resultDir, "seekminer.out")
        utils.file_truncate(outfile)
        print('SeekMiner run query {}'.format(queryName))
        seekMinerCmd = f'time {seekMinerBin} -x {dataDir}/{dsetFile} -i {dataDir}/{geneMapFile} ' \
                       f'-d {dataDir}/{dbDir} -p {dataDir}/{prepDir} ' \
                       f'-P {dataDir}/{platDir} -Q {dataDir}/quant2 ' \
                       f'-u {dataDir}/{sinfoDir} -n 1000 -b 200  ' \
                       f'-V CV -I LOI -z z_score -m -M -O ' \
                       f'-R {dataDir}/{dsetSizeFile} ' \
                       f'-q {queryFile} -o {resultDir} '
        utils.file_appendline(outfile, seekMinerCmd)
        if args.verbose:
            seekMinerCmd += " |& tee -a {}".format(outfile)
            print(seekMinerCmd, flush=True)
        else:
            seekMinerCmd += " &>> {}".format(outfile)
        os.system(seekMinerCmd)

    # Create filelists in preparation for running SeekEvaluator
    # create plot1a filelist
    for queryName in queryNames:
        print('Create filelists for plot1a {}'.format(queryName))
        goldDir = os.path.join(cfg.goldStdPath, queryName)
        resultDir = os.path.join(cfg.outputPath, queryName)
        outputDir = resultDir
        filelist1a.makePlot1aFilelists(goldDir, resultDir, outputDir,
                                       cfg.geneFile)

    # create plot1c filelist
    print('Create filelists for plot1c')
    queryNamesStr = ','.join(map(str, queryNames))
    filelist1c.makePlot1cFilelists(cfg.goldStdPath, cfg.outputPath, cfg.outputPath,
                                   cfg.geneFile, queryNamesStr)

    # Run Seek Evaluator
    # Run SeekEvaluator for plot 1a - aggregating for each query size
    for queryName in queryNames:
        # get the different query sizes from the fileslist names
        print('SeekEvaluator for query {}'.format(queryName))
        resultDir = os.path.join(cfg.outputPath, queryName)
        filepattern = os.path.join(resultDir, r'filelist_q*.query')
        for queryfilelist in glob.iglob(filepattern):
            goldfilelist = re.sub('query$', 'gold', queryfilelist)
            gscorefilelist = re.sub('query$', 'gscore', queryfilelist)
            excludefilelist = re.sub('query$', 'exclude', queryfilelist)
            includefilelist = re.sub('query$', 'include', queryfilelist)
            resultfile = re.sub('filelist_q', 'result_qgroup.', queryfilelist)
            resultfile = re.sub(r'\.query$', '', resultfile)
            seekEvalCmd = \
                f'{seekEvaluatorBin} -S {goldfilelist} -G {gscorefilelist} ' \
                f'-Q {queryfilelist} -X {excludefilelist} -Y {includefilelist} ' \
                f'-i {cfg.geneMap} -d {resultDir} --pr --x_per {cfg.recallPct} ' \
                f'-M -B -f &> {resultfile}'
            # print(seekEvalCmd)
            os.system(seekEvalCmd)

    # Run SeekEvaluator for plot 1c to get aggregate precision at different recall depth
    recall_depths = [.01, .02, .04, .06, .08, .1, .2, .4, .6, .8, 1]
    for depth in recall_depths:
        print('SeekEvaluator for depth {}'.format(depth))
        resultfile = os.path.join(cfg.outputPath,
                                  'result_recall_curve.{}'.format(int(depth*100)))
        seekEvalCmd = f'{seekEvaluatorBin} ' \
            f'-S {cfg.outputPath}/filelist_seek1c.gold ' \
            f'-G {cfg.outputPath}/filelist_seek1c.gscore ' \
            f'-Q {cfg.outputPath}/filelist_seek1c.query ' \
            f'-X {cfg.outputPath}/filelist_seek1c.query ' \
            f'-Y {cfg.outputPath}/filelist_seek1c.include ' \
            f'-i {cfg.geneMap} -M -B -f --pr --x_per {depth} ' \
            f'-d {cfg.outputPath} &> {resultfile}'
        # print(seekEvalCmd)
        os.system(seekEvalCmd)

    # Create the plots and result csv files
    for queryName in queryNames:
        inputDir = os.path.join(cfg.outputPath, queryName)
        plot1a.plot_group(inputDir, cfg.outputPath, int(cfg.recallPct*100))

    plot1c.plot_curve(cfg.outputPath, cfg.outputPath)

    # Compare the results to known good
    compare_errors = 0
    filepattern = os.path.join(cfg.goldStdPath, r'*.csv')
    for goldcsvfile in glob.iglob(filepattern):
        filename = os.path.basename(goldcsvfile)
        testcsvfile = os.path.join(cfg.outputPath, filename)
        pctdiff = cmpfiles.get_pct_error(goldcsvfile, testcsvfile,
                                         skiprows=1, skipcols=3)
        print('Compare: {}, pct_diff {}'.format(filename, pctdiff))
        if pctdiff > cfg.maxPctDiff:
            compare_errors += 1
            print('Failed: pct diff exceeded {}'.format(testcsvfile))
    if compare_errors > 0:
        sys.exit(-1)
    sys.exit(0)
