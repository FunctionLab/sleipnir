import sys
import os
import argparse
currPath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(currPath)
import seekUtils as sutils


if __name__=="__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--config', '-c', default=None, type=str, required=True,
                           help='config file')
    argParser.add_argument('--queryFile', '-q', default=None, type=str, required=True,
                           help='file with list of queries')
    argParser.add_argument('--outputDir', '-o', default=None, type=str, required=True,
                           help='output directory to put query results')
    argParser.add_argument('--concurrency', '-m', default=1, type=int,
                           help='number of concurrent process to split queries across')
    argParser.add_argument('--binDir', '-b', default="/data/seek/bin", type=str,
                           help='location of executable binaries')
    args = argParser.parse_args()

    if not os.path.exists(args.config):
        print("Error: -c config file doesn't exist")
        sys.exit(-1)

    if not os.path.exists(args.binDir):
        print("Error: -b binary executable directory doesn't exist")
        sys.exit(-1)

    if not os.path.exists(args.queryFile):
        print("Error: -q query file doesn't exist")
        sys.exit(-1)

    if not os.path.exists(args.outputDir):
        print("Error: -o output directory doesn't exist")
        sys.exit(-1)

    cfg = sutils.loadConfig(args.config);
    cfg.binDir = args.binDir
    sutils.runSeekMiner(cfg, args.queryFile, args.outputDir, concurrency=args.concurrency)
