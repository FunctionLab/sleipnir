'''
Transform a file that has the gold standard genes for a set of queries, one per line.
To a set of files where each file has the gold standard genes for one query.
Number the output file by the line number in the input file starting at 0.
'''
import os
import sys
import argparse
sys.path.append(os.path.dirname(__file__))
from utils import read_genes, write_goldstd


if __name__ == "__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--goldstd-file', '-g', type=str, required=True,
                           help='file with lists of gold standard genes, one list per line')
    argParser.add_argument('--result-dir', '-o', type=str, required=True,
                           help='directory where result files will be written')
    args = argParser.parse_args()

    gold = read_genes(args.goldstd_file)
    write_goldstd(gold, args.result_dir)
