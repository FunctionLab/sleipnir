#!/usr/bin/python

'''
To reproduce the plots from the Seek paper, plot supplemental 1a requires
grouping queries by query size (such as all queries with 2 genes, with 4 genes etc).
This program iterates through the output files from SeekMiner and groups them by
query size and creates a set of filelists for each query size to be used by SeekEvaluator.
The filelists expected by SeekEvaluator are query files, gscore files, goldStd files
include files, exclude files.
The steps include:
1. Iterate through a directory of query files
2. Determine the query size (num genes in search string)
3. Append files corresponding to the query file to the appropriate file lists
'''

import os
import sys
import re
import glob
import argparse
from pathlib import Path
import utils


def makePlot1aFilelists(golddir, resultdir, outputdir, genefile):
    golddir = Path(golddir).resolve().as_posix()
    resultdir = Path(resultdir).resolve().as_posix()
    outputdir = Path(outputdir).resolve().as_posix()
    genefile = Path(genefile).resolve().as_posix()

    # Note: query and gold standard files will be from the known good results dir
    #  in this case the 'golddir'
    # gscore files will be from the 'resultdir'

    filecount = 0
    filepattern = os.path.join(golddir, r'[0-9]*.query')
    for queryfile in glob.iglob(filepattern):
        filecount += 1
        # print(queryfile)
        # read file to find query size
        lines = utils.file_read(queryfile)
        assert len(lines) == 1
        num_genes = len(lines[0].split())
        outbasename = os.path.join(outputdir, 'filelist_q{}'.format(num_genes))

        # append the files to the appropriate file lists
        # query file list
        queryfilelist = outbasename + '.query'
        utils.file_appendline(queryfilelist, queryfile)

        # gold standard file list
        #  - gold standard file is the genes that should be correlated with the query genes
        goldfilelist = outbasename + '.gold'
        utils.file_appendline(goldfilelist, re.sub('query$', 'gold', queryfile))

        # query results 'gscore' file list
        #  - gscore are the resulting gene correlation scores from the SeekMiner query result
        gscorefilelist = outbasename + '.gscore'
        qbasename = os.path.basename(queryfile)
        gscorefile = os.path.join(resultdir, qbasename.replace('query', 'gscore'))
        utils.file_appendline(gscorefilelist, gscorefile)

        # include file list - genes to include in the results
        includefile = outbasename + '.include'
        utils.file_appendline(includefile, genefile)

        # excluse file list
        # exclude all genes in the query file
        excludefile = outbasename + '.exclude'
        utils.file_appendline(excludefile, queryfile)

    if filecount == 0:
        print("No matching files found")
        exit(-1)

    print("Num files processed: {}".format(filecount))


if __name__ == "__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--resultdir', '-r', type=str, required=True,
                           help='directory that contains the result files from SeekMiner')
    argParser.add_argument('--golddir', '-q', type=str, required=True,
                           help='directory that contains the query and gold standard files')
    argParser.add_argument('--outputdir', '-o', type=str, required=True,
                           help='directory to write the filelists to')
    argParser.add_argument('--genefile', '-g', type=str, required=True,
                           help='list of all genes for the include file')
    args = argParser.parse_args()

    makePlot1aFilelists(args.golddir, args.resultdir, args.outputdir, args.genefile)
