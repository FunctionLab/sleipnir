#!/usr/bin/python

'''
To reproduce the plot supplemental 1c from the Seek paper requires
grouping all query sizes together for all gene group sizes.
This program iterates through the output files from SeekMiner and groups them all
together making filelists to be used by SeekEvaluator.
The filelists expected by SeekEvaluator are query files, gscore files, goldStd files
include files, exclude files
'''

import os
import sys
import re
import glob
import argparse
#sys.path.append(os.path.dirname(__file__))
import utils


def makePlot1cFilelists(golddir, resultdir, outputdir, genefile, queryNames):
    golddir = utils.makeAbsolutePath(golddir)
    resultdir = utils.makeAbsolutePath(resultdir)
    outputdir = utils.makeAbsolutePath(outputdir)
    genefile = utils.makeAbsolutePath(genefile)

    queryfilelist = os.path.join(outputdir, 'filelist_seek1c.query')
    goldfilelist = os.path.join(outputdir, 'filelist_seek1c.gold')
    gscorefilelist = os.path.join(outputdir, 'filelist_seek1c.gscore')
    includefilelist = os.path.join(outputdir, 'filelist_seek1c.include')
    excludefilelist = os.path.join(outputdir, 'filelist_seek1c.exclude')

    filecount = 0
    queryNames = [name.rstrip() for name in queryNames.split(',')]
    for queryName in queryNames:
        queryDir = os.path.join(golddir, queryName)
        gscoreDir = os.path.join(resultdir, queryName)
        filepattern = os.path.join(queryDir, r'[0-9]*.query')
        for queryfile in glob.iglob(filepattern):
            filecount += 1
            # Add a line to each filelist
            utils.file_appendline(queryfilelist, queryfile)
            utils.file_appendline(goldfilelist, re.sub('query$', 'gold', queryfile))
            utils.file_appendline(includefilelist, genefile)
            utils.file_appendline(excludefilelist, queryfile)
            # gscore file is in the result dir rather than gold dir
            qbasename = os.path.basename(queryfile)
            gscorefile = os.path.join(gscoreDir, qbasename.replace('query', 'gscore'))
            utils.file_appendline(gscorefilelist, gscorefile)

    if filecount == 0:
        print("No matching files found")
        exit(-1)

    print("Num files processed: {}".format(filecount))


if __name__ == "__main__":
    # Iterate through a directory of query results (from SeekMiner)
    # Write out filelists (containing lists of filenames in the query groups)
    # The list files are for input to SeekEvaluator (query, gscore, goldStd, include, exclude)
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--resultdir', '-r', type=str, required=True,
                           help='directory that contains the test result files from SeekMiner')
    argParser.add_argument('--golddir', '-q', type=str, required=True,
                           help='directory that contains the query and gold standard files')
    argParser.add_argument('--outputdir', '-o', type=str, required=True,
                           help='directory to write the filelists to')
    argParser.add_argument('--queryNames', '-n', type=str, required=True,
                           help='list of query names corresponding to subdirectory ' \
                               'names in the results path')
    argParser.add_argument('--genefile', '-g', type=str, required=True,
                           help='list of all genes for the include file')
    args = argParser.parse_args()

    makePlot1cFilelists(args.golddir, args.resultdir, args.outputdir,
                        args.genefile, args.queryNames)
