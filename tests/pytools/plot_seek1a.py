'''
Script to plot the Seek results to reproduce the graphs 1a and 1c in the Seek paper
'''
import os
import sys
import glob
import argparse
import matplotlib.pyplot as plt
sys.path.append(os.path.dirname(__file__))
from utils import file_appendline, file_read
from struct_dict import StructDict

# Filename pattern output from SeekEvaluator that we will read in to plot
filebase = "result_qgroup."


def save_vals(results, filename):
    qsizes = list(results.keys())
    qsizes.sort()
    with open(filename, 'w') as fp:
        fp.write("qsize min max quart1 quart2 quart3\n")
        for i in qsizes:
            vals = results[i]
            outline = "{} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(
                i, vals.min, vals.max, vals.quart1, vals.quart2, vals.quart3
            )
            fp.write(outline)


def plot_group(inputdir, outputdir, recall_pct):
    filecount = 0
    group = os.path.basename(inputdir)
    filepattern = os.path.join(inputdir, filebase + r'*')
    res = StructDict()
    for resultfile in glob.iglob(filepattern):
        # read in the results files per query size
        # get the query size from the end part of the file name
        qsize = os.path.basename(resultfile).split('.')[-1]
        qsize = int(qsize)
        lines = file_read(resultfile)
        # 1st line is min/max
        # 2nd line is 1st, 2nd, 3rd quartile
        minmax = [float(x) for x in lines[0].split()]
        quartiles = [float(x) for x in lines[1].split()]
        vals = StructDict(
                {'min': minmax[0],
                 'max': minmax[1],
                 'quart1': quartiles[0],
                 'quart2': quartiles[1],
                 'quart3': quartiles[2]}
                )
        res[qsize] = vals
        filecount += 1

    # create an array of result values per query size
    qsizes = list(res.keys())
    qsizes.sort()
    # print(qsizes)
    quartile_1 = [res[qsize]['quart1'] for qsize in qsizes]
    quartile_2 = [res[qsize]['quart2'] for qsize in qsizes]
    quartile_3 = [res[qsize]['quart3'] for qsize in qsizes]
    # print(res)
    # print(quartile_1)
    # print(quartile_2)
    # print(quartile_3)
    print("Plot1a dir {} processed {} files".format(group, filecount))

    # plot the array
    plt.plot(qsizes, quartile_1, color="gray")
    plt.plot(qsizes, quartile_2, color="black")
    plt.plot(qsizes, quartile_3, color="gray")
    plt.ylim(ymin=0)
    plt.grid(which='major', linewidth='0.5', color='darkgray')
    plt.grid(which='minor', linewidth='0.5', color='lightgray')
    plt.fill_between(qsizes, quartile_1, quartile_3, color="gray")
    plt.title('Processes {} genes'.format(group))
    plt.xlabel('Query size')
    plt.ylabel('Fold precision at {}% over random'.format(recall_pct))
    # plt.draw()
    # plt.pause(1)
    filename = "{}r{}.pdf".format(group, recall_pct)
    outputfile = os.path.join(outputdir, filename)
    plt.savefig(outputfile)
    resfile = os.path.join(outputdir, filename.replace('pdf', 'csv'))
    save_vals(res, resfile)


if __name__ == "__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--input-dir', '-i', type=str, required=True,
                           help='input directory where the SeekEvaluator results are')
    argParser.add_argument('--result-dir', '-o', type=str, required=True,
                           help='directory where result files will be written')
    argParser.add_argument('--recall-pct', '-p', type=float, required=True,
                           help='Recall depth percent as specified in the SeekMiner run')
    args = argParser.parse_args()

    if args.recall_pct < 1:
        args.recall_pct *= 100

    plot_group(args.input_dir, args.result_dir, int(args.recall_pct))
