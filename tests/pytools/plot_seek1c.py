
'''
Script to plot the Seek query results to reproduce the graph 1c in the Seek paper
'''
import os
import sys
import glob
import argparse
import matplotlib.pyplot as plt
# sys.path.append(os.path.dirname(__file__))
from utils import file_read
from struct_dict import StructDict

# Filename pattern output from SeekEvaluator that we will read in to plot
filebase = "result_recall_curve."


def save_vals(results, filename):
    recall_depths = list(results.keys())
    recall_depths.sort()
    with open(filename, 'w') as fp:
        fp.write("recall_pct min max quart1 quart2 quart3\n")
        for depth in recall_depths:
            vals = results[depth]
            outline = "{} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(
                depth, vals.min, vals.max, vals.quart1, vals.quart2, vals.quart3
            )
            fp.write(outline)


def plot_curve(inputdir, outputdir):
    filecount = 0
    res = StructDict()
    filepattern = os.path.join(inputdir, filebase + r'*')
    for resultfile in glob.iglob(filepattern):
        # read in the results files per query size
        # get the query size from the end part of the file name
        recall_pct = os.path.basename(resultfile).split('.')[-1]
        recall_pct = int(recall_pct)
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
        res[recall_pct] = vals
        filecount += 1

    # create an array of result values per query size
    recall_depths = list(res.keys())
    recall_depths.sort()
    # print(recall_depths)
    quartile_1 = [res[depth]['quart1'] for depth in recall_depths]
    quartile_2 = [res[depth]['quart2'] for depth in recall_depths]
    quartile_3 = [res[depth]['quart3'] for depth in recall_depths]
    # print(res)
    # print(quartile_1)
    # print(quartile_2)
    # print(quartile_3)
    print("Plot1c processed {} files".format(filecount))

    # plot the array
    plt.plot(recall_depths, quartile_1, color="gray")
    plt.plot(recall_depths, quartile_2, color="black")
    plt.plot(recall_depths, quartile_3, color="gray")
    plt.ylim(ymin=0)
    plt.xscale('log')
    plt.grid(which='major', linewidth='0.5', color='darkgray')
    plt.grid(which='minor', linewidth='0.5', color='lightgray')
    plt.fill_between(recall_depths, quartile_1, quartile_3, color="gray")
    plt.title('Seek Search Performance')
    plt.xlabel('Recall Percent')
    plt.ylabel('Fold precision over random')
    # plt.draw()
    # plt.pause(1)
    filename = "precision_vs_depth.pdf"
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

    args = argParser.parse_args()

    plot_curve(args.input_dir, args.result_dir)
