'''
Given two files with ranked lists of genes or datasets, determine how closely the
rankings of the files agree using Spearmans correlation.
FileA and FileB can have multiple lines, each line will be correlated with the
corresponding line in the other file and a list of correlation values will be
returned with one value per line.
Note that FileB is compared against FileA. So FileA ranking is the standard
and for each value in FileA a rank is looked up in FileB. If the value doesn't
appear in FileB then the rank will be 0 for B.
'''
import sys
import argparse
import scipy.stats as stats
import utils


def files_rank_correlation(fileA, fileB, remove_substr=None):
    # read in the file data lines
    A_lines = utils.file_read(fileA)
    B_lines = utils.file_read(fileB)

    assert len(A_lines) == len(B_lines)

    correlations = []
    # for each line calculate the correlation between the two files
    for i in range(len(A_lines)):
        if remove_substr is not None:
            A_lines[i] = A_lines[i].replace(remove_substr, "")
            B_lines[i] = B_lines[i].replace(remove_substr, "")
        # assign an order rank to the values in the line
        A_dict_rank = {val: j for j, val in enumerate(A_lines[i].split())}
        B_dict_rank = {val: j for j, val in enumerate(B_lines[i].split())}
        A_rank = []
        B_rank = []
        # for each value in line A get the rank order in list A and B
        keys = list(A_dict_rank.keys())
        max_rank = len(keys)
        for key in keys:
            A_rank.append(A_dict_rank[key])
            B_rank.append(B_dict_rank.get(key, max_rank))
        # check the correlation between the lists
        if A_rank == B_rank:
            # lists are identical
            corr = 1
        else:
            # calculate the spearman's correlation coefficient between the two rankings
            corr, _ = stats.spearmanr(A_rank, B_rank)
        # print(A_rank)
        # print(B_rank)
        correlations.append(corr)

    return correlations


if __name__ == "__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--file-a', '-a', type=str, required=True,
                           help='comparison file A')
    argParser.add_argument('--file-b', '-b', type=str, required=True,
                           help='comparison file B')
    argParser.add_argument('--expected-correlation', '-e', type=float, default=0.95,
                           help='How correlated file A and B are expected to be.'
                           'default=0.95')
    args = argParser.parse_args()
    correlations = files_rank_correlation(args.file_a, args.file_b)
    print(correlations)
    for corr in correlations:
        if corr < args.expected_correlation:
            print('Correlation below threshold, {}'.format(corr))
            sys.exit(-1)
    sys.exit(0)
