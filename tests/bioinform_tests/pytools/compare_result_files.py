import sys
import numpy as np
import argparse

np.set_printoptions(suppress=True, precision=4, floatmode='maxprec_equal')


def smape(A, F):
    '''
    Symmetric mean absolute percentage error (SMAPE).
    A and F are numpy arrays of values to be compared.
    Find the percent difference between the two arrays for each data point,
    where the percent diff is calculated against their average value. Then return
    the average percent difference.
    '''
    return 100/A.size * np.sum(2 * np.abs(F - A) / (np.abs(A) + np.abs(F)))


def smpe(A, F):
    '''
    Calculate the percent difference of set of values by calculating the percent
    difference relative to each value and multiplying them together to make it
    a symmetric percentage. Then take the square root and normalize by the number
    of elements.
    (symmetric - i.e. rather than calculating the percent difference vs one or
    the other it is against both). Essentially this is the percent difference
    against one times the percent diff against the other multiplied together.
    '''
    return 100/A.size * np.sum(((F - A)**2 / (F * A))**0.5)


def get_pct_error(fileA, fileB, skiprows=0, skipcols=0, delim=' ', func='smape'):
    A = np.loadtxt(fileA, delimiter=delim, skiprows=skiprows)
    B = np.loadtxt(fileB, delimiter=delim, skiprows=skiprows)
    if skipcols > 0:
        # for i in range(skipcols):
        column_list = list(range(skipcols))
        A = np.delete(A, column_list, 1)
        B = np.delete(B, column_list, 1)
    # print(A)
    # print(B)
    if func == 'smape':
        pcterr = smape(A, B)
    elif func == 'smpe':
        pcterr = smpe(A, B)
    return pcterr


if __name__ == "__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--file-a', '-a', type=str, required=True,
                           help='comparison file A')
    argParser.add_argument('--file-b', '-b', type=str, required=True,
                           help='comparison file B')
    argParser.add_argument('--pct-tolerance', '-p', type=float, default=10.0,
                           help='How close file A and B are expected to be in'
                           ' percent difference, default=10.0')
    args = argParser.parse_args()
    pctdiff = get_pct_error(args.file_a, args.file_b, skiprows=1, skipcols=3)
    print(pctdiff)
    if pctdiff > args.pct_tolerance:
        sys.exit(-1)
    sys.exit(0)
