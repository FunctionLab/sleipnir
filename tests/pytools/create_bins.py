'''
Program to create the quantization bins for gene correlation values.
Take in a value range such as 0,1 and the number of bins.
example: create_bins.py --range 0,1 --num-bins 255
Note: You will generally want to specify one less bins than the int type can hold
to allow for a N/A or missing value.
'''
import sys
import argparse


def create_bins(bin_range, num_bins):
    # print("Creating {} bins in range {}".format(num_bins, bin_range))
    gap_size = (bin_range[1] - bin_range[0]) / float(num_bins-2)
    sep = ""
    for i in range(num_bins-1):
        print("{}{:.2f}".format(sep, i*gap_size), end ="")
        sep = ", "
    print("")


if __name__ == "__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--range', '-r', default="0,1", type=str,
                           help='numeric range of bins, default 0,1')
    argParser.add_argument('--num-bins', '-n', default=255, type=int,
                           help='number of bins')
    args = argParser.parse_args()

    # num bins must be > 2
    if args.num_bins <= 2:
        print("Invocation Error: --num-bins must be > 2")
        sys.exit(-1)

    bin_range = [int(s) for s in args.range.split(',')]
    create_bins(bin_range, args.num_bins)
