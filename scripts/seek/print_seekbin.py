import sys
import argparse
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', '-f', default=None,
                        help="binary file to parse")
    parser.add_argument('--outfile', '-o', default='binfile.txt',
                        help="name of output file to write results")
    args = parser.parse_args()

    if args.file is None:
        print('Must supply name of file to parse, using: -f <binfile>')
        sys.exit(-1)

    vals = []
    with open(args.file, 'rb') as fp:
        # The first 8 byte (long int) is the number of elements stored
        headerVals = np.fromfile(fp, count=1, dtype=np.ulonglong)
        numVals = headerVals[0]
        # The remaining are 4 byte float values, numVal of them
        vals = np.fromfile(fp, dtype=np.float32)
        assert len(vals) == numVals

    with open(args.outfile, 'w') as fp:
        for val in vals:
            fp.write(f'{val:.06f}\n')
