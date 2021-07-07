import os
import argparse

"""
Program to rename experiment.pcl to experiment.platform.pcl
Input is a file with list of the final names to create
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--filelist', '-l', default=None,
                        help="Text file with datasets to copy")
    args = parser.parse_args()

    if args.filelist is None:
        raise Exception("Must specify -l <filelist> to rename")

    with open(args.filelist) as fp:
        while line := fp.readline():
            line = line.rstrip()
            name, plat, ext = line.split('.', 3)
            src = f'{name}.{ext}'
            dst = line
            if not os.path.exists(src):
                print(f'Skipping missing file: {src}')
                continue
            print(f'mv {src} {dst}')
            os.rename(src, dst)
