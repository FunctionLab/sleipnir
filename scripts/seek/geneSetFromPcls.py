"""
Program that parses pcl or pcl.bin files and extracts the set of all genes contained.
If a gene-info file is provided it will limit the set to genes of type protein-coding
and rRNA.
"""
import os
import sys
import argparse
import subprocess


def getGenesFromPcl(filename, gene_set):
    with open(filename) as fp:
        while line := fp.readline():
            gene = line[:line.index("\t")]
            gene = gene.strip()
            gene_set.add(gene)


def getGenesfromPclBin(filename, gene_set, pcl2BinCmd):
    process = subprocess.Popen([pcl2BinCmd, '-E', '-i', filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # communicate blocks until the process completes
    out, err = process.communicate()
    result = out.decode('utf-8')
    lines = result.split('\n')
    for line in lines:
        if len(line.strip()) == 0:
            continue
        idx, gene = line.split()
        gene = gene.strip()
        gene_set.add(gene)


def readGeneInfoFile(filename):
    categoryMap = {}
    with open(filename) as fp:
        while line := fp.readline():
            line = line.strip()
            try:
                id, sid, category = line.split('\t')
            except Exception as err:
                print(f'error spliting line {line}')
                raise err
            categoryMap[id] = category
    return categoryMap


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', '-d', default=None, action='append', required=True,
                        help="directory of .pcl or .pcl.bin files to parse, can specify -d multiple times for multiple directories")
    parser.add_argument('--geneInfoFile', '-g', default=None,
                        help="Trim set of genes to protein-encoding and rRNa from this list")
    parser.add_argument('--geneLimitList', '-l', default=None,
                        help="Limit set of genes to ones in this list")
    parser.add_argument('--pcl2bin', '-p', default=None, required=True,
                        help="path to the pcl2bin program")
    parser.add_argument('--outfile', '-o', default=None,
                        help="path write the results")
    args = parser.parse_args()

    if args.dir is None:
        raise Exception("Must specify -d <directory> to parse")

    gene_set = set()

    # check all dirs exist
    for dir in args.dir:
        if not os.path.exists(dir):
            raise FileNotFoundError(f"Directory {dir} not found")

    for dir in args.dir:
        for file in os.listdir(dir):
            fullpath = os.path.join(dir, file)
            if os.path.exists(fullpath):
                if os.path.isfile(fullpath):
                    file_ext = os.path.splitext(fullpath)[1]
                    if file_ext == '.bin':
                        getGenesfromPclBin(fullpath, gene_set, args.pcl2bin)
                        pass
                    elif file_ext == '.pcl':
                        getGenesFromPcl(fullpath, gene_set)
                    else:
                        print(f'Error: file extension not recognized {file_ext}, {fullpath}')
            else:
                print(f'Error: file not found: {fullpath}')

    gene_list = list(gene_set)
    gene_list.sort()

    # load geneInfo file if specified
    if args.geneInfoFile is not None:
        categoryMap = readGeneInfoFile(args.geneInfoFile)
        filteredList = []
        for gene in gene_list:
            category = categoryMap.get(gene, 'Missing')
            # print(f'{gene} {category}')
            if category in ('protein-coding', 'rRNA'):
                filteredList.append(gene)
        gene_list = filteredList

    if args.geneLimitList is not None:
        with open(args.geneLimitList) as fp:
            limitGenes = fp.read().splitlines()
        filteredList = set(limitGenes).intersection(gene_list)
        gene_list = list(filteredList)
        gene_list.sort(key=int)

    # outpult the results
    fp = sys.stdout
    if args.outfile is not None:
        fp = open(args.outfile, 'w')
    for idx, gene in enumerate(gene_list):
        print(f'{idx+1}\t{gene}', file=fp)
    if args.outfile is not None:
        fp.close()
