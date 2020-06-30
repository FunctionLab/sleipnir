'''
Utility helper functions for common use by other scripts
'''
import os
import re
import toml
import json
import pathlib
from struct_dict import StructDict, recurseCreateStructDict


def file_read(filename):
    with open(filename, 'r') as fp:
        lines = fp.readlines()
    return lines


def file_appendline(filename, data):
    with open(filename, 'a') as fp:
        fp.write(data+'\n')


def file_truncate(filename):
    if os.path.exists(filename):
        os.truncate(filename, 0)


def checkAndMakePath(path):
    if not os.path.exists(path):
        os.mkdir(path)


def load_config_file(filename):
    file_suffix = pathlib.Path(filename).suffix
    if file_suffix == '.json':
        # load json
        with open(filename) as fp:
            cfg_dict = json.load(fp)
        # to write out config
        # with open('t1.json', 'w+') as fd:
        #   json.dump(cfg, fd, indent=2)
    elif file_suffix == '.toml':
        # load toml
        cfg_dict = toml.load(filename)
        # to write out config
        # with open('t1.toml', 'w+') as fd:
        #  toml.dump(cfg, fd)
    else:
        raise ValueError("Config file format not recognized (expecting .json or .toml)")
    cfg_struct = recurseCreateStructDict(cfg_dict)
    return cfg_struct


def read_genes(gene_file):
    '''
    Read in file with sets of genes, one set per line.
    Return a list of lists containing the gene sets.
    '''
    groups = []
    with open(gene_file, 'r') as fp:
        groups = [line.rstrip("\n") for line in fp]
    return groups


def write_goldstd(goldstd, outdir):
    '''
    Take as input a file with a list of gold standard genes one per line
    corresponding to the query genes one per line in a query file. Convert
    this to a set of files numbered 0 to numlines where each output file
    has one line worth of gold standard genes corresponding to one query.
    '''
    for i, g in enumerate(goldstd):
        fw = open("%s/%d.gold" % (outdir, i), "w")
        fw.write(g + "\n")
        fw.close()


def parse_gmt(filename, min_gene_count, max_gene_count):
    """
    Parse a file that has a correlated gene group per line.
    Only keep groups that meet the min, max gene count criteria.
    """
    groups = []
    with open(filename, 'r') as fp:
        for line in fp:
            vals = line.rstrip(os.linesep).split("\t")
            group_id = vals[0]
            desc = vals[1]
            genes = vals[2:]
            mtch = re.search(r'(.*)\((\d+)\)', desc)
            if mtch is None:
                raise ValueError("Description missing gene count, {}".format(desc))
            num_genes = int(mtch.group(2))
            desc = mtch.group(1).rstrip()
            assert len(genes) == num_genes
            if num_genes >= min_gene_count and num_genes <= max_gene_count:
                group = StructDict({'id': group_id, 'desc': desc, 'genes': genes, 'size': num_genes, 'line': line})
                groups.append(group)
    return groups
