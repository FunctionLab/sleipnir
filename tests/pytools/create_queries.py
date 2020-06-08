"""
This scirpt takes in a .gmt file of correlated sets of genes and creates 3 output
files. 1) base.query.txt, 2) base.goldStd.txt, 3) base.group.txt. The query.txt
files has the query string genes, the goldStd has the corresponding correlated
genes that should be found by the query, and the group.txt identifies the .gmt
group that each query is based on. The query_pct values indicate what percentage
of the group genes to use as a query. You can specify multiple query_pct values
as a comma separated list or with mulitple -q args.
"""
import os
import re
import math
import random
import argparse
from struct_dict import StructDict
from utils import parse_gmt


def create_queries(groups, query_count, max_query_genes, outBaseFilename,
                   query_pct=None, max_groups=None):
    """
    Create queries out of groups of biologically informative gene sets.
    Each group is a set of genes that are correlated with each other.
    One query line will be made for each query_count in the list, the query_count
    indicates what how many of the gene group to use as the query, the rest
    will be the target result genes.
    """
    queryFH = open(outBaseFilename + '.query.txt', 'w')
    goldStdFH = open(outBaseFilename + '.goldStd.txt', 'w')
    groupFH = open(outBaseFilename + '.group.txt', 'w')
    num_groups = 0
    for group in groups:
        if max_groups is not None:
            if num_groups >= max_groups:
                break
        num_groups += 1
        # if query_pct specified, translate that into the number of genes for each query
        if query_pct is not None and len(query_pct) > 0:
            num_query_genes = [(group.size * qpct / 100.0) for qpct in query_pct]
            if max(num_query_genes) > max_query_genes:
                # Using the max(query_pct), calculate a group.size that will keep the
                # number of query genes below max_query_genes
                max_qpct = max(query_pct)
                gsize = max_query_genes * 100.0 / max_qpct
                num_query_genes = [(gsize * qpct / 100.0) for qpct in query_pct]
            # Append the counts derived from pct to the count array
            query_count.extend(num_query_genes)
        group.queries = []
        for qnum in query_count:
            query = StructDict()
            if qnum < 10:
                qnum = math.ceil(qnum)
            else:
                qnum = math.floor(qnum)
            assert qnum <= max_query_genes
            assert qnum <= group.size - 1
            query.qnum = qnum
            # create a random list of index values for the query genes
            query.q_indicies = random.sample(range(group.size), qnum)
            # create the inverse list for the results genes
            query.res_indicies = list(set(range(group.size)) - set(query.q_indicies))
            assert len(query.q_indicies) == query.qnum
            assert len(query.res_indicies) == group.size - query.qnum
            assert set(query.q_indicies).isdisjoint(set(query.res_indicies))
            assert len(group.genes) == group.size
            assert max(query.q_indicies) < group.size
            assert max(query.res_indicies) < group.size
            # select the query genes and result genes
            query.genes = [group.genes[i] for i in query.q_indicies]
            query.results = [group.genes[i] for i in query.res_indicies]
            # write the query genes to a file
            queryFH.write("    ".join(query.genes) + "\n")
            # write the result genes to a file
            goldStdFH.write("    ".join(query.results) + "\n")
            # write the group to a file
            groupFH.write("\t".join([group.id, group.desc]) + "\n")
            # groupFH.write("\t".join(group.genes) + "\n")
            group.queries.append(query)
        # end for query percents
    # end for groups
    queryFH.close()
    goldStdFH.close()
    groupFH.close()


if __name__ == "__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--min-gene-count', '-min', default="20", type=int,
                           help='min number of genes in a correlated group, default 20')
    argParser.add_argument('--max-gene-count', '-max', default="200", type=int,
                           help='max number of genes in a correlated group, default 200')
    argParser.add_argument('--query-pct', '-q', default=[], action='append',
                           help='percent of genes to keep as query, can specify multiple')
    argParser.add_argument('--query-count', '-c', default=[], action='append',
                           help='number of genes to use as query, can specify multiple')
    argParser.add_argument('--input-file', '-i', type=str, required=True,
                           help='input file with list of correlated genes')
    argParser.add_argument('--output-file-base', '-o', type=str, required=True,
                           help='base filename for output files of query and gold standard genes')
    argParser.add_argument('--max-query-genes', default="50", type=int,
                           help='max number of genes allowed per query (seekcentral limit), default 50')
    argParser.add_argument('--max-groups', default=None, type=int,
                           help='max number of groups or queries ')
    args = argParser.parse_args()

    # Get the percentage of genes from a group to use for the query from the args
    query_pct = [int(y) for x in args.query_pct for y in x.split(',')]
    query_count = [int(y) for x in args.query_count for y in x.split(',')]
    print("min genes {}, max genes {}, query_pct {}, query_count {}".
        format(args.min_gene_count, args.max_gene_count, query_pct, query_count))

    groups = parse_gmt(args.input_file, args.min_gene_count, args.max_gene_count)
    create_queries(groups, query_count, args.max_query_genes, args.output_file_base, query_pct=query_pct, max_groups=args.max_groups)
    # print(groups)
