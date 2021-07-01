"""
P-value steps
1) Run 10,000 random queries 0-9999 (use the existing ones for now)
x - 2) From python open the .gscore files and read in the values
3) Collect the gene scores per gene (i.e. in an array, maybe nparray)
4) For starters - look at the distribution, are they clustered or evenly spread out
5) Make a binning scheme that evenly distributes the scores among bins. May need to keep a separate binning scheme for each gene
6) Write out the binning scheme and bin counts for each gene. (what format?)
7) Read in the binning scheme and bin counts to prepare to serve p-values
8) To get p-value, walk the binning scheme until the target gene score is > the bin boundary, then add up the bin counts in those bins and divide by total count in all bins.

Notes:
Load the array either as 
data[idx] = row
or read into data[:][colIdx]

Sort each gene's scores
Divide the number of scores by number of bins (i.e. 10000/1000 = 10 scores per bin)
Walk the scores, 10 at a time, the bin boundary will be set to halfway between the scores on either side
Store the bin boundaries as an array
Write the array out (HDF5? or other?)
Read in the array
To get the p-value for a gene, walk the bin boundaries to see which bin the gene falls into
Then divide that bin number by total bins to get the p-value number.

"""
import os
import sys
import time
import numpy as np
import glob
import argparse

numBins = 1000
randomDir = "/data/seek/Seek/random_all"

"""
The general strategy will be to group the results of a set of random queries into
a set of buckets or 'bins'. We will use 1000 bins by default to get a precision
of .001. The scores will be equally distributed to bins, so that the p-value can
be calculated by finding the bin that a score would belong in and then dividing
the bin offset (i.e. index) by total number of bins. So for example a score that
would go in the 100th bin would have a p-value of 100/1000 or 0.1. To implement
this we must calculate and store the bin boundaries so that we can later determine
which bin a gene score would fall into.
"""

def calcGeneScoreBoundaries(randomTrialData):
    """
    Calculate the bin boundaries given a set of gene score values from random queries.
    Arguments:
        data: A 2D numpy array, each row represents a random query trial and 
              cols are the genes (in gene_map.txt order) with their scores.

    Returns: None - instead writes out a file of the bin boundaries
    """
    # The randomTrialData rows are random trials and columns are genes
    # Transpose to be each row is a gene and columns are random trials
    data = randomTrialData.transpose()
    numGenes, numTrials = np.shape(data)

    # Sort each gene (row) by ascending score
    startTime = time.time()
    for idx in range(numGenes):
        data[idx] = np.sort(data[idx])
        # to sort by descending order
        # data[idx][::-1].sort()
        # if using descending, then change the for loop iteration below when creating boundaries
        # https://stackoverflow.com/questions/26984414/efficiently-sorting-a-numpy-array-in-descending-order
    print('Sort time: {}s'.format(time.time() - startTime))

    # Calculate the bin boundaries
    itemsPerBin = int(numTrials / numBins)
    geneScoreBounds = np.empty((numGenes, numBins+1), dtype=np.float16)
    for geneIdx in range(numGenes):
        scores = data[geneIdx]
        boundaryScores = []
        for idx in range(1, len(scores)+1, itemsPerBin):
            boundaryScores.append(np.float16(scores[-idx]))
        boundaryScores.append(np.float16(scores[0]))
        assert len(boundaryScores) == numBins+1
        bounds = np.array(boundaryScores, dtype=np.float16)
        geneScoreBounds[geneIdx] = bounds

    # Write array out to a file
    np.save('/tmp/pval-gscore-bins.npy', geneScoreBounds)
    return


def calcGeneRankBoundaries(randomTrialData):
    """
    Calculate the bin boundaries based on gene rank from a set of random queries.
    Arguments:
        data: A 2D numpy array, each row represents a random query trial and 
              cols are the genes (in gene_map.txt order) with their scores.

    Returns: None - instead writes out a file of the bin boundaries
    """

    # First determine rank that each score would give a gene in the random query trials
    numTrials, numGenes = np.shape(randomTrialData)
    rankData = np.empty((numTrials, numGenes))
    for idx in range(numTrials):
        # np.argsort will give the index each element would be if the list were sorted
        rankData[idx] = np.argsort(randomTrialData[idx])
    
    # Next transpose the data and sort each genes rank scores in acending order
    data = rankData.transpose()
    for idx in range(numGenes):
        data[idx] = np.sort(data[idx])

    # Calculate the bin boundaries
    itemsPerBin = int(numTrials / numBins)
    geneRankBounds = np.empty((numGenes, numBins+1), dtype=np.float16)
    for geneIdx in range(numGenes):
        ranks = data[geneIdx]
        boundaryRanks = []
        for idx in range(0, len(ranks), itemsPerBin):
            boundaryRanks.append(np.float16(ranks[idx]))
        boundaryRanks.append(np.float16(ranks[-1]))
        assert len(boundaryRanks) == numBins+1
        bounds = np.array(boundaryRanks, dtype=np.float16)
        geneRankBounds[geneIdx] = bounds

    # Write array out to a file
    np.save('/tmp/pval-grank-bins.npy', geneRankBounds)
    return


def createBoundsArray():
    # A .gscore file is the score for each gene in the gene order specified in the gene_map.txt file
    # get list of all gscore files
    fileList = glob.glob(os.path.join(randomDir, '*.gscore'))

    # get number of gene scores entries in each file
    numGenes = 0
    with open(fileList[0], 'rb') as f:
        # The first 8 byte long int is the number of elements stored
        vals = np.fromfile(f, count=1, dtype=np.ulonglong)
        numGenes = vals[0]

    print(f'Num files: {len(fileList)}, Num genes: {numGenes}')
    data = np.empty((len(fileList), numGenes))

    for idx, file in enumerate(fileList):
        with open(file, 'rb') as f:
            # The first 8 byte long int is the number of elements (gene scores) stored
            vals = np.fromfile(f, count=1, dtype=np.ulonglong)
            numElems = vals[0]
            # The remaining are 4 byte float values, i.e. the gene scores for each gene
            row = np.fromfile(f, dtype=np.float32)
        assert len(row) == numElems
        # A row contains the random gene scores for one query, one score per gene
        data[idx] = row
        # dataT[:, idx] = row  # This performance was twice as slow as the row based one above

    # At this point the data rows are random trials and columns are genes
    numTrials, numRows = np.shape(data)
    assert numRows == numGenes

    calcGeneScoreBoundaries(data)
    calcGeneRankBoundaries(data)
    return


def getGeneScorePval(boundaries, geneId, score):
    row = boundaries[geneId]
    # import pdb; pdb.set_trace()
    for idx, bound in enumerate(row):
        if score > bound:
            pval = idx / numBins
            return pval
    return 1.0

def getGeneRankPval(boundaries, geneId, rank):
    row = boundaries[geneId]
    # import pdb; pdb.set_trace()
    for idx, bound in enumerate(row):
        if rank < bound:
            pval = idx / numBins
            return pval
    return 1.0

def calcScorePvalues(scores):
    geneScoreBounds = np.load('/tmp/pval-gscore-bins.npy')
    for gScore in scores:
        geneId, score = gScore
        pval = getGeneScorePval(geneScoreBounds, geneId, score)
        print(f'gene {geneId}, score {score}, pval {pval}')


def calcRankPvalues(ranks):
    geneRankBounds = np.load('/tmp/pval-grank-bins.npy')
    for gRank in ranks:
        geneId, rank = gRank
        pval = getGeneRankPval(geneRankBounds, geneId, rank)
        print(f'gene {geneId}, rank {rank}, pval {pval}')  


if __name__ == "__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--genes', '-g', default=None, type=str,
                           help='list of geneIds for the corresponding scores')
    argParser.add_argument('--scores', '-s', default=None, type=str,
                           help='list of scores to get pval for')
    args = argParser.parse_args()
    # gene_scores = [[1123, 0.49]]
    # calcPvalues(gene_scores)
    gene_ranks = [[]]
    # createBoundsArray()



#if np.array_equal(data, dataT):
#    print('they are equal')

# If I read in all the data, then that is an array where
#  rows are trials and colums genes
# So read in all the rows and transpose the array
# Or read in all rows and make a pandas and pull off a column
# Then I need to process per gene and sort the random score values and create the bins
#
