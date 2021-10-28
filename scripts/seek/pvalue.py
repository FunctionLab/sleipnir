"""
P-value steps
1) Run 10,000 random queries 0-9999 (use the existing ones for now)
2) From python open the .gscore files and read in the values
3) Collect the gene scores per gene (i.e. in an nparray)
4) Look at the distribution, are they clustered or evenly spread out
5) Make a binning scheme that evenly distributes the scores among bins. May need to keep a separate binning scheme for each gene
6) Write out the binning scheme and bin counts for each gene. (what format?)
7) Read in the binning scheme and bin counts to prepare to serve p-values
8) To get p-value, walk the binning scheme until the target gene score is > the bin boundary, then divide num unwalked (remaining) bins by total bins.

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
Then (1 - bindIdx / totalBins) to get the p-value number.

"""
import os
import sys
import time
import numpy as np
import glob
import argparse

"""
The general strategy will be to group the results of a set of random queries into
a set of buckets or 'bins'. We will use 1000 bins by default to get a precision
of .001. The scores will be equally distributed to bins, so that the p-value can
be calculated by finding the bin that a score would belong in and then calculating
1 - binIndex / totalBins. So for example a score that
would go in the 100th bin would have a p-value of 1 - 100/1000 or 0.9. To implement
this we must calculate and store the bin boundaries so that we can later determine
which bin a gene score would fall into.
"""
class PValue():
    def __init__(self, numBins, randomDir):
        self.numBins = numBins
        self.randomDir = randomDir
        self.saveDir = '/tmp'
        self.numGenes = None
        self.geneScoreBounds = None
        self.geneRankBounds = None


    def calcGeneScoreBoundaries(self, randomTrialData):
        """
        Calculate the bin boundaries given a set of gene score values from random queries.
        Arguments:
            data: A 2D numpy array, each row represents a random query trial and
                cols are the genes (in gene_map.txt order) with their scores.

        Returns:
            Gene score boundaries with early boundaries having highest scores,
            i.e. boundaries are reverse sorted.
            None - instead writes out a file of the bin boundaries
        """
        # The randomTrialData: rows are random trials and columns are genes
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
        itemsPerBin = int(numTrials / self.numBins)
        geneScoreBounds = np.empty((numGenes, self.numBins+1), dtype=np.float16)
        for geneIdx in range(numGenes):
            scores = data[geneIdx]
            boundaryScores = []
            # Create the boundaries starting from the highest score,
            #  use negative indexing to do this, so iterate from 1 to len+1
            for idx in range(1, len(scores)+1, itemsPerBin):
                # Insert scores starting from the end of the array (neg indexing)
                boundaryScores.append(np.float16(scores[-idx]))
            boundaryScores.append(np.float16(scores[0]))
            assert len(boundaryScores) == self.numBins+1
            bounds = np.array(boundaryScores, dtype=np.float16)
            geneScoreBounds[geneIdx] = bounds

        return geneScoreBounds


    def calcGeneRankBoundaries(self, randomTrialData):
        """
        Calculate the bin boundaries based on gene rank from a set of random queries.
        Arguments:
            data: A 2D numpy array, each row represents a random query trial and
                cols are the genes (in gene_map.txt order) with their scores.

        Returns:
            Gene rank boundaries with early boundaries having highest rank (lowest value),
            None - instead writes out a file of the bin boundaries
        """

        # First determine rank that each score would give a gene in the random query trials
        numTrials, numGenes = np.shape(randomTrialData)
        rankData = np.empty((numTrials, numGenes))
        for idx in range(numTrials):
            # np.argsort will give the indicies that would sort the array
            # TODO - this isn't working the way you think, use scipy.stats.rankdata instead
            rankData[idx] = np.argsort(randomTrialData[idx])

        # Next transpose the data and sort each genes rank scores in acending order
        data = rankData.transpose()
        for idx in range(numGenes):
            data[idx] = np.sort(data[idx])

        # Calculate the bin boundaries
        itemsPerBin = int(numTrials / self.numBins)
        geneRankBounds = np.empty((numGenes, self.numBins+1), dtype=np.float16)
        for geneIdx in range(numGenes):
            ranks = data[geneIdx]
            boundaryRanks = []
            for idx in range(0, len(ranks), itemsPerBin):
                boundaryRanks.append(np.float16(ranks[idx]))
            boundaryRanks.append(np.float16(ranks[-1]))
            assert len(boundaryRanks) == self.numBins+1
            bounds = np.array(boundaryRanks, dtype=np.float16)
            geneRankBounds[geneIdx] = bounds

        return geneRankBounds


    def createBoundsArray(self):
        # A .gscore file is the score for each gene in the gene order specified in the gene_map.txt file
        # get list of all gscore files
        fileList = glob.glob(os.path.join(self.randomDir, '*.gscore'))

        # get number of gene scores entries in each file
        numGenes = 0
        with open(fileList[0], 'rb') as fp:
            # The first 8 byte long int is the number of elements stored
            vals = np.fromfile(fp, count=1, dtype=np.ulonglong)
            numGenes = vals[0]

        self.numGenes = numGenes

        print(f'Num files: {len(fileList)}, Num genes: {numGenes}')
        data = np.empty((len(fileList), numGenes))

        for idx, file in enumerate(fileList):
            with open(file, 'rb') as fp:
                # The first 8 byte long int is the number of elements (gene scores) stored
                vals = np.fromfile(fp, count=1, dtype=np.ulonglong)
                numElems = vals[0]
                # The remaining are 4 byte float values, i.e. the gene scores for each gene
                row = np.fromfile(fp, dtype=np.float32)
            assert len(row) == numElems
            # A row contains the random gene scores for one query, one score per gene
            data[idx] = row
            # dataT[:, idx] = row  # This performance was twice as slow as the row based one above

        # At this point the data rows are random trials and columns are genes
        numTrials, numCols = np.shape(data)
        assert numCols == numGenes

        self.geneScoreBounds = self.calcGeneScoreBoundaries(data)
        self.geneRankBounds = self.calcGeneRankBoundaries(data)
        # Write arrays out to a file
        np.save(os.path.join(self.saveDir, 'pval-gscore-bins.npy'), self.geneScoreBounds)
        np.save(os.path.join(self.saveDir, 'pval-grank-bins.npy'), self.geneRankBounds)
        np.save(os.path.join(self.saveDir, 'pval-num-genes.npy'), self.numGenes)
        return

    def loadBounds(self):
        self.geneScoreBounds = np.load(os.path.join(self.saveDir, 'pval-gscore-bins.npy'))
        self.geneRankBounds = np.load(os.path.join(self.saveDir, 'pval-grank-bins.npy'))
        numGenesArr = np.load(os.path.join(self.saveDir, 'pval-num-genes.npy'))
        self.numGenes = numGenesArr[0]

    def getGeneScorePval(self, boundaries, geneId, score):
        # TODO if the score is zero this current implementation would return 0.5
        #  Seems like we want to multiply by 2x what is calculated here?
        row = boundaries[geneId]
        # Note the boundaries are descending sorted (highest to lowest)
        if score >= 0:
            for idx, bound in enumerate(row):
                if score >= bound:
                    pval = idx / self.numBins
                    return pval
        else:
            idx = 0
            for bound in reversed(row):
                if score <= bound:
                    pval = idx / self.numBins
                    return pval
                idx += 1
        return 1.0

    def getGeneRankPval(self, boundaries, geneId, rank):
        row = boundaries[geneId]
        assert self.numGenes is not None
        # Note the boundaries are ascending sorted (lowest to highest)
        if rank < self.numGenes / 2:
            for idx, bound in enumerate(row):
                if rank <= bound:
                    pval = idx / self.numBins
                    return pval
        else:
            idx = 0
            for bound in reversed(row):
                if rank >= bound:
                    pval = idx / self.numBins
                    return pval
                idx += 1
        return 1.0

    def calcScorePvalues(self, scores):
        assert self.geneScoreBounds is not None
        for gScore in scores:
            geneId, score = gScore
            pval = self.getGeneScorePval(self.geneScoreBounds, geneId, score)
            print(f'gene {geneId}, score {score}, pval {pval}')


    def calcRankPvalues(self, ranks):
        assert self.geneRankBounds is not None
        for gRank in ranks:
            geneId, rank = gRank
            pval = self.getGeneRankPval(self.geneRankBounds, geneId, rank)
            print(f'gene {geneId}, rank {rank}, pval {pval}')



# End PValue Class

if __name__ == "__main__":
    argParser = argparse.ArgumentParser()
    argParser.add_argument('--genes', '-g', default=None, type=str,
                           help='list of geneIds for the corresponding scores')
    argParser.add_argument('--scores', '-s', default=None, type=str,
                           help='list of scores to get pval for')
    args = argParser.parse_args()

    numBins = 1000
    randomDir = "/data/seek/Seek/random_all"

    pvalClass = PValue(numBins, randomDir)

    # gene_scores = [[1123, 0.49]]
    # pvalClass.calcPvalues(gene_scores)
    gene_ranks = [[]]
    # pvalClass.createBoundsArray()



#if np.array_equal(data, dataT):
#    print('they are equal')

# If I read in all the data, then that is an array where
#  rows are trials and colums genes
# So read in all the rows and transpose the array
# Or read in all rows and make a pandas and pull off a column
# Then I need to process per gene and sort the random score values and create the bins
#
