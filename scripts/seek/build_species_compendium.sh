#!/bin/bash

########
### Script to run the various steps to build a species Seek compendium
#
# Required inputs: pcl files, quant2 file refine_bio_metadata.json file
#
# Additionally, the refinebio_platforms will be downloaded from
# https://api.refine.bio/v1/platforms/
#
# Edit the first 5 lines below for the species name, NCBI location and
#   base diretory information
#
# Typical species built:
#   yeast : saccharomyces_cerevisiae
#   mouse : mus_musculus
#   worm : caenorhabditis_elegans
#   fly : drosophila_melanogaster
#   zebrafish : danio_rerio
########

SPECIES_NAME='saccharomyces_cerevisiae'
SHORT_NAME='yeast'
NCBI_GENE_INFO='https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Fungi/Saccharomyces_cerevisiae.gene_info.gz'
OUTPUT_BASE_DIR='/data/seek/tmp'
OUTPUT_DIR="${OUTPUT_BASE_DIR}/${SHORT_NAME}"

NEW_PCL_DATA='/Genomics/function/pentacon/akwong/projects/modSeek_rpc/refine_bio_nonhuman/'
SCRIPTS_DIR=$(realpath ~/src/github/FunctionLab/sleipnir/scripts/seek)
SEEK_BIN=$(realpath ~/src/github/FunctionLab/sleipnir/Debug)

# Activate conda environment
if [ -z $CONDA_DEFAULT_ENV ] || [ $CONDA_DEFAULT_ENV != "genomics" ]; then
  source ~/.bashrc
  conda activate genomics
fi

# 1. make the species dir
if [ -d ${OUTPUT_DIR} ]; then
    echo "Error: species directory '${OUTPUT_DIR}' already exists!"
    exit -1
fi
mkdir ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

# check if required files exist
if [ ! -f ../refine_bio_meta.json ]; then
    echo "Expecting file refine_bio_meta.json "
    echo  "in output base directory ${OUTPUT_BASE_DIR}"
    exit -1
fi
if [ ! -f ../quant2 ]; then
    echo "Expecting file quant2 "
    echo  "in output base directory ${OUTPUT_BASE_DIR}"
    exit -1
fi
if [ ! -f ../refine_bio_platforms.json ]; then
    # Download the platform metadata from https://api.refine.bio/v1/platforms/
    wget https://api.refine.bio/v1/platforms/ -O ../refine_bio_platforms.json
fi

cp ../quant2 .

# 2. Get the gene info
wget ${NCBI_GENE_INFO} -O ${SPECIES_NAME}.gene_info.gz
gunzip ${SPECIES_NAME}.gene_info.gz

# 2.1 Get the gene types
cut -f 2,10 ${SPECIES_NAME}.gene_info > gene_types_map.txt

# 2.2 Get the coding gene list
grep -e protein-coding -e rRNA gene_types_map.txt | cut -f 1 > coding_gene_list.txt


# 3. Parse refine_bio info and make the dsetPlatMap file
# Note: first command below outputs file dset_map.txt
python ${SCRIPTS_DIR}/parseRefineBioMetadata.py  \
    -f ../refine_bio_meta.json \
    -p ../refine_bio_platforms.json \
    -s ${SPECIES_NAME} \
    --make-dset-map -o .
cut -f 3,4 dset_map.txt > dsetPlatMap.txt


# 4. Copy the PCL files and rename them to includ platform in file name
mkdir pcl
python ${SCRIPTS_DIR}/copyDatasetPcls.py \
    -m dsetPlatMap.txt -i ${NEW_PCL_DATA}/${SPECIES_NAME}/ -o pcl/
# 4.1 Copy the dataset list and revised dsetPlatMap file from pcl dir to species dir
cp pcl/datasets.txt .
cp pcl/dsetPlatMap.txt .


# 5. Get list of genes from PCL files and limit to those also in the coding_gene_list
python ${SCRIPTS_DIR}/geneSetFromPcls.py -p ${SEEK_BIN}/PCL2Bin \
    -d pcl -l coding_gene_list.txt -o gene_map.txt


# 6. Get the corresponding symbol names for the genes, make the symbol map file
cut -f 2,3 ${SPECIES_NAME}.gene_info > gene_symbols_all.txt
# Note FNR==NR is only true while processing the first file. And the last part
# of the command (after next) is only done while processing the second file.
awk 'BEGIN { FS = "[ \t]+" }; FNR==NR {a[$1]=toupper($0); next}; $2 in a {print a[$2]}' \
    gene_symbols_all.txt gene_map.txt > gene_symbols.txt


# 7. Build the species compendium
time python ${SCRIPTS_DIR}/seekCreateDB.py --all --dab-use-gene-set \
    -b ${SEEK_BIN} -i ./ -o ./ -p ./pcl -m 30 2>&1 | tee out.txt
