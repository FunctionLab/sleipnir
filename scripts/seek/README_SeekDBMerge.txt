Steps to merge new datasets into and existing seek database

The overall approach is to have two databases, a large one and a small one. SeekServer and SeekMiner can query across multiple databases. To enable multi-database mode a file is passed in with the config settings of the secondary databases. See SeekServer additional DB file example below.

When new datasets are to be added, they are merged into the small database to reduce the merge time. When the small database exceeds 10% of the large database size it is merged into the large database and an empty small database is started.

A number of helper scripts are added under sleipnir/scripts/seek directory to help with the incremental merge process.

To merge in new datasets, do the following:
1) Make a directory with the PCL files and a dataset_list.txt file of the new datasets
2) Activate a conda environment with python 3.x enabled (use conda_environment.yml to create the conda env)
3) Run the incremental merge script, it will:
  - Create a tmp new database out of the new PCL files.
  - Merge the new database into the small datbase, including metadata files.
  - Verify that the merge was successful (combined db files match inputs)
Note: the small and large dbs are referenced to make sure there is no overlap in the datasets and to use the same number of DB files.

Example Usage:
  ```conda activate genomics```
  ```python seekIncrementalMerge.py -p new/pcl/ -dn new/dataset.description.txt -ds small_db/dataset.description.txt -dl large_db/dataset.description.txt -s small_db/ -l large_db/ -b $sleipnir/Debug -o ./merged```

4) Move the original small DB directory to small_db.prev.num
   Rename the new merged db directory to small_db


How to configure Seek to use multiple DBs:
SeekServer Additional DB Settings File Example:
START
SINFO_DIR       /data/db_small/sinfo
GVAR_DIR        NA
PCL_DIR         data/db_small/pcl
PREP_DIR        data/db_small/prep
PLATFORM_DIR    data/db_small/plat
DB_DIR          data/db_small/db
DSET_MAP_FILE   data/db_small/dataset.map
GENE_MAP_FILE   data/db_small/gene_map.txt
QUANT_FILE      data/db_small/quant2
NUMBER_OF_DB    1000
END