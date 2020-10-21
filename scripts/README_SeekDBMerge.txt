Steps to merge new datasets into and existing seek database

The overall approach is to have two databases, a large one and a small one. SeekServer and SeekMiner can query across multiple databases. To enable multi-database mode a file is passed in with the config settings of the secondary databases. See SeekServer additional DB file example below.

When new datasets are to be added, they are merged into the small database to reduce the merge time. When the small database exceeds 10% of the large database size it is merged into the large database and an empty small database is started.

A number of helper scripts are added under sleipnir/scripts directory to help with the incremental merge process.

To merge in new datasets, do the following:
1) Make a directory with the PCL files and a dataset_list.txt file
2) Activate a conda environment with python 3.x enabled
3) Run the incremental merge script to create a new database out of the new PCL files. When it concludes it will indicate the command to run to then merge the new database into the small database.

  ```python $sleipnir/scripts/seekIncrementalMerge.py -p new/pcl/ -dn new/dataset.description.txt -ds small_db/dataset.description.txt -dl large_db/dataset.description.txt -s small_db/ -l large_db/ -b $sleipnir/Debug -o ./merged```

At the conclusion of that script it will specify the command to merge the new database into the small database, similar to:

  ```bash $sleipnir/scripts/seekCombineDB.sh -d small_db/ -n incr_dset_20201021_144003/ -o ./merged -b $sleipnir/Debug```

You can then run a verify script to check that the *.db files have been merged correctly:

  ```bash $sleipnir/scripts/seekVerifyMergeDB.sh -d small_db/db -n incr_dset_20201021_144003/db -c merged/db -b $sleipnir/Debug```

4) Move the original small DB directory to small_db.prev.num
   Rename the new merged db directory to small_db



SeekServer Additional DB Settings File Example:
START
SINFO_DIR       /data/db_small/sinfo
GVAR_DIR        NA
PREP_DIR        data/db_small/prep
PLATFORM_DIR    data/db_small/plat
DB_DIR          data/db_small/db
DSET_MAP_FILE   data/db_small/dataset.map
GENE_MAP_FILE   data/db_small/gene_map.txt
QUANT_FILE      data/db_small/quant2
NUMBER_OF_DB    1000
END