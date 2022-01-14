## **Overview:**
Sleipnir is a C++ library enabling efficient analysis, integration, mining, and machine learning over genomic data. This includes a particular focus on microarrays, since they make up the bulk of available data for many organisms, but Sleipnir can also integrate a wide variety of other data types, from pairwise physical interactions to sequence similarity or shared transcription factor binding sites.



## **Docs:**
Main documentation: https://functionlab.github.io/sleipnir-docs/

The Sleipnir wiki and bug reporting system are at: (*TBD*)

The file README.developer has notes for Sleipnir developers.

Sleipnir also includes the code to compile SEEK (the human coexpression search engine). See the link http://seek.princeton.edu/installation.jsp for information on its installation.

## **Code:**
The latest version of Sleipnir software can be obtained
by issuing the following command:

git clone https://github.com/FunctionLab/sleipnir.git

## **Build:**
1. Install g++, cmake

2. Install libraries
    - *On Mac:*
        - <code> brew install libsvm</code>
        - <code> brew install libomp</code>
        - <code> brew install thrift</code>
        - <code> brew install gsl</code>
        - <code> brew install boost</code>
    - *On CentOS Linux:*
        - <code> sudo yum install libsvm</code>
        - <code> sudo yum install libgomp</code>
        - <code> sudo yum install thrift-devel</code>
        - <code> sudo yum install gsl</code>
        - <code> sudo yum install boost</code>
    - *On Ubuntu Linux:*
        - <code> apt-get update</code>
        - <code> apt-get install build-essential</code>
        - <code> apt-get install libsvm-dev</code>
        - <code> apt-get install libomp-dev</code>
        - <code> apt-get install libthrift-dev</code>
        - <code> apt-get install libgsl-dev</code>
        - <code> apt-get install libboost-dev</code>
        - <code> apt-get install libboost-graph-dev</code>
        - <code> apt-get install libboost-regex-dev</code>
        - <code> apt-get install libreadline-dev</code>

3. Clone repository
    - <code> git clone https://github.com/FunctionLab/sleipnir.git </code>
    - <code> cd sleipnir </code>
    - <code> git submodule init </code>
    - <code> git submodule update </code>

4. Prep make files with cmake
    - <code> mkdir Debug </code>
    - <code> cd Debug/ </code>
    - <code> cmake -DCMAKE_BUILD_TYPE=Debug .. </code>
    - Alternately replace *'Debug'* with *'Release'* in all the above commands to make the release build

5. Build the code
    - (On Mac) - Edit sleipnir/src/libsvm.h
        - Replace: #include <libsvm/svm.h>
        - With: #include <svm.h>
    - <code> cd Debug/ </code>
    - <code>make </code>
        - In case of errors:
            - <code> make clean </code>
            - <code> make VERBOSE=1 </code>

## **Tests:**
0. One-time prep: create the conda environment (by default this will create the 'genomics' conda env)
    - <code>conda env create --file scripts/seek/conda_environment.yml</code>

1. Run the c++ unit tests
    - <code>Debug/tests/unit_tests</code>

2. Test the scripts for building and merging SEEK database compendiums
    - <code>conda activate genomics</code>
    - <code>python -m pytest -s -v scripts/seek/tests</code>

3. Run the SEEK system tests (test SeekMiner and SeekRPC)
    - <code>conda activate genomics</code>
    - <code>python -m pytest -s -v tests/</code>


4. Run Seek DB tests (test that the database gives expected bio-informative results). These tests can only be run where the full SEEK database is installed.
    - ```cd tests/bioinform_tests```
    - PREP: Install and init Git LFS (Large File Storage)
        - On Mac: ```brew install git-lfs```
        - On Centos: ```yum install git-lfs```
        - On Ubuntu: ```apt-get install git-lfs```
        - Initialize git-lfs: <code>git lfs install </code>
        - Refresh the gold standard tgz files (should be multipe MB in size)
            - <code> rm gold_standard_results/* </code>
            - <code> git restore gold_standard_results/* </code>
    - Run the tests:
        (The bioinform test has an option for different lengths of test, i.e. how many queries are run)
        - <code>bash run_paramtest.sh -v -s <path_to_seek_db> -b <path_to_seek_binaries> </code>
        - <code>bash run_querysize.sh -v -s <path_to_seek_db> -b <path_to_seek_binaries> </code>
        - <code>bash run_bioinform.sh -v -s <path_to_seek_db> -b <path_to_seek_binaries> -t [tiny,short,medium,long]</code>


