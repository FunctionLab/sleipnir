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
    - *On Linux:*

3. Clone repository
    - <code> git clone https://github.com/FunctionLab/sleipnir.git </code>
    - <code> cd sleipnir </code>
    - <code> git submodule init </code>
    - <code> git submodule update </code>

4. Prep make files with cmake
    - <code> mkdir Debug </code>
    - <code> cd Debug </code>
    - <code> cmake -DCMAKE_BUILD_TYPE=Debug .. </code>

5. Build the code
    - (On Mac) - Edit sleipnir/src/libsvm.h
        - Replace: #include <libsvm/svm.h>
        - With: #include <svm.h>
    - <code>make </code>
        - In case of errors:
            - <code> make clean </code>
            - <code> make VERBOSE=1 </code>

## **Tests:**
1. Run the c++ unit tests
    - <code>Debug/tests/unit_tests </code>

2. Run Seek system tests
    - <code>cd scripts/seek </code>
    - Init conda environment (first time only)
        - <code>conda env create --file conda_environment.yml </code>
    - <code>conda activate genomics </code>
    - <code>python -m pytest -s -v tests/ </code>

3. Run Seek DB tests (test that the database gives expected bio-informative results)
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
        - <code>bash run_paramtest.sh -v -s <path_to_seek_db> -b <path_to_seek_binaries> </code>
        - <code>bash run_querysize.sh -v -s <path_to_seek_db> -b <path_to_seek_binaries> </code>
        - <code>bash run_bioinform.sh -v -s <path_to_seek_db> -b <path_to_seek_binaries> </code>


