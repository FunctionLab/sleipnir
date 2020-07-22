/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library, the included executable tools, or any related
* code in your work, please cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/
#include "stdafx.h"
#include "cmdline.h"

bool oneGenePerFile = true;  // TODO - make this an input param

int main(int iArgs, char **aszArgs) {
    static const size_t c_iBuffer = 1024;
#ifdef WIN32
    pthread_win32_process_attach_np( );
#endif // WIN32
    gengetopt_args_info sArgs;
    ifstream ifsm;
    istream *pistm;
    vector <string> vecstrLine, vecstrGenes, vecstrDatasets;
    char acBuffer[c_iBuffer];
#ifndef NO_SMILE
    CBayesNetSmile BNSmile;
#endif
    size_t i;
    map <string, size_t> mapstriZeros;

    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }
    CMeta Meta(sArgs.verbosity_arg);

    // Open and read the gene_map file and create a vector of gene names
    // (i.e. numeric Entrez names)
    if (sArgs.input_arg) {
        ifsm.open(sArgs.input_arg);
        pistm = &ifsm;
    } else
        pistm = &cin;
    while (!pistm->eof()) {
        pistm->getline(acBuffer, c_iBuffer - 1);
        acBuffer[c_iBuffer - 1] = 0;
        vecstrLine.clear();
        CMeta::Tokenize(acBuffer, vecstrLine);
        if (vecstrLine.size() < 2) {
            cerr << "Ignoring line: " << acBuffer << endl;
            continue;
        }
        if (!(i = atoi(vecstrLine[0].c_str()))) {
            cerr << "Illegal gene ID: " << vecstrLine[0] << " for " << vecstrLine[1] << endl;
            return 1;
        }
        i--;
        if (vecstrGenes.size() <= i)
            vecstrGenes.resize(i + 1);
        vecstrGenes[i] = vecstrLine[1];
    }
    if (sArgs.input_arg)
        ifsm.close();

    // Read and process the zeros file
    // TODO - What does it do?
    if (sArgs.zeros_arg) {
        ifstream ifsm_zero;
        vector <string> vecstrLine;
        char acLine[1024];

        ifsm_zero.open(sArgs.zeros_arg);
        if (!ifsm_zero.is_open()) {
            cerr << "Couldn't open: " << sArgs.zeros_arg << endl;
            return 1;
        }
        while (!ifsm_zero.eof()) {
            ifsm_zero.getline(acLine, ARRAYSIZE(acLine) - 1);
            acLine[ARRAYSIZE(acLine) - 1] = 0;
            vecstrLine.clear();
            CMeta::Tokenize(acLine, vecstrLine);
            if (vecstrLine.empty())
                continue;
            mapstriZeros[vecstrLine[0]] = atoi(vecstrLine[1].c_str());
        }
    }

    size_t numFiles = min((size_t) sArgs.files_arg, vecstrGenes.size());
    if (sArgs.one_gene_per_file_flag == 1) {
        numFiles = vecstrGenes.size();
    }

    bool useNibble = false;
    if (sArgs.use_nibble_flag == 1) {
        useNibble = true;
    }

    CDatabase DB(useNibble);
    DB.SetMemmap(!!sArgs.memmap_flag);
    DB.SetBuffer(!!sArgs.buffer_flag);
    DB.SetBlockOut(sArgs.block_files_arg);
    DB.SetBlockIn(sArgs.block_datasets_arg);

    if (sArgs.network_arg) {
        if (sArgs.dataset_arg) {
            cerr << "Confused. Only network OR dataset list." << endl;
            return 1;
        }
#ifdef NO_SMILE
        cerr << "network option is disabled because SMILE was not used in building Data2DB."<< endl;
        return 1;
#endif

#ifndef NO_SMILE
        if (!BNSmile.Open(sArgs.network_arg)) {
            cerr << "Could not open: " << sArgs.network_arg << endl;
            return 1;
        }
        if (!DB.Open(vecstrGenes, sArgs.dir_in_arg, &BNSmile, sArgs.dir_out_arg,
                     numFiles, mapstriZeros)) {
            cerr << "Could not open data" << endl;
            return 1;
        }
#endif
    } else if (sArgs.dataset_arg) {
        // Read dataset_list file and create a vector of the datasets
        // Note: the list entries are typically of form dataset.platform and
        //  are the prefix of the file names for the dab fles.
        ifsm.open(sArgs.dataset_arg);
        while (!pistm->eof()) {
            pistm->getline(acBuffer, c_iBuffer - 1);
            if (acBuffer[0] == 0)
                break;
            acBuffer[c_iBuffer - 1] = 0;
            //If line contains multiple columns,
            //use the first column, which is the dataset column
            vector <string> tok;
            CMeta::Tokenize(acBuffer, tok, " \t");
            vecstrDatasets.push_back(tok[0]);
        }
        vecstrDatasets.resize(vecstrDatasets.size());
        ifsm.close();

        // Create the CDatabase, note the last arg is the min of num_files or num_genes
        // Rename to CreateDBFiles
        if (!DB.Open(vecstrGenes, vecstrDatasets, sArgs.dir_in_arg,
                              sArgs.dir_out_arg, numFiles, mapstriZeros)) {
            cerr << "Could not open data" << endl;
            return 1;
        }

    } else {
        cerr << "Must give a network or a dataset list." << endl;
        return 1;

    }

#ifdef WIN32
    pthread_win32_process_detach_np( );
#endif // WIN32
    return 0;
}
