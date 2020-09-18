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
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>


float logit(float x) {
    if (x == 0.0)//min float number
        return -15;
    //return -std::numeric_limits<float>::max();
    //return std::numeric_limits<float>::min();
    if (x == 1.0)//max float number
        return 15;
    //return std::numeric_limits<float>::max();

    return log(x / (1 - x));
}

float sigmoid(float x) {
    float exp_value;
    float return_value;

    /*** Exponential calculation ***/
    exp_value = exp((double) -x);

    /*** Final sigmoid value ***/
    return_value = 1 / (1 + exp_value);

    return return_value;
}

float equal_prior(float x, size_t p, size_t n) {
    float logit_value;
    float return_value;

    logit_value = logit(x) + log(((float) p) / n); // - log(n/p) ... subtract prior ... in effect give equal prior to
    return_value = sigmoid(logit_value);

    return return_value;
}

inline bool exists(const std::string &name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}


int main(int iArgs, char **aszArgs) {

    //	cout << sigmoid(logit(0.999999)) << endl;
    //	cout << sigmoid(logit(0.9)) << endl;
    //	cout << sigmoid(logit(0.1)) << endl;
    //	cout << sigmoid(logit(1)) << endl;
    //	cout << sigmoid(logit(0)) << endl;


    gengetopt_args_info sArgs{};
    int iRet = 0;
    size_t i, j, k, l;
    float d;
    DIR *dp;
    struct dirent *ep;
    CDat DatOut, DatCur;
    vector <size_t> veciGenesCur;

    // context weight filename
    std::string weight_filename;
    std::vector <string> input_files;
    std::string dab_dir;

    // store prior information
    vector <size_t> vpos;
    vector <size_t> vneg;

    // read context weights
    std::string line;
    std::vector <std::string> terms;
    std::map<std::string, double> term2weight;
    double weight_sum = 0.0;


    if (cmdline_parser(iArgs, aszArgs, &sArgs)) {
        cmdline_parser_print_help();
        return 1;
    }
    CMeta Meta(sArgs.verbosity_arg);

    //int type = sArgs.type_arg ; // type of combiner

    // check if directory valid
    if (sArgs.directory_arg) {
        dp = opendir(sArgs.directory_arg);
        if (dp != NULL) {
            (void) closedir(dp);
            dab_dir = sArgs.directory_arg;
        } else {
            cerr << "Couldn't open the directory: " << sArgs.directory_arg << '\n';
            return 1;
        }
    }

    if (sArgs.weights_arg) { // only if weight file given..
        weight_filename = sArgs.weights_arg;
        std::ifstream weight_file(weight_filename.c_str());
        cout << "Using weight file: " << weight_filename << endl;
        // read context weights
        if (weight_file.is_open()) {
            while (weight_file.good()) {
                std::string term;
                double weight;

                weight_file >> term >> weight;
                //weight = abs(weight);
                //cout << term << "\t" << weight << endl;
                if (term.length() < 1) {
                    continue;
                }
                if (weight < 0.0) {//ignore weights less than zero
                    continue;
                }
                cout << term << endl;
                terms.push_back(term);
                term2weight[term] = weight;
                weight_sum += weight;
            }
            weight_file.close();
        }

        cout << "The number of contexts with weight: " << terms.size() << endl;

        // read only networks/contexts with weights
        std::vector <string> subset_terms;
        std::string dabname;
        for (j = 0; j < terms.size(); ++j) {
            dabname = dab_dir + "/" + terms[j] + ".dab";
            if (exists(dabname)) {
                cout << dabname;
                input_files.push_back(dabname);
                subset_terms.push_back(terms[j]);
            }
            //cout << pos << endl;
        }
        terms = subset_terms;
    } else {

        dp = opendir(dab_dir.c_str());
        if (dp != nullptr) {
            while ((ep = readdir(dp))) {
                // skip . .. files and temp files with ~
                if (ep->d_name[0] == '.' || ep->d_name[strlen(ep->d_name) - 1] == '~')
                    continue;
                if (std::string(ep->d_name).substr(strlen(ep->d_name) - 4, 4) == string(".dab"))
                    input_files.push_back((string) sArgs.directory_arg + "/" + ep->d_name);
            }
            (void) closedir(dp);
        }

        //TODO: temporary hacking based on implementation ... should improve
        // sort so that input file and term array are corresponding with index

        std::sort(input_files.begin(), input_files.end());
        for (i = 0; i < input_files.size(); ++i) {
            terms.push_back(input_files[i]);
            term2weight[input_files[i]] = 1;
        }

        //std::sort(terms.begin(),terms.end());
        weight_sum = input_files.size();
    }


    // now iterate dat/dab networks
    for (i = 0; i < input_files.size(); ++i) {

        // open dat/dab network
        if (!DatCur.Open(input_files[i].c_str())) {
            cerr << "Couldn't open: " << input_files[i] << endl;
            return 1;
        }
        cerr << "opened: " << input_files[i] << endl;


        if (sArgs.prior_arg) {
            for (j = 0; j < DatCur.GetGenes(); ++j)
                for (k = (j + 1); k < DatCur.GetGenes(); ++k)
                    DatCur.Set(j, k, equal_prior(DatCur.Get(j, k), vpos[i], vneg[i]));

        }

        if (sArgs.logit_flag) {
            for (j = 0; j < DatCur.GetGenes(); ++j)
                for (k = (j + 1); k < DatCur.GetGenes(); ++k)
                    DatCur.Set(j, k, logit(DatCur.Get(j, k)));

        } else {
        }

        if (sArgs.znormalize_flag) {
            DatCur.Normalize(CDat::ENormalizeZScore);
        } else {
        }

        cerr << term2weight[terms[i]] << endl;


        // if open first network, we will just add edge weights to this CDat
        if (i == 0) {
            DatOut.Open(DatCur);
            for (j = 0; j < DatOut.GetGenes(); ++j)
                //cerr << "set " << j << endl;

                for (k = (j + 1); k < DatOut.GetGenes(); ++k) {
                    //                DatOut.Set( j, k, DatOut.Get( j, k ) );
                    DatOut.Set(j, k, DatOut.Get(j, k) * term2weight[terms[i]]);
                }
            continue;
        }

        //cerr << "map flag" << endl;	  
        if (sArgs.map_flag) {
            // Get gene index match
            veciGenesCur.clear();
            veciGenesCur.resize(DatOut.GetGenes());
            for (l = 0; l < DatOut.GetGenes(); l++) {
                veciGenesCur[l] = DatCur.GetGene(DatOut.GetGene(l));
                if (veciGenesCur[l] == -1) {
                    cerr << "ERROR: missing gene" << input_files[l] << endl;
                    return 1;
                }
            }
        }

        cerr << "add a network" << endl;

        // now add edges to Dat
        for (j = 0; j < DatOut.GetGenes(); ++j)
            for (k = (j + 1); k < DatOut.GetGenes(); ++k) {

                if (sArgs.map_flag) {
                    // we are assuming a fully connected network
                    if (CMeta::IsNaN(d = DatCur.Get(veciGenesCur[j], veciGenesCur[k]))) {
                        cerr << d << endl;
                        cerr << veciGenesCur[j] << endl;
                        cerr << veciGenesCur[k] << endl;
                        cerr << "ERROR: missing values" << input_files[i] << endl;
                        return 1;
                    }
                } else {
                    if (CMeta::IsNaN(d = DatCur.Get(j, k))) {
                        cerr << "ERROR: missing values" << input_files[i] << endl;
                        return 1;
                    }
                }

                //cerr << d << endl;
                //cerr << terms.size() << endl;
                //for(int q = 0 ; q < terms.size() ; q++)
                //  cerr << terms[q] << endl;
                //cerr << terms[i] << endl;
                //cerr << term2weight[terms[i]] << endl;

                DatOut.Set(j, k, DatOut.Get(j, k) + d *
                                                    term2weight[terms[i]]); //weight edge based on context weight TODO: right now this implementation depends on the order of terms vector and niput_file vector to be correspondingi
                //                DatOut.Set( j, k, DatOut.Get( j, k ) +  d ) ;//


            }
    }

    // now convert sum to mean
    for (j = 0; j < DatOut.GetGenes(); ++j)
        for (k = (j + 1); k < DatOut.GetGenes(); ++k) {
            DatOut.Set(j, k, DatOut.Get(j, k) / weight_sum);

            if (sArgs.logit_flag) {// transform back to probability values
                DatOut.Set(j, k, sigmoid(DatOut.Get(j, k)));
            }

        }
    //DatOut.Set( j, k, DatOut.Get( j, k ) / input_files.size( ) );

    DatOut.Save(sArgs.output_arg);
    return iRet;
}
