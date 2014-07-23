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
#include <algorithm>


template<class tType>
struct SCompareRank {                                                       
    const vector<tType>&    m_vecData;
    SCompareRank( const vector<tType>& vecData ) : m_vecData(vecData) { }
    bool operator()( size_t iOne, size_t iTwo ) const {
        return ( m_vecData[ iOne ] < m_vecData[ iTwo ] ); }
};


double WilcoxonRankSum( const vector<float> vecdValues, const vector<float> vecAnswers ) {
    std::vector<size_t> veciIndices;
    std::vector<float> vecdRanks;
    size_t              iIndex, iCount, iPos, iNeg, i, j; 
    double  dSum, d;

    veciIndices.resize( vecdValues.size( ) ); 
    for( i = 0; i < vecdValues.size( ); ++i )
        veciIndices[ i ] = i; 
    // Sort ranks by values
    std::sort( veciIndices.begin( ), veciIndices.end( ), SCompareRank<float>( vecdValues ) ); 
    vecdRanks.resize( veciIndices.size( ) ); 
    for( i = 0; i < vecdRanks.size( ); ++i ) {
        iIndex = veciIndices[ i ]; 
        // Handle ties
        if( !i || ( vecdValues[ iIndex ] != vecdValues[ veciIndices[ i - 1 ] ] ) ) {
            for( iCount = 0,j = i; j < veciIndices.size( ); ++j ) {
                if( vecdValues[ veciIndices[ j ] ] != vecdValues[ iIndex ] )
                    break;
                iCount++; }
            d = i + ( iCount - 1 ) / 2.0f; }
        vecdRanks[ iIndex ] = d; }

    for ( i = 0, dSum = 0, iPos = 0, iNeg = 0; i < vecAnswers.size(); i++ ) {
        if ( vecAnswers[ i ] > 0 ) {
            iPos += 1;
            dSum += vecdRanks[ i ];
        }
        else {
            iNeg += 1;
        }
    }

    dSum -= ( iPos * ( iPos - 1 ) ) / 2;  
    return ( dSum / iPos / iNeg ); 
}


int main( int iArgs, char** aszArgs ) {

	gengetopt_args_info	sArgs;
    int					iRet;
    size_t				i, j, k, l, iGene, iExp, iAnswer;
    float d, fScore, fSum, fAnswer;
    DIR* dp;
    struct dirent* ep;
    CDat					DatSum, DatCur, DatRef, Data, Answers;
    
    vector<size_t>				veciGenesCur, veciAnnotIdx, veciGenesetIdx;
    vector<int>					veciAnswers;	
    vector<float>				vecdValues, vecfAnswers;

    CPCL PCL, AnnotPCL, GenesetsPCL;
    std::vector<string> input_files;
    std::string dab_dir;
    CFullMatrix<float> GeneDeg = CFullMatrix<float>();
    vector<float> RefGeneDeg = vector<float>();
    CFullMatrix<float> MatSum = CFullMatrix<float>();
    vector<string> features;

    if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
        cmdline_parser_print_help( );
        return 1; }
        CMeta Meta( sArgs.verbosity_arg );

        // check if directory valid
        if(sArgs.directory_arg){
            dp = opendir (sArgs.directory_arg);
            if (dp != NULL){
                (void) closedir (dp);	    
                dab_dir = sArgs.directory_arg;
            }
        else{
            cerr << "Couldn't open the directory: " << sArgs.directory_arg << '\n';
            return 1;
        }
    }

    dp = opendir (dab_dir.c_str());
    if (dp != NULL){
        while (ep = readdir (dp)){
      // skip . .. files and temp files with ~
        if (ep->d_name[0] == '.' || ep->d_name[strlen(ep->d_name)-1] == '~') 
            continue;
        if (std::string(ep->d_name).substr(strlen(ep->d_name)-4,4).compare(string(".dab")) == 0)
            input_files.push_back((string)sArgs.directory_arg + "/" + ep->d_name);          
        }
        (void) closedir (dp);           
    }

    if ( sArgs.genesets_arg ) {
        GenesetsPCL.Open( sArgs.genesets_arg, 0 );
    }

    if ( sArgs.refnet_arg ) {
        DatRef.Open( sArgs.refnet_arg );
        RefGeneDeg.resize( DatRef.GetGenes() );
        // Calculate gene degrees
        for ( j = 0; j < DatRef.GetGenes( ); ++ j ) {
            fSum = 0;
            for ( k = 0; k < DatRef.GetGenes( ); ++ k ) {
                if( j == k ) continue;
                fSum += DatRef.Get( j, k );
            }
            RefGeneDeg[j] = fSum / ( DatRef.GetGenes()-1 );
        }
    }

    if ( sArgs.annot_arg ) {
        AnnotPCL.Open( sArgs.annot_arg, 0 );
    }
    else if ( sArgs.enorm_flag || sArgs.gnorm_flag || sArgs.nnorm_flag ) {
        for( i = 0; i < input_files.size( ); ++i ) {
            // open dat/dab network
            if( !DatCur.Open( input_files[ i ].c_str() ) ) {
                cerr << "Couldn't open: " << input_files[ i ] << endl;
                return 1; 
            }	    

            cout << "calculating degree - opened: " << input_files[i] << endl;
            if ( i == 0 ) {
                GeneDeg.Initialize( DatCur.GetGenes(), input_files.size() );
            }

            // Calculate gene degrees
            for ( j = 0; j < DatCur.GetGenes( ); ++ j ) {
                fSum = 0;
                for ( k = 0; k < DatCur.GetGenes( ); ++ k ) {
                    if( j == k ) continue;
                    fSum += DatCur.Get( j, k );
                }
                GeneDeg.Set( j, i, fSum / (DatCur.GetGenes()-1) );
            }
        }
    }


    if ( ! sArgs.annot_given ) {
    // now iterate dat/dab networks
    for( i = 0; i < input_files.size( ); ++i ) {
        // open dat/dab network
        if( !DatCur.Open( input_files[ i ].c_str() ) ) {
            cerr << "Couldn't open: " << input_files[ i ] << endl;
            return 1; 
        }	    
        cerr << "opened: " << input_files[ i ] << endl;

        // if open first network, we will just add edge weights to this CDat
        if( i == 0 ){
            DatSum.Open( DatCur.GetGeneNames() );  
            MatSum.Initialize( DatCur.GetGenes(), DatCur.GetGenes() );
        }

        // now add edges to Dat
        for( j = 0; j < DatSum.GetGenes( ); ++j ) {
            for( k = ( j + 1 ); k < DatSum.GetGenes( ); ++k ) {
                if( CMeta::IsNaN( d = DatCur.Get(  j, k ) ) ){
                    cerr << "ERROR: missing values" << input_files[ i ] << endl;
                    return 1;
                }

                if ( sArgs.gnorm_flag || sArgs.nnorm_flag ) {
                    MatSum.Set( j, k, ( i == 0 ? 0 : MatSum.Get( j, k ) ) + d / GeneDeg.Get( j, i ) ); 
                    MatSum.Set( k, j, ( i == 0 ? 0 : MatSum.Get( k, j ) ) + d / GeneDeg.Get( k, i ) ); 
                    continue;
                }
                else if ( sArgs.enorm_flag ) { 
                    d /= sqrt ( GeneDeg.Get( j, i ) * GeneDeg.Get( k, i ) );
                }
 
                DatSum.Set( j, k, ( i == 0 ? 0 : DatSum.Get( j, k ) ) + d) ; 
            }
	    }
    }
    }

    for( i = 0; i < input_files.size(); ++i ) {
        if( !DatCur.Open( input_files[ i ].c_str() ) ) {
            cerr << "Couldn't open: " << input_files[ i ] << endl;
            return 1; 
        }	    
        cerr <<  input_files[ i ] << endl;
        if ( i == 0 ) {
            if ( sArgs.genesets_arg ) {
                if ( sArgs.gene_flag ) 
                    PCL.Open( DatCur.GetGeneNames(), input_files, features );
                else
                    PCL.Open( GenesetsPCL.GetExperimentNames(), input_files, features );
            }
            else
                PCL.Open( DatCur.GetGeneNames(), input_files, features );
        }

        string netName = CMeta::Basename( input_files[ i ].c_str() );
        netName = CMeta::Deextension( netName );

        cerr << netName << endl;

        if ( sArgs.annot_arg && !i ) {
            veciAnnotIdx.resize( DatCur.GetGenes() );
            for ( j = 0; j < DatCur.GetGenes(); j++ ) {
            size_t idx = AnnotPCL.GetGene( DatCur.GetGene( j ) );
            veciAnnotIdx[ j ] = idx;
            }
        }
        if ( sArgs.genesets_arg && !i ) {
            veciGenesetIdx.resize( DatCur.GetGenes() );
            for ( j = 0; j < DatCur.GetGenes(); j++ ) {
            size_t idx = GenesetsPCL.GetGene( DatCur.GetGene( j ) );
            veciGenesetIdx[ j ] = idx;
            }
        }

        size_t iGenes, iGenesMax;
        iGenesMax = 1;
        if ( sArgs.genesets_arg ) 
            iGenesMax = GenesetsPCL.GetExperiments();
        for( iGenes = sArgs.geneset_idx_arg; iGenes < std::max((size_t)sArgs.geneset_idx_arg, iGenesMax); iGenes++ ) {
            float fSetScore = 0;
            size_t iSet = 0;

            for( j = 0; j < DatCur.GetGenes(); ++j ) {
                fScore = 0;
                float fGeneScore = 0;
                size_t iGeneScores = 0; 

                if ( sArgs.genesets_arg && ( veciGenesetIdx[ j ] == -1 || GenesetsPCL.Get( veciGenesetIdx[ j ], iGenes ) != 1 ) ) continue;

                if ( sArgs.annot_arg ) {

                    iExp = AnnotPCL.GetExperiment( netName );
                    string strName = DatCur.GetGene ( j );

                    vecdValues.clear();
                    vecfAnswers.clear();

                    for( k = 0; k < DatCur.GetGenes( ); ++k ) {
                        if ( j == k ) continue;
                        d = DatCur.Get( j, k );

                        if ( ( iGene = veciAnnotIdx[ k ] ) == -1 )
                            continue;

                        if ( ( fAnswer = AnnotPCL.Get( iGene, iExp ) ) == 0 )
                            continue;

                        vecfAnswers.push_back( fAnswer );
                        vecdValues.push_back( d );
                    }

                    fScore = WilcoxonRankSum( vecdValues, vecfAnswers );
                    PCL.Set( j, i, fScore);

                    continue;
                }

                for( k = 0; k < DatCur.GetGenes( ); ++k ) {
                    if ( j == k ) continue;
                    if ( sArgs.genesets_arg && k <= i && !sArgs.gene_flag ) continue;
                    if ( sArgs.genesets_arg && ( veciGenesetIdx[ k ] == -1 || GenesetsPCL.Get( veciGenesetIdx[ k ], iGenes ) != 1 ) ) continue;

                    d = DatCur.Get( j, k );

                    // Correct for background in network (i) 
                    if ( sArgs.gnorm_flag ) {
                        d /= GeneDeg.Get( j, i );
                    }
                    else if ( sArgs.nnorm_flag ) {
                        d /= GeneDeg.Get( k, i );
                    }
                    else if ( sArgs.enorm_flag ) {
                        d /= sqrt ( GeneDeg.Get( j, i ) * GeneDeg.Get( k, i ) );
                    }

                    if ( sArgs.backg_flag ) {
                        float fDiv;
                        // Divide by average edge weight across all networks
                        if ( sArgs.refnet_arg ) {
                            float refd;
                            refd = DatRef.Get( j, k );
                            if ( sArgs.gnorm_flag ) 
                                fDiv = RefGeneDeg[ j ]; 
                            else if ( sArgs.enorm_flag )
                                fDiv = sqrt( RefGeneDeg[ j ] * RefGeneDeg[ k ] ); 
                            else
                                fDiv = RefGeneDeg[ k ]; 

                            refd /= fDiv;

                            d -= fDiv;
                            if ( d < 0 ) d = 0;
                        }
                        else if ( sArgs.gnorm_flag ) {
                            fDiv = ( MatSum.Get( j, k ) / input_files.size() );
                            if ( d ) d /= fDiv;
                        }
                        else if ( sArgs.nnorm_flag ) {
                            fDiv = ( MatSum.Get( k, j ) / input_files.size() );
                            if ( d ) d /= fDiv;
                        }
                        else {
                            fDiv = ( DatSum.Get( j, k ) / input_files.size() );
                            if ( d ) d /= fDiv;
                        }
                    }

                    if ( sArgs.genesets_arg ) {
                        fSetScore += d;
                        iSet += 1;
                        fGeneScore += d;
                        iGeneScores += 1;
                    }
                    // Use connectivity score to weight edge, or use score itself
                    else if ( sArgs.log_weight_flag ) {
                        fScore += DatCur.Get( j, k ) * log2( d );
                    }
                    else if ( sArgs.weight_flag ) {
                        fScore += DatCur.Get( j, k ) * d;
                    }
                    else {
                        fScore += d;
                    }
                }            


                if ( !sArgs.genesets_arg ) {
                    fScore /= ( DatCur.GetGenes() - 1 );
                    PCL.Set( j, i, fScore);
                }
                else if( sArgs.genesets_arg && sArgs.gene_flag ) {
                    fGeneScore /= iGeneScores;
                    PCL.Set( j, i, fGeneScore );
                }
            }      

            if ( sArgs.genesets_arg && !sArgs.gene_flag )
                PCL.Set( iGenes, i, fSetScore / (float)iSet );
	    }
    }

    PCL.Save( sArgs.output_arg );

    return iRet; 
}
