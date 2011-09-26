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
#include <vector>

static const char	c_acDab[]	= ".dab";

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	ifstream			ifsm;
	CBayesNetSmile		BNSmile;
	CDat				Dat, DatLookup;
	istream*			pistm;
	vector<size_t>		veciGenes;
	CPCL				PCLLookup;
	vector<string>		vecstrNodes, vecstrGenes, vecstrDummy;
	size_t				i, j, k, iOne, iTwo, iGene, nStart, nEnd;
	float				dPrior;
	CDataMatrix			MatCPT;
	unsigned char		b;
	
	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg );

	if( !BNSmile.Open( sArgs.network_arg ) ) {
		cerr << "Could not open: " << sArgs.network_arg << endl;
		return 1; }
	if( !Dat.Open( sArgs.input_arg, !!sArgs.memmap_flag ) ) {
		cerr << "Could not open: " << sArgs.input_arg << endl;
		return 1; }
	if( sArgs.lookup_arg ) {
		ifsm.open( sArgs.lookup_arg );
		pistm = &ifsm; }
	else
		pistm = &cin;
	
	if( !DatLookup.Open( *pistm, CDat::EFormatText, 1 ) ) {
		cerr << "Could not open: " << ( sArgs.lookup_arg ? sArgs.lookup_arg : "stdin" ) << endl;
		return 1; }
	ifsm.close( );
	
	for( i = 0; i < DatLookup.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < DatLookup.GetGenes( ); ++j )
			if( DatLookup.Get( i, j ) == 1 )
				vecstrGenes.push_back( DatLookup.GetGene( i ) + " - " + DatLookup.GetGene( j ) );
	BNSmile.GetNodes( vecstrNodes );
	vecstrNodes[ 0 ] = sArgs.input_arg;
			
	///// Set for start idx and end idx	  
	nStart = 0;
	nEnd = 0;
	if( sArgs.start_arg > -1){
	  nStart = sArgs.start_arg + 1;
	  
	  if( nStart == 1 )
	    vecstrNodes.erase( vecstrNodes.begin() + 1 );
	  else if( nStart <= vecstrNodes.size() )
	    vecstrNodes.erase( vecstrNodes.begin() + 1, vecstrNodes.begin() + nStart );
	}	
	if( sArgs.end_arg > -1 ){
	  nEnd = (sArgs.end_arg+1) - nStart + 1;
	  
	  if( (nEnd+1) < vecstrNodes.size() )
	    vecstrNodes.erase( (vecstrNodes.begin() + nEnd + 1), vecstrNodes.end() );
	}
	
	PCLLookup.Open( vecstrGenes, vecstrNodes, vecstrDummy );
	
	veciGenes.resize( DatLookup.GetGenes( ) );
	for( i = 0; i < veciGenes.size( ); ++i )
		veciGenes[ i ] = Dat.GetGene( DatLookup.GetGene( i ) );
	for( iGene = i = 0; i < DatLookup.GetGenes( ); ++i ) {
		iOne = veciGenes[ i ];
		for( j = ( i + 1 ); j < DatLookup.GetGenes( ); ++j ) {
			iTwo = veciGenes[ j ];
			if( DatLookup.Get( i, j ) != 1 )
				continue;
			if( ( iOne != -1 ) && ( iTwo != -1 ) )
				PCLLookup.Set( iGene, 0, Dat.Get( iOne, iTwo ) );
			iGene++; } }

	BNSmile.GetCPT( 0, MatCPT );
	dPrior = 1 - MatCPT.Get( 0, 0 );
	for( i = 1; i < vecstrNodes.size( ); ++i ) {
		CDataPair	DatCur;
		string		strIn;
		size_t		iValue;
		
		cerr << "Dataset: " << vecstrNodes[ i ] << endl;
		strIn = (string)sArgs.directory_arg + "/" + vecstrNodes[ i ] + c_acDab;				
		if( !DatCur.Open( strIn.c_str( ), false, !!sArgs.memmap_flag ) ) {
			cerr << "Could not open: " << strIn << endl;
			return 1; }
		for( j = 0; j < veciGenes.size( ); ++j )
			veciGenes[ j ] = DatCur.GetGene( DatLookup.GetGene( j ) );
		for( iGene = j = 0; j < DatLookup.GetGenes( ); ++j ) {
			iOne = veciGenes[ j ];
			for( k = ( j + 1 ); k < DatLookup.GetGenes( ); ++k ) {
				if( DatLookup.Get( j, k ) != 1 )
					continue;
				iTwo = veciGenes[ k ];
				iValue = -1;
				if( ( iOne != -1 ) && ( iTwo != -1 ) )
					iValue = DatCur.Quantize( DatCur.Get( iOne, iTwo ) );
				if( ( iValue == -1 ) && ( ( b = BNSmile.GetDefault( i ) ) != (unsigned char)-1 ) )
					iValue = b;
				if( iValue != -1 )
					PCLLookup.Set( iGene, i, 1 - BNSmile.Evaluate( i, (unsigned char)iValue ) - dPrior );
				iGene++; } } }
	PCLLookup.Save( cout );

	return 0; }
