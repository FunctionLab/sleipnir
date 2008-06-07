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

int open_genes( const char* szFile, CGenes& Genes ) {
	ifstream	ifsm;

	if( szFile ) {
		ifsm.open( szFile );
		if( !Genes.Open( ifsm ) ) {
			cerr << "Could not open: " << szFile << endl;
			return 1; } }

	return 0; }

int main( int iArgs, char** aszArgs ) {
	static const size_t	c_iBuf	= 1024;
	char				szBuf[ c_iBuf ];
	gengetopt_args_info	sArgs;
	ifstream			ifsm;
	CDat				Dat, DatNew;
	CDat*				pDat;
	float				d, dCutoff;
	CGenome				Genome;
	CGenes				GenesIn( Genome ), GenesQr( Genome );
	int					iRet;
	size_t				i, j;
	vector<float>		vecdColors, vecdBorders;

	if( cmdline_parser2( iArgs, aszArgs, &sArgs, 0, 1, 0 ) && ( sArgs.config_arg &&
		cmdline_parser_configfile( sArgs.config_arg, &sArgs, 0, 0, 1 ) ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta = CMeta( sArgs.verbosity_arg );

	if( sArgs.features_arg ) {
		ifsm.open( sArgs.features_arg );
		if( !Genome.Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.features_arg << endl;
			return 1; }
		ifsm.close( ); }

	if( sArgs.colors_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.colors_arg );
		if( !ifsm.is_open( ) ) {
			cerr << "Could not open: " << sArgs.colors_arg << endl;
			return 1; }
		while( ifsm.peek( ) != EOF ) {
			ifsm.getline( szBuf, c_iBuf - 1 );
			vecdColors.push_back( (float)atof( szBuf ) ); }
		ifsm.close( ); }

	if( sArgs.borders_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.borders_arg );
		if( !ifsm.is_open( ) ) {
			cerr << "Could not open: " << sArgs.borders_arg << endl;
			return 1; }
		while( ifsm.peek( ) != EOF ) {
			ifsm.getline( szBuf, c_iBuf - 1 );
			vecdBorders.push_back( (float)atof( szBuf ) ); }
		ifsm.close( ); }

	if( iRet = open_genes( sArgs.genes_arg, GenesIn ) )
		return iRet;
	if( iRet = open_genes( sArgs.geneq_arg, GenesQr ) )
		return iRet;

	if( sArgs.input_arg ) {
		if( !Dat.Open( sArgs.input_arg, sArgs.memmap_flag && !( sArgs.normalize_flag ||
			sArgs.genes_arg || sArgs.geneq_arg || sArgs.knowns_arg ) ) ) {
			cerr << "Couldn't open: " << sArgs.input_arg << endl;
			return 1; } }
	else if( !Dat.Open( cin, CDat::EFormatText ) ) {
		cerr << "Couldn't open input" << endl;
		return 1; }
	pDat = &Dat;

	dCutoff = (float)( sArgs.cutoff_given ? sArgs.cutoff_arg : HUGE_VAL );
	if( GenesIn.GetGenes( ) ) {
		vector<size_t>	veciGenes;

		DatNew.Open( GenesIn.GetGeneNames( ) );
		veciGenes.resize( DatNew.GetGenes( ) );
		for( i = 0; i < veciGenes.size( ); ++i )
			veciGenes[ i ] = Dat.GetGene( DatNew.GetGene( i ) );
		for( i = 0; i < veciGenes.size( ); ++i ) {
			if( veciGenes[ i ] == -1 )
				continue;
			for( j = ( i + 1 ); j < veciGenes.size( ); ++j )
				if( veciGenes[ j ] != -1 )
					DatNew.Set( i, j, Dat.Get( veciGenes[ i ], veciGenes[ j ] ) ); }
		pDat = &DatNew; }
	if( sArgs.normalize_flag )
		pDat->Normalize( );
	if( GenesQr.GetGenes( ) ) {
		if( sArgs.cutoff_given )
			for( i = 0; i < pDat->GetGenes( ); ++i )
				for( j = ( i + 1 ); j < pDat->GetGenes( ); ++j )
					if( !CMeta::IsNaN( d = pDat->Get( i, j ) ) && ( d < sArgs.cutoff_arg ) )
						pDat->Set( i, j, CMeta::GetNaN( ) );
		if( !strcmp( sArgs.format_arg, "correl" ) ) {
			CMeasurePearson	MeasurePearson;
			float*			adCentroid;
			float*			adCur;
			size_t			iCur;
			vector<size_t>	veciCounts;
			vector<float>	vecdScores;

			veciCounts.resize( pDat->GetGenes( ) );
			adCentroid = new float[ pDat->GetGenes( ) ];
			for( i = 0; i < GenesQr.GetGenes( ); ++i ) {
				if( ( iCur = pDat->GetGene( GenesQr.GetGene( i ).GetName( ) ) ) == -1 )
					continue;
				for( j = 0; j < pDat->GetGenes( ); ++j )
					if( !CMeta::IsNaN( d = pDat->Get( iCur, j ) ) ) {
						adCentroid[ j ] += d;
						veciCounts[ j ]++; } }
			for( i = 0; i < pDat->GetGenes( ); ++i )
				adCentroid[ i ] /= veciCounts[ i ];

			vecdScores.resize( pDat->GetGenes( ) );
			adCur = new float[ pDat->GetGenes( ) ];
			for( i = 0; i < pDat->GetGenes( ); ++i ) {
				for( j = 0; j < pDat->GetGenes( ); ++j )
					adCur[ j ] = pDat->Get( i, j );
				vecdScores[ i ] = (float)MeasurePearson.Measure( adCentroid, pDat->GetGenes( ), adCur,
					pDat->GetGenes( ), IMeasure::EMapNone, NULL, NULL ); }
			delete[] adCur;
			delete[] adCentroid;
			for( i = 0; i < vecdScores.size( ); ++i )
				cout << pDat->GetGene( i ) << '\t' << vecdScores[ i ] << endl; }
		else {
			dCutoff = 0;
			if( vecdColors.empty( ) ) {
				vecdColors.resize( pDat->GetGenes( ) );
				fill( vecdColors.begin( ), vecdColors.end( ), 0.5f );
				for( i = 0; i < GenesQr.GetGenes( ); ++i )
					if( ( j = pDat->GetGene( GenesQr.GetGene( i ).GetName( ) ) ) != -1 )
						vecdColors[ j ] = 1; }
			pDat->FilterGenes( GenesQr, sArgs.hefalmp_flag ? CDat::EFilterHefalmp : CDat::EFilterPixie,
				sArgs.neighbors_arg, (float)sArgs.edges_arg ); } }
	if( sArgs.knowns_arg ) {
		CDat			DatKnowns;
		vector<size_t>	veciKnowns;
		size_t			iOne, iTwo;

		if( !DatKnowns.Open( sArgs.knowns_arg, !!sArgs.memmap_flag ) ) {
			cerr << "Could not open: " << sArgs.knowns_arg << endl;
			return 1; }
		veciKnowns.resize( pDat->GetGenes( ) );
		for( i = 0; i < veciKnowns.size( ); ++i )
			veciKnowns[ i ] = DatKnowns.GetGene( pDat->GetGene( i ) );
		for( i = 0; i < pDat->GetGenes( ); ++i )
			if( ( iOne = veciKnowns[ i ] ) != -1 )
				for( j = ( i + 1 ); j < pDat->GetGenes( ); ++j )
					if( ( ( iTwo = veciKnowns[ j ] ) != -1 ) &&
						!CMeta::IsNaN( d = DatKnowns.Get( iOne, iTwo ) ) && ( d > 0 ) )
						pDat->Set( i, j, CMeta::GetNaN( ) ); }

	if( !strcmp( sArgs.format_arg, "dot" ) )
		pDat->SaveDOT( cout, dCutoff, &Genome, false, true, vecdColors.empty( ) ? NULL : &vecdColors,
			vecdBorders.empty( ) ? NULL : &vecdBorders );
	else if( !strcmp( sArgs.format_arg, "gdf" ) )
		pDat->SaveGDF( cout, dCutoff );
	else if( !strcmp( sArgs.format_arg, "net" ) )
		pDat->SaveNET( cout, dCutoff );
	else if( !strcmp( sArgs.format_arg, "matisse" ) )
		pDat->SaveMATISSE( cout, dCutoff, &Genome );
	else if( !strcmp( sArgs.format_arg, "list" ) ) {
		vector<bool>					vecfQuery;
		map<size_t, float>				mapGenes;
		map<size_t, float>::iterator	iterGene;
		size_t							iGene;

		vecfQuery.resize( pDat->GetGenes( ) );
		for( i = 0; i < vecfQuery.size( ); ++i )
			vecfQuery[ i ] = GenesQr.IsGene( pDat->GetGene( i ) );
		for( i = 0; i < pDat->GetGenes( ); ++i )
			for( j = ( i + 1 ); j < pDat->GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = pDat->Get( i, j ) ) &&
					( CMeta::IsNaN( dCutoff ) || ( d > dCutoff ) ) &&
					( vecfQuery[ i ] != vecfQuery[ j ] ) ) {
					iGene = vecfQuery[ i ] ? j : i;
					if( ( iterGene = mapGenes.find( iGene ) ) == mapGenes.end( ) )
						mapGenes[ iGene ] = d;
					else
						iterGene->second += d; }
		for( iterGene = mapGenes.begin( ); iterGene != mapGenes.end( ); ++iterGene )
			cout << pDat->GetGene( iterGene->first ) << '\t' << iterGene->second << endl; }
	else if( !strcmp( sArgs.format_arg, "dat" ) )
		pDat->Save( cout, CDat::EFormatText );

	return 0; }
