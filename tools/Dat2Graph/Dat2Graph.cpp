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
	CDat				Dat;
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
	CMeta::Startup( sArgs.verbosity_arg );

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

	dCutoff = (float)( sArgs.cutoff_given ? sArgs.cutoff_arg : HUGE_VAL );
	if( GenesIn.GetGenes( ) )
		Dat.FilterGenes( GenesIn, CDat::EFilterInclude );
	if( GenesQr.GetGenes( ) ) {
		if( sArgs.cutoff_given )
			for( i = 0; i < Dat.GetGenes( ); ++i )
				for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
					if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) && ( d < sArgs.cutoff_arg ) )
						Dat.Set( i, j, CMeta::GetNaN( ) );
		if( !strcmp( sArgs.format_arg, "correl" ) ) {
			CMeasurePearson	MeasurePearson;
			float*			adCentroid;
			float*			adCur;
			size_t			iCur;
			vector<size_t>	veciCounts;
			vector<float>	vecdScores;

			veciCounts.resize( Dat.GetGenes( ) );
			adCentroid = new float[ Dat.GetGenes( ) ];
			for( i = 0; i < GenesQr.GetGenes( ); ++i ) {
				if( ( iCur = Dat.GetGene( GenesQr.GetGene( i ).GetName( ) ) ) == -1 )
					continue;
				for( j = 0; j < Dat.GetGenes( ); ++j )
					if( !CMeta::IsNaN( d = Dat.Get( iCur, j ) ) ) {
						adCentroid[ j ] += d;
						veciCounts[ j ]++; } }
			for( i = 0; i < Dat.GetGenes( ); ++i )
				adCentroid[ i ] /= veciCounts[ i ];

			vecdScores.resize( Dat.GetGenes( ) );
			adCur = new float[ Dat.GetGenes( ) ];
			for( i = 0; i < Dat.GetGenes( ); ++i ) {
				for( j = 0; j < Dat.GetGenes( ); ++j )
					adCur[ j ] = Dat.Get( i, j );
				vecdScores[ i ] = (float)MeasurePearson.Measure( adCentroid, Dat.GetGenes( ), adCur,
					Dat.GetGenes( ), IMeasure::EMapNone, NULL, NULL ); }
			delete[] adCur;
			delete[] adCentroid;
			for( i = 0; i < vecdScores.size( ); ++i )
				cout << Dat.GetGene( i ) << '\t' << vecdScores[ i ] << endl; }
		else {
			dCutoff = 0;
			if( vecdColors.empty( ) ) {
				vecdColors.resize( Dat.GetGenes( ) );
				fill( vecdColors.begin( ), vecdColors.end( ), 0.5f );
				for( i = 0; i < GenesQr.GetGenes( ); ++i )
					if( ( j = Dat.GetGene( GenesQr.GetGene( i ).GetName( ) ) ) != -1 )
						vecdColors[ j ] = 1; }
			Dat.FilterGenes( GenesQr, CDat::EFilterPixie, sArgs.neighbors_arg,
				!!strcmp( sArgs.format_arg, "list" ) ); } }
	if( sArgs.knowns_arg ) {
		CDat			DatKnowns;
		vector<size_t>	veciKnowns;
		size_t			iOne, iTwo;

		if( !DatKnowns.Open( sArgs.knowns_arg, !!sArgs.memmap_flag ) ) {
			cerr << "Could not open: " << sArgs.knowns_arg << endl;
			return 1; }
		veciKnowns.resize( Dat.GetGenes( ) );
		for( i = 0; i < veciKnowns.size( ); ++i )
			veciKnowns[ i ] = DatKnowns.GetGene( Dat.GetGene( i ) );
		for( i = 0; i < Dat.GetGenes( ); ++i )
			if( ( iOne = veciKnowns[ i ] ) != -1 )
				for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
					if( ( ( iTwo = veciKnowns[ j ] ) != -1 ) &&
						!CMeta::IsNaN( d = DatKnowns.Get( iOne, iTwo ) ) && ( d > 0 ) )
						Dat.Set( i, j, CMeta::GetNaN( ) ); }
	if( sArgs.normalize_flag )
		Dat.Normalize( );

	if( !strcmp( sArgs.format_arg, "dot" ) )
		Dat.SaveDOT( cout, dCutoff, &Genome, false, true, vecdColors.empty( ) ? NULL : &vecdColors,
			vecdBorders.empty( ) ? NULL : &vecdBorders );
	else if( !strcmp( sArgs.format_arg, "gdf" ) )
		Dat.SaveGDF( cout, dCutoff );
	else if( !strcmp( sArgs.format_arg, "net" ) )
		Dat.SaveNET( cout, dCutoff );
	else if( !strcmp( sArgs.format_arg, "matisse" ) )
		Dat.SaveMATISSE( cout, dCutoff, &Genome );
	else if( !strcmp( sArgs.format_arg, "list" ) ) {
		vector<bool>					vecfQuery;
		map<size_t, float>				mapGenes;
		map<size_t, float>::iterator	iterGene;
		size_t							iGene;

		vecfQuery.resize( Dat.GetGenes( ) );
		for( i = 0; i < vecfQuery.size( ); ++i )
			vecfQuery[ i ] = GenesQr.IsGene( Dat.GetGene( i ) );
		for( i = 0; i < Dat.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) &&
					( CMeta::IsNaN( dCutoff ) || ( d > dCutoff ) ) &&
					( vecfQuery[ i ] != vecfQuery[ j ] ) ) {
					iGene = vecfQuery[ i ] ? j : i;
					if( ( iterGene = mapGenes.find( iGene ) ) == mapGenes.end( ) )
						mapGenes[ iGene ] = d;
					else
						iterGene->second += d; }
		for( iterGene = mapGenes.begin( ); iterGene != mapGenes.end( ); ++iterGene )
			cout << Dat.GetGene( iterGene->first ) << '\t' << iterGene->second << endl; }
	else if( !strcmp( sArgs.format_arg, "dat" ) )
		Dat.Save( cout, CDat::EFormatText );

	CMeta::Shutdown( );
	return 0; }
