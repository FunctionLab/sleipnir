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

int read_genes( const char*, const char*, CGenome&, vector<string>&, vector<CGenes*>&, CDataPair& );
int write_posteriors( const string&, const CBayesNetSmile&, ostream& );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info			sArgs;
	CDatasetCompact				Data;
	size_t						i, j, iFunction;
	CGenome						Genome;
	ifstream					ifsm;
	CGenes						GenesIn( Genome ), GenesEx( Genome );
	vector<CGenes*>				vecpPositives;
	vector<float>				vecdQuants;
	vector<string>				vecstrPositives, vecstrIDs, vecstrNodes;
	int							iRet;
	ofstream					ofsmPosteriors;
	CBayesNetSmile				BNGlobal;
	IMeasure*					pMeasure;
	CMeasurePearson				Pearson;
	CMeasureEuclidean			Euclidean;
	CMeasureKendallsTau			KendallsTau;
	CMeasureKolmogorovSmirnov	KolmSmir;
	CMeasureSpearman			Spearman( true );
	CMeasureNegate				EuclideanNeg( &Euclidean, false );
	CMeasurePearNorm			PearNorm;
	IMeasure*					apMeasures[]	= { &Pearson, &EuclideanNeg, &KendallsTau,
		&KolmSmir, &Spearman, &PearNorm, NULL };

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) || !sArgs.inputs_num ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta = CMeta( sArgs.verbosity_arg, sArgs.random_arg );
#ifdef SMILEXML_LIB
	EnableXdslFormat( );
#endif
	pMeasure = NULL;
	for( i = 0; apMeasures[ i ]; ++i )
		if( !strcmp( apMeasures[ i ]->GetName( ), sArgs.distance_arg ) ) {
			pMeasure = apMeasures[ i ];
			break; }
	if( !pMeasure ) {
		cmdline_parser_print_help( );
		return 1; }

	if( sArgs.genes_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.genes_arg );
		if( !GenesIn.Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.genes_arg << endl;
			return 1; }
		ifsm.close( ); }
	if( sArgs.genex_arg ) {
		ifsm.clear( );
		ifsm.open( sArgs.genex_arg );
		if( !GenesEx.Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.genex_arg << endl;
			return 1; }
		ifsm.close( ); }

	if( sArgs.bins_arg ) {
		static const size_t	c_iBuffer	= 1024;
		char			acBuffer[ c_iBuffer ];
		vector<string>	vecstrQuants;

		ifsm.clear( );
		ifsm.open( sArgs.bins_arg );
		ifsm.getline( acBuffer, c_iBuffer - 1 );
		if( !( ifsm.is_open( ) && acBuffer[ 0 ] ) ) {
			cerr << "Could not open: " << sArgs.bins_arg << endl;
			return 1; }
		ifsm.close( );

		CMeta::Tokenize( acBuffer, vecstrQuants );
		vecdQuants.resize( vecstrQuants.size( ) );
		for( i = 0; i < vecdQuants.size( ); ++i )
			vecdQuants[ i ] = (float)atof( vecstrQuants[ i ].c_str( ) ); }
	else {
		static const float	c_adQuants[]	= {-1, 0, 1, 2, 3};

		vecdQuants.resize( ARRAYSIZE(c_adQuants) );
		copy( c_adQuants, c_adQuants + ARRAYSIZE(c_adQuants), vecdQuants.begin( ) ); }

	{
		static const float	c_adQuantsBinary[]	= {0.5, 1};
		CDataPair		Answers;
		vector<string>	vecstrPCLs, vecstrGenes;

		vecstrPCLs.resize( sArgs.inputs_num );
		copy( sArgs.inputs, sArgs.inputs + sArgs.inputs_num, vecstrPCLs.begin( ) );
		vecstrIDs.resize( sArgs.inputs_num );
		for( i = 0; i < sArgs.inputs_num; ++i )
			vecstrIDs[ i ] = CMeta::Basename( sArgs.inputs[ i ] );
		if( !BNGlobal.Open( vecstrIDs, vecdQuants.size( ) ) ) {
			cerr << "Failed to create a Bayesian network" << endl;
			return 1; }

		if( iRet = read_genes( sArgs.related_arg, sArgs.unrelated_arg, Genome, vecstrPositives, vecpPositives,
			Answers ) )
			return iRet;
		Answers.SetQuants( c_adQuantsBinary, ARRAYSIZE(c_adQuantsBinary) );

		if( !Data.Open( GenesIn, GenesEx, Answers, vecstrPCLs, sArgs.skip_arg, pMeasure, vecdQuants ) )
			return 1;
	}

	cerr << "Learning global network...";
	BNGlobal.Learn( &Data, 1, !!sArgs.zero_flag );
	if( !BNGlobal.Save( sArgs.global_arg ) ) {
		cerr << "Could not save: " << sArgs.global_arg << endl;
		return 1; }
	cerr << "done" << endl;

	ofsmPosteriors.open( sArgs.trusts_arg );
	if( !ofsmPosteriors.is_open( ) ) {
		cerr << "Could not open: " << sArgs.trusts_arg << endl;
		return 1; }
	BNGlobal.GetNodes( vecstrNodes );
	for( i = 1; i < vecstrNodes.size( ); ++i )
		ofsmPosteriors << '\t' << vecstrNodes[ i ];
	ofsmPosteriors << endl;

	for( iFunction = 0; iFunction < vecpPositives.size( ); ++iFunction ) {
		CBayesNetSmile	BNFunction;
		CDataMask		DataMask;
		string			strFile;
		CDat			DatOut;
		ofstream		ofsm;

		cerr << "Learning function: " << vecstrPositives[ iFunction ] << endl;

		DataMask.Attach( &Data );
		DataMask.FilterGenes( *vecpPositives[ iFunction ], CDat::EFilterTerm );
		BNFunction.Open( vecstrIDs, vecdQuants.size( ) );
		BNFunction.SetDefault( BNGlobal );
		BNFunction.Learn( &DataMask, 1, !!sArgs.zero_flag );
		strFile = (string)sArgs.output_arg + '/' + vecstrPositives[ iFunction ] + '.' +
			( sArgs.xdsl_flag ? "x" : "" ) + "dsl";
		if( !BNFunction.Save( strFile.c_str( ) ) ) {
			cerr << "Could not save: " << strFile << endl;
			return 1; }
		if( iRet = write_posteriors( vecstrPositives[ iFunction ], BNFunction, ofsmPosteriors ) )
			return iRet;

		strFile = (string)sArgs.predictions_arg + '/' + vecstrPositives[ iFunction ] + ( sArgs.dab_flag ?
			".dab" : ".dat" );
		DatOut.Open( DataMask.GetGeneNames( ) );
		if( !BNFunction.Evaluate( &DataMask, DatOut, !!sArgs.zero_flag ) ) {
			cerr << "Could not evaluate: " << strFile << endl;
			return 1; }
		DatOut.Invert( );
		if( sArgs.cutoff_arg > 0 )
			for( i = 0; i < DatOut.GetGenes( ); ++i )
				for( j = ( i + 1 ); j < DatOut.GetGenes( ); ++j )
					if( DatOut.Get( i, j ) < sArgs.cutoff_arg )
						DatOut.Set( i, j, CMeta::GetNaN( ) );
		DatOut.Save( strFile.c_str( ) ); }

	for( i = 0; i < vecpPositives.size( ); ++i )
		delete vecpPositives[ i ];

	return 0; }

int read_genes( const char* szPositives, const char* szNegatives, CGenome& Genome, vector<string>& vecstrPositives,
	vector<CGenes*>& vecpPositives, CDataPair& Answers ) {
	CGenes*			pGenes;
	string			strDir, strFile, strBase;
	CDat			DatNegatives;
	ifstream		ifsm;
	set<string>		setstrGenes;
	size_t			i;

	cerr << "Reading related gene pairs...";
	strDir = szPositives;
#ifdef _MSC_VER
	bool			fOK;
	HANDLE			hFind;
	WIN32_FIND_DATA	sFind;
	for( fOK = ( hFind = FindFirstFile( ( strDir + "\\*" ).c_str( ), &sFind ) ) !=
		INVALID_HANDLE_VALUE; fOK; fOK = !!FindNextFile( hFind, &sFind ) ) {
		if( sFind.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY )
			continue;
		strFile = strDir + '/' + ( strBase = sFind.cFileName );
#else // _MSC_VER
	DIR*			pDir;
	struct dirent*	pFind;
	struct stat		sStat;
	pDir = opendir( strDir.c_str( ) );
	while( pFind = readdir( pDir ) ) {
		strFile = strDir + '/' + ( strBase = pFind->d_name );
		stat( strFile.c_str( ), &sStat );
		if( S_ISDIR( sStat.st_mode ) )
			continue;
#endif // _MSC_VER

		ifsm.clear( );
		ifsm.open( strFile.c_str( ) );
		pGenes = new CGenes( Genome );
		if( !pGenes->Open( ifsm ) ) {
			cerr << "Could not open: " << strFile << endl;
			return 1; }
		ifsm.close( );
		vecstrPositives.push_back( CMeta::Deextension( strBase ) );
		vecpPositives.push_back( pGenes ); }
	cerr << "  done" << endl;

	cerr << "Reading unrelated gene pairs...";
	ifsm.clear( );
	ifsm.open( szNegatives );
	if( !DatNegatives.Open( ifsm, CDat::EFormatText, 0 ) ) {
		cerr << "Could not open: " << szNegatives << endl;
		return 1; }
	ifsm.close( );
	for( i = 0; i < DatNegatives.GetGenes( ); ++i )
		Genome.AddGene( DatNegatives.GetGene( i ) );
	cerr << "  done" << endl;

	return ( Answers.Open( DatNegatives, vecpPositives, Genome, true ) ? 0 : 1 ); }

int write_posteriors( const string& strFunction, const CBayesNetSmile& BNFunction, ostream& ostm ) {
	size_t					i, iNode;
	vector<unsigned char>	vecbDatum;
	vector<float>			vecdOut;
	float					dPrior;
	unsigned char			bValue;
	vector<string>			vecstrNodes;

	BNFunction.GetNodes( vecstrNodes );
	vecbDatum.resize( vecstrNodes.size( ) );
	for( i = 0; i < vecbDatum.size( ); ++i )
		vecbDatum[ i ] = 0;

	BNFunction.Evaluate( vecbDatum, vecdOut, false );
	dPrior = vecdOut[ vecdOut.size( ) - 1 ];
	ostm << strFunction;
	for( iNode = 1; iNode < vecstrNodes.size( ); ++iNode ) {
		float	d;

		vecbDatum[ iNode - 1 ] = 0;
		vecdOut.clear( );
		for( bValue = 0; bValue < BNFunction.GetValues( iNode ); ++bValue ) {
			vecbDatum[ iNode ] = bValue + 1;
			BNFunction.Evaluate( vecbDatum, vecdOut, false ); }

		d = 0;
		for( i = 0; i < vecdOut.size( ); ++i )
			d += fabs( dPrior - vecdOut[ i ] );
		ostm << '\t' << ( d / BNFunction.GetValues( iNode ) ); }
	ostm << endl;

	return 0; }
