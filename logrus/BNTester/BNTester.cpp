#include "stdafx.h"
#include "cmdline.h"

static int Genes( const char*, CGenes& );

int main( int iArgs, char** aszArgs ) {
	IBayesNet*				pBN;
	CDatasetCompact			Data;
	CDatasetCompactMap		DataMap;
	CDataset				DataFull;
	IDataset*				pData;
	CDat					Dat;
	CGenome					Genome;
	CGenes					GenesIn( Genome ), GenesEx( Genome ), GenesOv( Genome );
	ifstream				ifsm;
	vector<vector<float> >	vecvecdResults;
	float					d;
	size_t					i, j, k;
	gengetopt_args_info		sArgs;
	int						iRet;
	vector<bool>			vecfGenes;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );
#if !( defined(_MSC_VER) && defined(_DEBUG) )
	EnableXdslFormat( );
#endif // !( defined(_MSC_VER) && defined(_DEBUG) )

	CBayesNetSmile	BNSmile( !!sArgs.group_flag );
#ifdef PNL_ENABLED
	CBayesNetPNL	BNPNL( !!sArgs.group_flag );
#endif // PNL_ENABLED
	CBayesNetFN		BNFN;

	if( sArgs.function_flag ) {
		if( !BNFN.Open( sArgs.input_arg ) ) {
			cerr << "Couldn't open: " << sArgs.input_arg << endl;
			return 1; }
		pBN = &BNFN; }
	else {
		if( !BNSmile.Open( sArgs.input_arg ) ) {
			cerr << "Couldn't open: " << sArgs.input_arg << endl;
			return 1; }
#ifdef PNL_ENABLED
		if( sArgs.pnl_flag ) {
			BNSmile.Convert( BNPNL );
			pBN = &BNPNL; }
		else
#endif // PNL_ENABLED
			pBN = &BNSmile; }

	if( ( iRet = Genes( sArgs.genes_arg, GenesIn ) ) ||
		( iRet = Genes( sArgs.genee_arg, GenesOv ) ) ||
		( iRet = Genes( sArgs.genex_arg, GenesEx ) ) )
		return iRet;
	if( pBN->IsContinuous( ) ) {
		if( !DataFull.Open( sArgs.datadir_arg, pBN ) ) {
			cerr << "Couldn't open: " << sArgs.datadir_arg << endl;
			return 1; }
		DataFull.FilterGenes( GenesIn, CDat::EFilterInclude );
		DataFull.FilterGenes( GenesEx, CDat::EFilterExclude );
		pData = &DataFull; }
	else if( sArgs.database_arg ) {
		CDatabase			Database;
		CBayesNetMinimal	BNMinimal;

		if( !BNMinimal.Open( BNSmile ) ) {
			cerr << "Couldn't open: " << sArgs.input_arg << endl;
			return 1; }

	else if( sArgs.dataset_arg ) {
		if( !DataMap.Open( sArgs.dataset_arg ) ) {
			cerr << "Couldn't open: " << sArgs.dataset_arg << endl;
			return 1; }
		if( sArgs.genes_arg && !DataMap.FilterGenes( sArgs.genes_arg, CDat::EFilterInclude ) ) {
			cerr << "Couldn't open: " << sArgs.genes_arg << endl;
			return 1; }
		if( sArgs.genex_arg && !DataMap.FilterGenes( sArgs.genex_arg, CDat::EFilterExclude ) ) {
			cerr << "Couldn't open: " << sArgs.genex_arg << endl;
			return 1; }
		if( !sArgs.everything_flag )
			DataMap.FilterAnswers( );
		pData = &DataMap; }
	else {
		if( !Data.Open( sArgs.datadir_arg, pBN, GenesIn, GenesEx ) ) {
			cerr << "Couldn't open: " << sArgs.datadir_arg << endl;
			return 1; }
		pData = &Data; }
	pData->FilterGenes( GenesOv, CDat::EFilterTerm );

	if( sArgs.output_arg )
		Dat.Open( sArgs.genes_arg ? GenesIn.GetGeneNames( ) : pData->GetGeneNames( ) );
	vecvecdResults.clear( );
	cerr << "Evaluating..." << endl;
	if( sArgs.output_arg )
		pBN->Evaluate( pData, Dat, !!sArgs.zero_flag );
	else
		pBN->Evaluate( pData, vecvecdResults, !!sArgs.zero_flag );

	if( sArgs.output_arg ) {
		cerr << "Saving..." << endl;
		for( i = 0; i < Dat.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) )
					Dat.Set( i, j, 1 - d );
		Dat.Save( sArgs.output_arg ); }
	else {
		cerr << "Storing..." << endl;
		for( k = i = 0; i < pData->GetGenes( ); ++i )
			for( j = ( i + 1 ); j < pData->GetGenes( ); ++j ) {
				if( !pData->IsExample( i, j ) )
					continue;
				d = vecvecdResults[ k++ ][ 0 ];
				if( !pBN->IsContinuous( ) )
					d = 1 - d;
				cout << pData->GetGene( i ) << '\t' << pData->GetGene( j ) << '\t' << d <<
					endl; }
		cout.flush( ); }

	CMeta::Shutdown( );
	return 0; }

static int Genes( const char* szGenes, CGenes& Genes ) {
	ifstream	ifsm;

	if( !szGenes )
		return 0;

	ifsm.open( szGenes );
	if( !Genes.Open( ifsm ) ) {
		cerr << "Couldn't open: " << szGenes << endl;
		return 1; }
	ifsm.close( );
	return 0; }
