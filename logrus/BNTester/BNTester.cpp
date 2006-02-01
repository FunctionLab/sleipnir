#include "stdafx.h"
#include "cmdline.h"

static int Genes( const char*, CGenes& );

int main( int iArgs, char** aszArgs ) {
	IBayesNet*				pBN;
	CDatasetCompact			Data;
	CDatasetCompactMap		DataMap;
	CDatasetCompact*		pData;
	CDat					Dat;
	CGenome					Genome;
	CGenes					GenesIn( Genome ), GenesEx( Genome ), GenesOv( Genome );
	ifstream				ifsm;
	ofstream				ofsm;
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
#if defined(_MSC_VER) && !defined(_DEBUG)
	EnableXdslFormat( );
#endif // defined(_MSC_VER) && !defined(_DEBUG)

	CBayesNetSmile	BNSmile( !!sArgs.group_flag );
	CBayesNetPNL	BNPNL( !!sArgs.group_flag );

	if( !BNSmile.Open( sArgs.input_arg ) ) {
		cerr << "Couldn't open: " << sArgs.input_arg << endl;
		return 1; }
	if( ( iRet = Genes( sArgs.genes_arg, GenesIn ) ) ||
		( iRet = Genes( sArgs.genee_arg, GenesOv ) ) ||
		( iRet = Genes( sArgs.genex_arg, GenesEx ) ) )
		return iRet;
	if( sArgs.dataset_arg ) {
		if( !DataMap.Open( sArgs.dataset_arg ) ) {
			cerr << "Couldn't open: " << sArgs.dataset_arg << endl;
			return 1; }
		if( sArgs.genes_arg && !DataMap.FilterGenes( sArgs.genes_arg, CDat::EFilterInclude ) ) {
			cerr << "Couldn't open: " << sArgs.genes_arg << endl;
			return 1; }
		if( sArgs.genex_arg && !DataMap.FilterGenes( sArgs.genex_arg, CDat::EFilterExclude ) ) {
			cerr << "Couldn't open: " << sArgs.genex_arg << endl;
			return 1; }
		DataMap.FilterAnswers( );
		pData = &DataMap; }
	else {
		if( !Data.Open( sArgs.datadir_arg, &BNSmile, GenesIn, GenesEx ) ) {
			cerr << "Couldn't open: " << sArgs.datadir_arg << endl;
			return 1; }
		pData = &Data; }
	pData->FilterGenes( GenesOv, CDat::EFilterInclude );

	if( sArgs.output_arg )
		Dat.Open( pData->GetGeneNames( ) );
	vecvecdResults.clear( );
	cerr << "Evaluating..." << endl;
	if( sArgs.pnl_flag ) {
		BNSmile.Convert( BNPNL );
		pBN = &BNPNL; }
	else
		pBN = &BNSmile;
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
		ofsm.open( sArgs.output_arg, ios_base::binary );
		Dat.Save( ofsm, true );
		ofsm.close( ); }
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
