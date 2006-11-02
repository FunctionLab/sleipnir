#include "stdafx.h"
#include "cmdline.h"

int Dat2Dab( const CGenes&, const gengetopt_args_info& );
int DabShuffle( const char*, bool );
int DabFinder( const char*, bool );
void DatFlip( CDat& );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	int					iRet;
	CGenome				Genome;
	CGenes				Genes( Genome );
	ifstream			ifsm;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( sArgs.genes_arg ) {
		ifsm.open( sArgs.genes_arg );
		if( !Genes.Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.genes_arg << endl;
			return 1; }
		ifsm.close( ); }
	iRet = Dat2Dab( Genes, sArgs );

	CMeta::Shutdown( );
	return iRet; }

int Dat2Dab( const CGenes& Genes, const gengetopt_args_info& sArgs ) {
	CDat		Dat;
	ifstream	ifsm;
	ofstream	ofsm;
	CPCL		PCL;

	if( sArgs.pcl_arg ) {
		ifsm.open( sArgs.pcl_arg );
		if( !PCL.Open( ifsm, sArgs.skip_arg ) ) {
			cerr << "Couldn't open: " << sArgs.pcl_arg << endl;
			return 1; }
		ifsm.close( );
		{
			CMeasurePearNorm	PearNorm;
//			CDatPCL				DatPCL( PCL, &PearNorm );

//			DatPCL.Save( cout, false );
		} }
	else if( !!sArgs.output_arg ) {
		Dat.Open( cin, false );
		if( sArgs.normalize_flag || sArgs.zscore_flag )
			Dat.Normalize( !!sArgs.normalize_flag );
		if( sArgs.flip_flag )
			DatFlip( Dat );
		if( Genes.GetGenes( ) )
			Dat.FilterGenes( Genes, CDat::EFilterInclude );
		ofsm.open( sArgs.output_arg, ios_base::binary );
		Dat.Save( ofsm, true );
		ofsm.close( ); }
	else {
		ifsm.open( sArgs.input_arg, ios_base::binary );
		Dat.Open( ifsm, true );
		if( sArgs.normalize_flag || sArgs.zscore_flag )
			Dat.Normalize( !!sArgs.normalize_flag );
		if( sArgs.flip_flag )
			DatFlip( Dat );
		if( Genes.GetGenes( ) )
			Dat.FilterGenes( Genes, CDat::EFilterInclude );
		ifsm.close( );
		Dat.Save( cout, false );
		cout.flush( ); }

	return 0; }

int DabShuffle( const char* szFile, bool fCreate ) {
	CDat		Dat;
	ifstream	ifsm;
	ofstream	ofsm;
	size_t		i, j;

	ifsm.open( szFile, ios_base::binary );
	Dat.Open( ifsm, true );
	ifsm.close( );

	srand( (int)time( NULL ) );
	for( i = 0; i < Dat.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
			Dat.Set( i, j, !!( rand( ) % 2 ) );

	ofsm.open( szFile, ios_base::binary );
	Dat.Save( ofsm, true );
	ofsm.close( );

	return 0; }

int DabFinder( const char* szFile, bool fCreate ) {
	CDat		Dat;
	ifstream	ifsm;
	size_t		i, j;
	float		d;

	ifsm.open( szFile, ios_base::binary );
	Dat.Open( ifsm, true );
	ifsm.close( );

	for( i = 0; i < Dat.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
			if( ( d = Dat.Get( i, j ) ) != -1 )
				cout << Dat.GetGene( i ) << '\t' << Dat.GetGene( j ) << '\t' << d << endl;

	return 0; }

void DatFlip( CDat& Dat ) {
	size_t	i, j;
	float	d;

	for( i = 0; i < Dat.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
			if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) )
				Dat.Set( i, j, 1 - d ); }
