#include "stdafx.h"

int main( int iArgs, char** aszArgs ) {
	CPCL				PCL;
	CDat				Dat;
	vector<size_t>		veciClusters;
	size_t				i, j, iClusters, iSize;
	double				dDiameter;
	vector<string>		vecstrGenes;
	ofstream			ofsm;
	const char*			szOutput;
	bool				fAutoc;
	CMeasureKendallsTau	KendallsTau;

	if( iArgs < 4 ) {
		cerr << "Usage: " << aszArgs[ 0 ] << " <diameter> <size> <output.dab> [autocorrelate]" <<
			endl;
		return 1; }
	dDiameter = atof( aszArgs[ 1 ] );
	iSize = atoi( aszArgs[ 2 ] );
	szOutput = aszArgs[ 3 ];
	fAutoc = ( iArgs > 4 ) ? !!atoi( aszArgs[ 4 ] ) : false;

	PCL.Open( cin, 0 );
	iClusters = CClustQTC::Cluster( PCL.Get( ), &KendallsTau, (float)dDiameter, iSize,
		fAutoc, veciClusters );

	vecstrGenes.reserve( PCL.GetGenes( ) );
	for( i = 0; i < PCL.GetGenes( ); ++i )
		vecstrGenes.push_back( PCL.GetGene( i ) );
	Dat.Open( vecstrGenes );
	for( i = 0; i < veciClusters.size( ); ++i ) {
		if( ( veciClusters[ i ] + 1 ) == iClusters )
			continue;
		for( j = ( i + 1 ); j < veciClusters.size( ); ++j ) {
			if( ( veciClusters[ j ] + 1 ) == iClusters )
				continue;
			Dat.Set( i, j, ( veciClusters[ i ] == veciClusters[ j ] ) ? 1.0f : 0.0f ); } }
	ofsm.open( szOutput, ios_base::binary );
	Dat.Save( ofsm, true );
	ofsm.close( );

	return 0; }
