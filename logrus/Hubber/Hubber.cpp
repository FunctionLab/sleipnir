#include "stdafx.h"
#include "cmdline.h"

struct SDatum {
	float			m_dHubbiness;
	float			m_dHubbinessStd;
	float			m_dCliquiness;
	float			m_dCliquinessStd;
	vector<size_t>	m_veciCliquish;
	vector<float>	m_vecdCliquish;
};

void count( const CDat&, SDatum&, const CGenes* );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CGenome				Genome;
	CDat				Dat;
	size_t				iGenes;
	SDatum				sDatum;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( sArgs.input_arg ) {
		if( !Dat.Open( sArgs.input_arg, !!sArgs.memmap_flag && !sArgs.normalize_flag ) ) {
			cerr << "Could not open: " << sArgs.input_arg << endl;
			return 1; } }
	else if( !Dat.Open( cin, CDat::EFormatText ) ) {
		cerr << "Could not open input" << endl;
		return 1; }
	if( sArgs.normalize_flag )
		Dat.Normalize( );

	cout << "name	size	hubbiness	hubbiness std.	cliquiness	cliquiness std." << endl;
	count( Dat, sDatum, NULL );
	cout << "total	" << Dat.GetGenes( ) << '\t' << sDatum.m_dHubbiness << '\t' <<
		sDatum.m_dHubbinessStd << '\t' << sDatum.m_dCliquiness << '\t' << sDatum.m_dCliquinessStd <<
		endl;

	for( iGenes = 0; iGenes < sArgs.inputs_num; ++iGenes ) {
		CGenes		Genes( Genome );
		ifstream	ifsm;

		if( !( iGenes % 25 ) )
			cerr << iGenes << '/' << sArgs.inputs_num << endl;
		ifsm.open( sArgs.inputs[ iGenes ] );
		if( !Genes.Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.inputs[ iGenes ] << endl;
			return 1; }
		ifsm.close( );
		count( Dat, sDatum, &Genes );
		cout << CMeta::Basename( sArgs.inputs[ iGenes ] ) << '\t' << Genes.GetGenes( ) << '\t' <<
			sDatum.m_dHubbiness << '\t' << sDatum.m_dHubbinessStd << '\t' << sDatum.m_dCliquiness <<
			'\t' << sDatum.m_dCliquinessStd << endl; }

	CMeta::Shutdown( );
	return 0; }

void count( const CDat& Dat, SDatum& sDatum, const CGenes* pGenes ) {
	size_t			i, j, iOne, iTwo;
	float			d;
	vector<float>	vecdHub, vecdClique;
	vector<size_t>	veciHub, veciClique, veciGenes;

	if( pGenes ) {
		veciGenes.resize( Dat.GetGenes( ) );
		for( i = 0; i < veciGenes.size( ); ++i )
			veciGenes[ i ] = pGenes->GetGene( Dat.GetGene( i ) ); }
	sDatum.m_vecdCliquish.resize( Dat.GetGenes( ) );
	sDatum.m_veciCliquish.resize( Dat.GetGenes( ) );
	veciHub.resize( pGenes ? pGenes->GetGenes( ) : Dat.GetGenes( ) );
	vecdHub.resize( veciHub.size( ) );
	veciClique.resize( veciHub.size( ) );
	vecdClique.resize( veciHub.size( ) );
	for( i = 0; i < Dat.GetGenes( ); ++i ) {
		iOne = pGenes ? veciGenes[ i ] : i;
		for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j ) {
			iTwo = pGenes ? veciGenes[ j ] : j;
			if( ( ( iOne == -1 ) && ( iTwo == -1 ) ) || CMeta::IsNaN( d = Dat.Get( i, j ) ) )
				continue;
			if( iOne != -1 ) {
				veciHub[ iOne ]++;
				vecdHub[ iOne ] += d;
				if( iTwo != -1 ) {
					veciClique[ iOne ]++;
					veciClique[ iTwo ]++;
					vecdClique[ iOne ] += d;
					vecdClique[ iTwo ] += d; } }
			if( iTwo != -1 ) {
				veciHub[ iTwo ]++;
				vecdHub[ iTwo ] += d; } } }
	for( i = 0; i < vecdHub.size( ); ++i ) {
		vecdHub[ i ] /= veciHub[ i ];
		vecdClique[ i ] /= veciClique[ i ]; }

	sDatum.m_dHubbiness = (float)CStatistics::Average( vecdHub );
	sDatum.m_dHubbinessStd = (float)sqrt( CStatistics::Variance( vecdHub, sDatum.m_dHubbiness ) );
	sDatum.m_dCliquiness = (float)CStatistics::Average( vecdClique );
	sDatum.m_dCliquinessStd = (float)sqrt( CStatistics::Variance( vecdClique, sDatum.m_dCliquiness ) ); }
