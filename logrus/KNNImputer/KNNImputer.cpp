#include "stdafx.h"
#include "cmdline.h"

struct SNeighbors {
	CFullMatrix<pair<size_t, float> >	m_MatNeighbors;
	vector<size_t>						m_veciMin;

	void Initialize( size_t iColumns, size_t iNeighbors ) {
		size_t	i, j;

		m_veciMin.resize( iColumns );
		m_MatNeighbors.Initialize( iNeighbors, iColumns );
		for( i = 0; i < m_MatNeighbors.GetRows( ); ++i )
			for( j = 0; j < m_MatNeighbors.GetColumns( ); ++j )
				m_MatNeighbors.Get( i, j ).second = -FLT_MAX; }

	bool Add( size_t iNeighbor, float dSim, const float* adValues ) {
		size_t	i, j;
		bool	fRet;

		for( fRet = false,i = 0; i < m_MatNeighbors.GetColumns( ); ++i )
			if( !CMeta::IsNaN( adValues[ i ] ) && ( dSim > m_MatNeighbors.Get( m_veciMin[ i ], i ).second ) ) {
				fRet = true;
				m_MatNeighbors.Get( m_veciMin[ i ], i ).first = iNeighbor;
				m_MatNeighbors.Get( m_veciMin[ i ], i ).second = dSim;

				for( m_veciMin[ i ] = 0,j = 1; j < m_MatNeighbors.GetRows( ); ++j )
					if( m_MatNeighbors.Get( j, i ).second < m_MatNeighbors.Get( m_veciMin[ i ], i ).second )
						m_veciMin[ i ] = j; }

		return fRet; }
};

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CPCL				PCL;
	CDat				Dat;
	vector<SNeighbors>	vecsNeighbors;
	size_t				i, j, k;
	float				d, dSum;
	const float*		ad;
	int					iRet;
	ofstream			ofsm;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( iRet = CPCL::Distance( sArgs.input_arg, sArgs.skip_arg, sArgs.distance_arg, false, false,
		!!sArgs.autocorrelate_flag, sArgs.genes_arg, CMeta::GetNaN( ), sArgs.limit_arg, PCL, Dat ) ) {
		cmdline_parser_print_help( );
		return iRet; }

	vecsNeighbors.resize( Dat.GetGenes( ) );
	for( i = 0; i < vecsNeighbors.size( ); ++i )
		vecsNeighbors[ i ].Initialize( PCL.GetExperiments( ), sArgs.neighbors_arg );
	for( i = 0; i < Dat.GetGenes( ); ++i ) {
		if( !( i % 100 ) )
			cerr << "Finding neighbors for gene " << i << '/' << Dat.GetGenes( ) << endl;
		// This won't work if only a subset of genes has been selected for the dat...
		ad = PCL.Get( i );
		for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j ) {
			d = Dat.Get( i, j );
			vecsNeighbors[ i ].Add( j, d, PCL.Get( j ) );
			vecsNeighbors[ j ].Add( i, d, ad ); } }

	for( i = 0; i < PCL.GetGenes( ); ++i )
		for( j = 0; j < PCL.GetExperiments( ); ++j ) {
			if( !CMeta::IsNaN( PCL.Get( i, j ) ) )
				continue;
			for( d = dSum = 0,k = 0; k < (size_t)sArgs.neighbors_arg; ++k ) {
				const pair<size_t, float>&	prNeighbor	= vecsNeighbors[ i ].m_MatNeighbors.Get( k, j );

				dSum += prNeighbor.second;
				d += PCL.Get( prNeighbor.first, j ) * prNeighbor.second; }
			PCL.Set( i, j, d / dSum ); }

	if( sArgs.output_arg ) {
		ofsm.open( sArgs.output_arg );
		PCL.Save( ofsm );
		ofsm.close( ); }
	else {
		PCL.Save( cout );
		cout.flush( ); }

	CMeta::Shutdown( );
	return 0; }
