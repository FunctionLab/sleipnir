#include "stdafx.h"
#include "cmdline.h"

struct SNeighbors {
	CFullMatrix<pair<size_t, float> >	m_MatNeighbors;
	vector<size_t>						m_veciMin;
	vector<size_t>						m_veciColumns;

	void Initialize( const float* adValues, size_t iValues, size_t iNeighbors ) {
		size_t	i, j;

		m_veciColumns.clear( );
		for( i = 0; i < iValues; ++i )
			if( CMeta::IsNaN( adValues[ i ] ) )
				m_veciColumns.push_back( i );

		m_veciMin.resize( m_veciColumns.size( ) );
		m_MatNeighbors.Initialize( iNeighbors, m_veciColumns.size( ) );
		for( i = 0; i < m_MatNeighbors.GetRows( ); ++i )
			for( j = 0; j < m_MatNeighbors.GetColumns( ); ++j )
				m_MatNeighbors.Get( i, j ).second = -FLT_MAX; }

	bool Add( size_t iNeighbor, float dSim, const float* adValues ) {
		size_t	i, j, iCol;
		bool	fRet;

		for( fRet = false,i = 0; i < m_veciColumns.size( ); ++i ) {
			iCol = m_veciColumns[ i ];
			if( !CMeta::IsNaN( adValues[ iCol ] ) && ( dSim > m_MatNeighbors.Get( m_veciMin[ i ], i ).second ) ) {
				fRet = true;
				m_MatNeighbors.Get( m_veciMin[ i ], i ).first = iNeighbor;
				m_MatNeighbors.Get( m_veciMin[ i ], i ).second = dSim;

				for( m_veciMin[ i ] = 0,j = 1; j < m_MatNeighbors.GetRows( ); ++j )
					if( m_MatNeighbors.Get( j, i ).second < m_MatNeighbors.Get( m_veciMin[ i ], i ).second )
						m_veciMin[ i ] = j; } }

		return fRet; }

	size_t GetColumn( size_t iColumn ) const {
		size_t	i;

		for( i = 0; i < m_veciColumns.size( ); ++i )
			if( m_veciColumns[ i ] == iColumn )
				return i;

		return -1; }
};

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CPCL				PCL;
	CDat				Dat;
	vector<SNeighbors>	vecsNeighbors;
	size_t				i, j, k, iCol, iMissing, iOne, iTwo;
	float				d, dSum;
	const float*		ad;
	int					iRet;
	ofstream			ofsm;
	vector<float>		vecdMissing;
	vector<size_t>		veciMissing;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg );

	if( iRet = CPCL::Distance( sArgs.input_arg, sArgs.skip_arg, sArgs.distance_arg, false, false,
		!!sArgs.autocorrelate_flag, sArgs.genes_arg, CMeta::GetNaN( ), sArgs.limit_arg, PCL, Dat ) ) {
		cmdline_parser_print_help( );
		return iRet; }

	vecdMissing.resize( PCL.GetGenes( ) );
	for( i = 0; i < vecdMissing.size( ); ++i ) {
		for( iMissing = j = 0; j < PCL.GetExperiments( ); ++j )
			if( CMeta::IsNaN( PCL.Get( i, j ) ) )
				iMissing++;
		if( ( vecdMissing[ i ] = (float)( PCL.GetExperiments( ) - iMissing ) / PCL.GetExperiments( ) ) >=
			sArgs.missing_arg )
			veciMissing.push_back( i ); }
	vecsNeighbors.resize( veciMissing.size( ) );
	for( i = 0; i < vecsNeighbors.size( ); ++i )
		vecsNeighbors[ i ].Initialize( PCL.Get( veciMissing[ i ] ), PCL.GetExperiments( ), sArgs.neighbors_arg );
	for( i = 0; i < veciMissing.size( ); ++i ) {
		if( !( i % 100 ) )
			cerr << "Finding neighbors for gene " << i << '/' << veciMissing.size( ) << endl;
		ad = PCL.Get( iOne = veciMissing[ i ] );
		for( j = ( i + 1 ); j < veciMissing.size( ); ++j ) {
			d = Dat.Get( iOne, iTwo = veciMissing[ j ] );
			vecsNeighbors[ i ].Add( iTwo, d, PCL.Get( iTwo ) );
			vecsNeighbors[ j ].Add( iOne, d, ad ); } }

	for( iOne = i = 0; i < PCL.GetGenes( ); ++i ) {
		if( vecdMissing[ i ] < sArgs.missing_arg ) {
			PCL.MaskGene( i );
			continue; }
		{
			const SNeighbors&	sGene	= vecsNeighbors[ iOne++ ];

			for( j = 0; j < PCL.GetExperiments( ); ++j ) {
				if( !CMeta::IsNaN( PCL.Get( i, j ) ) )
					continue;

				iCol = sGene.GetColumn( j );
				for( d = dSum = 0,k = 0; k < (size_t)sArgs.neighbors_arg; ++k ) {
					const pair<size_t, float>&	prNeighbor	= sGene.m_MatNeighbors.Get( k, iCol );

					dSum += prNeighbor.second;
					d += PCL.Get( prNeighbor.first, j ) * prNeighbor.second; }
			PCL.Set( i, j, d / dSum ); }
		} }

	if( sArgs.output_arg ) {
		ofsm.open( sArgs.output_arg );
		PCL.Save( ofsm );
		ofsm.close( ); }
	else {
		PCL.Save( cout );
		cout.flush( ); }

	CMeta::Shutdown( );
	return 0; }
