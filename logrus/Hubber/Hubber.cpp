#include "stdafx.h"
#include "cmdline.h"

static const char	c_szDab[]	= ".dab";

struct SDatum {
	float						m_dHubbiness;
	float						m_dHubbinessStd;
	float						m_dCliquiness;
	float						m_dCliquinessStd;
	vector<pair<size_t,float> >	m_vecprSpecific;
};

void hubs( const CDat&, vector<float>& );
void cliques( const CDat&, const vector<float>&, bool, SDatum&, const CGenes* );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CGenome				Genome;
	CDat				Dat;
	size_t				i, iGenes;
	vector<float>		vecdHub;
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

	hubs( Dat, vecdHub );
	if( sArgs.genes_arg == -1 ) {
		cout << "Function";
		for( i = 0; i < Dat.GetGenes( ); ++i )
			cout << '\t' << Dat.GetGene( i ); }
	else {
		cliques( Dat, vecdHub, true, sDatum, NULL );
		cout << "name	size	hubbiness	hubbiness std.	cliquiness	cliquiness std." << endl;
		cout << "total	" << Dat.GetGenes( ) << '\t' << sDatum.m_dHubbiness << '\t' <<
			sDatum.m_dHubbinessStd << '\t' << sDatum.m_dCliquiness << '\t' << sDatum.m_dCliquinessStd; }
	cout << endl;

	for( iGenes = 0; iGenes < sArgs.inputs_num; ++iGenes ) {
		CGenes		Genes( Genome );
		ifstream	ifsm;
		size_t		i;

		if( !( iGenes % 25 ) )
			cerr << iGenes << '/' << sArgs.inputs_num << endl;
		ifsm.open( sArgs.inputs[ iGenes ] );
		if( !Genes.Open( ifsm ) ) {
			cerr << "Could not open: " << sArgs.inputs[ iGenes ] << endl;
			return 1; }
		ifsm.close( );
		cliques( Dat, vecdHub, sArgs.genes_arg != -1, sDatum, &Genes );
		cout << CMeta::Basename( sArgs.inputs[ iGenes ] );
		if( sArgs.genes_arg == -1 )
			for( i = 0; i < sDatum.m_vecprSpecific.size( ); ++i )
				cout << '\t' << ( sDatum.m_vecprSpecific[ i ].second *
					( Genes.IsGene( Dat.GetGene( sDatum.m_vecprSpecific[ i ].first ) ) ? -1 : 1 ) );
		else {
			cout << '\t' << Genes.GetGenes( ) << '\t' << sDatum.m_dHubbiness << '\t' <<
				sDatum.m_dHubbinessStd << '\t' << sDatum.m_dCliquiness << '\t' <<
				sDatum.m_dCliquinessStd;
			for( i = 0; i < min( (size_t)sArgs.genes_arg, sDatum.m_vecprSpecific.size( ) ); ++i )
				cout << '\t' << Dat.GetGene( sDatum.m_vecprSpecific[ i ].first ) << '|' <<
					sDatum.m_vecprSpecific[ i ].second << '|' <<
					( Genes.IsGene( Dat.GetGene( sDatum.m_vecprSpecific[ i ].first ) ) ? 1 : 0 ); }
		cout << endl;

		if( sArgs.clip_arg ) {
			CDat			DatOut;
			size_t			j;
			vector<bool>	vecfGenes;

			vecfGenes.resize( Dat.GetGenes( ) );
			for( i = 0; i < vecfGenes.size( ); ++i )
				vecfGenes[ i ] = Genes.IsGene( Dat.GetGene( i ) );
			DatOut.Open( Dat.GetGeneNames( ) );
			for( i = 0; i < DatOut.GetGenes( ); ++i )
				for( j = ( i + 1 ); j < DatOut.GetGenes( ); ++j )
					if( vecfGenes[ i ] || vecfGenes[ j ] )
						DatOut.Set( i, j, Dat.Get( i, j ) );
			DatOut.Save( ( (string)sArgs.clip_arg + '/' + CMeta::Deextension( CMeta::Basename(
				sArgs.inputs[ iGenes ] ) ) + c_szDab ).c_str( ) ); } }

	CMeta::Shutdown( );
	return 0; }

void hubs( const CDat& Dat, vector<float>& vecdHub ) {
	size_t			i, j;
	float			d;
	vector<size_t>	veciHub;

	vecdHub.resize( Dat.GetGenes( ) );
	veciHub.resize( vecdHub.size( ) );
	for( i = 0; i < vecdHub.size( ); ++i )
		vecdHub[ i ] = 0;
	for( i = 0; i < Dat.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j ) {
			if( CMeta::IsNaN( d = Dat.Get( i, j ) ) )
				continue;
			veciHub[ i ]++;
			vecdHub[ i ] += d;
			veciHub[ j ]++;
			vecdHub[ j ] += d; }
	for( i = 0; i < vecdHub.size( ); ++i )
		vecdHub[ i ] /= veciHub[ i ]; }

struct SSorter {

	bool operator()( const pair<size_t,float>& prOne, const pair<size_t,float>& prTwo ) const {

		return ( prOne.second > prTwo.second ); }
};

void cliques( const CDat& Dat, const vector<float>& vecdHub, bool fSort, SDatum& sDatum,
	const CGenes* pGenes ) {
	size_t			i, j;
	float			d;
	vector<float>	vecdClique;
	vector<size_t>	veciClique;
	vector<bool>	vecfOutside;

	vecfOutside.resize( Dat.GetGenes( ) );
	if( pGenes )
		for( i = 0; i < vecfOutside.size( ); ++i )
			if( !pGenes->IsGene( Dat.GetGene( i ) ) )
				vecfOutside[ i ] = true;
	veciClique.resize( Dat.GetGenes( ) );
	vecdClique.resize( Dat.GetGenes( ) );
	for( i = 0; i < Dat.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j ) {
			if( CMeta::IsNaN( d = Dat.Get( i, j ) ) )
				continue;
			if( !vecfOutside[ i ] ) {
				veciClique[ j ]++;
				vecdClique[ j ] += d; }
			if( !vecfOutside[ j ] ) {
				veciClique[ i ]++;
				vecdClique[ i ] += d; } }
	for( sDatum.m_dCliquiness = sDatum.m_dCliquinessStd = 0,i = 0; i < vecdClique.size( ); ++i ) {
		vecdClique[ i ] /= veciClique[ i ];
		if( !vecfOutside[ i ] ) {
			sDatum.m_dCliquiness += vecdClique[ i ];
			sDatum.m_dCliquinessStd += vecdClique[ i ] * vecdClique[ i ]; } }
	i = pGenes ? pGenes->GetGenes( ) : Dat.GetGenes( );
	sDatum.m_dCliquiness /= i;
	sDatum.m_dCliquinessStd = (float)sqrt( ( sDatum.m_dCliquinessStd / i ) -
		( sDatum.m_dCliquiness * sDatum.m_dCliquiness ) );

	sDatum.m_dHubbiness = sDatum.m_dHubbinessStd = 0;
	for( i = 0; i < vecdClique.size( ); ++i )
		sDatum.m_dHubbiness += vecdClique[ i ];
	sDatum.m_dHubbiness = ( sDatum.m_dHubbiness - sDatum.m_dCliquiness ) / vecdClique.size( );

	sDatum.m_vecprSpecific.resize( Dat.GetGenes( ) );
	for( i = 0; i < sDatum.m_vecprSpecific.size( ); ++i ) {
		sDatum.m_vecprSpecific[ i ].first = i;
		sDatum.m_vecprSpecific[ i ].second = vecdClique[ i ] / vecdHub[ i ]; }
	if( fSort )
		sort( sDatum.m_vecprSpecific.begin( ), sDatum.m_vecprSpecific.end( ), SSorter( ) ); }
