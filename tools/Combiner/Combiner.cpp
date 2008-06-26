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

static int MainDABs( const gengetopt_args_info& );
static int MainDATs( const gengetopt_args_info& );
static int MainPCLs( const gengetopt_args_info& );
static int MainModules( const gengetopt_args_info& );

static const TPFnCombiner	c_apfnCombiners[]	= { MainPCLs, MainDATs, MainDABs, MainModules, NULL };
static const char*			c_aszCombiners[]	= { "pcl", "dat", "dab", "module", NULL };
static const char			c_szMean[]			= "mean";
static const char			c_szGMean[]			= "gmean";
static const char			c_szHMean[]			= "hmean";
static const char			c_szMax[]			= "max";
static const char			c_szMin[]			= "min";
static const char			c_szVote[]			= "vote";

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	int					iRet;
	size_t				i;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta = CMeta( sArgs.verbosity_arg );

	for( i = 0; c_aszCombiners[ i ]; ++i )
		if( !strcmp( c_aszCombiners[ i ], sArgs.type_arg ) ) {
			iRet = c_apfnCombiners[ i ]( sArgs );
			break; }

	return iRet; }

int MainPCLs( const gengetopt_args_info& sArgs ) {
	CPCL						PCL, PCLNew;
	size_t						i, j, iArg, iExp, iGene;
	vector<string>				vecstrGenes, vecstrExps, vecstrFeatures;
	set<string>					setstrGenes;
	set<string>::const_iterator	iterGenes;
	ifstream					ifsm;
	ofstream					ofsm;

	for( iArg = 0; iArg < sArgs.inputs_num; ++iArg ) {
		ifsm.clear( );
		ifsm.open( sArgs.inputs[ iArg ] );
		PCL.Open( ifsm, sArgs.skip_arg );
		if( !iArg )
			for( i = 0; i < PCL.GetFeatures( ); ++i )
				vecstrFeatures.push_back( PCL.GetFeature( i ) );
		for( i = 0; i < PCL.GetExperiments( ); ++i )
			vecstrExps.push_back( PCL.GetExperiment( i ) );
		for( i = 0; i < PCL.GetGenes( ); ++i )
			setstrGenes.insert( PCL.GetGene( i ) );
		ifsm.close( ); }
	vecstrGenes.resize( setstrGenes.size( ) );
	copy( setstrGenes.begin( ), setstrGenes.end( ), vecstrGenes.begin( ) );

	PCLNew.Open( vecstrGenes, vecstrExps, vecstrFeatures );
	iExp = 0;
	for( iArg = 0; iArg < sArgs.inputs_num; ++iArg ) {
		cerr << "Processing " << sArgs.inputs[ iArg ] << "..." << endl;
		ifsm.clear( );
		ifsm.open( sArgs.inputs[ iArg ] );
		PCL.Open( ifsm, sArgs.skip_arg );
		for( i = 0; i < PCLNew.GetGenes( ); ++i )
			if( ( iGene = PCL.GetGene( vecstrGenes[ i ] ) ) != -1 ) {
				if( !iArg )
					for( j = 1; j < PCLNew.GetFeatures( ); ++j )
						PCLNew.SetFeature( i, j, PCL.GetFeature( iGene, j ) );
				for( j = 0; j < PCL.GetExperiments( ); ++j )
					PCLNew.Set( i, iExp + j, PCL.Get( iGene, j ) ); }
		iExp += PCL.GetExperiments( );
		ifsm.close( ); }

	if( sArgs.output_arg ) {
		ofsm.open( sArgs.output_arg );
		PCLNew.Save( ofsm );
		ofsm.close( ); }
	else {
		PCLNew.Save( cout );
		cout.flush( ); }

	return 0; }

int MainDATs( const gengetopt_args_info& sArgs ) {
	CDataset					Dataset;
	CDat						DatOut, DatCur;
	CHalfMatrix<unsigned char>	HMatCounts;
	size_t						i, j, k, iOne, iTwo;
	vector<size_t>				veciGenes;
	float						d;
	vector<string>				vecstrFiles;

	if( !sArgs.inputs_num )
		return 1;

	vecstrFiles.resize( sArgs.inputs_num );
	copy( sArgs.inputs, sArgs.inputs + sArgs.inputs_num, vecstrFiles.begin( ) );
	if( !Dataset.OpenGenes( vecstrFiles ) ) {
		cerr << "Couldn't open: " << vecstrFiles[ 0 ];
		for( i = 1; i < vecstrFiles.size( ); ++i )
			cerr << ", " << vecstrFiles[ i ];
		cerr << endl;
		return 1; }

	DatOut.Open( Dataset.GetGeneNames( ), false, sArgs.memmap_flag ? sArgs.output_arg : NULL );
	if( !strcmp( c_szMax, sArgs.method_arg ) )
		d = -FLT_MAX;
	else if( !strcmp( c_szMin, sArgs.method_arg ) )
		d = FLT_MAX;
	else
		d = 0;
	for( i = 0; i < DatOut.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < DatOut.GetGenes( ); ++j )
			DatOut.Set( i, j, d );
	if( !d ) {
		HMatCounts.Initialize( DatOut.GetGenes( ) );
		HMatCounts.Clear( ); }
	veciGenes.resize( DatOut.GetGenes( ) );
	for( i = 0; i < sArgs.inputs_num; ++i ) {
		if( !DatCur.Open( sArgs.inputs[ i ], !!sArgs.memmap_flag && !sArgs.normalize_flag ) ) {
			cerr << "Couldn't open: " << sArgs.inputs[ i ] << endl;
			return 1; }
		cerr << "Opened: " << sArgs.inputs[ i ] << endl;
		if( sArgs.normalize_flag )
			DatCur.Normalize( CDat::ENormalizeZScore );
		for( j = 0; j < veciGenes.size( ); ++j )
			veciGenes[ j ] = DatCur.GetGene( DatOut.GetGene( j ) );
		for( j = 0; j < veciGenes.size( ); ++j ) {
			if( ( iOne = veciGenes[ j ] ) == -1 )
				continue;
			for( k = ( j + 1 ); k < veciGenes.size( ); ++k ) {
				if( ( ( iTwo = veciGenes[ k ] ) == -1 ) || CMeta::IsNaN( d = DatCur.Get( iOne, iTwo ) ) )
					continue;
				if( !strcmp( c_szMean, sArgs.method_arg ) ) {
					DatOut.Get( j, k ) += d;
					HMatCounts.Get( j, k )++; }
				else if( !strcmp( c_szGMean, sArgs.method_arg ) ) {
					DatOut.Get( j, k ) *= d;
					HMatCounts.Get( j, k )++; }
				else if( !strcmp( c_szHMean, sArgs.method_arg ) ) {
					DatOut.Get( j, k ) += 1 / d;
					HMatCounts.Get( j, k )++; }
				else if( !strcmp( c_szMax, sArgs.method_arg ) ) {
					if( d > DatOut.Get( j, k ) )
						DatOut.Set( j, k, d ); }
				else if( !strcmp( c_szMin, sArgs.method_arg ) ) {
					if( d < DatOut.Get( j, k ) )
						DatOut.Set( j, k, d ); } } } }
	for( i = 0; i < DatOut.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < DatOut.GetGenes( ); ++j )
			if( !strcmp( c_szMean, sArgs.method_arg ) )
				DatOut.Set( i, j, ( k = HMatCounts.Get( i, j ) ) ? ( DatOut.Get( i, j ) / k ) :
					CMeta::GetNaN( ) );
			else if( !strcmp( c_szGMean, sArgs.method_arg ) )
				DatOut.Set( i, j, ( k = HMatCounts.Get( i, j ) ) ?
					(float)pow( (double)DatOut.Get( i, j ), 1.0 / k ) : CMeta::GetNaN( ) );
			else if( !strcmp( c_szHMean, sArgs.method_arg ) )
				DatOut.Set( i, j, ( k = HMatCounts.Get( i, j ) ) ? ( k / DatOut.Get( i, j ) ) :
					CMeta::GetNaN( ) );
			else if( !strcmp( c_szMax, sArgs.method_arg ) ) {
				if( DatOut.Get( i, j ) == -FLT_MAX )
					DatOut.Set( i, j, CMeta::GetNaN( ) ); }
			else if( !strcmp( c_szMin, sArgs.method_arg ) ) {
				if( DatOut.Get( i, j ) == FLT_MAX )
					DatOut.Set( i, j, CMeta::GetNaN( ) ); }

	if( !sArgs.memmap_flag )
		DatOut.Save( sArgs.output_arg );

	return 0; }

static int MainDABs( const gengetopt_args_info& sArgs ) {
	CDatasetCompact	Dataset;
	size_t			i;
	vector<string>	vecstrFiles;
	ofstream		ofsm;

	if( !sArgs.inputs_num )
		return 1;

	vecstrFiles.resize( sArgs.inputs_num );
	copy( sArgs.inputs, sArgs.inputs + sArgs.inputs_num, vecstrFiles.begin( ) );
	if( !Dataset.Open( vecstrFiles ) ) {
		cerr << "Couldn't open: " << vecstrFiles[ 0 ];
		for( i = 1; i < vecstrFiles.size( ); ++i )
			cerr << ", " << vecstrFiles[ i ];
		cerr << endl;
		return 1; }

	if( sArgs.output_arg ) {
		ofsm.open( sArgs.output_arg );
		Dataset.Save( ofsm, true );
		ofsm.close( ); }
	else {
		Dataset.Save( cout, false );
		cout.flush( ); }

	return 0; }

struct SSimpleome {
	map<string, size_t>	m_mapstriGenes;
	map<size_t, string>	m_mapistrGenes;

	size_t Get( const string& strGene ) {
		map<string, size_t>::const_iterator	iterGene;
		size_t								iRet;

		if( ( iterGene = m_mapstriGenes.find( strGene ) ) != m_mapstriGenes.end( ) )
			return iterGene->second;

		m_mapstriGenes[ strGene ] = iRet = m_mapstriGenes.size( );
		m_mapistrGenes[ iRet ] = strGene;
		return iRet; }

	const string& Get( size_t iGene ) const {

		return m_mapistrGenes.find( iGene )->second; }
};

struct SModule {
	float				m_dSpecificity;
	set<size_t>			m_setiGenes;
	set<const SModule*>	m_setpsChildren;

	static float Open( const char* szFile, SSimpleome& sSimpleome, vector<SModule*>& vecpsModules ) {
		static const size_t	c_iBuffer	= 131072;
		ifstream		ifsm;
		char			acBuffer[ c_iBuffer ];
		vector<string>	vecstrLine;
		float			dRet;
		size_t			i, iModule;
		SModule*		psModule;

		ifsm.open( szFile );
		if( !ifsm.is_open( ) ) {
			cerr << "Could not open: " << szFile << endl;
			return CMeta::GetNaN( ); }
		for( iModule = 0; !ifsm.eof( ); ++iModule ) {
			ifsm.getline( acBuffer, c_iBuffer - 1 );
			if( !acBuffer[ 0 ] )
				iModule--; }
		ifsm.close( );
		if( !iModule ) {
			cerr << "No modules found in: " << szFile << endl;
			return CMeta::GetNaN( ); }
		cerr << "Found " << --iModule << " modules in: " << szFile << endl;
		vecpsModules.resize( iModule );

		dRet = CMeta::GetNaN( );
		ifsm.clear( );
		ifsm.open( szFile );
		for( iModule = 0; !ifsm.eof( ); ++iModule ) {
			ifsm.getline( acBuffer, c_iBuffer - 1 );
			acBuffer[ c_iBuffer - 1 ] = 0;
			if( CMeta::IsNaN( dRet ) ) {
				dRet = (float)atof( acBuffer );
				iModule--;
				continue; }
			vecstrLine.clear( );
			CMeta::Tokenize( acBuffer, vecstrLine );
			if( vecstrLine.size( ) < 2 ) {
				iModule--;
				continue; }
			vecpsModules[ iModule ] = psModule = new SModule( );
			psModule->m_dSpecificity = (float)atof( vecstrLine[ 0 ].c_str( ) );
			for( i = 1; i < vecstrLine.size( ); ++i )
				psModule->m_setiGenes.insert( sSimpleome.Get( vecstrLine[ i ] ) ); }

		return dRet; }

	float Jaccard( const SModule& sModule ) const {
		set<size_t>::const_iterator	iterGene;
		size_t						iUnion, iIntersection;

		iIntersection = 0;
		iUnion = sModule.m_setiGenes.size( );
		for( iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene )
			if( sModule.m_setiGenes.find( *iterGene ) != sModule.m_setiGenes.end( ) )
				iIntersection++;
			else
				iUnion++;

		return ( (float)iIntersection / iUnion ); }

	size_t Intersection( const SModule& sModule ) const {
		size_t						iRet;
		set<size_t>::const_iterator	iterGene;

		for( iRet = 0,iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene )
			if( sModule.m_setiGenes.find( *iterGene ) != sModule.m_setiGenes.end( ) )
				iRet++;

		return iRet; }

	void Merge( const SModule& sModule ) {
		set<size_t>::const_iterator	iterGene;

		for( iterGene = sModule.m_setiGenes.begin( ); iterGene != sModule.m_setiGenes.end( ); ++iterGene )
			m_setiGenes.insert( *iterGene );
		m_dSpecificity = ( m_dSpecificity + sModule.m_dSpecificity ) / 2; }

	bool IsChild( const SModule* psModule ) const {

		return ( m_setpsChildren.find( psModule ) != m_setpsChildren.end( ) ); }

	void Save( ostream& ostm, float dCutoff, const SSimpleome& sSimpleome ) const {
		set<const SModule*>::const_iterator	iterChild;
		set<size_t>::const_iterator			iterGene;

		ostm << this << '\t' << dCutoff << '\t' << m_dSpecificity << '\t';
		for( iterChild = m_setpsChildren.begin( ); iterChild != m_setpsChildren.end( ); ++iterChild ) {
			if( iterChild != m_setpsChildren.begin( ) )
				ostm << '|';
			ostm << *iterChild; }
		for( iterGene = m_setiGenes.begin( ); iterGene != m_setiGenes.end( ); ++iterGene )
			ostm << '\t' << sSimpleome.Get( *iterGene );
		ostm << endl; }
};

struct SSorterModules {
	const vector<float>&	m_vecdModules;

	SSorterModules( const vector<float>& vecdModules ) : m_vecdModules(vecdModules) { }

	bool operator()( size_t iOne, size_t iTwo ) const {

		return ( m_vecdModules[ iOne ] > m_vecdModules[ iTwo ] ); }
};

int MainModules( const gengetopt_args_info& sArgs ) {
	vector<float>				vecdModules;
	vector<vector<SModule*> >	vecvecpsModules;
	vector<size_t>				veciIndices;
	size_t						i, j, iOuter, iCutoffOne, iCutoffTwo, iModuleOne, iModuleTwo;
	SSimpleome					sSimpleome;
	bool						fDone;
	float						d;
	ofstream					ofsm;
	ostream*					postm;

	vecdModules.resize( sArgs.inputs_num );
	vecvecpsModules.resize( sArgs.inputs_num );
	veciIndices.resize( sArgs.inputs_num );
	for( i = 0; i < vecvecpsModules.size( ); ++i ) {
		veciIndices[ i ] = i;
		if( CMeta::IsNaN( vecdModules[ i ] = SModule::Open( sArgs.inputs[ i ], sSimpleome,
			vecvecpsModules[ i ] ) ) )
			return 1; }
	sort( veciIndices.begin( ), veciIndices.end( ), SSorterModules( vecdModules ) );

	for( iOuter = 0,fDone = false; !fDone; ++iOuter ) {
		fDone = true;
		cerr << "Outer loop: " << iOuter << endl;
		for( iCutoffOne = 0; iCutoffOne < veciIndices.size( ); ++iCutoffOne ) {
			vector<SModule*>&	vecpsModulesOne	= vecvecpsModules[ veciIndices[ iCutoffOne ] ];

			cerr << "Merging cutoff: " << vecdModules[ veciIndices[ iCutoffOne ] ] << endl;
			for( iModuleOne = 0; iModuleOne < vecpsModulesOne.size( ); ++iModuleOne ) {
				SModule*	psOne	= vecpsModulesOne[ iModuleOne ];

				if( !psOne )
					continue;
				for( iModuleTwo = ( iModuleOne + 1 ); iModuleTwo < vecpsModulesOne.size( ); ++iModuleTwo ) {
					SModule*	psTwo	= vecpsModulesOne[ iModuleTwo ];
					
					if( !psTwo )
						continue;
					if( ( d = psOne->Jaccard( *psTwo ) ) >= sArgs.jaccard_arg ) {
						cerr << "Merging @" << d << ' ' << iCutoffOne << ':' << iModuleOne << " (" <<
							psOne->m_setiGenes.size( ) << ") " << iCutoffOne << ':' << iModuleTwo << " (" <<
							psTwo->m_setiGenes.size( ) << ')' << endl;
						psOne->Merge( *psTwo );
						delete vecpsModulesOne[ iModuleTwo ];
						vecpsModulesOne[ iModuleTwo ] = NULL;
						iModuleTwo--;
						fDone = false; } } }
			for( iCutoffTwo = ( iCutoffOne + 1 ); iCutoffTwo < veciIndices.size( ); ++iCutoffTwo ) {
				vector<SModule*>&	vecpsModulesTwo	= vecvecpsModules[ veciIndices[ iCutoffTwo ] ];

				for( iModuleOne = 0; iModuleOne < vecpsModulesOne.size( ); ++iModuleOne ) {
					SModule*	psOne	= vecpsModulesOne[ iModuleOne ];

					if( !psOne )
						continue;
					for( iModuleTwo = 0; iModuleTwo < vecpsModulesTwo.size( ); ++iModuleTwo ) {
						SModule*	psTwo	= vecpsModulesTwo[ iModuleTwo ];
						
						if( !psTwo )
							continue;
						if( ( d = psOne->Jaccard( *psTwo ) ) >= sArgs.jaccard_arg ) {
							cerr << "Merging @" << d << ' ' << iCutoffOne << ':' << iModuleOne << " (" <<
								psOne->m_setiGenes.size( ) << ") " << iCutoffTwo << ':' << iModuleTwo <<
								" (" << psTwo->m_setiGenes.size( ) << ')' << endl;
							psOne->Merge( *psTwo );
							delete vecpsModulesTwo[ iModuleTwo ];
							vecpsModulesTwo[ iModuleTwo ] = NULL;
							iModuleTwo--;
							fDone = false; } } } } } }

	for( iCutoffOne = 0; iCutoffOne < veciIndices.size( ); ++iCutoffOne ) {
		vector<SModule*>&	vecpsModulesOne	= vecvecpsModules[ veciIndices[ iCutoffOne ] ];

		for( iCutoffTwo = ( iCutoffOne + 1 ); iCutoffTwo < veciIndices.size( ); ++iCutoffTwo ) {
			vector<SModule*>&	vecpsModulesTwo	= vecvecpsModules[ veciIndices[ iCutoffTwo ] ];

			for( iModuleOne = 0; iModuleOne < vecpsModulesOne.size( ); ++iModuleOne ) {
				SModule*	psOne	= vecpsModulesOne[ iModuleOne ];

				if( !psOne )
					continue;
				for( iModuleTwo = 0; iModuleTwo < vecpsModulesTwo.size( ); ++iModuleTwo ) {
					SModule*	psTwo	= vecpsModulesTwo[ iModuleTwo ];
					
					if( !psTwo )
						continue;
					if( !psTwo->IsChild( psOne ) && ( ( d = ( (float)psOne->Intersection( *psTwo ) /
						psOne->m_setiGenes.size( ) ) ) >= sArgs.intersection_arg ) ) {
						cerr << vecdModules[ veciIndices[ iCutoffOne ] ] << ':' << iModuleOne <<
							" child of " << vecdModules[ veciIndices[ iCutoffTwo ] ] << ':' << iModuleTwo <<
							" (" << psOne->m_setiGenes.size( ) << ", " << psTwo->m_setiGenes.size( ) << ", " <<
							d << ')' << endl;
						psTwo->m_setpsChildren.insert( psOne ); } } } } }

	if( sArgs.output_arg ) {
		ofsm.open( sArgs.output_arg );
		postm = &ofsm; }
	else
		postm = &cout;
	for( iCutoffOne = 0; iCutoffOne < veciIndices.size( ); ++iCutoffOne ) {
		vector<SModule*>&	vecpsModulesOne	= vecvecpsModules[ veciIndices[ iCutoffOne ] ];

		for( iModuleOne = 0; iModuleOne < vecpsModulesOne.size( ); ++iModuleOne ) {
			SModule*	psOne	= vecpsModulesOne[ iModuleOne ];

			if( psOne )
				psOne->Save( *postm, vecdModules[ veciIndices[ iCutoffOne ] ], sSimpleome ); } }
	if( sArgs.output_arg )
		ofsm.close( );

	for( i = 0; i < vecvecpsModules.size( ); ++i )
		for( j = 0; j < vecvecpsModules[ i ].size( ); ++j )
			if( vecvecpsModules[ i ][ j ] )
				delete vecvecpsModules[ i ][ j ];

	return 0; }
