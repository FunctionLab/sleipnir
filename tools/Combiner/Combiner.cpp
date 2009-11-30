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
static int MainRevDATs( const gengetopt_args_info& );
static int MainPCLs( const gengetopt_args_info& );
static int MainModules( const gengetopt_args_info& );

static const TPFnCombiner	c_apfnCombiners[]	= { MainPCLs, MainDATs, MainDABs, MainModules, MainRevDATs, NULL };
static const char*			c_aszCombiners[]	= { "pcl", "dat", "dab", "module", "revdat", NULL };
static const char			c_szMean[]			= "mean";
static const char			c_szGMean[]			= "gmean";
static const char			c_szHMean[]			= "hmean";
static const char			c_szMax[]			= "max";
static const char			c_szMin[]			= "min";
static const char			c_szSum[]			= "sum";
static const char			c_szVote[]			= "vote";

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	int					iRet;
	size_t				i;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg );

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
		if( !PCL.Open( ifsm, sArgs.skip_arg ) ) {
			cerr << "Could not open: " << sArgs.inputs[ iArg ] << endl;
			return 1; }
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

enum EMethod {
	EMethodMean,
	EMethodSum,
	EMethodGMean,
	EMethodHMean,
	EMethodMax,
	EMethodMin
};

int MainDATs( const gengetopt_args_info& sArgs ) {
	CDataset				Dataset;
	CDat					DatOut, DatCur;
	CHalfMatrix<float>		MatCounts;
	size_t					i, j, k, iOne, iTwo, iA, iB;
	vector<vector<size_t> >	vecveciGenes;
	float					d, dWeight1, dWeight2;
	vector<string>			vecstrFiles, vecstrTerms;
	CPCL					PCLWeights( false );
	CGenome					Genome;
	CGenes					GenesIn( Genome );
	vector<CGenes*>			vecpTerms;
	vector<set<size_t> >	vecsetiGenes;
	EMethod					eMethod;

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
	if( sArgs.weights_arg && !PCLWeights.Open( sArgs.weights_arg, 0 ) ) {
		cerr << "Could not open: " << sArgs.weights_arg << endl;
		return 1; }
	if( sArgs.genes_arg && !GenesIn.Open( sArgs.genes_arg ) ) {
		cerr << "Could not open: " << sArgs.genes_arg << endl;
		return 1; }
	if( sArgs.terms_arg && !CGenes::Open( sArgs.terms_arg, Genome, vecstrTerms, vecpTerms ) ) {
		cerr << "Could not open: " << sArgs.terms_arg << endl;
		return 1; }

	if( !strcmp( c_szSum, sArgs.method_arg ) )
		eMethod = EMethodSum;
	else if( !strcmp( c_szGMean, sArgs.method_arg ) )
		eMethod = EMethodGMean;
	else if( !strcmp( c_szHMean, sArgs.method_arg ) )
		eMethod = EMethodHMean;
	else if( !strcmp( c_szMax, sArgs.method_arg ) )
		eMethod = EMethodMax;
	else if( !strcmp( c_szMin, sArgs.method_arg ) )
		eMethod = EMethodMin;
	else
		eMethod = EMethodMean;

	DatOut.Open( vecstrTerms.empty( ) ? Dataset.GetGeneNames( ) : vecstrTerms, false, sArgs.memmap_flag ? sArgs.output_arg : NULL );
	if( !strcmp( c_szMax, sArgs.method_arg ) )
		d = -FLT_MAX;
	else if( !strcmp( c_szMin, sArgs.method_arg ) )
		d = FLT_MAX;
	else if( !strcmp( c_szGMean, sArgs.method_arg ) )
		d = 1;
	else
		d = 0;
	for( i = 0; i < DatOut.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < DatOut.GetGenes( ); ++j )
			DatOut.Set( i, j, d );
	if( fabs( d ) < 2 ) {
		MatCounts.Initialize( DatOut.GetGenes( ) );
		MatCounts.Clear( ); }
	vecsetiGenes.resize( DatOut.GetGenes( ) );
	for( i = 0; i < sArgs.inputs_num; ++i ) {
		if( !DatCur.Open( sArgs.inputs[ i ], !!sArgs.memmap_flag && !sArgs.normalize_flag && !GenesIn.GetGenes( ) ) ) {
			cerr << "Couldn't open: " << sArgs.inputs[ i ] << endl;
			return 1; }
		if( PCLWeights.GetGenes( ) ) {
			if( ( j = PCLWeights.GetGene( CMeta::Deextension( CMeta::Basename( sArgs.inputs[ i ] ) ) ) ) == -1 ) {
				cerr << "Ignoring unweighted graph: " << sArgs.inputs[ i ] << endl;
				continue; }
			dWeight1 = PCLWeights.Get( j, 0 );
			dWeight2 = sArgs.reweight_flag ? 1 : dWeight1; }
		else
			dWeight1 = dWeight2 = 1;
		cerr << "Opened: " << sArgs.inputs[ i ] << endl;
		if( sArgs.normalize_flag )
			DatCur.Normalize( CDat::ENormalizeZScore );
		if( GenesIn.GetGenes( ) )
			DatCur.FilterGenes( GenesIn, CDat::EFilterInclude );
		vecveciGenes.resize( DatCur.GetGenes( ) );
		for( j = 0; j < vecveciGenes.size( ); ++j )
			vecveciGenes[j].clear( );
		for( j = 0; j < DatOut.GetGenes( ); ++j ) {
			vecsetiGenes[j].clear( );
			if( vecstrTerms.empty( ) ) {
				if( ( iOne = DatCur.GetGene( DatOut.GetGene( j ) ) ) != -1 )
					vecveciGenes[iOne].push_back( j ); }
			else
				for( k = 0; k < vecpTerms[j]->GetGenes( ); ++k )
					if( ( iOne = DatCur.GetGene( vecpTerms[j]->GetGene( k ).GetName( ) ) ) != -1 ) {
						vecveciGenes[iOne].push_back( j );
						vecsetiGenes[j].insert( iOne ); } }
		for( j = 0; j < DatCur.GetGenes( ); ++j )
			for( k = ( j + 1 ); k < DatCur.GetGenes( ); ++k ) {
				if( CMeta::IsNaN( d = DatCur.Get( j, k ) ) )
					continue;
				for( iA = 0; iA < vecveciGenes[j].size( ); ++iA ) {
					iOne = vecveciGenes[j][iA];
					if( vecsetiGenes[iOne].find( k ) != vecsetiGenes[iOne].end( ) )
						continue;
					for( iB = 0; iB < vecveciGenes[k].size( ); ++iB ) {
						iTwo = vecveciGenes[k][iB];
						if( vecsetiGenes[iTwo].find( j ) != vecsetiGenes[iTwo].end( ) )
							continue;
						switch( eMethod ) {
							case EMethodGMean:
								DatOut.Get( iOne, iTwo ) *= pow( d, dWeight1 );
								MatCounts.Get( iOne, iTwo ) += dWeight2;
								break;

							case EMethodHMean:
								DatOut.Get( iOne, iTwo ) += dWeight1 / d;
								MatCounts.Get( iOne, iTwo ) += dWeight2;
								break;

							case EMethodMax:
								if( d > DatOut.Get( iOne, iTwo ) )
									DatOut.Set( iOne, iTwo, d );
								break;

							case EMethodMin:
								if( d < DatOut.Get( iOne, iTwo ) )
									DatOut.Set( iOne, iTwo, d );
								break;

							default:
								DatOut.Get( iOne, iTwo ) += dWeight1 * d;
								MatCounts.Get( iOne, iTwo ) += dWeight2; } } } } }
	for( i = 0; i < DatOut.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < DatOut.GetGenes( ); ++j )
			switch( eMethod ) {
				case EMethodMean:
					DatOut.Set( i, j, ( d = MatCounts.Get( i, j ) ) ? ( DatOut.Get( i, j ) / d ) :
						CMeta::GetNaN( ) );
					break;

				case EMethodGMean:
					DatOut.Set( i, j, ( d = MatCounts.Get( i, j ) ) ?
						(float)pow( (double)DatOut.Get( i, j ), 1.0 / d ) : CMeta::GetNaN( ) );
					break;

				case EMethodHMean:
					DatOut.Set( i, j, ( d = MatCounts.Get( i, j ) ) ? ( d / DatOut.Get( i, j ) ) :
						CMeta::GetNaN( ) );
					break;

				case EMethodMax:
					if( DatOut.Get( i, j ) == -FLT_MAX )
						DatOut.Set( i, j, CMeta::GetNaN( ) );
					break;

				case EMethodMin:
					if( DatOut.Get( i, j ) == FLT_MAX )
						DatOut.Set( i, j, CMeta::GetNaN( ) ); }

	if( !sArgs.memmap_flag )
		DatOut.Save( sArgs.output_arg );

	for( i = 0; i < vecpTerms.size( ); ++i )
		delete vecpTerms[i];
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

int MainRevDATs( const gengetopt_args_info& sArgs ) {
	vector<string>			vecstrTerms, vecstrGenes;
	vector<CGenes*>			vecpTerms;
	CGenome					Genome;
	CDat					DatOut;
	CHalfMatrix<uint32_t>	MatCounts;
	size_t					i, j, k, iOne, iTwo, iArg, iA, iB;
	float					d;

	if( !sArgs.terms_arg ) {
		cerr << "Terms argument required" << endl;
		return 1; }
	if( !CGenes::Open( sArgs.terms_arg, Genome, vecstrTerms, vecpTerms ) ) {
		cerr << "Could not open: " << sArgs.terms_arg << endl;
		return 1; }

	{
		set<string>		setstrGenes;

		for( i = 0; i < vecpTerms.size( ); ++i )
			for( j = 0; j < vecpTerms[i]->GetGenes( ); ++j )
				setstrGenes.insert( vecpTerms[i]->GetGene( j ).GetName( ) );
		vecstrGenes.resize( setstrGenes.size( ) );
		copy( setstrGenes.begin( ), setstrGenes.end( ), vecstrGenes.begin( ) );
	}

	DatOut.Open( vecstrGenes, false );
	for( i = 0; i < DatOut.GetGenes( ); ++i )
		memset( DatOut.Get( i ), 0, ( DatOut.GetGenes( ) - i - 1 ) * sizeof(*DatOut.Get( i )) );
	MatCounts.Initialize( DatOut.GetGenes( ) );
	MatCounts.Clear( );
	for( iArg = 0; iArg < sArgs.inputs_num; ++iArg ) {
		CDat					DatIn;
		vector<vector<size_t> >	vecveciGenes;

		if( !DatIn.Open( sArgs.inputs[iArg] ) ) {
			cerr << "Could not open: " << sArgs.inputs[iArg] << endl;
			return 1; }
		vecveciGenes.resize( DatIn.GetGenes( ) );
		for( i = 0; i < DatIn.GetGenes( ); ++i ) {
			for( j = 0; j < vecstrTerms.size( ); ++j )
				if( DatIn.GetGene( i ) == vecstrTerms[j] )
					break;
			if( j < vecstrTerms.size( ) ) {
				for( k = 0; k < vecpTerms[j]->GetGenes( ); ++k )
					if( ( iOne = DatOut.GetGene( vecpTerms[j]->GetGene( k ).GetName( ) ) ) != -1 )
						vecveciGenes[i].push_back( iOne ); }
			else
				cerr << "Unrecognized gene: " << DatIn.GetGene( i ) << endl; }
		for( i = 0; i < DatIn.GetGenes( ); ++i ) {
			if( vecveciGenes[i].empty( ) )
				continue;
			for( j = ( i + 1 ); j < DatIn.GetGenes( ); ++j ) {
				if( CMeta::IsNaN( d = DatIn.Get( i, j ) ) )
					continue;
				for( iA = 0; iA < vecveciGenes[i].size( ); ++iA ) {
					iOne = vecveciGenes[i][iA];
					for( iB = 0; iB < vecveciGenes[j].size( ); ++iB ) {
						iTwo = vecveciGenes[j][iB];
						MatCounts.Get( iOne, iTwo )++;
						DatOut.Get( iOne, iTwo ) += d; } } } } }

	for( i = 0; i < DatOut.GetGenes( ); ++i )
		for( j = ( i + 1 ); j < DatOut.GetGenes( ); ++j )
			if( k = MatCounts.Get( i, j ) )
				DatOut.Get( i, j ) /= k;
			else
				DatOut.Set( i, j, CMeta::GetNaN( ) );
	DatOut.Save( sArgs.output_arg );

	for( i = 0; i < vecpTerms.size( ); ++i )
		delete vecpTerms[i];
	return 0; }
