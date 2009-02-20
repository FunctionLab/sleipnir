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
#include "cmdline.h"

enum EFile {
	EFilePCL,
	EFileWIG,
	EFileError
};

int main_postprocess( const gengetopt_args_info&, CCoalesceMotifLibrary& );
int main_test( const gengetopt_args_info&, const CPCL&, const CFASTA&, CCoalesceMotifLibrary& );
int main_test2( const gengetopt_args_info&, const CPCL&, const CFASTA&, CCoalesceMotifLibrary& );
bool recluster( const gengetopt_args_info&, size_t, CCoalesceMotifLibrary&, const CHierarchy&,
	const vector<CCoalesceCluster>&, const vector<string>&, const CPCL&, size_t& );
bool trim( const gengetopt_args_info&, const CPCL&, vector<CCoalesceCluster>& );

EFile open_pclwig( const char* szFile, size_t iSkip, CFASTA& FASTA, CPCL& PCL ) {

	if( FASTA.Open( szFile ) )
		return EFilePCL;
	if( PCL.Open( szFile, iSkip ) )
		return EFileWIG;

	cerr << "Could not open: " << szFile << endl;
	return EFileError; }

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info			sArgs;
	CFASTA						FASTA, FASTANucleosomes, FASTAConservation;
	CPCL						PCL;
	CCoalesce					Coalesce;
	vector<CCoalesceCluster>	vecClusters;
	size_t						i, j;
	set<string>					setstrTypes;

#ifdef WIN32
	pthread_win32_process_attach_np( );
#endif // WIN32

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg, sArgs.random_arg );

	CCoalesceMotifLibrary		Motifs( sArgs.k_arg );
	Motifs.SetPenaltyGap( (float)sArgs.penalty_gap_arg );
	Motifs.SetPenaltyMismatch( (float)sArgs.penalty_mismatch_arg );

	if( sArgs.postprocess_arg )
		return main_postprocess( sArgs, Motifs );

	if( sArgs.sequences_arg ) {
		vector<string>	vecstrTypes;

		CMeta::Tokenize( sArgs.sequences_arg, vecstrTypes, "," );
		for( i = 0; i < vecstrTypes.size( ); ++i )
			setstrTypes.insert( vecstrTypes[ i ] ); }

	if( !PCL.Open( sArgs.input_arg, sArgs.skip_arg ) ) {
		cerr << "Could not open: " << ( sArgs.input_arg ? sArgs.input_arg : "stdin" ) << endl;
		return 1; }
	if( sArgs.fasta_arg && !FASTA.Open( sArgs.fasta_arg, setstrTypes ) ) {
		cerr << "Could not open: " << sArgs.fasta_arg << endl;
		return 1; }

/*
size_t					iGene	= 4;
uint32_t				iOne;
float					d;
vector<SFASTASequence>	vecsSequences;
SCoalesceModifiers		sMods;
SCoalesceModifierCache	sCache( sMods );
vector<float>			vecdScores;
vector<size_t>			veciLengths;
CCoalesceGeneScores		GeneScores;
float*					ad;

GeneScores.SetGenes( PCL.GetGenes( ) );
iOne = Motifs.Open( "C:1 G:3 A:3 T:3 G:3 A:3 G:3 (A:1)|(C:1)" );
FASTA.Get( FASTA.GetGene( PCL.GetGene( iGene ) ), vecsSequences );
for( i = 0; i < vecsSequences.size( ); ++i )
	if( vecsSequences[ i ].m_strType == "5" ) {
		d = Motifs.GetMatch( vecsSequences[ i ].m_vecstrSequences[ 0 ], iOne, 0, sCache );
		cerr << ( d / vecsSequences[ i ].m_vecstrSequences[ 0 ].length( ) ) << endl;

		GeneScores.Add( iGene, Motifs, vecsSequences[ i ], sCache, iOne, vecdScores, veciLengths );
		ad = GeneScores.Get( GeneScores.GetType( vecsSequences[ i ].m_strType ),
			CCoalesceSequencerBase::ESubsequenceTotal, iGene );
		cerr << ad[ iOne ] << endl;
		return 0; }
//*/
//return main_test2( sArgs, PCL, FASTA, Motifs );

	if( sArgs.datasets_arg ) {
		static const size_t	c_iBuffer	= 131072;
		ifstream			ifsm;
		char				acBuffer[ c_iBuffer ];

		ifsm.open( sArgs.datasets_arg );
		if( !ifsm.is_open( ) ) {
			cerr << "Could not open: " << sArgs.datasets_arg << endl;
			return 1; }
		while( !ifsm.eof( ) ) {
			vector<string>	vecstrLine;
			set<size_t>		setiConditions;

			ifsm.getline( acBuffer, c_iBuffer - 1 );
			acBuffer[ c_iBuffer - 1 ] = 0;
			CMeta::Tokenize( acBuffer, vecstrLine );
			for( i = 0; i < vecstrLine.size( ); ++i )
				for( j = 0; j < PCL.GetExperiments( ); ++j )
					if( vecstrLine[ i ] == PCL.GetExperiment( j ) ) {
						setiConditions.insert( j );
						break; }
			Coalesce.AddDataset( setiConditions ); } }

	Coalesce.SetMotifs( Motifs );
	Coalesce.SetProbabilityGene( (float)sArgs.prob_gene_arg );
	Coalesce.SetPValueCondition( (float)sArgs.pvalue_cond_arg );
	Coalesce.SetPValueMotif( (float)sArgs.pvalue_motif_arg );
	Coalesce.SetPValueCorrelation( (float)sArgs.pvalue_correl_arg );
	Coalesce.SetNumberCorrelation( sArgs.number_correl_arg );
	Coalesce.SetPValueMerge( (float)sArgs.pvalue_merge_arg );
	Coalesce.SetCutoffMerge( (float)sArgs.cutoff_merge_arg );
	Coalesce.SetBasesPerMatch( sArgs.bases_arg );
	Coalesce.SetSizeMinimum( sArgs.size_minimum_arg );
	Coalesce.SetSizeMerge( sArgs.size_merge_arg );
	Coalesce.SetSizeMaximum( sArgs.size_maximum_arg );
	Coalesce.SetThreads( sArgs.threads_arg );
	if( sArgs.output_arg )
		Coalesce.SetDirectoryIntermediate( sArgs.output_arg );
	if( sArgs.cache_arg )
		Coalesce.SetSequenceCache( sArgs.cache_arg );
	if( sArgs.nucleosomes_arg ) {
		if( !FASTANucleosomes.Open( sArgs.nucleosomes_arg ) ) {
			cerr << "Could not open: " << sArgs.nucleosomes_arg << endl;
			return 1; }
		Coalesce.AddWiggle( FASTANucleosomes ); }
	if( sArgs.conservation_arg ) {
		if( !FASTAConservation.Open( sArgs.conservation_arg ) ) {
			cerr << "Could not open: " << sArgs.conservation_arg << endl;
			return 1; }
		Coalesce.AddWiggle( FASTAConservation ); }
	if( !Coalesce.Cluster( PCL, FASTA, vecClusters ) ) {
		cerr << "Clustering failed" << endl;
		return 1; }

	for( i = 0; i < vecClusters.size( ); ++i )
		vecClusters[ i ].Save( cout, i, PCL, &Motifs );

#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; }

int main_postprocess( const gengetopt_args_info& sArgs, CCoalesceMotifLibrary& Motifs ) {
	string						strDir, strFile, strBase;
	vector<CCoalesceCluster>	vecClustersFrom, vecClustersTo;
	bool						fFailed;
	CDistanceMatrix				MatSim;
	size_t						i, j, iDatasets;
	CHierarchy*					pHier;
	vector<string>				vecstrClusters;
	CPCL						PCL;
	ifstream					ifsm;

	if( sArgs.input_arg ) {
		if( !PCL.Open( sArgs.input_arg, sArgs.skip_arg ) ) {
			cerr << "Could not open: " << sArgs.input_arg << endl;
			return 1; } }
	else if( !PCL.Open( cin, sArgs.skip_arg ) ) {
		cerr << "Could not open: stdin" << endl;
		return 1; }
	if( sArgs.known_motifs_arg ) {
		ifsm.open( sArgs.known_motifs_arg );
		if( !Motifs.OpenKnown( ifsm ) ) {
			cerr << "Could not open: " << sArgs.known_motifs_arg << endl;
			return 1; }
		ifsm.close( ); }
	fFailed = false;
	strDir = sArgs.postprocess_arg;
	FOR_EACH_DIRECTORY_FILE(strDir, strFile)
		if( !CMeta::IsExtension( strFile, CPCL::GetExtension( ) ) )
			continue;
		strBase = CMeta::Deextension( strFile );
		strFile = strDir + '/' + strFile;

		if( !fFailed )
			vecClustersFrom.push_back( CCoalesceCluster( ) );
		if( fFailed = ( ( iDatasets = vecClustersFrom.back( ).Open( strFile, sArgs.skip_arg, PCL,
			&Motifs ) ) == -1 ) ) {
			cerr << "Could not open: " << strFile << endl;
			continue; }
		if( !( vecstrClusters.size( ) % 50 ) )
			cerr << "Opened cluster " << vecstrClusters.size( ) << "..." << endl;
		vecstrClusters.push_back( strBase ); }
	if( fFailed )
		vecClustersFrom.pop_back( );
	if( vecClustersFrom.empty( ) ) {
		char	szBuffer[ 16 ];

		fFailed = false;
		ifsm.clear( );
		ifsm.open( sArgs.postprocess_arg );
		while( !ifsm.eof( ) ) {
			if( !fFailed )
				vecClustersFrom.push_back( CCoalesceCluster( ) );
			if( fFailed = ( ( iDatasets = vecClustersFrom.back( ).Open( ifsm, PCL, &Motifs ) ) == -1 ) )
				continue;
			if( !( vecstrClusters.size( ) % 50 ) )
				cerr << "Opened cluster " << vecstrClusters.size( ) << "..." << endl;
			sprintf_s( szBuffer, "c%06d", vecstrClusters.size( ) );
			vecstrClusters.push_back( szBuffer ); }
		if( fFailed )
			vecClustersFrom.pop_back( ); }

	cerr << "Trimming clusters..." << endl;
	if( !trim( sArgs, PCL, vecClustersFrom ) )
		return 1;

	cerr << "Calculating cluster similarities..." << endl;
	MatSim.Initialize( vecClustersFrom.size( ) );
	for( i = 0; i < MatSim.GetSize( ); ++i )
		for( j = ( i + 1 ); j < MatSim.GetSize( ); ++j ) {
			if( sArgs.verbosity_arg >= 7 )
				cerr << "Comparing clusters:	" << i << '\t' << j << endl;
			MatSim.Set( i, j, vecClustersFrom[ i ].GetSimilarity( vecClustersFrom[ j ], PCL.GetGenes( ),
				iDatasets ) ); }

	i = 0;
	return ( ( ( pHier = CClustHierarchical::Cluster( MatSim ) ) &&
		recluster( sArgs, MatSim.GetSize( ) * ( MatSim.GetSize( ) - 1 ) / 2, Motifs, *pHier, vecClustersFrom,
		vecstrClusters, PCL, i ) ) ? 0 : 1 ); }

bool recluster( const gengetopt_args_info& sArgs, size_t iPairs, CCoalesceMotifLibrary& Motifs,
	const CHierarchy& Hier, const vector<CCoalesceCluster>& vecClustersFrom,
	const vector<string>& vecstrClustersFrom, const CPCL& PCL, size_t& iID ) {

	if( Hier.IsGene( ) || ( Hier.GetSimilarity( ) >= sArgs.cutoff_postprocess_arg ) ) {
		CCoalesceCluster	Cluster;

		cerr << "Creating output cluster " << iID << endl;
		if( !Cluster.Open( Hier, vecClustersFrom, vecstrClustersFrom,
			(float)sArgs.fraction_postprocess_arg, (float)sArgs.cutoff_merge_arg, sArgs.max_motifs_arg,
			&Motifs ) )
			return false;
		if( Cluster.GetGenes( ).size( ) < (size_t)sArgs.size_minimum_arg ) {
			cerr << "Cluster too small: " << Cluster.GetGenes( ).size( ) << endl;
			return true; }

		Cluster.RemoveMotifs( Motifs, (float)sArgs.min_zscore_arg );
		if( !Cluster.LabelMotifs( Motifs, (float)sArgs.penalty_gap_arg, (float)sArgs.penalty_mismatch_arg,
			(float)sArgs.known_cutoff_arg ) )
			return false;
		if( sArgs.output_arg )
			Cluster.Save( sArgs.output_arg, iID, PCL, &Motifs );
		Cluster.Save( cout, iID, PCL, &Motifs, (float)sArgs.min_info_arg, (float)sArgs.penalty_gap_arg,
			(float)sArgs.penalty_mismatch_arg, !!sArgs.remove_rcs_flag );
		iID++;
		return true; }

	return ( recluster( sArgs, iPairs, Motifs, Hier.Get( false ), vecClustersFrom, vecstrClustersFrom, PCL,
		iID ) && recluster( sArgs, iPairs, Motifs, Hier.Get( true ), vecClustersFrom, vecstrClustersFrom,
		PCL, iID ) ); }

bool trim( const gengetopt_args_info& sArgs, const CPCL& PCL, vector<CCoalesceCluster>& vecClusters ) {
	vector<size_t>				veciCounts, veciGenes, veciRemove;
	CHalfMatrix<size_t>			MatCounts;
	CDistanceMatrix				MatScores;
	size_t						i, j, iOne, iCluster;
	float						dAve, dStd, d, dCutoff;
	vector<float>				vecdScores;

	veciCounts.resize( PCL.GetGenes( ) );
	MatCounts.Initialize( PCL.GetGenes( ) );
	MatCounts.Clear( );
	for( iCluster = 0; iCluster < vecClusters.size( ); ++iCluster ) {
		const CCoalesceCluster&	Cluster	= vecClusters[ iCluster ];

		veciGenes.resize( Cluster.GetGenes( ).size( ) );
		copy( Cluster.GetGenes( ).begin( ), Cluster.GetGenes( ).end( ), veciGenes.begin( ) );
		for( i = 0; i < veciGenes.size( ); ++i ) {
			veciCounts[ iOne = veciGenes[ i ] ]++;
			for( j = ( i + 1 ); j < veciGenes.size( ); ++j )
				MatCounts.Get( iOne, veciGenes[ j ] )++; } }
	MatScores.Initialize( MatCounts.GetSize( ) );
	for( i = 0; i < MatScores.GetSize( ); ++i )
		for( j = ( i + 1 ); j < MatScores.GetSize( ); ++j ) {
			iOne = MatCounts.Get( i, j );
			iCluster = veciCounts[ i ] + veciCounts[ j ] - iOne;
			MatScores.Set( i, j, iCluster ? ( (float)iOne / iCluster ) : CMeta::GetNaN( ) ); }

	dAve = dStd = 0;
	for( iOne = i = 0; i < MatScores.GetSize( ); ++i )
		for( j = ( i + 1 ); j < MatScores.GetSize( ); ++j )
			if( !CMeta::IsNaN( d = MatScores.Get( i, j ) ) ) {
				iOne++;
				dAve += d;
				dStd += d * d; }
	dAve /= iOne;
	dStd = sqrt( ( dStd / iOne ) - ( dAve * dAve ) );
	dCutoff = dAve + ( (float)sArgs.cutoff_trim_arg * dStd );
	if( sArgs.verbosity_arg >= 6 )
		cerr << "Cutoff " << dCutoff << " (" << dAve << ',' << dStd << ')' << endl;

	for( iCluster = 0; iCluster < vecClusters.size( ); ++iCluster ) {
		CCoalesceCluster&	Cluster	= vecClusters[ iCluster ];

		if( sArgs.verbosity_arg >= 6 )
			cerr << "Trimming cluster " << iCluster << endl;
		veciGenes.resize( Cluster.GetGenes( ).size( ) );
		copy( Cluster.GetGenes( ).begin( ), Cluster.GetGenes( ).end( ), veciGenes.begin( ) );
		vecdScores.resize( veciGenes.size( ) );
		for( i = 0; i < vecdScores.size( ); ++i ) {
			iOne = veciGenes[ i ];
			for( j = ( i + 1 ); j < vecdScores.size( ); ++j ) {
				d = MatScores.Get( iOne, veciGenes[ j ] );
				vecdScores[ i ] += d;
				vecdScores[ j ] += d; } }
		for( i = 0; i < vecdScores.size( ); ++i )
			vecdScores[ i ] /= veciGenes.size( ) - 1;
		if( sArgs.verbosity_arg >= 7 )
			for( i = 0; i < vecdScores.size( ); ++i )
				cerr << "Gene " << i << " (" << PCL.GetGene( veciGenes[ i ] ) << ") " << vecdScores[ i ] <<
					endl;

		veciRemove.clear( );
		for( i = 0; i < vecdScores.size( ); ++i )
			if( vecdScores[ i ] < dCutoff ) {
				if( sArgs.verbosity_arg >= 6 )
					cerr << "Removing " << PCL.GetGene( veciGenes[ i ] ) << " at " << vecdScores[ i ] << endl;
				veciRemove.push_back( veciGenes[ i ] ); }
		Cluster.RemoveGenes( veciRemove ); }

	return true; }

struct STestFASTA {
	string			m_strSequence;
	const CFASTA*	m_pFASTA;
	size_t			m_iGene;
};

void* ThreadTestFASTA( void* pData ) {
	STestFASTA*				psData;
	vector<SFASTASequence>	vecsSequences;

	psData = (STestFASTA*)pData;
	psData->m_pFASTA->Get( psData->m_iGene, vecsSequences );
	if( !( vecsSequences.empty( ) || vecsSequences[ 0 ].m_vecstrSequences.empty( ) ) )
		psData->m_strSequence = vecsSequences[ 0 ].m_vecstrSequences[ 0 ];

	return NULL; }

struct SThreadCombineMotif {
	size_t							m_iOffset;
	size_t							m_iStep;
	const std::vector<size_t>*		m_pveciPCL2FASTA;
	CCoalesceGeneScores*			m_pGeneScores;
	const CCoalesceMotifLibrary*	m_pMotifs;
	uint32_t						m_iMotif;
	const CFASTA*					m_pFASTA;
	const SCoalesceModifiers*		m_psModifiers;
};

void* ThreadCombineMotif( void* pData ) {
	SThreadCombineMotif*	psData;
	size_t					i, j;
	vector<float>			vecdScores;
	vector<size_t>			veciLengths;

	psData = (SThreadCombineMotif*)pData;
	SCoalesceModifierCache	sModifiers( *psData->m_psModifiers );

	for( i = psData->m_iOffset; i < psData->m_pveciPCL2FASTA->size( ); i += psData->m_iStep )
		if( (*psData->m_pveciPCL2FASTA)[ i ] != -1 ) {
			vector<SFASTASequence>	vecsSequences;

			if( psData->m_pFASTA->Get( (*psData->m_pveciPCL2FASTA)[ i ], vecsSequences ) ) {
				sModifiers.Get( i );
				for( j = 0; j < vecsSequences.size( ); ++j )
					if( !psData->m_pGeneScores->Add( i, *psData->m_pMotifs, vecsSequences[ j ],
						sModifiers, psData->m_iMotif, vecdScores, veciLengths ) )
						break; } }

	return NULL; }

int main_test( const gengetopt_args_info& sArgs, const CPCL& PCL, const CFASTA& FASTA, CCoalesceMotifLibrary& Motifs ) {
	CCoalesceGeneScores			GeneScores;
	vector<size_t>				veciPCL2FASTA;
	SCoalesceModifiers			sMods;
	vector<vector<float> >		vecvecdCounts;
	vector<size_t>				veciLengths;
	SCoalesceModifierCache		sModifiers( sMods );
	size_t						i, j, k;
	uint32_t					iMotifOne, iMotifTwo;
	vector<uint32_t>			veciMotifs;
	vector<float>				vecdScores;
	float						d;
	CCoalesceCluster			Cluster, Pot;
	CCoalesceGroupHistograms	HistsCluster( 12, 1.0f / sArgs.bases_arg );
	CCoalesceGroupHistograms	HistsPot( 12, 1.0f / sArgs.bases_arg );
	CCoalesceSequencerBase::ESubsequence	eSubsequence;
	unsigned short				s;
	vector<pthread_t>			vecpthdThreads;
	vector<SThreadCombineMotif>	vecsThreads;

/*
	vector<STestFASTA>			vecsThreads;

	vecpthdThreads.resize( sArgs.threads_arg );
	vecsThreads.resize( vecpthdThreads.size( ) );
	for( i = 0; i < 100; ++i ) {
		iFASTA = rand( ) % FASTA.GetGenes( );
		cerr << i << '\t' << FASTA.GetGene( iFASTA ) << endl;
		for( j = 0; j < vecpthdThreads.size( ); ++j ) {
			vecsThreads[ j ].m_iGene = iFASTA;
			vecsThreads[ j ].m_pFASTA = &FASTA;
			vecsThreads[ j ].m_strSequence = "";
			pthread_create( &vecpthdThreads[ j ], NULL, ThreadTestFASTA, &vecsThreads[ j ] ); }
		for( j = 0; j < vecpthdThreads.size( ); ++j ) {
			pthread_join( vecpthdThreads[ j ], NULL );
			if( j && ( vecsThreads[ j ].m_strSequence != vecsThreads[ 0 ].m_strSequence ) )
				cerr << "Thread " << j << " failed:" << endl << vecsThreads[ 0 ].m_strSequence << endl <<
					vecsThreads[ j ].m_strSequence << endl; }
		cerr << vecsThreads[ 0 ].m_strSequence << endl; }
	return 0;
//*/

	veciPCL2FASTA.resize( PCL.GetGenes( ) );
	for( i = 0; i < veciPCL2FASTA.size( ); ++i )
		veciPCL2FASTA[ i ] = FASTA.GetGene( PCL.GetGene( i ) );
	sMods.Initialize( PCL );

	GeneScores.SetGenes( PCL.GetGenes( ) );
	for( i = 0; i < veciPCL2FASTA.size( ); ++i )
		if( veciPCL2FASTA[ i ] != -1 ) {
			vector<SFASTASequence>	vecsSequences;

			if( FASTA.Get( veciPCL2FASTA[ i ], vecsSequences ) ) {
				sModifiers.Get( i );
				for( j = 0; j < vecsSequences.size( ); ++j )
					if( !GeneScores.Add( i, Motifs, vecsSequences[ j ], sModifiers, vecvecdCounts,
						veciLengths ) )
						return 1; } }

	Pot.SetGenes( PCL.GetGenes( ) );
	Pot.CalculateHistograms( GeneScores, HistsPot, NULL );
	for( i = 0; i < PCL.GetGenes( ); ++i )
		if( ( (float)rand( ) / RAND_MAX ) < ( 1.0 / 100 ) )
			Cluster.Add( i, Pot );

	veciMotifs.resize( 20 );
	vecpthdThreads.resize( sArgs.threads_arg );
	vecsThreads.resize( vecpthdThreads.size( ) );
	for( i = 0; i < 100; ++i ) {
		Cluster.CalculateHistograms( GeneScores, HistsCluster, &HistsPot );
		Cluster.Snapshot( GeneScores, HistsCluster );
		Pot.Snapshot( GeneScores, HistsPot );
		cerr << "Cluster	" << Cluster.GetGenes( ).size( ) << endl << "Pot	" <<
			Pot.GetGenes( ).size( ) << endl;

		for( j = 0; j < veciMotifs.size( ); ++j ) {
			iMotifOne = rand( ) % Motifs.GetMotifs( );
			iMotifTwo = rand( ) % Motifs.GetMotifs( );
			veciMotifs[ j ] = Motifs.Merge( iMotifOne, iMotifTwo, FLT_MAX, false );
			cerr << "Merged:" << endl << Motifs.GetMotif( iMotifOne ) << endl <<
				Motifs.GetMotif( iMotifTwo ) << endl << Motifs.GetMotif( veciMotifs[ j ] ) << endl;

			for( k = 0; k < vecsThreads.size( ); ++k ) {
				vecsThreads[ k ].m_iOffset = k;
				vecsThreads[ k ].m_iStep = vecsThreads.size( );
				vecsThreads[ k ].m_pveciPCL2FASTA = &veciPCL2FASTA;
				vecsThreads[ k ].m_pGeneScores = &GeneScores;
				vecsThreads[ k ].m_pMotifs = &Motifs;
				vecsThreads[ k ].m_iMotif = veciMotifs[ j ];
				vecsThreads[ k ].m_pFASTA = &FASTA;
				vecsThreads[ k ].m_psModifiers = &sMods;
				if( pthread_create( &vecpthdThreads[ k ], NULL, ThreadCombineMotif, &vecsThreads[ k ] ) ) {
					cerr << "CCoalesceImpl::CombineMotifs( " << sArgs.threads_arg <<
						" ) could not combine motif: " << Motifs.GetMotif( veciMotifs[ j ] ).c_str( ) << endl;
					return 1; } }
			for( k = 0; k < vecpthdThreads.size( ); ++k )
				pthread_join( vecpthdThreads[ k ], NULL );
			for( k = 0; k < PCL.GetGenes( ); ++k )
				if( veciPCL2FASTA[ k ] != -1 ) {
					if( Cluster.IsGene( k ) )
						HistsCluster.Add( GeneScores, k, false, veciMotifs[ j ] );
					else
						HistsPot.Add( GeneScores, k, false, veciMotifs[ j ] ); } }
/*
		for( j = 0; j < veciPCL2FASTA.size( ); ++j )
			if( ( iFASTA = veciPCL2FASTA[ j ] ) != -1 ) {
				vector<SFASTASequence>	vecsSequences;
				float					dFive;

				if( FASTA.Get( iFASTA, vecsSequences ) ) {
					sModifiers.Get( j );
					dFive = 0;
					for( k = 0; k < vecsSequences.size( ); ++k ) {
						if( !GeneScores.Add( j, Motifs, vecsSequences[ k ],
							sModifiers, iMotifThree, vecdScores, veciLengths ) ) {
							cerr << "FAILURE!" << endl;
							return 1; }
						if( vecsSequences[ k ].m_strType == "5" )
							dFive = Motifs.GetMatch( vecsSequences[ k ].m_vecstrSequences[ 0 ],
								iMotifThree, 0, sModifiers ) / vecsSequences[ k ].m_vecstrSequences[ 0 ].length( ); }

					cerr << "Gene	" << j << '\t' << PCL.GetGene( j ) << endl;
					if( dFive != ( d = GeneScores.Get( GeneScores.GetType( "5" ), CCoalesceSequencerBase::ESubsequenceTotal, j, iMotifThree ) ) ) {
						cerr << "C	" << dFive << '\t' << d << endl;
						return 1; }
					for( k = 0; k < GeneScores.GetTypes( ); ++k )
						for( eSubsequence = CCoalesceSequencerBase::ESubsequenceBegin;
							(size_t)eSubsequence < GeneScores.GetSubsequences( k );
							eSubsequence = (CCoalesceSequencerBase::ESubsequence)( (size_t)eSubsequence + 1 ) ) {
							const float*	ad	= GeneScores.Get( k, eSubsequence, j );

							d = GeneScores.Get( k, eSubsequence, j, iMotifThree );
							if( ad ) {
								if( ad[ iMotifThree ] != d ) {
									cerr << "A	" << GeneScores.GetType( k ) << '\t' << eSubsequence << '\t' << ad[ iMotifThree ] <<
										'\t' << GeneScores.Get( k, eSubsequence, j, iMotifThree ) << endl;
									return 1; } }
							else if( d != 0 ) {
								cerr << "B	" << GeneScores.GetType( k ) << '\t' << eSubsequence << '\t' <<
									'\t' << GeneScores.Get( k, eSubsequence, j, iMotifThree ) << endl;
								return 1; } }
					if( Cluster.IsGene( j ) )
						HistsCluster.Add( GeneScores, j, false, iMotifThree );
					else
						HistsPot.Add( GeneScores, j, false, iMotifThree ); } }
//*/

/*
		for( j = 0; j < PCL.GetGenes( ); ++j )
			if( veciPCL2FASTA[ j ] != -1 ) {
				vector<SFASTASequence>	vecsSequences;

				FASTA.Get( veciPCL2FASTA[ j ], vecsSequences );
				for( k = 0; k < vecsSequences.size( ); ++k ) {
					const SFASTASequence&	sSequence	= vecsSequences[ k ];
					float*					ad;

					if( sSequence.m_vecstrSequences.size( ) != 1 )
						continue;
					d = Motifs.GetMatch( sSequence.m_vecstrSequences[ 0 ], iMotifThree, 0,
						sModifiers ) / sSequence.m_vecstrSequences[ 0 ].length( );
					if( !( ad = GeneScores.Get( GeneScores.GetType( sSequence.m_strType ),
						CCoalesceSequencerBase::ESubsequenceTotal, j ) ) || ( d != ad[ iMotifThree ] ) ) {
						cerr << j << '\t' << PCL.GetGene( j ) << '\t' << sSequence.m_strType << '\t' <<
							( ad ? ad[ iMotifThree ] : -1 ) << '\t' << d << endl;
						return 1; } } }
//*/

		for( j = 0; j < PCL.GetGenes( ); ++j )
			for( k = 0; k < GeneScores.GetTypes( ); ++k )
				for( eSubsequence = CCoalesceSequencerBase::ESubsequenceBegin;
					(size_t)eSubsequence < GeneScores.GetSubsequences( k );
					eSubsequence = (CCoalesceSequencerBase::ESubsequence)( (size_t)eSubsequence + 1 ) ) {
					if( !GeneScores.Get( k, eSubsequence, j ) )
						continue;
					for( size_t m = 0; m < veciMotifs.size( ); ++m ) {
						d = GeneScores.Get( k, eSubsequence, j, veciMotifs[ m ] );
						s = Cluster.IsGene( j ) ?
							HistsCluster.Get( k, eSubsequence ).Get( veciMotifs[ m ], d ) :
							HistsPot.Get( k, eSubsequence ).Get( veciMotifs[ m ], d );
						if( !s ) {
							cerr << "D	" << j << '\t' << PCL.GetGene( j ) << '\t' <<
								GeneScores.GetType( k ) << '\t' << eSubsequence << '\t' << d << '\t' <<
								Cluster.IsGene( j ) << endl;
							cerr << Motifs.GetMotif( veciMotifs[ m ] ) << endl;
							cerr << HistsCluster.Get( k, eSubsequence ).Save( veciMotifs[ m ] ) << endl;
							cerr << HistsPot.Get( k, eSubsequence ).Save( veciMotifs[ m ] ) << endl;
							return 1; } } }

		for( j = 0; j < PCL.GetGenes( ); ++j ) {
			if( ( (float)rand( ) / RAND_MAX ) > ( 1.0 / 100 ) )
				continue;
			if( Cluster.IsGene( j ) )
				Pot.Add( j, Cluster );
			else
				Cluster.Add( j, Pot ); } }

	return 0; }

int main_test2( const gengetopt_args_info& sArgs, const CPCL& PCL, const CFASTA& FASTA, CCoalesceMotifLibrary& Motifs ) {
	CCoalesceGeneScores			GeneScores;
	size_t						i, j, iIter, iOrig;
	uint32_t					iOne, iTwo, iThree;
	SCoalesceModifiers			sMods;
	vector<vector<float> >		vecvecdCounts, vecvecdOrig;
	vector<size_t>				veciLengths;
	vector<float>				vecdScores;
	SCoalesceModifierCache		sModifiers( sMods );
	vector<size_t>				veciPCL2FASTA;
	float						d;

	veciPCL2FASTA.resize( PCL.GetGenes( ) );
	for( i = 0; i < veciPCL2FASTA.size( ); ++i )
		veciPCL2FASTA[ i ] = FASTA.GetGene( PCL.GetGene( i ) );
	sMods.Initialize( PCL );

	GeneScores.SetGenes( PCL.GetGenes( ) );
	for( i = 0; i < veciPCL2FASTA.size( ); ++i )
		if( veciPCL2FASTA[ i ] != -1 ) {
			vector<SFASTASequence>	vecsSequences;

			if( FASTA.Get( veciPCL2FASTA[ i ], vecsSequences ) ) {
				sModifiers.Get( i );
				for( j = 0; j < vecsSequences.size( ); ++j )
					if( !GeneScores.Add( i, Motifs, vecsSequences[ j ], sModifiers, vecvecdCounts,
						veciLengths ) )
						return 1; } }

	vecvecdOrig.resize( PCL.GetGenes( ) );
	for( i = 0; i < vecvecdOrig.size( ); ++i ) {
		vecvecdOrig[ i ].resize( iOrig = Motifs.GetMotifs( ) );
		for( j = 0; j < vecvecdOrig[ i ].size( ); ++j )
			vecvecdOrig[ i ][ j ] = GeneScores.Get( 0, CCoalesceSequencerBase::ESubsequenceTotal, i,
				(uint32_t)j ); }

	for( iIter = 0; iIter < 1000; ++iIter ) {
		iOne = rand( ) % Motifs.GetMotifs( );
		iTwo = rand( ) % Motifs.GetMotifs( );
		iThree = Motifs.Merge( iOne, iTwo, FLT_MAX, false );

		for( i = 0; i < veciPCL2FASTA.size( ); ++i )
			if( veciPCL2FASTA[ i ] != -1 ) {
				vector<SFASTASequence>	vecsSequences;

				if( FASTA.Get( veciPCL2FASTA[ i ], vecsSequences ) ) {
					sModifiers.Get( i );
					for( j = 0; j < vecsSequences.size( ); ++j )
						if( !GeneScores.Add( i, Motifs, vecsSequences[ j ], sModifiers, iThree, vecdScores,
							veciLengths ) ) {
							cerr << "FAILURE" << endl;
							return 1; } } }

		for( i = 0; i < PCL.GetGenes( ); ++i )
			for( j = 0; j < vecvecdOrig[ i ].size( ); ++j )
				if( vecvecdOrig[ i ][ j ] != ( d = GeneScores.Get( 0,
					CCoalesceSequencerBase::ESubsequenceTotal, i, (uint32_t)j ) ) ) {
					cerr << i << '\t' << PCL.GetGene( i ) << '\t' << j << '\t' << vecvecdOrig[ i ][ j ] <<
						'\t' << d << endl;
					cerr << Motifs.GetMotif( (uint32_t)j ) << endl;
					return 1; } }

	return 0; }
