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

static const char	c_acAlphabetExons[]	= "ACGT";
static const char	c_acAlphabetTotal[]	= "ACGTacgt";

struct STF {
	void Initialize( double dAverage, double dStdev, size_t iGenes, double dPGene, size_t iConditions,
		double dPCondition, size_t iMin, size_t iMax ) {
		size_t	i, iLength;

		m_strBS = "";
		iLength = ( rand( ) % ( iMax - iMin ) ) + iMin;
		for( i = 0; i < iLength; ++i )
			m_strBS += c_acAlphabetExons[ rand( ) % ( ARRAYSIZE(c_acAlphabetExons) - 1 ) ];
		m_dActivity = (float)( ( CStatistics::SampleNormalStandard( ) * dStdev ) + dAverage );
		if( ( (float)rand( ) / RAND_MAX ) < 0.5 )
			m_dActivity *= -1;

		m_veciGenes.reserve( (size_t)( iGenes * dPGene + 1 ) );
		for( i = 0; i < iGenes; ++i )
			if( ( (float)rand( ) / RAND_MAX ) < dPGene )
				m_veciGenes.push_back( i );
		m_veciConditions.reserve( (size_t)( iConditions * dPCondition + 1 ) );
		for( i = 0; i < iConditions; ++i )
			if( ( (float)rand( ) / RAND_MAX ) < dPCondition )
				m_veciConditions.push_back( i ); }

	void Save( ostream& ostm ) const {
		size_t	i;

		ostm << m_strBS << endl;
		ostm << m_dActivity << endl;
		for( i = 0; i < m_veciGenes.size( ); ++i )
			ostm << ( i ? "\t" : "" ) << m_veciGenes[ i ];
		ostm << endl;
		for( i = 0; i < m_veciConditions.size( ); ++i )
			ostm << ( i ? "\t" : "" ) << m_veciConditions[ i ];
		ostm << endl; }

	string			m_strBS;
	float			m_dActivity;
	vector<size_t>	m_veciGenes;
	vector<size_t>	m_veciConditions;
};

int main( int iArgs, char** aszArgs ) {
	static const size_t					c_iBuffer	= 16;
	gengetopt_args_info					sArgs;
	CPCL								PCL;
	vector<string>						vecstrGenes, vecstrExperiments, vecstrFeatures;
	size_t								i, j, k, iCopies;
	char								acBuffer[ c_iBuffer ];
	vector<STF>							vecsTFs;
	vector<CHMM>						vecHMMs;
	CFASTA								FASTA;
	ofstream							ofsm;
	map<string, size_t>					mapstriTypes;
	map<string, size_t>::const_iterator	iterType;
	string								strSequence;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg, sArgs.random_arg );

	vecsTFs.resize( sArgs.tfs_arg );
	for( i = 0; i < vecsTFs.size( ); ++i ) {
		vecsTFs[ i ].Initialize( sArgs.tf_mean_arg, sArgs.tf_stdev_arg, sArgs.genes_arg, sArgs.tf_gene_arg,
			sArgs.conditions_arg, sArgs.tf_condition_arg, sArgs.tf_min_arg, sArgs.tf_max_arg );
		vecsTFs[ i ].Save( cout ); }

	vecstrGenes.resize( sArgs.genes_arg );
	for( i = 0; i < vecstrGenes.size( ); ++i ) {
		sprintf_s( acBuffer, "G%07d", i );
		vecstrGenes[ i ] = acBuffer; }
	vecstrExperiments.resize( sArgs.conditions_arg );
	for( i = 0; i < vecstrExperiments.size( ); ++i ) {
		sprintf_s( acBuffer, "E%07d", i );
		vecstrExperiments[ i ] = acBuffer; }
	vecstrFeatures.push_back( "GID" );
	vecstrFeatures.push_back( "TFS" );
	vecstrFeatures.push_back( "GWEIGHT" );
	PCL.Open( vecstrGenes, vecstrExperiments, vecstrFeatures );

	for( i = 0; i < PCL.GetGenes( ); ++i ) {
		PCL.SetFeature( i, 2, "1" );
		for( j = 0; j < PCL.GetExperiments( ); ++j )
			PCL.Set( i, j, (float)( ( CStatistics::SampleNormalStandard( ) * sArgs.stdev_arg ) +
				sArgs.mean_arg ) ); }
	for( i = 0; i < vecsTFs.size( ); ++i )
		for( j = 0; j < vecsTFs[ i ].m_veciGenes.size( ); ++j ) {
			size_t			iGene	= vecsTFs[ i ].m_veciGenes[ j ];
			ostringstream	ossm;

			ossm << PCL.GetFeature( iGene, 1 );
			if( !ossm.str( ).empty( ) )
				ossm << ", ";
			ossm << i << '@' << vecsTFs[ i ].m_dActivity;
			PCL.SetFeature( iGene, 1, ossm.str( ) );
			for( k = 0; k < vecsTFs[ i ].m_veciConditions.size( ); ++k )
				PCL.Get( iGene, vecsTFs[ i ].m_veciConditions[ k ] ) +=
					(float)( ( CStatistics::SampleNormalStandard( ) * sArgs.stdev_arg ) +
					vecsTFs[ i ].m_dActivity ); }
	if( sArgs.output_pcl_arg ) {
		ofsm.open( sArgs.output_pcl_arg );
		if( !ofsm.is_open( ) ) {
			cerr << "Could not open: " << sArgs.output_pcl_arg << endl;
			return 1; }
		PCL.Save( ofsm );
		ofsm.close( ); }

	if( sArgs.fasta_arg ) {
		set<string>::const_iterator	iterType;

		if( !FASTA.Open( sArgs.fasta_arg ) ) {
			cerr << "Could not open: " << sArgs.fasta_arg << endl;
			return 1; }
		vecHMMs.resize( FASTA.GetTypes( ).size( ) );
		for( i = 0; i < vecHMMs.size( ); ++i )
			vecHMMs[ i ].Open( sArgs.degree_arg, c_acAlphabetTotal );
		for( iterType = FASTA.GetTypes( ).begin( ); iterType != FASTA.GetTypes( ).end( ); ++iterType )
			mapstriTypes[ *iterType ] = mapstriTypes.size( );
		for( i = 0; i < FASTA.GetGenes( ); ++i ) {
			vector<SFASTASequence>	vecsSequences;

			if( FASTA.Get( i, vecsSequences ) )
				for( j = 0; j < vecsSequences.size( ); ++j ) {
					const SFASTASequence&	sSequence	= vecsSequences[ j ];
					string					strCur, strSequence;

					for( k = 0; k < sSequence.m_vecstrSequences.size( ); ++k ) {
						strCur = sSequence.m_vecstrSequences[ k ];
						if( sSequence.m_fIntronFirst && !( k % 2 ) )
							transform( strCur.begin( ), strCur.end( ), strCur.begin( ), ::tolower );
						strSequence += strCur; }
					if( !vecHMMs[ mapstriTypes[ sSequence.m_strType ] ].Add( strSequence ) ) {
						cerr << "Could not add sequence: " << strSequence << endl;
						return 1; } } } }
	else {
		vecHMMs.resize( 1 );
		vecHMMs[ 0 ].Open( sArgs.degree_arg, c_acAlphabetExons );
		vecHMMs[ 0 ].SetUniform( );
		mapstriTypes[ "" ] = 0; }

	if( sArgs.output_fasta_arg ) {
		vector<string>			vecstrTypes;
		vector<vector<size_t> >	vecveciTFs;

		CMeta::Tokenize( sArgs.tf_types_arg, vecstrTypes, "," );
		vecveciTFs.resize( PCL.GetGenes( ) );
		for( i = 0; i < vecsTFs.size( ); ++i )
			for( j = 0; j < vecsTFs[ i ].m_veciGenes.size( ); ++j )
				vecveciTFs[ vecsTFs[ i ].m_veciGenes[ j ] ].push_back( i );
		ofsm.clear( );
		ofsm.open( sArgs.output_fasta_arg );
		if( !ofsm.is_open( ) ) {
			cerr << "Could not open: " << sArgs.output_fasta_arg << endl;
			return 1; }
		for( iterType = mapstriTypes.begin( ); iterType != mapstriTypes.end( ); ++iterType ) {
			const CHMM&	HMM			= vecHMMs[ iterType->second ];
			string		strHeader;
			bool		fTFs;

			fTFs = ( find( vecstrTypes.begin( ), vecstrTypes.end( ), iterType->first ) !=
				vecstrTypes.end( ) );
			if( !iterType->first.empty( ) )
				strHeader = strHeader + '\t' + iterType->first;
			for( i = 0; i < PCL.GetGenes( ); ++i ) {
				ofsm << "> " << PCL.GetGene( i ) << strHeader << endl;
				j = ( sArgs.seq_max_arg <= sArgs.seq_min_arg ) ? sArgs.seq_min_arg :
					( ( rand( ) % ( sArgs.seq_max_arg - sArgs.seq_min_arg ) ) + sArgs.seq_min_arg );
				strSequence = HMM.Get( j );
				if( fTFs )
					for( j = 0; j < vecveciTFs[ i ].size( ); ++j ) {
						const STF&	sTF	= vecsTFs[ vecveciTFs[ i ][ j ] ];

						iCopies = ( rand( ) % ( sArgs.tf_copx_arg - sArgs.tf_copm_arg ) ) + sArgs.tf_copm_arg;
						for( k = 0; k < iCopies; ++k ) {
							string	strBS	= sTF.m_strBS;
							size_t	iTarget	= rand( ) % ( strSequence.length( ) - strBS.length( ) );

							if( islower( strSequence[ iTarget ] ) )
								transform( strBS.begin( ), strBS.end( ), strBS.begin( ), ::tolower );
							strSequence.replace( iTarget, strBS.length( ), strBS ); } }
				for( j = 0; j < strSequence.length( ); j += sArgs.wrap_arg )
					ofsm << strSequence.substr( j, sArgs.wrap_arg ) << endl; } }
		ofsm.close( ); }

	return 0; }
