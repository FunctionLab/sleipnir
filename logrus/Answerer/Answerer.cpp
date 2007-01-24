#include "stdafx.h"
#include "cmdline.h"

int read_genes( const char*, CGenome&, vector<CGenes*>& );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CGenome				Genome;
	CDat				Dat;
	vector<CGenes*>		vecpPositives, vecpNegatives;
	size_t				i, j, iPositives, iNegatives;
	int					iRet;
	float				d, dProb;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta::Startup( sArgs.verbosity_arg, sArgs.random_arg );

	if( sArgs.genome_arg ) {
		ifstream	ifsm;

		ifsm.open( sArgs.genome_arg );
		if( !Genome.Open( ifsm ) ) {
			cerr << "Couldn't open: " << sArgs.genome_arg << endl;
			return 1; } }
	if( ( sArgs.positives_arg && ( iRet = read_genes( sArgs.positives_arg, Genome, vecpPositives ) ) ) ||
		( sArgs.negatives_arg && ( iRet = read_genes( sArgs.negatives_arg, Genome, vecpNegatives ) ) ) )
		return iRet;

	if( sArgs.input_arg ) {
		CDat	DatPositives;

		if( !DatPositives.Open( sArgs.input_arg, true ) ) {
			cerr << "Could not open: " << sArgs.input_arg << endl;
			return 1; }
		if( !Dat.Open( DatPositives, vecpNegatives, Genome, true ) ) {
			cerr << "Could not open " << sArgs.input_arg << " with negatives" << endl;
			return 1; } }
	else
		Dat.Open( vecpPositives, vecpNegatives, (float)sArgs.overlap_arg, Genome );
	if( sArgs.interactions_given ) {
		for( iPositives = i = 0; i < Dat.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
				if( !CMeta::IsNaN( d = Dat.Get( i, j ) ) && d )
					iPositives++;
		iNegatives = (size_t)( ( (float)( iPositives * ( Dat.GetGenes( ) - 1 ) ) /
			sArgs.interactions_arg ) + 0.5f );
		i = ( Dat.GetGenes( ) * ( Dat.GetGenes( ) - 1 ) / 2 ) - iPositives;
		if( iNegatives > i ) {
			cerr << "Impossibly large negative interaction requirement: " << iNegatives << endl;
			return 1; }
		dProb = (float)iNegatives / i;
		for( i = 0; i < Dat.GetGenes( ); ++i )
			for( j = ( i + 1 ); j < Dat.GetGenes( ); ++j )
				if( ( CMeta::IsNaN( d = Dat.Get( i, j ) ) || !d ) &&
					( ( rand( ) / RAND_MAX ) < dProb ) )
					Dat.Set( i, j, 0 ); }
	if( sArgs.test_arg ) {
		set<string>					setstrGenes;
		set<string>::const_iterator	iterGenes;
		ostream&					ostm	= sArgs.output_arg ? cout : cerr;

		for( i = 0; i < vecpPositives.size( ); ++i ) {
			size_t	iGenes;

			if( iGenes = vecpPositives[ i ]->GetGenes( ) )
				setstrGenes.insert( vecpPositives[ i ]->GetGene( rand( ) % iGenes ).GetName( ) ); }
		while( setstrGenes.size( ) < ( sArgs.test_arg * Dat.GetGenes( ) ) )
			setstrGenes.insert( Dat.GetGene( rand( ) % Dat.GetGenes( ) ) );
		for( iterGenes = setstrGenes.begin( ); iterGenes != setstrGenes.end( ); ++iterGenes )
			ostm << *iterGenes << endl; }

	Dat.Save( sArgs.output_arg );

	for( i = 0; i < vecpPositives.size( ); ++i )
		delete vecpPositives[ i ];
	for( i = 0; i < vecpNegatives.size( ); ++i )
		delete vecpNegatives[ i ];
	CMeta::Shutdown( );
	return 0; }

int read_genes( const char* szDir, CGenome& Genome, vector<CGenes*>& vecpGenes ) {
	CGenes*		pGenes;
	string		strDir, strFile;

	strDir = szDir;
#ifdef _MSC_VER
	bool			fOK;
	HANDLE			hFind;
	WIN32_FIND_DATA	sFind;
	for( fOK = ( hFind = FindFirstFile( ( strDir + "\\*" ).c_str( ), &sFind ) ) !=
		INVALID_HANDLE_VALUE; fOK; fOK = !!FindNextFile( hFind, &sFind ) ) {
		if( sFind.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY )
			continue;
		strFile = strDir + '/' + sFind.cFileName;
#else // _MSC_VER
	DIR*			pDir;
	struct dirent*	pFind;
	struct stat		sStat;
	pDir = opendir( szDir );
	while( pFind = readdir( pDir ) ) {
		strFile = strDir + '/' + pFind->d_name;
		stat( strFile.c_str( ), &sStat );
		if( S_ISDIR( sStat.st_mode ) )
			continue;
#endif // _MSC_VER
		ifstream	ifsm;

		ifsm.open( strFile.c_str( ) );
		pGenes = new CGenes( Genome );
		if( !pGenes->Open( ifsm ) ) {
			cerr << "Couldn't open: " << strFile << endl;
			return 1; }
		vecpGenes.push_back( pGenes ); }

	return 0; }
