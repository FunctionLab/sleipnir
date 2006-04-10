#include "stdafx.h"
#include "cmdline.h"

int read_genes( const char*, CGenome&, vector<CGenes*>& );

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info	sArgs;
	CGenome				Genome;
	CDat				Dat;
	vector<CGenes*>		vecpPositives, vecpNegatives;
	size_t				i;
	int					iRet;

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
	if( ( iRet = read_genes( sArgs.positives_arg, Genome, vecpPositives ) ) ||
		( iRet = read_genes( sArgs.negatives_arg, Genome, vecpNegatives ) ) )
		return iRet;

	Dat.Open( vecpPositives, vecpNegatives, Genome );
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

	if( sArgs.output_arg ) {
		ofstream	ofsm;

		ofsm.open( sArgs.output_arg );
		Dat.Save( ofsm, true ); }
	else
		Dat.Save( cout, false );

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
