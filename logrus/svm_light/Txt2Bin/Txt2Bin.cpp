#include "stdafx.h"

int Bin2Txt( const char* );
int Txt2Bin( const char*, unsigned int, istream& );

int main( int iArgs, char** aszArgs ) {
	unsigned int	iWords;
	const char*		szFile;

	if( iArgs != 3 ) {
		cerr << "Usage: " << aszArgs[ 0 ] << " <words> <file>" << endl;
		return 1; }

	iWords = atoi( aszArgs[ 1 ] );
	szFile = aszArgs[ 2 ];

	return ( iWords ? Txt2Bin( szFile, iWords, cin ) : Bin2Txt( szFile ) ); }

int Txt2Bin( const char* szOutput, unsigned int iWords, istream& istm ) {
	char			szToken[ c_iLineLength ];
	ofstream		ofsmOutput;
	float			d;
	char*			pc;
	float*			adWords;
	string			strComment;
	size_t			i, iDocs;

	adWords = new float[ iWords ];
	ofsmOutput.open( szOutput, ios_base::out | ios_base::binary );
	iDocs = 0;
	ofsmOutput.write( (char*)&iWords, sizeof(iWords) );
	ofsmOutput.write( (char*)&iDocs, sizeof(iDocs) );
	istm >> szToken;
	while( istm.peek( ) != EOF ) {
		if( !szToken[ 0 ] )
			break;
		iDocs++;
		d = (float)atof( szToken );
		ofsmOutput.write( (char*)&d, sizeof(d) );
		strComment = "";
		memset( adWords, 0, sizeof(*adWords) * iWords );
		for( istm >> szToken; szToken[ 0 ]; istm >> szToken ) {
			if( szToken[ 0 ] == '#' ) {
				strComment = szToken + 1;
				istm.getline( szToken, c_iLineLength - 1 );
				strComment += szToken;
				continue; }
			if( !( pc = strchr( szToken, ':' ) ) )
				break;
			*pc = 0;
			adWords[ atoi( szToken ) - 1 ] = (float)atof( pc + 1 ); }
		ofsmOutput.write( (char*)adWords, sizeof(*adWords) * iWords );
		i = strComment.length( );
		ofsmOutput.write( (char*)&i, sizeof(i) );
		ofsmOutput.write( strComment.c_str( ), (streamsize)strComment.length( ) ); }
		
	delete[] adWords;

	ofsmOutput.seekp( sizeof(iWords) );
	ofsmOutput.write( (char*)&iDocs, sizeof(iDocs) );
	ofsmOutput.close( );

	return 0; }

int Bin2Txt( const char* szInput ) {
	char			szComment[ c_iLineLength ];
	ifstream		ifsmInput;
	unsigned int	iDoc, iWord, iWords, iDocs, i;
	float			d;
	float*			adWords;

	ifsmInput.open( szInput, ios_base::binary );
	ifsmInput.read( (char*)&iWords, sizeof(iWords) );
	ifsmInput.read( (char*)&iDocs, sizeof(iDocs) );
	adWords = new float[ iWords ];
	for( iDoc = 0; iDoc < iDocs; ++iDoc ) {
		ifsmInput.read( (char*)&d, sizeof(d) );
		if( d == 1 )
			cout << "+1";
		else if( d == -1 )
			cout << "-1";
		else
			cout << d;

		ifsmInput.read( (char*)adWords, sizeof(*adWords) * iWords );
		for( iWord = 0; iWord < iWords; ++iWord )
			if( adWords[ iWord ] )
				cout << '\t' << ( iWord + 1 ) << ':' << adWords[ iWord ];

		ifsmInput.read( (char*)&i, sizeof(i) );
		if( i ) {
			ifsmInput.read( szComment, i );
			szComment[ i ] = 0;
			cout << "\t#" << szComment; }

		cout << endl; }
	delete[] adWords;
	cout.flush( );

	return 0; }
