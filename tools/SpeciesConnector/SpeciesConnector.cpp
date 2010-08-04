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

static const char	c_acDab[]	= ".dab";

struct SSorter {
	const vector<size_t>&	m_veciSizes;

	SSorter( const vector<size_t>& veciSizes ) : m_veciSizes(veciSizes) { }

	bool operator()( size_t iOne, size_t iTwo ) const {

		return ( m_veciSizes[ iOne ] > m_veciSizes[ iTwo ] ); }
};

float find_value( size_t iOne, size_t iTwo, const CDat& Dat ) {

	return ( ( ( iOne == -1 ) || ( iTwo == -1 ) ) ? CMeta::GetNaN( ) : Dat.Get( iOne, iTwo ) ); }

size_t sample ( vector<float> flag ){
	size_t						i;
	float						UnifSamp;
		
	UnifSamp = ( float )rand( ) / RAND_MAX; 
	for( i = 0; i < flag.size( ); ++i ){
		if( UnifSamp <= flag[ i ] )
			break;}

	return i;}

int main( int iArgs, char** aszArgs ) {
	size_t						NGibbs, burn;
	gengetopt_args_info			sArgs;
	CDataPair					DatOne, DatTwo;
	size_t						i, j, k, iDatOne, iDatTwo, iGeneOne, iGeneTwo, count1, count2;
	size_t						iCountJoint, iValueOne, iValueTwo;
	vector<size_t>				veciGenesOne, veciGenesTwo, veciDefaults, veciSizes;
	vector<vector<size_t> >		vecveciJoint;
	vector<string>				vecstrTInputs, vecstrInputs, vecstrxInputs, vecstrlInputs, vecstrlxInputs, vecstrGenes;
	float						dValueOne, dValueTwo;
	map<string, size_t>			mapZeros;
	CBayesNetSmile				BNIn;
				
	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg, sArgs.random_arg );

	if( sArgs.inputs_num%2 !=0 ){
		cerr << "Number of inputs should be even." <<  endl;
		return 1;}
	
	NGibbs = burn = sArgs.gibbs_arg;

	vecstrTInputs.resize( sArgs.inputs_num );
	copy( sArgs.inputs, sArgs.inputs + sArgs.inputs_num, vecstrTInputs.begin( ) );

	vecstrInputs.resize( sArgs.inputs_num/2 );
	vecstrxInputs.resize( sArgs.inputs_num/2 );
	vecstrlInputs.resize( sArgs.inputs_num/2 );
	vecstrlxInputs.resize( sArgs.inputs_num/2 );
	
	for( count1 = count2 = 0, i = 0; i < vecstrTInputs.size( ); ++i ){
		if( strcmp( &vecstrTInputs[ i ][ vecstrTInputs[ i ].rfind( "." ) + 1 ] , "dab" ) == 0 ){
			vecstrInputs[ count1 ] = vecstrTInputs[ i ];
			vecstrlInputs[ count1 ].resize( vecstrTInputs[ i ].rfind( "." ) - vecstrTInputs[ i ].rfind( "/" ) - 1 );
			vecstrlInputs[ count1 ] = vecstrTInputs[ i ].substr( vecstrTInputs[ i ].rfind( "/" ) + 1, vecstrTInputs[ i ].rfind( "." ) - vecstrTInputs[ i ].rfind( "/" ) - 1 ); 
			count1++;}
		else if( strcmp( &vecstrTInputs[ i ][ vecstrTInputs[ i ].rfind( "." ) + 1 ], "xdsl" ) == 0 ){ 
			vecstrxInputs[ count2 ] = vecstrTInputs[ i ];
			vecstrlxInputs[ count2 ].resize( vecstrTInputs[ i ].rfind( "." ) - vecstrTInputs[ i ].rfind( "/" ) - 1 );
			vecstrlxInputs[ count2 ] = vecstrTInputs[ i ].substr( vecstrTInputs[ i ].rfind( "/" ) + 1, vecstrTInputs[ i ].rfind( "." ) - vecstrTInputs[ i ].rfind( "/" ) - 1 ); 
			count2++;}
		else{
			cerr << "Input file types should be xdsl and dab." <<  endl;
			return 1;}}

	vector<size_t>				vecIndInputs;
	vecIndInputs.resize( vecstrlInputs.size( ) );
	for( i = 0; i < vecstrlInputs.size( ); ++i ){
		for( j = 0; j < vecstrlxInputs.size( ); ++j ){
			if( !vecstrlInputs[ i ].compare( vecstrlxInputs[ j ] ) ){
				vecIndInputs[ i ] = j;
				break;}}}

	{	vector<size_t>	veciIndices, veciSizes;

		veciSizes.resize( vecstrInputs.size( ) );
		veciIndices.resize( vecstrInputs.size( ) );
		for( i = 0; i < vecstrInputs.size( ); ++i ) {
			CDat		Dat;
			ifstream	ifsm;

			veciIndices[ i ] = i;
			ifsm.open( vecstrInputs[ i ].c_str( ) );
			Dat.OpenGenes( ifsm, true );
			veciSizes[ i ] = Dat.GetGenes( ); }
		sort( veciIndices.begin( ), veciIndices.end( ), SSorter( veciSizes ) );
		CMeta::Permute( vecstrInputs, veciIndices ); }
	
	{	CDataset	Data;

		Data.OpenGenes( vecstrInputs );

		vecstrGenes.resize( Data.GetGenes( ) );
		copy( Data.GetGeneNames( ).begin( ), Data.GetGeneNames( ).end( ), vecstrGenes.begin( ) );}

	veciSizes.resize( vecstrInputs.size( ) );
	for( i = 0; i < veciSizes.size( ); ++i ) {
		DatOne.OpenQuants( vecstrInputs[ i ].c_str( ) );
		veciSizes[ i ] = DatOne.GetValues( );}

	if( sArgs.zeros_arg ) {
		ifstream		ifsm;
		vector<string>	vecstrZeros;
		char			acLine[ 1024 ];

		ifsm.open( sArgs.zeros_arg );
		if( !ifsm.is_open( ) ) {
			cerr << "Could not open: " << sArgs.zeros_arg << endl;
			return 1; }
		while( !ifsm.eof( ) ) {
			ifsm.getline( acLine, ARRAYSIZE(acLine) - 1 );
			acLine[ ARRAYSIZE(acLine) - 1 ] = 0;
			vecstrZeros.clear( );
			CMeta::Tokenize( acLine, vecstrZeros );
			if( vecstrZeros.empty( ) )
				continue;
			mapZeros[ vecstrZeros[ 0 ] ] = atoi( vecstrZeros[ 1 ].c_str( ) ); } }

	vector<vector<float> >					vecPrior;
	vector<vector<vector<float> > >			vecDataGSpZero;
	vector<vector<vector<float> > >			vecDataGSpOne;
	vector<vector<string> >					vecvecSpDat;
	CDataMatrix		MatCPT;
	vector<string>	vecstrFiles;

	vecPrior.resize( sArgs.inputs_num/2 );
	vecDataGSpZero.resize( sArgs.inputs_num/2 );
	vecDataGSpOne.resize( sArgs.inputs_num/2 );
	vecvecSpDat.resize( sArgs.inputs_num/2 );
	for( i = 0; i < vecPrior.size(  ); ++i ){
		if( !BNIn.Open( vecstrxInputs[ vecIndInputs[ i ] ].c_str( ) ) ) {
			cerr << "Couldn't open: " << endl;
			return 1;}
		BNIn.GetCPT( 0, MatCPT );
		vecPrior[ i ].resize( MatCPT.GetRows( ) );
		vecPrior[ i ][ 0 ] = MatCPT.Get( 0, 0 );
		vecPrior[ i ][ 1 ] = MatCPT.Get( 1, 0 );

		BNIn.GetNodes( vecstrFiles );
		vecstrFiles.erase( vecstrFiles.begin( ) );
		vecDataGSpZero[ i ].resize( vecstrFiles.size( ) );
		vecDataGSpOne[ i ].resize( vecstrFiles.size( ) );
		vecvecSpDat[ i ].resize( vecstrFiles.size( ) );
		for( j = 0; j < vecstrFiles.size( ); ++j ){
			vecvecSpDat[ i ][ j ] = vecstrFiles[ j ];
			BNIn.GetCPT( j+1, MatCPT );
			vecDataGSpZero[ i ][ j ].resize( MatCPT.GetRows( ) );
			vecDataGSpOne[ i ][ j ].resize( MatCPT.GetRows( ) );
			for( k = 0; k < MatCPT.GetRows( ); ++k ){
				vecDataGSpZero[ i ][ j ][ k ] = MatCPT.Get( k , 0 );
				vecDataGSpOne[ i ][ j ][ k ] = MatCPT.Get( k , 1 );}}
		vecstrFiles.clear( );}

	vector<vector<vector<vector<float> > > >		vec4OSpGSp;
	vec4OSpGSp.resize( sArgs.inputs_num/2 );

	for( iDatOne = 0; iDatOne < vecstrInputs.size( ); ++iDatOne ) {
			vec4OSpGSp[ iDatOne ].resize( sArgs.inputs_num/2 );
			for( iDatTwo = ( iDatOne + 1 ); iDatTwo < vecstrInputs.size( ); ++iDatTwo ) {
				vec4OSpGSp[ iDatOne ][ iDatTwo ].resize( veciSizes[ iDatOne ] );
				for( i = 0; i < veciSizes[ iDatOne ]; ++i )
					vec4OSpGSp[ iDatOne ][ iDatTwo ][ i ].resize( veciSizes[ iDatTwo ] );}}
	
	for( iDatOne = 0; iDatOne < vecstrInputs.size( ); ++iDatOne ) {
			vec4OSpGSp[ iDatOne ].resize( sArgs.inputs_num/2 );
			for( iDatTwo = 0; iDatTwo < iDatOne; ++iDatTwo ) {
				vec4OSpGSp[ iDatOne ][ iDatTwo ].resize( veciSizes[ iDatOne ] );
				for( i = 0; i < veciSizes[ iDatOne ]; ++i )
					vec4OSpGSp[ iDatOne ][ iDatTwo ][ i ].resize( veciSizes[ iDatTwo ] );}}
	
	vector<float>		veciValueOnep;
	vector<float>		veciValueTwop;
	
	for( iDatOne = 0; iDatOne < vecstrInputs.size( ); ++iDatOne ) {
			if( !( DatOne.Open( vecstrInputs[ iDatOne ].c_str( ), false, !!sArgs.memmap_flag ) ||
				DatOne.Open( vecstrInputs[ iDatOne ].c_str( ), true, !!sArgs.memmap_flag ) ) ) {
					cerr << "Could not open: " << vecstrInputs[ iDatOne ] << endl;
					return 1; }

			if( sArgs.genex_arg && !DatOne.FilterGenes( sArgs.genex_arg, CDat::EFilterExclude ) ){
				cerr << "Couldn't open: " << sArgs.genex_arg << endl;
				return 1; }

			veciGenesOne.resize( vecstrGenes.size( ) );
			for( i = 0; i < vecstrGenes.size( ); ++i )
				veciGenesOne[ i ] = DatOne.GetGene( vecstrGenes[ i ] );

			if( ( iDatOne + 1 ) == vecstrInputs.size( ) )
				break;

			for( iDatTwo = ( iDatOne + 1 ); iDatTwo < vecstrInputs.size( ); ++iDatTwo ) {
				if( !( DatTwo.Open( vecstrInputs[ iDatTwo ].c_str( ), false, !!sArgs.memmap_flag ) ||
					DatTwo.Open( vecstrInputs[ iDatTwo ].c_str( ), true, !!sArgs.memmap_flag ) ) ) {
						cerr << "Could not open: " << vecstrInputs[ iDatTwo ] << endl;
						return 1; }

				if( sArgs.genex_arg && !DatTwo.FilterGenes( sArgs.genex_arg, CDat::EFilterExclude ) ){
					cerr << "Couldn't open: " << sArgs.genex_arg << endl;
					return 1; }

				veciGenesTwo.resize( vecstrGenes.size( ) );
				for( i = 0; i < veciGenesTwo.size( ); ++i )
					veciGenesTwo[ i ] = DatTwo.GetGene( vecstrGenes[ i ] );
				
				vecveciJoint.resize( veciSizes[ iDatOne ] );
				for( i = 0; i < vecveciJoint.size( ); ++i ) {
					vecveciJoint[ i ].resize( veciSizes[ iDatTwo ] );
					fill( vecveciJoint[ i ].begin( ), vecveciJoint[ i ].end( ), 0 ); }
				
				for( i = iCountJoint = 0; i < vecstrGenes.size( ); ++i ) {
					iGeneOne = veciGenesOne[ i ];
					iGeneTwo = veciGenesTwo[ i ];
					for( j = ( i + 1 ); j < vecstrGenes.size( ); ++j ) {
						dValueTwo = find_value( iGeneTwo, veciGenesTwo[ j ], DatTwo );
						dValueOne = find_value( iGeneOne, veciGenesOne[ j ], DatOne );
						if( ( !CMeta::IsNaN( dValueTwo = find_value( iGeneTwo, veciGenesTwo[ j ], DatTwo ) ) ) && 
							( !CMeta::IsNaN( dValueOne = find_value( iGeneOne, veciGenesOne[ j ], DatOne ) ) ) ){
								iValueTwo = DatTwo.Quantize( dValueTwo );
								iValueOne = DatOne.Quantize( dValueOne );
								iCountJoint++;
								vecveciJoint[ iValueOne ][ iValueTwo ]++;}}}

				size_t				Flag716 = 0;
				for( i = 0; i < vecveciJoint.size( ); ++i )
					for( j = 0; j < vecveciJoint[ i ].size( ); ++j )
						Flag716 += vecveciJoint[ i ][ j ];

				if( Flag716 < 10 ){
					cerr << "Total number of common edges between the answer files of " << vecstrlInputs[ iDatOne ].c_str( ) << " and " 
						<< vecstrlInputs[ iDatTwo ].c_str( ) << " should be much greater than 10" << endl;
					return 1;}

				Flag716 = Flag716/10;

				for( i = 0; i < vecveciJoint.size( ); ++i ){
					for( j = 0; j < vecveciJoint[ i ].size( ); ++j ){
						vecveciJoint[ i ][ j ] += Flag716;
						iCountJoint += Flag716;}}

				vector<vector<float> >		vecveciTJoint;
				vecveciTJoint.resize( vecveciJoint.size( ) );
				for( i = 0; i < vecveciTJoint.size( ); ++i ){
					vecveciTJoint[ i ].resize( vecveciJoint[ i ].size( ) );
					for( j = 0; j < vecveciTJoint.size( ); ++j ){ 
						vecveciTJoint[ i ][ j ] = (float)vecveciJoint[ i ][ j ] / iCountJoint;}}				
									 							
				veciValueOnep.resize( veciSizes[ iDatOne ] );
				veciValueTwop.resize( veciSizes[ iDatTwo ] );
				fill( veciValueOnep.begin( ), veciValueOnep.end( ), 0.0f );
				fill( veciValueTwop.begin( ), veciValueTwop.end( ), 0.0f );
				
				for( i = 0; i < vecveciJoint.size( ); ++i ){
					for( j = 0; j < vecveciJoint[ i ].size( ); ++j ){
						veciValueOnep[ i ] += vecveciTJoint[ i ][ j ];
						veciValueTwop[ j ] += vecveciTJoint[ i ][ j ];}} 

				for( i = 0; i < vecveciJoint.size( ); ++i ){
					for( j = 0; j < vecveciJoint[ i ].size(  ); ++j ){
						vec4OSpGSp[ iDatOne ][ iDatTwo ][ i ][ j ] = vecveciTJoint[ i ][ j ] / veciValueOnep[ i ];
						vec4OSpGSp[ iDatTwo ][ iDatOne ][ i ][ j ] = vecveciTJoint[ j ][ i ] / veciValueTwop[ i ];}}}}

	// Gibbs Magic
	vector<float>				vecMax;
	vecMax.resize( sArgs.inputs_num/2 );
	for( i = 0; i < vecMax.size( ); ++i ){
		float Max = 0;
		for( j = 0; j < vecDataGSpZero[ i ].size(  ); ++j ){
			float		temp1;
			float		temp2;
			temp1 = *max_element( vecDataGSpZero[ i ][ j ].begin( ), vecDataGSpZero[ i ][ j ].end( ) );
			temp2 = *max_element( vecDataGSpOne[ i ][ j ].begin( ), vecDataGSpOne[ i ][ j ].end( ) );
			Max += ( temp1 >= temp2 ) ? log ( temp1 ) : log( temp2 );}
				
		for( j = 0; j < sArgs.inputs_num/2; ++j ){
			if( j != i ){
				vector<float>	temp3;
				temp3.resize( vec4OSpGSp[ i ][ j ].size( ) );
				for( k = 0; k < vec4OSpGSp[ i ][ j ].size( ) ; ++k )
					temp3[ k ] = *max_element( vec4OSpGSp[ i ][ j ][ k ].begin( ), vec4OSpGSp[ i ][ j ][ k ].end( ) );
				Max += log ( *max_element( temp3.begin( ), temp3.end( ) ) );}}
		vecMax[ i ] = Max;}

	// Initialization
	vector<size_t>					vecTSamples;
	vecTSamples.resize( sArgs.inputs_num/2 );
	for( i = 0; i < vecTSamples.size( ); ++i  ){
		vector<float>			temp4;
		temp4.resize( vecPrior[ i ].size( ) );
		temp4[ 0 ] = vecPrior[ i ][ 0 ];
		for( j = 1; j < vecPrior[ i ].size( ); ++j )
			temp4[ j ] = temp4[ j - 1 ] + vecPrior[ i ][ j ];

		vecTSamples[ i ] = sample( temp4 );}
	
	// Sampling
	CDat		DatOut00, DatOut01, DatOut10; 
	DatOut00.Open( vecstrlInputs );
	DatOut01.Open( vecstrlInputs );
	DatOut10.Open( vecstrlInputs );
	
	for( i = 0; i < NGibbs; ++i ){ 
		vector<vector<size_t> >		vecvecDSp;
		vecvecDSp.resize( sArgs.inputs_num/2 );
		for( j = 0; j < sArgs.inputs_num/2; ++j ){
//			const vector<float>&	vecdCur	= ( vecTSamples[j] ? vecDataGSpOne : vecDataGSpZero );

			vecvecDSp[ j ].resize( vecDataGSpZero[ j ].size(  ) );
			if( !vecTSamples[ j ] ){
				for( k = 0; k < vecDataGSpZero[ j ].size( ); ++k ){
					vector<float>		temp5;
					size_t				l;
					temp5.resize( vecDataGSpZero[ j ][ k ].size( ) );
					temp5[ 0 ] = vecDataGSpZero[ j ][ k ][ 0 ];
					for( l = 1; l < vecDataGSpZero[ j ][ k ].size( ); ++l )
						temp5[ l ] = temp5[ l - 1 ] + vecDataGSpZero[ j ][ k ][ l ];
					vecvecDSp[ j ][ k ] = sample( temp5 );}}
			else{
				for( k = 0; k < vecDataGSpOne[ j ].size(  ); ++k ){
					vector<float>		temp5;
					size_t				l;
					temp5.resize( vecDataGSpOne[ j ][ k ].size( ) );
					temp5[ 0 ] = vecDataGSpOne[ j ][ k ][ 0 ];
					for( l = 1; l < vecDataGSpOne[ j ][ k ].size( ); ++l )
						temp5[ l ] = temp5[ l - 1 ] + vecDataGSpOne[ j ][ k ][ l ];
					vecvecDSp[ j ][ k ] = sample( temp5 );}}

			vector<float>			temp6;
			temp6.resize( vecPrior[ j ].size( ) );
			temp6[ 0 ] = vecPrior[ j ][ 0 ];
			for( k = 1; k < vecPrior[ j ].size( ); ++k )
				temp6[ k ] = temp6[ k - 1 ] + vecPrior[ j ][ k ];
			
			float		sum;
			do{ 
				vecTSamples[ j ] = sample( temp6 );
				sum = 0;
				for( k = 0; k < vecvecDSp[ j ].size( ); ++k ){
// same fix here for multi vectors
					if( !vecTSamples[ j ] )
						sum += log ( vecDataGSpZero[ j ][ k ][ vecvecDSp[ j ][ k ] ] );
					else
						sum += log ( vecDataGSpOne[ j ][ k ][ vecvecDSp[ j ][ k ] ] );}

				for( k = 0; k < sArgs.inputs_num/2; ++k )
					if( k != j )
						sum += log ( vec4OSpGSp[ j ][ k ][ vecTSamples[ j ] ][ vecTSamples[ k ] ] );
			} while ( log ( (float)rand( ) / RAND_MAX ) > ( sum - vecMax[ j ] ) );
		}
		 
		if( i >= ( NGibbs - burn ) ){
			for( j = 0; j < sArgs.inputs_num/2; ++j ){
				for( k = j+1; k < sArgs.inputs_num/2; ++k ){
					if( !vecTSamples[ j ] && !vecTSamples[ k ] ){
						if( CMeta::IsNaN( DatOut00.Get( j, k ) ) )
							DatOut00.Set( j, k, 1 );
						else
							DatOut00.Set( j, k, ++DatOut00.Get( j, k ) );}
					if( !vecTSamples[ j ] && vecTSamples[ k ] ){
						if( CMeta::IsNaN( DatOut01.Get( j, k ) ) )
							DatOut01.Set( j, k, 1 );
						else 
							DatOut01.Set( j, k, ++DatOut01.Get( j, k ) );}
					if( vecTSamples[ j ] && !vecTSamples[ k ] ){
						if( CMeta::IsNaN( DatOut10.Get( j, k ) ) )
							DatOut10.Set( j, k, 1 );
						else
							DatOut10.Set( j, k, ++DatOut10.Get( j, k ) );}}}}	
	}	

	_mkdir( sArgs.directory_arg );

	vector<vector<vector<vector<float> > > >		vec4FSpGSp;
	vec4FSpGSp.resize( sArgs.inputs_num/2 );
	
	for( iDatOne = 0; iDatOne < vecstrInputs.size( ); ++iDatOne ) {
			vec4FSpGSp[ iDatOne ].resize( sArgs.inputs_num/2 );
			for( iDatTwo = ( iDatOne + 1 ); iDatTwo < vecstrInputs.size( ); ++iDatTwo ) {
				vec4FSpGSp[ iDatOne ][ iDatTwo ].resize( veciSizes[ iDatOne ] );
				for( i = 0; i < veciSizes[ iDatOne ]; ++i )
					vec4FSpGSp[ iDatOne ][ iDatTwo ][ i ].resize( veciSizes[ iDatTwo ] );}}

	for( j = 0; j < sArgs.inputs_num/2; ++j ){
		for( k = j+1; k < sArgs.inputs_num/2; ++k ){
			double			DJValue;
			if( !CMeta::IsNaN( DatOut00.Get( j, k ) ) ){
				DJValue = DatOut00.Get( j, k );
				DatOut00.Set( j, k, ( (float) DJValue / burn ) );
				vec4FSpGSp[ j ][ k ][ 0 ][ 0 ] = DatOut00.Get( j, k );}
			if( !CMeta::IsNaN( DatOut01.Get( j, k ) ) ){
				DJValue = DatOut01.Get( j, k );
				DatOut01.Set( j, k, ( (float) DJValue / burn ) );
				vec4FSpGSp[ j ][ k ][ 0 ][ 1 ] = DatOut01.Get( j, k );}
			if( !CMeta::IsNaN( DatOut10.Get( j, k ) ) ){
				DJValue = DatOut10.Get( j, k );
				DatOut10.Set( j, k, ( (float) DJValue / burn ) );
				vec4FSpGSp[ j ][ k ][ 1 ][ 0 ] = DatOut10.Get( j, k );}

			float		First3Sum = vec4FSpGSp[ j ][ k ][ 0 ][ 0 ] + vec4FSpGSp[ j ][ k ][ 0 ][ 1 ] + vec4FSpGSp[ j ][ k ][ 1 ][ 0 ];

			vec4FSpGSp[ j ][ k ][ 1 ][ 1 ] = ( First3Sum == 1 ) ? 0.0f : ( 1.0f - First3Sum );}}
	
	DatOut00.Save( ( ( string ) sArgs.directory_arg + "/Learned00.dab" ).c_str( ) );
	DatOut01.Save( ( ( string ) sArgs.directory_arg + "/Learned01.dab" ).c_str( ) );
	DatOut10.Save( ( ( string ) sArgs.directory_arg + "/Learned10.dab" ).c_str( ) );
	
    for( iDatOne = 0; iDatOne < vecstrInputs.size( ); ++iDatOne ) {
        vec4FSpGSp[ iDatOne ].resize( sArgs.inputs_num/2 );
        for( iDatTwo = 0; iDatTwo < iDatOne; ++iDatTwo ) {
            vec4FSpGSp[ iDatOne ][ iDatTwo ].resize( veciSizes[ iDatOne ] );
            for( i = 0; i < veciSizes[ iDatOne ]; ++i )
                vec4FSpGSp[ iDatOne ][ iDatTwo ][ i ].resize( veciSizes[ iDatTwo ] );}}

    for( iDatOne = 0; iDatOne < vec4FSpGSp.size( ); ++iDatOne ) {
        for( iDatTwo = ( iDatOne + 1 ); iDatTwo < vec4FSpGSp[ iDatOne ].size( ); ++iDatTwo ) {
            vector<vector<float> >			vecvecFTJoint;
            float							MaxFJoint = 0;
            size_t							ZeroFCount = 0;
            vecvecFTJoint.resize( vec4FSpGSp[ iDatOne ][ iDatTwo ].size( ) );
            for( j = 0; j < vec4FSpGSp[ iDatOne ][ iDatTwo ].size( ); ++j ){
                vecvecFTJoint[ j ].resize( vec4FSpGSp[ iDatOne ][ iDatTwo ][ j ].size( ) );
                fill( vecvecFTJoint[ j ].begin( ), vecvecFTJoint[ j ].end( ), 0.0f );
                for( k = 0; k < vec4FSpGSp[ iDatOne ][ iDatTwo ][ j ].size( ); ++k ){
                    vecvecFTJoint[ j ][ k ] = vec4FSpGSp[ iDatOne ][ iDatTwo ][ j ][ k ];
                    if( !vecvecFTJoint[ j ][ k ] )
                        ZeroFCount++;
                    if( vecvecFTJoint[ j ][ k ] > MaxFJoint )
                        MaxFJoint = vecvecFTJoint[ j ][ k ];}}

            if( ZeroFCount ){
                float							BiasF;
                size_t							BiasFFlag = 0;

                BiasF = ( float ) MaxFJoint / ZeroFCount / veciSizes[ iDatTwo ] / veciSizes[ iDatOne ] / 100;

                for( i = 0; i < vecvecFTJoint.size( ); ++i )
                    for( j = 0; j < vecvecFTJoint[ i ].size( ); ++j )
                        if( !vecvecFTJoint[ i ][ j ] )
                            vecvecFTJoint[ i ][ j ] = BiasF;

                for( i = 0; i < vecvecFTJoint.size( ); ++i ){
                    for( j = 0; j < vecvecFTJoint[ i ].size( ); ++j ){
                        if( !( vecvecFTJoint[ i ][ j ] - MaxFJoint ) ){
                            vecvecFTJoint[ i ][ j ] -= ( BiasF * ZeroFCount );
                            BiasFFlag = 1;
                            break;}}
                    if( BiasFFlag )
                        break;}}
        
            vector<float>		vecFValueOnep;
            vector<float>		vecFValueTwop;

            vecFValueOnep.resize( veciSizes[ iDatOne ] );
            vecFValueTwop.resize( veciSizes[ iDatTwo ] );
            fill( vecFValueOnep.begin( ), vecFValueOnep.end( ), 0.0f );
            fill( vecFValueTwop.begin( ), vecFValueTwop.end( ), 0.0f );

            for( i = 0; i < vecvecFTJoint.size( ); ++i ){
                for( j = 0; j < vecvecFTJoint[ i ].size( ); ++j ){
                    vecFValueOnep[ i ] += vecvecFTJoint[ i ][ j ];
                    vecFValueTwop[ j ] += vecvecFTJoint[ i ][ j ];}} 

            for( i = 0; i < vecvecFTJoint.size( ); ++i ){
                for( j = 0; j < vecvecFTJoint[ i ].size(  ); ++j ){
                    vec4FSpGSp[ iDatOne ][ iDatTwo ][ i ][ j ] = vecvecFTJoint[ i ][ j ] / vecFValueOnep[ i ];
                    vec4FSpGSp[ iDatTwo ][ iDatOne ][ i ][ j ] = vecvecFTJoint[ j ][ i ] / vecFValueTwop[ i ];}}}}
    
    vector<vector<string> >		vecvecstrInputs;
    size_t						countstrInputs;
    vecvecstrInputs.resize( vecvecSpDat.size( ) );
    countstrInputs = 0;
    for( i = 0; i < vecvecSpDat.size( ); ++i ){
        vecvecstrInputs[ i ].resize( vecvecSpDat[ i ].size( ) ); 
        for( j = 0; j < vecvecSpDat[ i ].size( ); ++j ){
            vecvecstrInputs[ i ][ j ] =  ( string )sArgs.ndirectory_arg + '/' + vecstrlInputs[ i ] + '/' + vecvecSpDat[ i ][ j ] + c_acDab;
            countstrInputs++;}}

    vector<string>				vecstrFInputs;
    vector<vector<size_t> >		InputMaps;				
    vecstrFInputs.resize( countstrInputs );
    countstrInputs = 0;
    InputMaps.resize( vecvecstrInputs.size( ) );
    for( i = 0; i < vecvecstrInputs.size( ); ++i ){
        InputMaps[ i ].resize( vecvecstrInputs[ i ].size( ) );
        for( j = 0; j < vecvecstrInputs[ i ].size( ); ++j ){
            vecstrFInputs[ countstrInputs ] = vecvecstrInputs[ i ][ j ].c_str( );
            InputMaps[ i ][ j ] = countstrInputs;
            countstrInputs++;}}
    
    vector<string>		vecstrFGenes;
    
    {	CDataset	DataF;

        DataF.OpenGenes( vecstrFInputs );
        vecstrFGenes.resize( DataF.GetGenes( ) );
        copy( DataF.GetGeneNames( ).begin( ), DataF.GetGeneNames( ).end( ), vecstrFGenes.begin( ) );}

    vector<CDataPair*>			DataF;
    size_t						iDatOneF, iDatTwoF, iGeneOneF, iGeneTwoF;
    
    DataF.resize( vecstrFInputs.size( ) );
    for( i = 0; i < vecstrFInputs.size( ); ++i ){
        DataF[ i ] = new CDataPair( );
        if( !( DataF[ i ]->Open( vecstrFInputs[ i ].c_str( ), false, !!sArgs.memmap_flag ) ||
            DataF[ i ]->Open( vecstrFInputs[ i ].c_str( ), true, !!sArgs.memmap_flag ) ) ){
                cerr << "Could not open:" << vecstrFInputs[ i ] << endl;
                return 1;}}
        
    vector<CDat*>			DatOutCS, DatOutCSH, DatOutB;	

    DatOutCS.resize( vecvecstrInputs.size( ) );
    for( i = 0; i < DatOutCS.size( ); ++i ){
        DatOutCS[ i ] = new CDat( );
        DatOutCS[ i ]->Open( vecstrFGenes );}

    DatOutCSH.resize( vecvecstrInputs.size( ) );
    for( i = 0; i < DatOutCSH.size( ); ++i ){
        DatOutCSH[ i ] = new CDat( );
        DatOutCSH[ i ]->Open( vecstrFGenes );}

	DatOutB.resize( vecvecstrInputs.size( ) );
    for( i = 0; i < DatOutB.size( ); ++i ){
        DatOutB[ i ] = new CDat( );
        DatOutB[ i ]->Open( vecstrFGenes );}

    vector<vector<vector<size_t> > >		vec3GeneIndex;
    vec3GeneIndex.resize( vecvecstrInputs.size( ) );
    for( i = 0; i < vecvecstrInputs.size( ); ++i ){
        vec3GeneIndex[ i ].resize( vecvecstrInputs[ i ].size( ) );
        for( j = 0; j < vecvecstrInputs[ i ].size( ); ++j ){
            vec3GeneIndex[ i ][ j ].resize( vecstrFGenes.size( ) );
            for( k = 0; k < vecstrFGenes.size( ); ++k ){
                vec3GeneIndex[ i ][ j ][ k ] = DataF[ InputMaps[ i ][ j ] ]->GetGene( vecstrFGenes[ k ] );}}}
            
    for( iDatOneF = 0; iDatOneF < vecstrFGenes.size( ); ++iDatOneF ){ 
        for( iDatTwoF = ( iDatOneF + 1 ); iDatTwoF < vecstrFGenes.size( ); ++iDatTwoF ){ 
            vector<float>		vecDataIntZero;
            vector<float>		vecDataIntOne;
            vecDataIntZero.resize( vecvecstrInputs.size( ) );
            fill( vecDataIntZero.begin( ), vecDataIntZero.end( ), 0.0f );
            vecDataIntOne.resize( vecvecstrInputs.size( ) );
            fill( vecDataIntOne.begin( ), vecDataIntOne.end( ), 0.0f );
            for( i = 0; i < vecvecstrInputs.size( ); ++i ){
                float		sumZero = 0;
                float		sumOne = 0;
                for( j = 0; j < vecvecstrInputs[ i ].size( ); ++j ){
                    iGeneOneF = vec3GeneIndex[ i ][ j ][ iDatOneF ];
                    iGeneTwoF = vec3GeneIndex[ i ][ j ][ iDatTwoF ];
                    float		DFValue;
                    DFValue = ( ( iGeneOneF == -1 ) || ( iGeneTwoF == -1 ) ) ? CMeta::GetNaN( ) : DataF[ InputMaps[ i ][ j ] ]->Get( iGeneOneF, iGeneTwoF );
                    if( !CMeta::IsNaN( DFValue ) ){
                        size_t		DatPFI = DataF[ InputMaps[ i ][ j ] ]->Quantize( DFValue );
                        sumZero += log( vecDataGSpZero[ i ][ j ][ DatPFI ] );
                        sumOne += log( vecDataGSpOne[ i ][ j ][ DatPFI ] );}}
                vecDataIntZero[ i ] = sumZero;
                vecDataIntOne[ i ] = sumOne;}

            for( i = 0; i < vecvecstrInputs.size( ); ++i ){
                float		sumOne = 0, sumZero = 0, Final = 0 ;
				float		sumOneB = 0, sumZeroB = 0, FinalB = 0 ;
				float		sumOneT = 0, sumZeroT = 0, FinalT = 0;
                for( k = 0; k < vecvecstrInputs.size( ); ++k ){
                    if( k != i ){
                        sumOne += log( exp( log( vec4FSpGSp[ i ][ k ][ 1 ][ 0 ] ) + vecDataIntZero[ k ] ) + exp( log( vec4FSpGSp[ i ][ k ][ 1 ][ 1 ] ) + vecDataIntOne[ k ] ) );
                        sumZero += log( exp( log( vec4FSpGSp[ i ][ k ][ 0 ][ 0 ] ) + vecDataIntZero[ k ] ) + exp( log( vec4FSpGSp[ i ][ k ][ 0 ][ 1 ] ) + vecDataIntOne[ k ] ) );}}
                sumOneT = sumOne + log( vecPrior[ i ][ 1 ] );
				sumOne += ( vecDataIntOne[ i ] + log( vecPrior[ i ][ 1 ] ) );
				sumOneB =  vecDataIntOne[ i ] + log( vecPrior[ i ][ 1 ] );
                sumZeroT = sumZero + log( vecPrior[ i ][ 0 ] );
				sumZero += ( vecDataIntZero[ i ] + log( vecPrior[ i ][ 0 ] ) );
				sumZeroB =  vecDataIntZero[ i ] + log( vecPrior[ i ][ 0 ] );
                Final = (float) ( 1 / ( 1 + exp ( sumZero - sumOne ) ) );
				FinalT = (float) ( 1 / ( 1 + exp ( sumZeroT - sumOneT ) ) );
				FinalB = (float) ( 1 / ( 1 + exp ( sumZeroB - sumOneB ) ) );
				DatOutCS[ i ]->Set( iDatOneF, iDatTwoF, Final );
				DatOutCSH[ i ]->Set( iDatOneF, iDatTwoF, FinalT );
				DatOutB[ i ]->Set( iDatOneF, iDatTwoF, FinalB );}}}
                
    for( i = 0; i < vecvecstrInputs.size( ); ++i ){
        DatOutCS[ i ]->Save( ( ( string ) sArgs.directory_arg + '/' + vecstrlInputs[ i ] + c_acDab ).c_str( ) );
		DatOutCSH[ i ]->Save( ( ( string ) sArgs.directory_arg + '/' + vecstrlInputs[ i ] + 'h' + c_acDab ).c_str( ) );
		DatOutB[ i ]->Save( ( ( string ) sArgs.directory_arg + '/' + vecstrlInputs[ i ] + 'b' + c_acDab ).c_str( ) );
        delete DataF[ i ];
    delete DatOutCS[ i ];
	delete DatOutCSH[ i ];
	delete DatOutB[ i ];}	
	
	return 0;}
