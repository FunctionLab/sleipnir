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

float find_value( size_t iOne, size_t iTwo, const CDat& Dat ) {

	return ( ( ( iOne == -1 ) || ( iTwo == -1 ) ) ? CMeta::GetNaN( ) : Dat.Get( iOne, iTwo ) ); }

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info			sArgs;
	CDataPair					DatOne, DatTwo;
	size_t						i, j, k, iDatOne, iDatTwo, iGeneOne, iGeneTwo, count1, count2;
	size_t						iCountJoint, iValueOne, iValueTwo;
	vector<size_t>				veciGenesOne, veciGenesTwo, veciSizes;
	vector<vector<size_t> >		vecveciJoint;
	vector<string>				vecstrTInputs, vecstrInputs, vecstrxInputs, vecstrlInputs, vecstrlxInputs, vecstrGenes, vecstrAllGenes;
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

	{	CDataset	Data;

		Data.OpenGenes( vecstrInputs );

		vecstrAllGenes.resize( Data.GetGenes( ) );
		copy( Data.GetGeneNames( ).begin( ), Data.GetGeneNames( ).end( ), vecstrAllGenes.begin( ) );}

	if( sArgs.genex_arg ) {

		CGenome				Genome;		
		CGenes				GenesEx( Genome );

		if( !GenesEx.Open( sArgs.genex_arg ) ) {
			cerr << "Could not open: " << sArgs.genex_arg << endl;
			return 1; } 

		vecstrGenes.resize( vecstrAllGenes.size( ) - GenesEx.GetGenes( ) );

		size_t		GenesExCount = 0;
		for( i = 0; i < vecstrAllGenes.size( ); ++i ){
			if( !GenesEx.IsGene( vecstrAllGenes[ i ] ) ){
				vecstrGenes[ GenesExCount ] = vecstrAllGenes[ i ];
				GenesExCount++;}}}
	else{
		vecstrGenes.resize( vecstrAllGenes.size( ) );
		for( i = 0; i < vecstrAllGenes.size( ); ++i )
			vecstrGenes[ i ] = vecstrAllGenes[ i ];}
	
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
			cerr << "Couldn't open: " << vecstrxInputs[ vecIndInputs[ i ] ] << endl;
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
						if( ( !CMeta::IsNaN( dValueTwo ) ) && ( !CMeta::IsNaN( dValueOne ) ) ){
								iValueTwo = DatTwo.Quantize( dValueTwo );
								iValueOne = DatOne.Quantize( dValueOne );
								iCountJoint++;
								vecveciJoint[ iValueOne ][ iValueTwo ]++;}
						if( ( CMeta::IsNaN( dValueTwo ) ) && ( !CMeta::IsNaN( dValueOne ) ) ){
								iValueTwo = 0;
								iValueOne = DatOne.Quantize( dValueOne );
								iCountJoint++;
								vecveciJoint[ iValueOne ][ iValueTwo ]++;}
						if( ( !CMeta::IsNaN( dValueTwo ) ) && ( CMeta::IsNaN( dValueOne ) ) ){
								iValueTwo = DatTwo.Quantize( dValueTwo );
								iValueOne = 0;
								iCountJoint++;
								vecveciJoint[ iValueOne ][ iValueTwo ]++;}
						if( ( CMeta::IsNaN( dValueTwo ) ) && ( CMeta::IsNaN( dValueOne ) ) ){
								iValueTwo = 0;
								iValueOne = 0;
								iCountJoint++;
								vecveciJoint[ iValueOne ][ iValueTwo ]++;}}}

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

	_mkdir( sArgs.directory_arg );
    
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
                        sumOne += log( exp( log( vec4OSpGSp[ i ][ k ][ 1 ][ 0 ] ) + vecDataIntZero[ k ] ) + exp( log( vec4OSpGSp[ i ][ k ][ 1 ][ 1 ] ) + vecDataIntOne[ k ] ) );
                        sumZero += log( exp( log( vec4OSpGSp[ i ][ k ][ 0 ][ 0 ] ) + vecDataIntZero[ k ] ) + exp( log( vec4OSpGSp[ i ][ k ][ 0 ][ 1 ] ) + vecDataIntOne[ k ] ) );}}
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