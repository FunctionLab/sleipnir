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

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info			sArgs;
	size_t						i, j, k, iDatOne, iDatTwo, iGeneOne, iGeneTwo, iZero;
	size_t						iCountJoint, iValueOne, iValueTwo;
	vector<vector<size_t> >		vecveciJoint;
	vector<string>				vecstrInputs, vecstrxInputs, vecstrGenes;
	map<string, size_t>			mapZeros;
	CBayesNetSmile				BNIn;
	vector<size_t>				veciGenesOneI, veciGenesTwoI;
	float						Threshold;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }
	CMeta Meta( sArgs.verbosity_arg, sArgs.random_arg );	

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
	
	_mkdir( sArgs.odirectory_arg );

	vecstrInputs.resize( sArgs.inputs_num );
	copy( sArgs.inputs, sArgs.inputs + sArgs.inputs_num, vecstrInputs.begin( ) );

	vecstrxInputs.resize( sArgs.inputs_num );

	for( i = 0; i < vecstrInputs.size( ); ++i ){
		if( strcmp( &vecstrInputs[ i ][ vecstrInputs[ i ].rfind( "." ) + 1 ], "xdsl" ) == 0 ){ 
			vecstrxInputs[ i ].resize( vecstrInputs[ i ].rfind( "." ) - vecstrInputs[ i ].rfind( "/" ) - 1 );
			vecstrxInputs[ i ] = vecstrInputs[ i ].substr( vecstrInputs[ i ].rfind( "/" ) + 1, vecstrInputs[ i ].rfind( "." ) - vecstrInputs[ i ].rfind( "/" ) - 1 );}
		else{
			cerr << "inputs file types should be xdsl." <<  endl;
			return 1;}}

	vector<vector<float> >					vecPrior;
	vector<vector<vector<float> > >			vecDataGSpZero;
	vector<vector<vector<float> > >			vecDataGSpOne;
	vector<vector<string> >					vecvecSpDat;
	CDataMatrix		MatCPT;
	vector<string>	vecstrFiles;

	vecPrior.resize( sArgs.inputs_num );
	vecDataGSpZero.resize( sArgs.inputs_num );
	vecDataGSpOne.resize( sArgs.inputs_num );
	vecvecSpDat.resize( sArgs.inputs_num );
	for( i = 0; i < vecPrior.size(  ); ++i ){
		if( !BNIn.Open( vecstrInputs[ i ].c_str( ) ) ) {
			cerr << "Couldn't open: " << vecstrInputs[ i ] << endl;
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

	if( sArgs.learn_flag ){

		_mkdir( sArgs.l0directory_arg );
		_mkdir( sArgs.l1directory_arg );
		_mkdir( sArgs.jdirectory_arg );

		Threshold = ( sArgs.threshold_arg == -1 ) ? 0.5f : float( sArgs.threshold_arg );
		
		vector<vector<string> >		vecvecstrInputs;
		size_t						countstrInputs;
		vecvecstrInputs.resize( vecvecSpDat.size( ) );
		countstrInputs = 0;
		for( i = 0; i < vecvecSpDat.size( ); ++i ){
			vecvecstrInputs[ i ].resize( vecvecSpDat[ i ].size( ) ); 
			for( j = 0; j < vecvecSpDat[ i ].size( ); ++j ){
				vecvecstrInputs[ i ][ j ] =  ( string )sArgs.ddirectory_arg + '/' + vecstrxInputs[ i ] + '/' + vecvecSpDat[ i ][ j ] + c_acDab;
				countstrInputs++;}}
		
		vector<string>				vecstrFDAInputs;
		size_t						countstrFDAInputs;
		vecstrFDAInputs.resize( countstrInputs + vecvecstrInputs.size( ) );
		countstrFDAInputs = 0;
		for( i = 0; i < vecvecstrInputs.size( ); ++i ){
			for( j = 0; j < vecvecstrInputs[ i ].size( ); ++j ){
				vecstrFDAInputs[ countstrFDAInputs ] = vecvecstrInputs[ i ][ j ].c_str( );
				countstrFDAInputs++;}
			vecstrFDAInputs[ countstrFDAInputs ] = ( string )sArgs.adirectory_arg + '/' + vecstrxInputs[ i ] + c_acDab;
			countstrFDAInputs++;}
		
		vector<string>		vecstrFGenes;
		{	CDataset	DataF;

		DataF.OpenGenes( vecstrFDAInputs );
		vecstrFGenes.resize( DataF.GetGenes( ) );
		copy( DataF.GetGeneNames( ).begin( ), DataF.GetGeneNames( ).end( ), vecstrFGenes.begin( ) );}
		
		if( sArgs.genelist_flag ){
			for( i = 0; i < vecstrFGenes.size( ); ++i )
				cout << vecstrFGenes[ i ] << endl;
			return 1;}

		vector<CDat*>			DatOutB;
		DatOutB.resize( vecvecstrInputs.size( ) );
		for( i = 0; i < DatOutB.size( ); ++i ){
			DatOutB[ i ] = new CDat( );
			DatOutB[ i ]->Open( vecstrFGenes );}

		vector<CDat*>		vecDataIntZero;
		vector<CDat*>		vecDataIntOne;
		vecDataIntZero.resize( vecvecstrInputs.size( ) );
		vecDataIntOne.resize( vecvecstrInputs.size( ) );
		for( i = 0; i < vecvecstrInputs.size( ); ++i ){
			vecDataIntZero[ i ] = new CDat( );
			vecDataIntOne[ i ] = new CDat( );
			vecDataIntZero[ i ]->Open( vecstrFGenes );
			vecDataIntOne[ i ]->Open( vecstrFGenes );}
		
		for( i = 0; i < vecvecstrInputs.size( ); ++i ){
			map<string,size_t>::const_iterator	iterZero;
			for( j = 0; j < vecvecstrInputs[ i ].size( ); ++j ){
				iZero = ( ( iterZero = mapZeros.find( vecvecstrInputs[ i ][ j ] ) ) == mapZeros.end( ) ) ? -1 : iterZero->second;

				CDataPair		DataF;
				if( !( DataF.Open( vecvecstrInputs[ i ][ j ].c_str( ), false, !!sArgs.memmap_flag ) ||
					DataF.Open( vecvecstrInputs[ i ][ j ].c_str( ), true, !!sArgs.memmap_flag ) ) ){
						cerr << "Could not open:" << vecvecstrInputs[ i ][ j ] << endl;
						return 1;}
				vector<size_t>		vecGeneIndex;
				vecGeneIndex.resize( vecstrFGenes.size( ) );
				for( k = 0; k < vecstrFGenes.size( ); ++k )
					vecGeneIndex[ k ] = DataF.GetGene( vecstrFGenes[ k ] );

				for( iDatOne = 0; iDatOne < vecstrFGenes.size( ); ++iDatOne ){ 
					for( iDatTwo = ( iDatOne + 1 ); iDatTwo < vecstrFGenes.size( ); ++iDatTwo ){ 
						iGeneOne = vecGeneIndex[ iDatOne ];
						iGeneTwo = vecGeneIndex[ iDatTwo ];
						float		DFValue;
						DFValue = ( ( iGeneOne == -1 ) || ( iGeneTwo == -1 ) ) ? CMeta::GetNaN( ) : DataF.Get( iGeneOne, iGeneTwo );
						if( !CMeta::IsNaN( DFValue ) ) {
							size_t		DatPFI = DataF.Quantize( DFValue );
							float		Zero, One;

							if( CMeta::IsNaN( Zero = vecDataIntZero[ i ]->Get( iDatOne, iDatTwo ) ) )
								vecDataIntZero[ i ]->Set( iDatOne, iDatTwo, log( vecDataGSpZero[ i ][ j ][ DatPFI ] ) );
							else
								vecDataIntZero[ i ]->Set( iDatOne, iDatTwo, Zero + log( vecDataGSpZero[ i ][ j ][ DatPFI ] ) );

							if( CMeta::IsNaN( One = vecDataIntOne[ i ]->Get( iDatOne, iDatTwo ) ) )
								vecDataIntOne[ i ]->Set( iDatOne, iDatTwo, log( vecDataGSpOne[ i ][ j ][ DatPFI ] ) );
							else
								vecDataIntOne[ i ]->Set( iDatOne, iDatTwo, One + log( vecDataGSpOne[ i ][ j ][ DatPFI ] ) );}

						if( CMeta::IsNaN( DFValue ) && ( iZero != -1 ) ){
							float		Zero, One;

							if( CMeta::IsNaN( Zero = vecDataIntZero[ i ]->Get( iDatOne, iDatTwo ) ) )
								vecDataIntZero[ i ]->Set( iDatOne, iDatTwo, log( vecDataGSpZero[ i ][ j ][ iZero ] ) );
							else
								vecDataIntZero[ i ]->Set( iDatOne, iDatTwo, Zero + log( vecDataGSpZero[ i ][ j ][ iZero ] ) );

							if( CMeta::IsNaN( One = vecDataIntOne[ i ]->Get( iDatOne, iDatTwo ) ) )
								vecDataIntOne[ i ]->Set( iDatOne, iDatTwo, log( vecDataGSpOne[ i ][ j ][ iZero ] ) );
							else
								vecDataIntOne[ i ]->Set( iDatOne, iDatTwo, One + log( vecDataGSpOne[ i ][ j ][ iZero ] ) );}}}}

			for( iDatOne = 0; iDatOne < vecstrFGenes.size( ); ++iDatOne ){ 
				for( iDatTwo = ( iDatOne + 1 ); iDatTwo < vecstrFGenes.size( ); ++iDatTwo ){ 
					float		sumOneB = 0, sumZeroB = 0, FinalB = 0 ;
					float		Zero, One;

					if( CMeta::IsNaN( Zero = vecDataIntZero[ i ]->Get( iDatOne, iDatTwo ) ) ){
						vecDataIntZero[ i ]->Set( iDatOne, iDatTwo, 0.0f);
						sumZeroB =  log( vecPrior[ i ][ 0 ] );}
					else 
						sumZeroB =  Zero + log( vecPrior[ i ][ 0 ] );

					if( CMeta::IsNaN( One = vecDataIntOne[ i ]->Get( iDatOne, iDatTwo ) ) ){
						vecDataIntOne[ i ]->Set( iDatOne, iDatTwo, 0.0f );
						sumOneB =  log( vecPrior[ i ][ 1 ] );}
					else
						sumOneB =  One + log( vecPrior[ i ][ 1 ] );

					FinalB = (float) ( 1 / ( 1 + exp ( sumZeroB - sumOneB ) ) );
					DatOutB[ i ]->Set( iDatOne, iDatTwo, FinalB );}}}

		for( i = 0; i < vecvecstrInputs.size( ); ++i ){
			DatOutB[ i ]->Save( ( ( string ) sArgs.odirectory_arg + '/' + vecstrxInputs[ i ] + 'b' + c_acDab ).c_str( ) );
			vecDataIntZero[ i ]->Save( ( ( string ) sArgs.l0directory_arg + '/' + vecstrxInputs[ i ] + c_acDab ).c_str( ) );	
			vecDataIntOne[ i ]->Save( ( ( string ) sArgs.l1directory_arg + '/' + vecstrxInputs[ i ] + c_acDab ).c_str( ) );}	

		vector<size_t>		MapstrGenes;	
		if( sArgs.genex_arg ) {
			CGenome				Genome;		
			CGenes				GenesEx( Genome );
			if( !GenesEx.Open( sArgs.genex_arg ) ) {
				cerr << "Could not open: " << sArgs.genex_arg << endl;
				return 1; } 
			vecstrGenes.resize( vecstrFGenes.size( ) - GenesEx.GetGenes( ) );
			MapstrGenes.resize( vecstrGenes.size( ) );
			size_t		GenesExCount = 0;
			for( i = 0; i < vecstrFGenes.size( ); ++i ){
				if( !GenesEx.IsGene( vecstrFGenes[ i ] ) ){
					vecstrGenes[ GenesExCount ] = vecstrFGenes[ i ];
					MapstrGenes[ GenesExCount ] = i;
					GenesExCount++;}}}
		else{
			vecstrGenes.resize( vecstrFGenes.size( ) );
			MapstrGenes.resize( vecstrFGenes.size( ) );
			for( i = 0; i < vecstrFGenes.size( ); ++i ){
				vecstrGenes[ i ] = vecstrFGenes[ i ];
				MapstrGenes[ i ] = i;}}
		
		CDat		DatOut00, DatOut01, DatOut10, DatOut11; 
		DatOut00.Open( vecstrxInputs );
		DatOut01.Open( vecstrxInputs );
		DatOut10.Open( vecstrxInputs );
		DatOut11.Open( vecstrxInputs );

		if( sArgs.uniformjoint_flag ){
			for( iDatOne = 0; iDatOne < vecstrInputs.size( ); ++iDatOne ){
				for( iDatTwo = ( iDatOne + 1 ); iDatTwo < vecstrInputs.size( ); ++iDatTwo ){
					for( i = 0; i < vecPrior[ iDatOne ].size( ); ++i ){
						for( j = 0; j < vecPrior[ iDatTwo ].size( ); ++j ){
							float	UnifProb = 1.0f / ( vecPrior[ iDatOne ].size( ) * vecPrior[ iDatTwo ].size( ) );
							DatOut00.Set( iDatOne, iDatTwo, UnifProb );
							DatOut01.Set( iDatOne, iDatTwo, UnifProb );
							DatOut10.Set( iDatOne, iDatTwo, UnifProb );
							DatOut11.Set( iDatOne, iDatTwo, UnifProb );}}}}}
					
		if( !sArgs.uniformjoint_flag ){

			for( iDatOne = 0; iDatOne < vecstrInputs.size( ); ++iDatOne ) {
				if( ( iDatOne + 1 ) == vecstrInputs.size( ) )
					break;

				veciGenesOneI.resize( vecstrGenes.size( ) );
				CDataset	DataInd;
				vector<string>				vecstrIndInputs;
				vecstrIndInputs.resize( vecvecstrInputs[ iDatOne ].size( ) + 1 );
				for( i = 0; i < vecvecstrInputs[ iDatOne ].size( ); ++i ){
					vecstrIndInputs[ i ] = vecvecstrInputs[ iDatOne ][ i ].c_str( );}
				vecstrIndInputs[ i ] = ( string )sArgs.adirectory_arg + '/' + vecstrxInputs[ iDatOne ] + c_acDab;
				DataInd.OpenGenes( vecstrIndInputs );
				for( i = 0; i < vecstrGenes.size( ); ++i )
					veciGenesOneI[ i ] = DataInd.GetGene( vecstrGenes[ i ] );

				for( iDatTwo = ( iDatOne + 1 ); iDatTwo < vecstrInputs.size( ); ++iDatTwo ) {
					veciGenesTwoI.resize( vecstrGenes.size( ) );
					CDataset	DataInd;
					vector<string>				vecstrIndInputs;
					vecstrIndInputs.resize( vecvecstrInputs[ iDatTwo ].size( ) + 1 );
					for( i = 0; i < vecvecstrInputs[ iDatTwo ].size( ); ++i ){
						vecstrIndInputs[ i ] = vecvecstrInputs[ iDatTwo ][ i ].c_str( );}
					vecstrIndInputs[ i ] = ( string )sArgs.adirectory_arg + '/' + vecstrxInputs[ iDatTwo ] + c_acDab;
					DataInd.OpenGenes( vecstrIndInputs );
					for( i = 0; i < vecstrGenes.size( ); ++i )
						veciGenesTwoI[ i ] = DataInd.GetGene( vecstrGenes[ i ] );

					vecveciJoint.resize( vecPrior[ iDatOne ].size( ) );
					for( i = 0; i < vecveciJoint.size( ); ++i ) {
						vecveciJoint[ i ].resize( vecPrior[ iDatTwo ].size( ) );
						fill( vecveciJoint[ i ].begin( ), vecveciJoint[ i ].end( ), 0 ); }

					for( i = iCountJoint = 0; i < vecstrGenes.size( ); ++i ) {
						if( ( veciGenesOneI[ i ] == -1 ) || ( veciGenesTwoI[ i ] == -1 ) )
							continue;

						for( j = ( i + 1 ); j < vecstrGenes.size( ); ++j ) {
							if( ( veciGenesOneI[ j ] == -1 ) || ( veciGenesTwoI[ j ] == -1 ) )
								continue;

							iValueOne = DatOutB[ iDatOne ]->Get( MapstrGenes[ i ], MapstrGenes[ j ] ) >= Threshold ? 1 : 0;
							iValueTwo = DatOutB[ iDatTwo ]->Get( MapstrGenes[ i ], MapstrGenes[ j ] ) >= Threshold ? 1 : 0;
							vecveciJoint[ iValueOne ][ iValueTwo ]++;
							iCountJoint++;}}

					for( i = 0; i < vecveciJoint.size( ); ++i ){
						for( j = 0; j < vecveciJoint.size( ); ++j ){
							vecveciJoint[ i ][ j ]++;
							iCountJoint++;}}					

					vector<vector<float> >		vecveciTJoint;
					vecveciTJoint.resize( vecveciJoint.size( ) );
					for( i = 0; i < vecveciTJoint.size( ); ++i ){
						vecveciTJoint[ i ].resize( vecveciJoint[ i ].size( ) );
						for( j = 0; j < vecveciTJoint.size( ); ++j ){ 
							vecveciTJoint[ i ][ j ] = ( (float)vecveciJoint[ i ][ j ] ) / iCountJoint;}}

					DatOut00.Set( iDatOne, iDatTwo, vecveciTJoint[ 0 ][ 0 ] );
					DatOut01.Set( iDatOne, iDatTwo, vecveciTJoint[ 0 ][ 1 ] );
					DatOut10.Set( iDatOne, iDatTwo, vecveciTJoint[ 1 ][ 0 ] );
					DatOut11.Set( iDatOne, iDatTwo, vecveciTJoint[ 1 ][ 1 ] );}}}
				
		DatOut00.Save( ( ( string ) sArgs.jdirectory_arg + "/Learned00" + c_acDab ).c_str( ) );
		DatOut01.Save( ( ( string ) sArgs.jdirectory_arg + "/Learned01" + c_acDab ).c_str( ) );	
		DatOut10.Save( ( ( string ) sArgs.jdirectory_arg + "/Learned10" + c_acDab ).c_str( ) );
		DatOut11.Save( ( ( string ) sArgs.jdirectory_arg + "/Learned11" + c_acDab ).c_str( ) );
		
		for( i = 0; i < vecvecstrInputs.size( ); ++i ){
			delete DatOutB[ i ];
			delete vecDataIntZero[ i ];	
			delete vecDataIntOne[ i ];}
		
		return 1;}
    
	else{

		vector<string>		vecDataIntZero;
		vector<string>		vecDataIntOne;
		vecDataIntZero.resize( vecstrxInputs.size( ) );
		vecDataIntOne.resize( vecstrxInputs.size( ) );
		for( i = 0; i < vecstrxInputs.size( ); ++i ){
			vecDataIntZero[ i ] = ( string )sArgs.l0directory_arg + '/' + vecstrxInputs[ i ] + c_acDab; 
			vecDataIntOne[ i ] = ( string )sArgs.l1directory_arg + '/' + vecstrxInputs[ i ] + c_acDab;}

		vector<CDataPair*>			DataZero;
		DataZero.resize( vecDataIntZero.size( ) );
		for( i = 0; i < vecDataIntZero.size( ); ++i ){
			DataZero[ i ] = new CDataPair( );
			if( !( DataZero[ i ]->Open( vecDataIntZero[ i ].c_str( ), false, !!sArgs.memmap_flag ) ||
				DataZero[ i ]->Open( vecDataIntZero[ i ].c_str( ), true, !!sArgs.memmap_flag ) ) ){
					cerr << "Could not open:" << vecDataIntZero[ i ] << endl;
					return 1;}}
		
		vector<CDataPair*>			DataOne;
		DataOne.resize( vecDataIntOne.size( ) );
		for( i = 0; i < vecDataIntOne.size( ); ++i ){
			DataOne[ i ] = new CDataPair( );
			if( !( DataOne[ i ]->Open( vecDataIntOne[ i ].c_str( ), false, !!sArgs.memmap_flag ) ||
				DataOne[ i ]->Open( vecDataIntOne[ i ].c_str( ), true, !!sArgs.memmap_flag ) ) ){
					cerr << "Could not open:" << vecDataIntOne[ i ] << endl;
					return 1;}}

		vector<string>		vecstrFGenes;
		{	CDataset	DataF;

		DataF.OpenGenes( vecDataIntZero );
		vecstrFGenes.resize( DataF.GetGenes( ) );
		copy( DataF.GetGeneNames( ).begin( ), DataF.GetGeneNames( ).end( ), vecstrFGenes.begin( ) );}

		vector<vector<vector<vector<float> > > >		vec4OSpGSp;
		size_t			JointDim = vecPrior[ 0 ].size( );

		vector<vector<string> >		VecLearnedNames;
		VecLearnedNames.resize( JointDim );
		VecLearnedNames[ 0 ].resize( JointDim );
		VecLearnedNames[ 1 ].resize( JointDim );
		VecLearnedNames[ 0 ][ 0 ] = "Learned00";
		VecLearnedNames[ 0 ][ 1 ] = "Learned01";
		VecLearnedNames[ 1 ][ 0 ] = "Learned10";
		VecLearnedNames[ 1 ][ 1 ] = "Learned11";

		vec4OSpGSp.resize( vecstrxInputs.size( ) );
		for( iDatOne = 0; iDatOne < vecstrxInputs.size( ); ++iDatOne ) {
			vec4OSpGSp[ iDatOne ].resize( vecstrxInputs.size( ) );
			for( iDatTwo = 0; iDatTwo < vecstrxInputs.size( ); ++iDatTwo ) {
				if( iDatTwo == iDatOne )
					continue;
				vec4OSpGSp[ iDatOne ][ iDatTwo ].resize( JointDim );
				for( i = 0; i < JointDim; ++i ){
					vec4OSpGSp[ iDatOne ][ iDatTwo ][ i ].resize( JointDim );}}}
	
	for( i = 0; i < JointDim; ++i ){
		for( j = 0; j < JointDim; ++j ){
			CDataPair		DataF;
			string			JointFile = ( string )sArgs.jdirectory_arg + VecLearnedNames[ i ][ j ] + c_acDab; 
			if( !( DataF.Open( JointFile.c_str( ), false, !!sArgs.memmap_flag ) ||
					DataF.Open( JointFile.c_str( ), true, !!sArgs.memmap_flag ) ) ){
						cerr << "Could not open:" << JointFile << endl;
						return 1;}

			vector<size_t>	veciSpecies;
			veciSpecies.resize( vecstrxInputs.size( ) );
			for( k = 0; k < vecstrxInputs.size( ); ++k )
				veciSpecies[ k ] = DataF.GetGene( vecstrxInputs[ k ] );
			
			for( iDatOne = 0; iDatOne < vecstrxInputs.size( ); ++iDatOne ) {
				for( iDatTwo = ( iDatOne + 1 ); iDatTwo < vecstrxInputs.size( ); ++iDatTwo ) {
					vec4OSpGSp[ iDatOne ][ iDatTwo ][ i ][ j ] = DataF.Get( veciSpecies[ iDatOne ], veciSpecies[ iDatTwo ] );
					vec4OSpGSp[ iDatTwo ][ iDatOne ][ j ][ i ] = DataF.Get( veciSpecies[ iDatOne ], veciSpecies[ iDatTwo ] );}}}}
	
    vector<float>		veciValueOnep;
	vector<float>		veciValueTwop;

	for( iDatOne = 0; iDatOne < vecstrxInputs.size( ); ++iDatOne ) {
		for( iDatTwo = ( iDatOne + 1 ); iDatTwo < vecstrxInputs.size( ); ++iDatTwo ) {
			
			veciValueOnep.resize( JointDim );
			veciValueTwop.resize( JointDim );

			fill( veciValueOnep.begin( ), veciValueOnep.end( ), 0.0f );
			fill( veciValueTwop.begin( ), veciValueTwop.end( ), 0.0f );

			for( i = 0; i < JointDim; ++i ){
				for( j = 0; j < JointDim; ++j ){
					veciValueOnep[ i ] += vec4OSpGSp[ iDatOne ][ iDatTwo ][ i ][ j ];
					veciValueTwop[ j ] += vec4OSpGSp[ iDatOne ][ iDatTwo ][ i ][ j ];}}

			for( i = 0; i < JointDim; ++i ){
				for( j = 0; j < JointDim; ++j ){
					vec4OSpGSp[ iDatOne ][ iDatTwo ][ i ][ j ] /= veciValueOnep[ i ];
					vec4OSpGSp[ iDatTwo ][ iDatOne ][ i ][ j ] /= veciValueTwop[ i ];}}}}
    
	vector<CDat*>			DatOutCS;
	
	DatOutCS.resize( vecstrxInputs.size( ) );
	for( i = 0; i < DatOutCS.size( ); ++i ){
		DatOutCS[ i ] = new CDat( );
		DatOutCS[ i ]->Open( vecstrFGenes );}

	vector<vector<size_t> >		veciGenesOne, veciGenesZero;
	veciGenesOne.resize( vecDataIntOne.size( ) );
	veciGenesZero.resize( vecDataIntZero.size( ) );
	for( i = 0; i < vecstrxInputs.size( ); ++i ){
		veciGenesOne[ i ].resize( vecstrFGenes.size( ) );
		veciGenesZero[ i ].resize( vecstrFGenes.size( ) );
		for( j = 0; j < vecstrFGenes.size( ); ++j ){
			veciGenesOne[ i ][ j ] = DataOne[ i ]->GetGene( vecstrFGenes[ j ] );
			veciGenesZero[ i ][ j ] = DataZero[ i ]->GetGene( vecstrFGenes[ j ] );}}

	if( !sArgs.holdout_flag ){
		for( i = 0; i < vecstrxInputs.size( ); ++i ){
			for( iDatOne = 0; iDatOne < vecstrFGenes.size( ); ++iDatOne ){ 
				for( iDatTwo = ( iDatOne + 1 ); iDatTwo < vecstrFGenes.size( ); ++iDatTwo ){ 
					float		sumOne = 0, sumZero = 0, Final = 0 ;
					for( k = 0; k < vecstrxInputs.size( ); ++k ){
						if( k != i ){
							sumOne += log( exp( log( vec4OSpGSp[ i ][ k ][ 1 ][ 0 ] ) + DataZero[ k ]->Get( veciGenesZero[ k ][ iDatOne ], veciGenesZero[ k ][ iDatTwo ] ) ) + exp( log( vec4OSpGSp[ i ][ k ][ 1 ][ 1 ] ) + DataOne[ k ]->Get( veciGenesOne[ k ][ iDatOne ], veciGenesOne[ k ][ iDatTwo ] ) ) );
							sumZero += log( exp( log( vec4OSpGSp[ i ][ k ][ 0 ][ 0 ] ) + DataZero[ k ]->Get( veciGenesZero[ k ][ iDatOne ], veciGenesZero[ k ][ iDatTwo ] ) ) + exp( log( vec4OSpGSp[ i ][ k ][ 0 ][ 1 ] ) + DataOne[ k ]->Get( veciGenesOne[ k ][ iDatOne ], veciGenesOne[ k ][ iDatTwo ] ) ) );}}
					sumOne += ( DataOne[ i ]->Get( veciGenesOne[ i ][ iDatOne ], veciGenesOne[ i ][ iDatTwo ] ) + log( vecPrior[ i ][ 1 ] ) );
					sumZero += ( DataZero[ i ]->Get( veciGenesZero[ i ][ iDatOne ], veciGenesZero[ i ][ iDatTwo ] ) + log( vecPrior[ i ][ 0 ] ) );
					Final = (float) ( 1 / ( 1 + exp ( sumZero - sumOne ) ) );
					DatOutCS[ i ]->Set( iDatOne, iDatTwo, Final );}}}}

	else{
		for( i = 0; i < vecstrxInputs.size( ); ++i ){
			for( iDatOne = 0; iDatOne < vecstrFGenes.size( ); ++iDatOne ){ 
				for( iDatTwo = ( iDatOne + 1 ); iDatTwo < vecstrFGenes.size( ); ++iDatTwo ){ 
					float		sumOne = 0, sumZero = 0, Final = 0 ;
					for( k = 0; k < vecstrxInputs.size( ); ++k ){
						if( k != i ){
							sumOne += log( exp( log( vec4OSpGSp[ i ][ k ][ 1 ][ 0 ] ) + DataZero[ k ]->Get( veciGenesZero[ k ][ iDatOne ], veciGenesZero[ k ][ iDatTwo ] ) ) + exp( log( vec4OSpGSp[ i ][ k ][ 1 ][ 1 ] ) + DataOne[ k ]->Get( veciGenesOne[ k ][ iDatOne ], veciGenesOne[ k ][ iDatTwo ] ) ) );
							sumZero += log( exp( log( vec4OSpGSp[ i ][ k ][ 0 ][ 0 ] ) + DataZero[ k ]->Get( veciGenesZero[ k ][ iDatOne ], veciGenesZero[ k ][ iDatTwo ] ) ) + exp( log( vec4OSpGSp[ i ][ k ][ 0 ][ 1 ] ) + DataOne[ k ]->Get( veciGenesOne[ k ][ iDatOne ], veciGenesOne[ k ][ iDatTwo ] ) ) );}}
					sumOne += log( vecPrior[ i ][ 1 ] );
					sumZero += log( vecPrior[ i ][ 0 ] );
					Final = (float) ( 1 / ( 1 + exp ( sumZero - sumOne ) ) );
					DatOutCS[ i ]->Set( iDatOne, iDatTwo, Final );}}}}
	
	for( i = 0; i < vecstrxInputs.size( ); ++i ){
		DatOutCS[ i ]->Save( ( ( string ) sArgs.odirectory_arg + '/' + vecstrxInputs[ i ] + c_acDab ).c_str( ) );
		delete DatOutCS[ i ];
		delete DataOne[ i ];
		delete DataZero[ i ];}}

	return 0;} 