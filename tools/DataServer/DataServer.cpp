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
#include "measure.h"
#include "DataServer.h"

const CDataServer::TPFNProcessor	CDataServer::c_apfnProcessors[]	=
	{&CDataServer::ProcessDatasetSearch, &CDataServer::ProcessDatasetMeasure};

struct SData {
    vector<size_t>* m_veciGenes; 
    vector<float>* m_vecfDataWeights;
    vector<float>* m_vecfScores; 
    vector<float>* m_vecfTotal;
    CDataServer* m_Server;
    //CPCL* m_QueryCorScores;
    //CPCL* m_QueryCorTotal;
    CFullMatrix<float>* m_QueryCorScores;
    CFullMatrix<float>* m_QueryCorTotal;
    
};

struct SWeight {
    size_t m_iData;
    //CPCL* m_CorPCL;
    CFullMatrix<float>* m_CorMat;
    float m_cutoff;
    float m_base;
    vector<float>* m_vecfDataWeights;
    CDataServer* m_Server;
    size_t m_iThread;
};

int main( int iArgs, char** aszArgs ) {
	gengetopt_args_info			sArgs;
    static const size_t c_iBuffer = 1024;
    char acBuffer[ c_iBuffer ];

	CServer						Server;
    CDatabase Database(false);
    CPCL VarPCL;
    vector<string> vecstrDatasets;
    ifstream ifsm;
    int iRet;
    vector<size_t> veciPCLGeneIdx, veciPCLDataIdx;
    size_t i;

    iRet = cmdline_parser2( iArgs, aszArgs, &sArgs, 0, 1, 0 );

	cerr << "Loading the database..." << endl;
	if( !Database.Open( sArgs.database_arg ) ) {
		cerr << "Could not open: " << sArgs.database_arg << endl;
		return 1; 
    }

    if( !VarPCL.Open( sArgs.variances_arg, 0 ) ) {
        cerr << "Could not open: " << sArgs.database_arg << endl;
        return 1;
    }

    if( sArgs.datasets_arg ) { 
        ifsm.clear( );
        ifsm.open(sArgs.datasets_arg);
        while(!ifsm.eof()){
            ifsm.getline(acBuffer, c_iBuffer -1); 
            if(acBuffer[0]==0)
                break;
            acBuffer[c_iBuffer-1] = 0; 
            vector<string> tok; 
            CMeta::Tokenize(acBuffer, tok, " \t");
            vecstrDatasets.push_back(tok[1]);
        }
        ifsm.close();
    }  

	cerr << "Quant file: " << sArgs.quant_arg << endl;

    // Map PCL genes to CDatabase genes
    veciPCLGeneIdx.resize( Database.GetGenes() );
    for ( i = 0; i < veciPCLGeneIdx.size(); i++ ) {
        veciPCLGeneIdx[i] = VarPCL.GetGene( Database.GetGene( i ) );
    } 

    veciPCLDataIdx.resize( vecstrDatasets.size() );
    for ( i = 0; i < veciPCLDataIdx.size(); i++ ) {
        veciPCLDataIdx[i] = VarPCL.GetExperiment( vecstrDatasets[i] );
    } 


	SDataServerData	sData( Database, vecstrDatasets, VarPCL, veciPCLGeneIdx, veciPCLDataIdx, sArgs.quant_arg, sArgs.threads_arg );
	CDataServer	DataServer( 0, "", sData );

    cerr << "Maximum number of threads: " << sArgs.threads_arg << endl;

	Server.Initialize( sArgs.port_arg, sArgs.timeout_arg, &DataServer );
#ifdef WIN32
	pthread_win32_process_attach_np( );
#endif // WIN32
	Server.Start( );
#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32

	return 0; }


CDataServer::CDataServer( SOCKET iSocket, const string& strConnection, const SDataServerData& sData ) :
	m_iSocket(iSocket), m_strConnection(strConnection), m_sData(sData) {

	if( m_strConnection.length( ) > 0 )
		cerr << "New connection from: " << m_strConnection << endl; }

CDataServer::~CDataServer( ) {

}

IServerClient* CDataServer::NewInstance( SOCKET iSocket, uint32_t iHost, uint16_t sPort ) {
	string	strConnection;
	char	acBuffer[ 16 ];
	in_addr	sAddr;

#pragma warning(disable : 4996)
	sprintf( acBuffer, "%hu", sPort );
#pragma warning(default : 4996)
	sAddr.s_addr = htonl( iHost );
	strConnection = (string)inet_ntoa( sAddr ) + ":" + acBuffer;
	return new CDataServer( iSocket, strConnection, m_sData ); }

void CDataServer::Destroy( ) {

	cerr << "Disconnected: " << m_strConnection << endl;

	delete this; }

bool CDataServer::ProcessMessage( const vector<unsigned char>& vecbMessage ) {
	size_t	i, iProcessed, iOffset;

	for( iOffset = 0; iOffset < vecbMessage.size( ); iOffset += ( iProcessed + 1 ) ) {
		cerr << "LOG	" << time( NULL ) << '\t' << m_strConnection << endl; //'\t' << hex;
		if( vecbMessage[ iOffset ] >= ARRAYSIZE(c_apfnProcessors) ) {
			cerr << m_strConnection << " unknown opcode: " << (int)vecbMessage[ iOffset ] << endl;
			return false; }
		else {
			cerr << m_strConnection << " opcode: " << (int)vecbMessage[ iOffset ] << endl;
		}
		if( ( iProcessed = (this->*c_apfnProcessors[ vecbMessage[ iOffset ] ])( vecbMessage,
			iOffset + 1 ) ) == -1 )
			return false; }

	return true; }


void* GetScoresThread( void* pData ) {
    SData* data = (SData*)pData;

    data->m_Server->GetScores( *(data->m_veciGenes), *(data->m_vecfDataWeights), 
        *(data->m_vecfScores), *(data->m_vecfTotal), data->m_QueryCorScores, data->m_QueryCorTotal );

    return NULL; 
}



void CDataServer::GetScores( const vector<size_t>& veciGenes, const vector<float>& vecfDataWeights,
        vector<float>& vecfScores, vector<float>& vecfTotal, 
        CFullMatrix<float>* AllScores = NULL, CFullMatrix<float>* AllCounts = NULL ) {

    size_t iGene, iTargetPCL, iDataPCL, iGenePCL, q, i, j, iOffset;
    vector<unsigned char> vecbData;
    float v, t;

    // Initialize data structures
    vecfScores.resize( GetDatabase().GetGenes() );
    fill( vecfScores.begin(), vecfScores.end(), 0 );
    vecfTotal.resize( GetDatabase().GetGenes() );
    fill( vecfTotal.begin(), vecfTotal.end(), 0 );

    if ( AllScores ) {
        AllScores->Clear();
        AllCounts->Clear();
    }

    // Iterate over query genes
    for ( q = 0; q < veciGenes.size(); q++ ) {
        iGene = veciGenes[q];
        iGenePCL = GetPCLGeneIdx()[ iGene ]; 

        if ( iGene == -1 || iGenePCL == -1 ) {
            cerr << "Missing gene: " << GetDatabase().GetGene( iGene ) << endl;
            continue;
        }
        cerr << q << ": " << GetDatabase().GetGene( iGene ) << endl;
 
        vecbData.clear();
        GetDatabase().Get( iGene, vecbData ); 

        // Iterate over all genes
        for ( i = 0; i < GetDatabase().GetGenes(); i++ ) {
            iOffset = i * GetDatabase().GetDatasets();
            
            for ( j = 0; j < GetDatabase().GetDatasets(); j++ ) {
                iDataPCL = GetPCLDataIdx()[j];
                if ( iDataPCL == -1 ) continue;

                v = (float)vecbData[ iOffset + j ];     
                if ( v == 0xFF ) continue;

                // Convert to z-score from bin value
                v = -5 + v*(10.0/255.0);
            
                // Ignore negative z-scores
                if ( v < 0 ) v = 0;

                // Store all scores
                if ( AllScores && AllCounts ) {
                    t = GetVarPCL().Get( iGenePCL, iDataPCL );
                    if ( CMeta::IsNaN( t ) ) continue;

                    AllScores->Set( i, j, AllScores->Get( i, j ) + ( v * t ) );
                    AllCounts->Set( i, j, AllCounts->Get( i, j ) + t );
                }
                // Summarize scores
                else {
                    t = GetVarPCL().Get( iGenePCL, iDataPCL ) * vecfDataWeights[ j ];
                    if ( CMeta::IsNaN( t ) ) continue;
                    v *= t;
                    vecfScores[i] += v;
                    vecfTotal[i] += t; 
                }
            }
        } 
    }
}
	
void* GetDataWeightsThread( void* pData ) {
    SWeight* data = (SWeight*)pData;

    data->m_Server->GetDataWeights( data->m_iData, *(data->m_CorMat), 
        data->m_cutoff, data->m_base, *(data->m_vecfDataWeights), data->m_iThread );

    return NULL; 
}


void CDataServer::GetDataWeights( size_t iData, CFullMatrix<float>& CorMat, float cutoff, float s, 
        vector<float>& vecfDataWeights, int iThread = -1 ) {
    size_t i, j, iDataPCL;
    float d1, d2, w1, w2, fSim;
    CMeasurePearson Pearson;
    vector<float> adOne, adTwo, adW1, adW2;

    cerr << "Calculating correlations for thread: " << iThread << endl;

    for( j = 0; j < GetDatabase().GetDatasets(); j++ ) {
        if ( iThread >= 0 && j % GetMaxThreads() != iThread )
            continue;

        adOne.clear();
        adTwo.clear();
        adW1.clear();
        adW2.clear();

        for( i = 0; i < GetDatabase().GetGenes(); i++ ) {
            d1 = CorMat.Get( i, j ); 
            d2 = CorMat.Get( i, iData );

            if ( CMeta::IsNaN( d1 ) || CMeta::IsNaN( d2 ) ) continue;
            if ( d1 < cutoff && d2 < cutoff ) continue;

            adOne.push_back( d1 );
            adTwo.push_back( d2 );

            w1 = pow ( s, d1 ) - 1;
            w2 = pow ( s, d2 ) - 1;

            adW1.push_back( w1 );
            adW2.push_back( w2 );
        }
        fSim = 0;

        if ( adOne.size() ) {
            fSim = (float) Pearson.Measure( &adOne[0], adOne.size(), &adTwo[0], adTwo.size(), 
                IMeasure::EMapNone, &adW1[0], &adW2[0] );
        }

        if ( fSim < 0 ) fSim = 0;

        vecfDataWeights[ j ] = fSim;
    }
}


size_t CDataServer::ProcessDatasetSearch( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	size_t		iStart, i, j, t;
	uint32_t	iGene, iDataset, iSize;
    float   iCut, iExp;
    vector<string> vecstrFeatures;
    vector<size_t> veciGenes;
    vector<float> vecfDataWeights, vecfScores, vecfTotal;
    vector<pthread_t> vecpthdThreads;
    vector<SData> vecsData;
    vector<SWeight> vecsWeight;
    CFullMatrix<float> QueryCorScores, QueryCorTotal;

    size_t iThreads = 1; 

	if( ( iOffset + sizeof(iDataset) ) > vecbMessage.size( ) )
		return -1;
	iStart = iOffset;
	iDataset = *(uint32_t*)&vecbMessage[ iOffset ];

    cerr << "Processing Search" << endl;
    cerr << "Dataset: " << GetDatasetNames()[iDataset] << endl; 
    
    iOffset += sizeof(iDataset);
	if( iOffset + sizeof(iCut) > vecbMessage.size( ) )
		return -1;
    iCut = *(float*)&vecbMessage[ iOffset ];

    iOffset += sizeof(iCut);
	if( iOffset + sizeof(iExp) > vecbMessage.size( ) )
		return -1;
    iExp = *(float*)&vecbMessage[ iOffset ];

    cerr << "Parameters: " << iCut << ", " << iExp << endl;

    for( iOffset += sizeof(iExp); ( iOffset + sizeof(iGene) ) <= vecbMessage.size();
        iOffset += sizeof(iGene) ) {
        iGene = *(uint32_t*)&vecbMessage[ iOffset ];
        veciGenes.push_back(iGene); 
    }

    iThreads = std::min( (int)veciGenes.size(), GetMaxThreads() );    

    cerr << "Setting threads to: " << iThreads << endl;


    // ------------------------------------------------------------------------/
    // Get correlations to query genes
    
    // Initialize dataset weights to equal weighting
    vecfDataWeights.resize( GetDatabase().GetDatasets() );
    fill( vecfDataWeights.begin(), vecfDataWeights.end(), 1 );

    QueryCorScores.Initialize( GetDatabase().GetGenes(), GetDatabase().GetDatasets() );
    QueryCorTotal.Initialize( GetDatabase().GetGenes(), GetDatabase().GetDatasets() );
    QueryCorScores.Clear();
    QueryCorTotal.Clear();

    vecsData.resize( iThreads/2 );
    vecpthdThreads.resize( iThreads/2 );
    for( i = 0; i < vecsData.size(); i++ ) {
        vecsData[ i ].m_vecfDataWeights = &vecfDataWeights;         
        vecsData[ i ].m_vecfScores = new vector<float>(); 
        vecsData[ i ].m_vecfTotal = new vector<float>(); 
        vecsData[ i ].m_veciGenes = new vector<size_t>(); 
        vecsData[ i ].m_Server = (CDataServer*)this;
        vecsData[ i ].m_QueryCorScores = new CFullMatrix<float>(); 
        vecsData[ i ].m_QueryCorTotal = new CFullMatrix<float>();

        vecsData[ i ].m_QueryCorScores->Initialize( GetDatabase().GetGenes(), GetDatabase().GetDatasets() );
        vecsData[ i ].m_QueryCorTotal->Initialize( GetDatabase().GetGenes(), GetDatabase().GetDatasets() );

        for( iGene = 0; iGene < veciGenes.size(); iGene++ ) {
            if ( iGene % (iThreads/2) == i ) {
                vecsData[ i ].m_veciGenes->push_back( veciGenes[ iGene ] );
            }
        }
        cerr << "Creating " << i << endl;
        pthread_create( &vecpthdThreads[ i ], NULL, GetScoresThread, &vecsData[ i ] ); 
    }
    for( i = 0; i < vecpthdThreads.size(); i++ ) {
	cerr << "Wait: " << i << endl;
        pthread_join( vecpthdThreads[ i ], NULL );
    }

    // Collect results
    for ( t = 0; t < vecsData.size(); t++ ) {
        for ( i = 0; i < GetDatabase().GetGenes(); i++ ) {
            for ( j = 0; j < GetDatabase().GetDatasets(); j++ ) {
                QueryCorScores.Set( i, j, 
                    QueryCorScores.Get( i, j ) + vecsData[ t ].m_QueryCorScores->Get( i, j ) ); 
                QueryCorTotal.Set( i, j, 
                    QueryCorTotal.Get( i, j ) + vecsData[ t ].m_QueryCorTotal->Get( i, j ) ); 
            }
        }
    }

    // Calculate average
    for ( i = 0; i < GetDatabase().GetGenes(); i++ ) {
        for ( j = 0; j < GetDatabase().GetDatasets(); j++ ) {
            QueryCorScores.Set( i, j, 
                QueryCorScores.Get( i, j ) / QueryCorTotal.Get( i, j ) ); 
        }
    }

    for ( t = 0; t < vecsData.size(); t++ ) {
        delete vecsData[ t ].m_vecfScores;
        delete vecsData[ t ].m_vecfTotal;
        delete vecsData[ t ].m_veciGenes;
        delete vecsData[ t ].m_QueryCorScores;
        delete vecsData[ t ].m_QueryCorTotal;
    }

    // ------------------------------------------------------------------------/
    // Calculate correlations
    // GetDataWeights( GetVarPCL().GetExperiment( GetDatasetNames()[ iDataset ] ), 
    //    QueryCorScores, iCut, iExp, vecfDataWeights );

    vecfDataWeights.resize( GetDatabase().GetDatasets() );
    fill( vecfDataWeights.begin(), vecfDataWeights.end(), 0 );

    vecsWeight.resize( iThreads );
    vecpthdThreads.resize( iThreads );
    for( i = 0; i < vecsWeight.size(); i++ ) {

        vecsWeight[ i ].m_iData = iDataset;
        vecsWeight[ i ].m_CorMat = &QueryCorScores;
        vecsWeight[ i ].m_cutoff = iCut;
        vecsWeight[ i ].m_base = iExp;
        vecsWeight[ i ].m_vecfDataWeights = &vecfDataWeights;
        vecsWeight[ i ].m_Server = this;
        vecsWeight[ i ].m_iThread = i;

        cerr << "Creating " << i << endl;
        pthread_create( &vecpthdThreads[ i ], NULL, GetDataWeightsThread, &vecsWeight[ i ] ); 
    }

    for( i = 0; i < vecpthdThreads.size(); i++ ) {
        pthread_join( vecpthdThreads[ i ], NULL );
    }


    // ------------------------------------------------------------------------/
    // Calculated weighted correlations
    vecsData.resize( iThreads );
    vecpthdThreads.resize( iThreads );
    for( i = 0; i < vecsData.size(); i++ ) {
        vecsData[ i ].m_vecfDataWeights = &vecfDataWeights;         
        vecsData[ i ].m_vecfScores = new vector<float>(); 
        vecsData[ i ].m_vecfTotal = new vector<float>(); 
        vecsData[ i ].m_veciGenes = new vector<size_t>(); 
        vecsData[ i ].m_Server = (CDataServer*)this;
        vecsData[ i ].m_QueryCorScores = NULL;
        vecsData[ i ].m_QueryCorTotal = NULL;

        for( iGene = 0; iGene < veciGenes.size(); iGene++ ) {
            if ( iGene % iThreads == i ) {
                vecsData[ i ].m_veciGenes->push_back( veciGenes[ iGene ] );
            }
        }
        cerr << "Creating " << i << endl;
        pthread_create( &vecpthdThreads[ i ], NULL, GetScoresThread, &vecsData[ i ] ); 
    }

    for( i = 0; i < vecpthdThreads.size(); i++ ) {
        pthread_join( vecpthdThreads[ i ], NULL );
    }

    vecfScores.resize( GetDatabase().GetGenes() );
    fill( vecfScores.begin(), vecfScores.end(), 0 );
    vecfTotal.resize( GetDatabase().GetGenes() );
    fill( vecfTotal.begin(), vecfTotal.end(), 0 );

    // Collect results
    for ( i = 0; i < vecsData.size(); i++ ) {
        for ( iGene = 0; iGene < GetDatabase().GetGenes(); iGene++ ) {
            vecfScores[ iGene ] += vecsData[ i ].m_vecfScores->at( iGene );
            vecfTotal[ iGene ] += vecsData[ i ].m_vecfTotal->at( iGene );
        }
    }

    // Calculate average
    for ( i = 0; i < vecfScores.size(); i++ ) {
        vecfScores[i] /= vecfTotal[i];
    }

    // Clean up
    for ( t = 0; t < vecsData.size(); t++ ) {
        delete vecsData[ t ].m_vecfScores;
        delete vecsData[ t ].m_vecfTotal;
        delete vecsData[ t ].m_veciGenes;
    }

    // Send scores
    iSize = (uint32_t)( sizeof(float) * ( GetDatabase().GetGenes() + GetDatabase().GetDatasets() ) ); 
    send( m_iSocket, (char*)&iSize, sizeof(iSize), 0 );
    send( m_iSocket, (char*)&vecfScores[0], sizeof(float)*vecfScores.size(), 0 );
    send( m_iSocket, (char*)&vecfDataWeights[0], sizeof(float)*vecfDataWeights.size(), 0 );
    return ( iOffset - iStart );
}




size_t CDataServer::ProcessDatasetMeasure( const vector<unsigned char>& vecbMessage, size_t iOffset ) {
	size_t iStart, iRet;
	uint32_t iLen, iSize;
	CPCL PCL;
	CDataPair Dat;
	string pclstr, dabstr;

	IMeasure::EMap eM = IMeasure::EMapNone;

	if( ( iOffset + sizeof(iLen) ) > vecbMessage.size( ) )
		return -1;
	iStart = iOffset;
	iLen = *(uint32_t*)&vecbMessage[ iOffset ];
	pclstr = string((char *)&vecbMessage[ iOffset + sizeof(iLen) ], (size_t)iLen);

    cerr << "Processing Measure" << endl;
    cerr << "Filename: " << pclstr << endl;
	iOffset += sizeof(iLen) + (size_t)iLen;

	iRet = CPCL::Distance( pclstr.c_str(), 0, NULL, "pearnorm", false, true, false, NULL, 
		CMeta::GetNaN(), -1, PCL, Dat, eM, false, 0, GetMaxThreads());

	if ( iRet < 0 )
		return -1;

	dabstr = CMeta::Deextension( pclstr ) + ".qdab";
	cerr << "Saving: " << dabstr << endl; 

	Dat.OpenQuants( GetQuantFile().c_str() );
	Dat.Quantize();
	Dat.Save( dabstr.c_str() ); 

    iSize = (uint32_t)( dabstr.length() ); 
    send( m_iSocket, (char*)&iSize, sizeof(iSize), 0 );
    send( m_iSocket, dabstr.c_str(), dabstr.length(), 0 );

	return (iOffset - iStart );
}
