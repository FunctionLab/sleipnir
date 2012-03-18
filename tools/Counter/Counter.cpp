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

class CRegularize;

static const char	c_acDab[]	= ".dab";
static const char	c_acDat[]	= ".dat";
static const char	c_acQDab[]	= ".qdab";
static const char	c_acQuant[]	= ".quant";
static const char	c_acTxt[]	= ".txt";

typedef CFullMatrix<size_t>	CCountMatrix;

struct SLearn {
    CCountMatrix*	    m_pMatCounts;
    const CGenes*	    m_pGenes;
    const CGenes*           m_pUbikGenes;
    const CDataPair*	    m_pAnswers;
    const CDatFilter*	    m_pDat;
    size_t	            m_iZero;
    CRegularize*	    m_pRegularize;
    size_t		    m_iDat;
    bool                    m_bInPos;
    bool                    m_bInNeg;
    bool                    m_bBridgePos;
    bool                    m_bBridgeNeg;
    bool                    m_bOutPos;
    bool                    m_bOutNeg;
};

struct SEvaluate {
    const CBayesNetMinimal*	m_pBN;
    const CDataPair*		m_pDat;
    const CGenes*		m_pGenes;
    CDat*			m_pYes;
    CDat*			m_pNo;
    size_t			m_iZero;
    size_t			m_iNode;
    const vector<size_t>*	m_pveciGenes;
    bool			m_fFirst;
    string			m_strName;
    bool			m_bLogratio;
};

struct SEvaluate2 {
    size_t				m_iBegin;
    size_t				m_iEnd;
    const CBayesNetMinimal*		m_pBN;
    const vector<CBayesNetMinimal*>*	m_pvecpBNs;
    const CDataPair*			m_pDat;
    CDat*				m_pYes;
    CDat*				m_pNo;
    size_t				m_iZero;
    const vector<size_t>*		m_pveciGenes;
    string				m_strName;
    const vector<size_t>*		m_pveciBNs;
    size_t				m_iNode;
    pthread_spinlock_t			m_sLock;
};

class CRegularize {
public:
    ~CRegularize( ) {

        Clear( );
    }

    bool Open( const char* szFile, CGenome& Genome, const CDat& Answers ) {
        vector<CGenes*>						vecpGenes;
        size_t								i, j, k;
        map<string, size_t>					mapstriGenes;
        map<string, size_t>::const_iterator	iterGene;

        Clear( );
        if( !Genome.Open( szFile, vecpGenes ) )
            return false;

        m_vecsetiGenes.resize( vecpGenes.size( ) );
        m_vecveciGenes.resize( Answers.GetGenes( ) );
        for( i = 0; i < m_vecsetiGenes.size( ); ++i ) {
            const CGenes&	Genes	= *vecpGenes[i];

            for( j = 0; j < Genes.GetGenes( ); ++j ) {
                const string&	strGene	= Genes.GetGene( j ).GetName( );

                if( ( iterGene = mapstriGenes.find( strGene ) ) == mapstriGenes.end( ) )
                    mapstriGenes[strGene] = k = Answers.GetGene( strGene );
                else
                    k = iterGene->second;
                if( k != -1 ) {
                    m_vecsetiGenes[i].insert( k );
                    m_vecveciGenes[k].push_back( i );
                }
            }
        }
        for( i = 0; i < vecpGenes.size( ); ++i )
            delete vecpGenes[i];

        return true;
    }

    void SaveCounts( ostream& ostm, const vector<string>& vecstrNames ) {
        size_t						iDataset, iValue, iGroup;
        const CFullMatrix<size_t>*	pMat;

        for( iDataset = 0; iDataset < m_vecpMatCounts.size( ); ++iDataset ) {
            if( !( pMat = m_vecpMatCounts[iDataset] ) )
                continue;
            ostm << vecstrNames[iDataset] << endl;
            for( iGroup = 0; iGroup < m_vecsetiGenes.size( ); ++iGroup ) {
                for( iValue = 0; iValue < pMat->GetRows( ); ++iValue )
                    ostm << ( iValue ? "\t" : "" ) << pMat->Get( iValue, iGroup );
                ostm << endl;
            }
        }
    }

    void Save( ostream& ostm, const vector<string>& vecstrNames ) {
        size_t						iDatasetOne, iDatasetTwo, iValueOne, iValueTwo, iGroup, iCountOne, iCountTwo;
        const CFullMatrix<size_t>*	pMatOne;
        const CFullMatrix<size_t>*	pMatTwo;
        vector<double>				vecdOne, vecdTwo, vecdGroups;
        vector<size_t>				veciCounts;
        CFullMatrix<double>			MatJoint;
        double						dOne, dTwo, dJoint, dMI;

        veciCounts.resize( m_vecsetiGenes.size( ) );
        fill( veciCounts.begin( ), veciCounts.end( ), 0 );
        for( iCountOne = iDatasetOne = 0; iDatasetOne < m_vecpMatCounts.size( ); ++iDatasetOne ) {
            if( !( pMatOne = m_vecpMatCounts[iDatasetOne] ) )
                continue;
            for( iValueOne = 0; iValueOne < pMatOne->GetRows( ); ++iValueOne )
                for( iGroup = 0; iGroup < pMatOne->GetColumns( ); ++iGroup ) {
                    iCountTwo = pMatOne->Get( iValueOne, iGroup );
                    veciCounts[iGroup] += iCountTwo;
                    iCountOne += iCountTwo;
                }
        }
        vecdGroups.resize( veciCounts.size( ) );
        for( iGroup = 0; iGroup < vecdGroups.size( ); ++iGroup )
            vecdGroups[iGroup] = (double)veciCounts[iGroup] / iCountOne;

        for( iDatasetOne = 0; iDatasetOne < m_vecpMatCounts.size( ); ++iDatasetOne ) {
            if( !( pMatOne = m_vecpMatCounts[iDatasetOne] ) )
                continue;
            /*
            for( iGroup = 0; iGroup < m_vecsetiGenes.size( ); ++iGroup ) {
            for( iCountOne = iValueOne = 0; iValueOne < pMatOne->GetRows( ); ++iValueOne )
            iCountOne += pMatOne->Get( iValueOne, iGroup );
            for( iValueOne = 0; iValueOne < pMatOne->GetRows( ); ++iValueOne )
            cerr << ( iValueOne ? "\t" : "" ) << ( (double)pMatOne->Get( iValueOne, iGroup ) / iCountOne );
            cerr << endl; }
            //*/
            vecdOne.resize( pMatOne->GetRows( ) );
            for( iDatasetTwo = iDatasetOne; iDatasetTwo < m_vecpMatCounts.size( ); ++iDatasetTwo ) {
                if( !( pMatTwo = m_vecpMatCounts[iDatasetTwo] ) )
                    continue;
                vecdTwo.resize( pMatTwo->GetRows( ) );
                MatJoint.Initialize( vecdOne.size( ), vecdTwo.size( ) );

                fill( vecdOne.begin( ), vecdOne.end( ), 0.0f );
                fill( vecdTwo.begin( ), vecdTwo.end( ), 0.0f );
                MatJoint.Clear( );
                for( iGroup = 0; iGroup < m_vecsetiGenes.size( ); ++iGroup ) {
                    for( iCountOne = iValueOne = 0; iValueOne < vecdOne.size( ); ++iValueOne )
                        iCountOne += pMatOne->Get( iValueOne, iGroup );
                    for( iCountTwo = iValueTwo = 0; iValueTwo < vecdTwo.size( ); ++iValueTwo )
                        iCountTwo += pMatTwo->Get( iValueTwo, iGroup );
                    for( iValueOne = 0; iValueOne < vecdOne.size( ); ++iValueOne ) {
                        dOne = (double)pMatOne->Get( iValueOne, iGroup ) / ( iCountOne ? iCountOne : 1 );
                        vecdOne[iValueOne] += dOne * vecdGroups[iGroup];
                        for( iValueTwo = 0; iValueTwo < vecdTwo.size( ); ++iValueTwo ) {
                            dTwo = (double)pMatTwo->Get( iValueTwo, iGroup ) / ( iCountTwo ? iCountTwo : 1 );
                            if( !iValueOne )
                                vecdTwo[iValueTwo] += dTwo * vecdGroups[iGroup];
                            MatJoint.Get( iValueOne, iValueTwo ) += dOne * dTwo * vecdGroups[iGroup];
                        }
                    }
                }
                /*
                				if( iDatasetOne == iDatasetTwo ) {
                					MatJoint.Clear( );
                					for( iValueOne = 0; iValueOne < vecdOne.size( ); ++iValueOne )
                						MatJoint.Set( iValueOne, iValueOne, vecdOne[iValueOne] ); }
                //*/
                /*
                cerr << "One: " << vecstrNames[iDatasetOne] << endl;
                for( iValueOne = 0; iValueOne < vecdOne.size( ); ++iValueOne )
                cerr << ( iValueOne ? "\t" : "" ) << vecdOne[iValueOne];
                cerr << endl;
                cerr << "Two: " << vecstrNames[iDatasetTwo] << endl;
                for( iValueTwo = 0; iValueTwo < vecdTwo.size( ); ++iValueTwo )
                cerr << ( iValueTwo ? "\t" : "" ) << vecdTwo[iValueTwo];
                cerr << endl;
                cerr << "Joint:" << endl;
                for( iValueOne = 0; iValueOne < vecdOne.size( ); ++iValueOne ) {
                for( iValueTwo = 0; iValueTwo < vecdTwo.size( ); ++iValueTwo )
                cerr << ( iValueTwo ? "\t" : "" ) << MatJoint.Get( iValueOne, iValueTwo );
                cerr << endl; }
                //*/

                for( dMI = 0,iValueOne = 0; iValueOne < vecdOne.size( ); ++iValueOne ) {
                    dOne = vecdOne[iValueOne];
                    for( iValueTwo = 0; iValueTwo < vecdTwo.size( ); ++iValueTwo )
                        if( dJoint = MatJoint.Get( iValueOne, iValueTwo ) )
                            dMI += dJoint * ( dJoint ? log( dJoint / dOne / vecdTwo[iValueTwo] ) : 0 );
                }
//cerr << "MI: " << dMI << endl;
                dMI = ( dMI < 0 ) ? 0 : ( dMI / log( 2.0f ) );

                ostm << vecstrNames[iDatasetOne] << '\t' << vecstrNames[iDatasetTwo] << '\t' << dMI << endl;
            }
        }
    }

    void Add( size_t iDataset, const CDatFilter& Dat, size_t iOne, size_t iTwo, size_t iValue ) {
        size_t					i;
        CFullMatrix<size_t>*	pMat;

        if( m_vecsetiGenes.empty( ) )
            return;

        while( m_vecpMatCounts.size( ) <= iDataset )
            m_vecpMatCounts.push_back( NULL );
        if( !( pMat = m_vecpMatCounts[iDataset] ) ) {
            m_vecpMatCounts[iDataset] = pMat = new CFullMatrix<size_t>( );
            pMat->Initialize( Dat.GetValues( ), m_vecsetiGenes.size( ) );
            pMat->Clear( );
        }
        for( i = 0; i < m_vecveciGenes[iOne].size( ); ++i )
            pMat->Get( iValue, m_vecveciGenes[iOne][i] )++;
        for( i = 0; i < m_vecveciGenes[iTwo].size( ); ++i )
            pMat->Get( iValue, m_vecveciGenes[iTwo][i] )++;
    }

private:
    vector<set<size_t> >			m_vecsetiGenes;
    vector<vector<size_t> >			m_vecveciGenes;
    vector<CFullMatrix<size_t>*>	m_vecpMatCounts;

    void Clear( ) {
        size_t					i;
        CFullMatrix<size_t>*	pMat;

        for( i = 0; i < m_vecpMatCounts.size( ); ++i )
            if( pMat = m_vecpMatCounts[i] )
                delete pMat;
        m_vecpMatCounts.clear( );
        m_vecsetiGenes.clear( );
        m_vecveciGenes.clear( );
    }
};

void* learn( void* );
void* evaluate( void* );
void* evaluate2( void* );
void* finalize( void* );
int main_count( const gengetopt_args_info&, const map<string, size_t>&, const CGenes&, const CGenes&,
                const CGenes&, const CGenes&, const CGenes& );
int main_xdsls( const gengetopt_args_info&, const map<string, size_t>&, const map<string, size_t>&,
                const vector<string>& );
int main_inference( const gengetopt_args_info&, const map<string, size_t>&, const map<string, size_t>& );
int main_inference2( const gengetopt_args_info&, const map<string, size_t>&, const map<string, size_t>& );

int main( int iArgs, char** aszArgs ) {
    gengetopt_args_info	sArgs;
    map<string, size_t>	mapstriZeros, mapstriDatasets;
    vector<string>		vecstrContexts;
    int					iRet;
    size_t				i;
    CGenome				Genome;
    CGenes				GenesIn( Genome ), GenesEx( Genome ), GenesEd( Genome ), GenesTm( Genome ), GenesUbik( Genome);

#ifdef WIN32
    pthread_win32_process_attach_np( );
#endif // WIN32
    if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
        cmdline_parser_print_help( );
        return 1;
    }
    CMeta Meta( sArgs.verbosity_arg );
    if( sArgs.pseudocounts_arg < 0 )
        sArgs.pseudocounts_arg = CMeta::GetNaN( );

    if( sArgs.zeros_arg ) {
        ifstream		ifsm;
        vector<string>	vecstrLine;
        char			acLine[ 1024 ];

        ifsm.open( sArgs.zeros_arg );
        if( !ifsm.is_open( ) ) {
            cerr << "Couldn't open: " << sArgs.zeros_arg << endl;
            return 1;
        }
        while( !ifsm.eof( ) ) {
            ifsm.getline( acLine, ARRAYSIZE(acLine) - 1 );
            acLine[ ARRAYSIZE(acLine) - 1 ] = 0;
            vecstrLine.clear( );
            CMeta::Tokenize( acLine, vecstrLine );
            if( vecstrLine.empty( ) )
                continue;
            mapstriZeros[ vecstrLine[ 0 ] ] = atoi( vecstrLine[ 1 ].c_str( ) );
        }
    }

    if( sArgs.datasets_arg ) {
        ifstream		ifsm;
        char			acLine[ 1024 ];
        vector<string>	vecstrLine, vecstrDatasets;

        ifsm.open( sArgs.datasets_arg );
        if( !ifsm.is_open( ) ) {
            cerr << "Could not open: " << sArgs.datasets_arg << endl;
            return 1;
        }
        while( !ifsm.eof( ) ) {
            ifsm.getline( acLine, ARRAYSIZE(acLine) - 1 );
            acLine[ ARRAYSIZE(acLine) - 1 ] = 0;
            if( !acLine[ 0 ] )
                continue;
            vecstrLine.clear( );
            CMeta::Tokenize( acLine, vecstrLine );
            if( vecstrLine.size( ) < 2 ) {
                cerr << "Illegal datasets line: " << acLine << endl;
                return 1;
            }
            if( ( iRet = atol( vecstrLine[ 0 ].c_str( ) ) ) <= 0 ) {
                cerr << "Dataset indices must be greater than zero: " << acLine << endl;
                return 1;
            }
            mapstriDatasets[ vecstrLine[ 1 ] ] = --iRet;
            if( vecstrDatasets.size( ) <= (size_t)iRet )
                vecstrDatasets.resize( iRet + 1 );
            vecstrDatasets[ iRet ] = vecstrLine[ 1 ];
        }
        ifsm.close( );
        for( i = 0; i < vecstrDatasets.size( ); ++i )
            if( vecstrDatasets[ i ].empty( ) ) {
                cerr << "Missing dataset: " << ( i + 1 ) << endl;
                return 1;
            }
    }

    if( sArgs.contexts_arg ) {
        ifstream		ifsm;
        char			acLine[ 2048 ];
        vector<string>	vecstrLine;

        ifsm.open( sArgs.contexts_arg );
        if( !ifsm.is_open( ) ) {
            cerr << "Could not open: " << sArgs.contexts_arg << endl;
            return 1;
        }
        while( !ifsm.eof( ) ) {
            ifsm.getline( acLine, ARRAYSIZE(acLine) - 1 );
            acLine[ ARRAYSIZE(acLine) - 1 ] = 0;
            if( !acLine[ 0 ] )
                continue;
            vecstrLine.clear( );
            CMeta::Tokenize( acLine, vecstrLine );
            if( vecstrLine.size( ) < 3 ) {
                cerr << "Illegal contexts line: " << acLine << endl;
                return 1;
            }
            if( ( atoi( vecstrLine[ 0 ].c_str( ) ) ) != ( vecstrContexts.size( ) + 1 ) ) {
                cerr << "Inconsistent context ID: " << acLine << endl;
                return 1;
            }
            vecstrContexts.push_back( vecstrLine[ 2 ] );
        }
        ifsm.close( );
    }

    if( sArgs.genes_arg ) {
        if( !GenesIn.Open( sArgs.genes_arg ) ) {
            cerr << "Could not open: " << sArgs.genes_arg << endl;
            return 1;
        }
    }
    if( sArgs.genex_arg ) {
        if( !GenesEx.Open( sArgs.genex_arg ) ) {
            cerr << "Could not open: " << sArgs.genex_arg << endl;
            return 1;
        }
    }
    if( sArgs.genet_arg ) {
        if( !GenesTm.Open( sArgs.genet_arg ) ) {
            cerr << "Could not open: " << sArgs.genet_arg << endl;
            return 1;
        }
    }
    if( sArgs.genee_arg ) {
        if( !GenesEd.Open( sArgs.genee_arg ) ) {
            cerr << "Could not open: " << sArgs.genee_arg << endl;
            return 1;
        }
    }
    if( sArgs.ubiqg_arg ) {
        if( !GenesUbik.Open( sArgs.ubiqg_arg ) ) {
            cerr << "Could not open: " << sArgs.ubiqg_arg << endl;
            return 1;
        }
    }

    if( sArgs.answers_arg )
        iRet = main_count( sArgs, mapstriZeros, GenesIn, GenesEx, GenesTm, GenesEd, GenesUbik );
    else if( sArgs.counts_arg )
        iRet = main_xdsls( sArgs, mapstriZeros, mapstriDatasets, vecstrContexts );
    else if( sArgs.networks_arg )
        iRet = sArgs.genewise_flag ? main_inference2( sArgs, mapstriZeros, mapstriDatasets ) :
               main_inference( sArgs, mapstriZeros, mapstriDatasets );

#ifdef WIN32
    pthread_win32_process_detach_np( );
#endif // WIN32
    return iRet;
}

int main_count( const gengetopt_args_info& sArgs, const map<string, size_t>& mapstriZeros,
                const CGenes& GenesIn, const CGenes& GenesEx, const CGenes& GenesTm, const CGenes& GenesEd, const CGenes& GenesUbik ) {
    size_t				i, j, k, m, iTerm, iThread;
    vector<vector<CCountMatrix*>* >	vecpvecpMats;
    vector<CCountMatrix*>		vecpMatRoots;
    vector<CGenes*>			vecpGenes;
    CDataPair				Answers, Dat;
    CDatFilter				Filter, FilterIn, FilterEx, FilterTm, FilterEd;
    CDatFilter*				pFilter;
    string	    			strFile;
    vector<pthread_t>	    		vecpthdThreads;
    vector<SLearn>			vecsData;
    map<string, size_t>::const_iterator	iterZero;
    CGenome			    	Genome;
    vector<string>			vecstrNames;
    CRegularize			        Regularize;
    if( !Answers.Open( sArgs.answers_arg, false, !!sArgs.memmap_flag ) ) {
        cerr << "Could not open: " << sArgs.answers_arg << endl;
        return 1;
    }
    if( !sArgs.inputs_num && sArgs.reggroups_arg && !Regularize.Open( sArgs.reggroups_arg, Genome, Answers ) ) {
        cerr << "Could not open: " << sArgs.reggroups_arg << endl;
        return 1;
    }

    vecpGenes.resize( sArgs.inputs_num );
    for( i = 0; i < vecpGenes.size( ); ++i ) {
        ifstream	ifsm;

        vecpGenes[ i ]  = new CGenes( Genome );
        ifsm.open( sArgs.inputs[ i ] );
        if( !vecpGenes[ i ]->Open( ifsm ) ) {
            cerr << "Couldn't open: " << sArgs.inputs[ i ] << endl;
            return 1;
        }
    }
    if( !vecpGenes.size( ) ) {
        vecpGenes.insert( vecpGenes.begin( ), new CGenes( Genome ) );
        vecpGenes[ 0 ]->Open( Answers.GetGeneNames( ) );
    }
    vecpMatRoots.resize( vecpGenes.size( ) );
    vecpthdThreads.resize( vecpMatRoots.size( ) );
    vecsData.resize( vecpthdThreads.size( ) );
    for( iTerm = 0; iTerm < vecpMatRoots.size( ); iTerm += iThread ) {
        cerr << "Learning root " << iTerm << '/' << vecpMatRoots.size( ) << endl;
        for( iThread = 0; ( ( sArgs.threads_arg == -1 ) || ( iThread < (size_t)sArgs.threads_arg ) ) &&
                ( ( iTerm + iThread ) < vecpMatRoots.size( ) ); ++iThread ) {
            i = iTerm + iThread;
            vecsData[ i ].m_pMatCounts = vecpMatRoots[ i ] = new CCountMatrix( );
            vecsData[ i ].m_pDat = NULL;
            vecsData[ i ].m_iDat = -1;
            vecsData[ i ].m_pGenes = vecpGenes[ i ];
	    vecsData[ i ].m_pUbikGenes = &GenesUbik;
            vecsData[ i ].m_pAnswers = &Answers;
            vecsData[ i ].m_iZero = -1;
            vecsData[ i ].m_pRegularize = &Regularize;
            vecsData[ i ].m_bInPos = sArgs.ctxtpos_flag;
            vecsData[ i ].m_bInNeg = sArgs.ctxtneg_flag;
            vecsData[ i ].m_bBridgePos = sArgs.bridgepos_flag;
            vecsData[ i ].m_bBridgeNeg = sArgs.bridgeneg_flag;
            vecsData[ i ].m_bOutPos = sArgs.outpos_flag;
            vecsData[ i ].m_bOutNeg = sArgs.outneg_flag;
            if( pthread_create( &vecpthdThreads[ i ], NULL, learn, &vecsData[ i ] ) ) {
                cerr << "Couldn't create root thread: " << sArgs.inputs[ i ] << endl;
                return 1;
            }
        }
        for( i = 0; i < iThread; ++i )
            pthread_join( vecpthdThreads[ iTerm + i ], NULL );
    }
    FOR_EACH_DIRECTORY_FILE((string)sArgs.directory_arg, strFile)
    string					strName;
    vector<CCountMatrix*>*	pvecpMatCounts;
    
    if( CMeta::IsExtension( strFile, c_acDab ) ) {
        i = strFile.rfind( '.' );
        strName = (string) sArgs.directory_arg + "/" + strFile.substr( 0, i ) + c_acDab;
    } else if( CMeta::IsExtension( strFile, c_acQDab ) ) {
        i = strFile.rfind( '.' );
        strName = (string) sArgs.directory_arg + "/" + strFile.substr( 0, i ) + c_acQDab;
    } else if( CMeta::IsExtension( strFile, c_acDat ) ) {
	i = strFile.rfind( '.' );
	strName = (string) sArgs.directory_arg + "/" + strFile.substr( 0, i ) + c_acDat;
    } else {
        continue;
    }

    if( !(Dat.Open( strName.c_str( ), false, !!sArgs.memmap_flag )) ) {
        cerr << "Couldn't open: " << strName << endl;
        return 1;
    }
    cerr << "Processing: " << strName << endl;
    strName = CMeta::Filename( CMeta::Deextension( CMeta::Basename( strName.c_str( ) ) ) );
    vecstrNames.push_back( strName );

    Filter.Attach( Dat );
    pFilter = &Filter;
    if( GenesIn.GetGenes( ) ) {
        FilterIn.Attach( *pFilter, GenesIn, CDat::EFilterInclude, &Answers );
        pFilter = &FilterIn;
    }
    if( GenesEx.GetGenes( ) ) {
        FilterEx.Attach( *pFilter, GenesEx, CDat::EFilterExclude, &Answers );
        pFilter = &FilterEx;
    }
    if( GenesTm.GetGenes( ) ) {
        FilterTm.Attach( *pFilter, GenesTm, CDat::EFilterTerm, &Answers );
        pFilter = &FilterTm;
    }
    if( GenesEd.GetGenes( ) ) {
        FilterEd.Attach( *pFilter, GenesEd, CDat::EFilterEdge, &Answers );
        pFilter = &FilterEd;
    }

    pvecpMatCounts = new vector<CCountMatrix*>( );
    pvecpMatCounts->resize( vecpGenes.size( ) );
    for( iTerm = 0; iTerm < vecpMatRoots.size( ); iTerm += iThread ) {
        cerr << "Learning term " << iTerm << '/' << vecpMatRoots.size( ) << endl;
        for( iThread = 0; ( ( sArgs.threads_arg == -1 ) || ( iThread < (size_t)sArgs.threads_arg ) ) &&
                ( ( iTerm + iThread ) < vecpMatRoots.size( ) ); ++iThread ) {
            i = iTerm + iThread;
            vecsData[ i ].m_pMatCounts = (*pvecpMatCounts)[ i ] = new CCountMatrix( );
            vecsData[ i ].m_pDat = pFilter;
            vecsData[ i ].m_iDat = vecstrNames.size( ) - 1;
            vecsData[ i ].m_pGenes = vecpGenes[ i ];
	    vecsData[ i ].m_pUbikGenes = &GenesUbik;
            vecsData[ i ].m_pAnswers = &Answers;
            vecsData[ i ].m_iZero = ( ( iterZero = mapstriZeros.find( strName ) ) == mapstriZeros.end( ) ) ? -1 : iterZero->second;
            vecsData[ i ].m_pRegularize = &Regularize;
            vecsData[ i ].m_bInPos = sArgs.ctxtpos_flag;
            vecsData[ i ].m_bInNeg = sArgs.ctxtneg_flag;
            vecsData[ i ].m_bBridgePos = sArgs.bridgepos_flag;
            vecsData[ i ].m_bBridgeNeg = sArgs.bridgeneg_flag;
            vecsData[ i ].m_bOutPos = sArgs.outpos_flag;
            vecsData[ i ].m_bOutNeg = sArgs.outneg_flag;
            if( pthread_create( &vecpthdThreads[ i ], NULL, learn, &vecsData[ i ] ) ) {
                cerr << "Couldn't create root thread: " << sArgs.inputs[ i ] << endl;
                return 1;
            }
        }
        for( i = 0; i < iThread; ++i )
            pthread_join( vecpthdThreads[ iTerm + i ], NULL );
    }
    vecpvecpMats.push_back( pvecpMatCounts );
}
#ifdef _MSC_VER
FindClose( hSearch );
#else // _MSC_VER
closedir( pDir );
#endif // _MSC_VER

for( i = 0; i < vecpMatRoots.size( ); ++i ) {
    ofstream	ofsm;

    ofsm.open( ( (string)sArgs.output_arg + '/' + ( sArgs.inputs ?
                 CMeta::Deextension( CMeta::Basename( sArgs.inputs[ i ] ) ) : "global" ) + c_acTxt ).c_str( ) );
    ofsm << ( sArgs.inputs ? CMeta::Deextension( CMeta::Basename( sArgs.inputs[ i ] ) ) : "global" ) <<
         '\t' << vecpvecpMats.size( ) << endl;
    for( j = 0; j < vecpMatRoots[ i ]->GetRows( ); ++j )
        ofsm << ( j ? "\t" : "" ) << vecpMatRoots[ i ]->Get( j, 0 );
    ofsm << endl;
    for( j = 0; j < vecpvecpMats.size( ); ++j ) {
        ofsm << vecstrNames[ j ] << endl;
        for( k = 0; k < (*vecpvecpMats[ j ])[ i ]->GetColumns( ); ++k ) {
            for( m = 0; m < (*vecpvecpMats[ j ])[ i ]->GetRows( ); ++m )
                ofsm << ( m ? "\t" : "" ) << (*vecpvecpMats[ j ])[ i ]->Get( m, k );
            ofsm << endl;
        }
    }
}
if( sArgs.reggroups_arg )
    Regularize.Save( cout, vecstrNames );

for( i = 0; i < vecpvecpMats.size( ); ++i ) {
    for( j = 0; j < vecpvecpMats[ i ]->size( ); ++j )
        delete (*vecpvecpMats[ i ])[ j ];
    delete vecpvecpMats[ i ];
}
for( i = 0; i < vecpMatRoots.size( ); ++i ) {
    delete vecpMatRoots[ i ];
    delete vecpGenes[ i ];
}

return 0;
}

void* learn( void* pData ) {
    SLearn*		psData;
    size_t		i, j, iAnswer, iVal, iOne, iTwo;
    vector<bool>	vecfGenes, vecfUbik;
    vector<size_t>	veciGenes;

    psData = (SLearn*)pData;
    if (psData->m_pUbikGenes->GetGenes( )) {
	vecfUbik.resize( psData->m_pAnswers->GetGenes( ) );
	for( i = 0; i < vecfUbik.size( ); ++i) {
	    vecfUbik[ i ] = psData->m_pUbikGenes->IsGene( psData->m_pAnswers->GetGene( i ) );
	}
    }
    vecfGenes.resize( psData->m_pAnswers->GetGenes( ) );
    for( i = 0; i < vecfGenes.size( ); ++i )
        vecfGenes[ i ] = psData->m_pGenes->IsGene( psData->m_pAnswers->GetGene( i ) );
    if( psData->m_pDat ) {
        psData->m_pMatCounts->Initialize( psData->m_pDat->GetValues( ), psData->m_pAnswers->GetValues( ) );
        veciGenes.resize( psData->m_pAnswers->GetGenes( ) );
        for( i = 0; i < veciGenes.size( ); ++i )
            veciGenes[ i ] = psData->m_pDat->GetGene( psData->m_pAnswers->GetGene( i ) );
    }
    else
        psData->m_pMatCounts->Initialize( psData->m_pAnswers->GetValues( ), 1 );
    psData->m_pMatCounts->Clear( );
    for( i = 0; i < psData->m_pAnswers->GetGenes( ); ++i ) {
        if( psData->m_pDat )
            iOne = veciGenes[ i ];
        for( j = ( i + 1 ); j < psData->m_pAnswers->GetGenes( ); ++j ) {
	    iAnswer = psData->m_pAnswers->Quantize( psData->m_pAnswers->Get( i, j ) );
            if( iAnswer == -1 ) {
                continue;
            }
	    
	    if ( CMeta::SkipEdge( !!iAnswer, i, j, vecfGenes, vecfUbik, psData->m_bInPos, psData->m_bInNeg, psData->m_bBridgePos, psData->m_bBridgeNeg, psData->m_bOutPos, psData->m_bOutNeg ) ) {
		continue;
	    }

            if( psData->m_pDat ) {
                iTwo = veciGenes[ j ];
                iVal = -1;
                if( ( iOne != -1 ) && ( iTwo != -1 ) )
                    iVal = psData->m_pDat->Quantize( psData->m_pDat->Get( iOne, iTwo ) );
                if( iVal == -1 )
                    iVal = psData->m_iZero;
                if( iVal == -1 )
                    continue;
                psData->m_pMatCounts->Get( iVal, iAnswer )++;
                psData->m_pRegularize->Add( psData->m_iDat, *psData->m_pDat, i, j, iVal );
            }
            else
                psData->m_pMatCounts->Get( iAnswer, 0 )++;
        }
    }

    return NULL;
}

int main_xdsls( const gengetopt_args_info& sArgs, const map<string, size_t>& mapstriZeros,
                const map<string, size_t>& mapstriDatasets, const vector<string>& vecstrContexts ) {
    static const size_t	c_iBuffer	= 1024;
    char						szBuffer[ c_iBuffer ];
    string						strFile;
    CBayesNetMinimal			BNDefault;
    CBayesNetMinimal*			pBN;
    vector<CBayesNetMinimal*>	vecpBNs;
    ifstream					ifsm;
    vector<string>				vecstrLine;
    ofstream					ofsm;
    uint32_t					iSize;
    size_t						i;
    vector<float>				vecdAlphas;
    vector<unsigned char>		vecbZeros;
    map<string, size_t>::const_iterator	iterZero, iterDataset;
    float						d;

    if( mapstriDatasets.empty( ) ) {
        cerr << "No datasets given" << endl;
        return 1;
    }
    if( sArgs.alphas_arg ) {
        vecdAlphas.resize( mapstriDatasets.size( ) );
        ifsm.clear( );
        ifsm.open( sArgs.alphas_arg );
        if( !ifsm.is_open( ) ) {
            cerr << "Could not open: " << sArgs.alphas_arg << endl;
            return 1;
        }
        while( !ifsm.eof( ) ) {
            ifsm.getline( szBuffer, c_iBuffer - 1 );
            szBuffer[ c_iBuffer - 1 ] = 0;
            if( !szBuffer[ 0 ] )
                continue;
            vecstrLine.clear( );
            CMeta::Tokenize( szBuffer, vecstrLine );
            if( vecstrLine.size( ) != 2 ) {
                cerr << "Illegal alphas line: " << szBuffer << endl;
                return 1;
            }
            if( ( iterDataset = mapstriDatasets.find( vecstrLine[ 0 ] ) ) == mapstriDatasets.end( ) )
                cerr << "Dataset in counts but not database: " << vecstrLine[ 0 ] << endl;
            else
                vecdAlphas[ iterDataset->second ] = (float)atof( vecstrLine[ 1 ].c_str( ) );
        }
        ifsm.close( );
    }

    vecbZeros.resize( mapstriDatasets.size( ) );
    fill( vecbZeros.begin( ), vecbZeros.end( ), 0xFF );
    for( iterZero = mapstriZeros.begin( ); iterZero != mapstriZeros.end( ); ++iterZero )
        if( ( iterDataset = mapstriDatasets.find( iterZero->first ) ) == mapstriDatasets.end( ) )
            cerr << "Unknown dataset in zeros file: " << iterZero->first << endl;
        else
            vecbZeros[ iterDataset->second ] = (unsigned char)iterZero->second;

    if( !BNDefault.OpenCounts( sArgs.default_arg, mapstriDatasets, vecbZeros, vecdAlphas,
                               sArgs.pseudocounts_arg ) ) {
        cerr << "Could not open default counts: " << ( sArgs.default_arg ? sArgs.default_arg : "not given" ) <<
             endl;
        return 1;
    }
    if( sArgs.regularize_flag ) {
        d = BNDefault.Regularize( vecdAlphas );
        BNDefault.OpenCounts( sArgs.default_arg, mapstriDatasets, vecbZeros, vecdAlphas, d );
    }

    for( i = 0; i < vecstrContexts.size( ); ++i ) {
        strFile = (string)sArgs.counts_arg + '/' + CMeta::Filename( vecstrContexts[ i ] ) + c_acTxt;
        cerr << "Processing: " << strFile << endl;
        pBN = new CBayesNetMinimal( );
        if( !pBN->OpenCounts( strFile.c_str( ), mapstriDatasets, vecbZeros, vecdAlphas, sArgs.pseudocounts_arg,
                              &BNDefault ) )
            return 1;
        if( sArgs.regularize_flag ) {
            d = pBN->Regularize( vecdAlphas );
            pBN->OpenCounts( strFile.c_str( ), mapstriDatasets, vecbZeros, vecdAlphas, d, &BNDefault );
        }
        vecpBNs.push_back( pBN );
    }

    cerr << "Created " << ( vecpBNs.size( ) + 1 ) << " Bayesian classifiers" << endl;

    if( sArgs.output_arg ) {
        if( sArgs.smile_flag ) {
            CBayesNetSmile						BNSmile;
            vector<string>						vecstrNames;
            map<string, size_t>::const_iterator	iterName;

            vecstrNames.resize( mapstriDatasets.size( ) );
            for( iterName = mapstriDatasets.begin( ); iterName != mapstriDatasets.end( ); ++iterName )
                vecstrNames[ iterName->second ] = iterName->first;

            BNSmile.Open( BNDefault, vecstrNames );
            BNSmile.Save( ( (string)sArgs.output_arg + '/' + BNDefault.GetID( ) +
                            ( sArgs.xdsl_flag ? ".x" : "." ) + "dsl" ).c_str( ) );
            for( i = 0; i < vecpBNs.size( ); ++i ) {
                BNSmile.Open( *vecpBNs[ i ], vecstrNames );
                BNSmile.Save( ( (string)sArgs.output_arg + '/' + vecpBNs[ i ]->GetID( ) +
                                ( sArgs.xdsl_flag ? ".x" : "." ) + "dsl" ).c_str( ) );
            }
        }
        else {
            ofsm.open( sArgs.output_arg, ios_base::binary );
            BNDefault.Save( ofsm );
            iSize = (uint32_t)vecpBNs.size( );
            ofsm.write( (const char*)&iSize, sizeof(iSize) );
            for( i = 0; i < vecpBNs.size( ); ++i )
                vecpBNs[ i ]->Save( ofsm );
            ofsm.close( );
        }
    }
    for( i = 0; i < vecpBNs.size( ); ++i )
        delete vecpBNs[ i ];

    return 0;
}

int main_inference( const gengetopt_args_info& sArgs, const map<string, size_t>& mapstriZeros,
                    const map<string, size_t>& mapstriDatasets ) {
    map<string, size_t>::const_iterator	iterDataset;
    vector<size_t>						veciGenes;
    size_t								i, j, iTerm, iThread;
    vector<CGenes*>						vecpGenes;
    vector<CDat*>						vecpYes, vecpNo;
    vector<string>						vecstrTmps;
    CGenome								Genome;
    vector<SEvaluate>					vecsData;
    CBayesNetMinimal					BNDefault;
    vector<CBayesNetMinimal*>			vecpBNs;
    ifstream							ifsm;
    uint32_t							iSize;
    map<string, size_t>::const_iterator	iterZero;
    vector<pthread_t>					vecpthdThreads;
    bool								fFirst;

    ifsm.open( sArgs.networks_arg, ios_base::binary );
    if( !BNDefault.Open( ifsm ) ) {
        cerr << "Could not open: " << sArgs.networks_arg << endl;
        return 1;
    }
    ifsm.read( (char*)&iSize, sizeof(iSize) );
    vecpBNs.resize( iSize );
    for( i = 0; i < vecpBNs.size( ); ++i ) {
        vecpBNs[ i ] = new CBayesNetMinimal( );
        if( !vecpBNs[ i ]->Open( ifsm ) ) {
            cerr << "Could not open: " << sArgs.networks_arg << endl;
            return 1;
        }
    }
    ifsm.close( );

    /*
    size_t k, m;
    for( i = 0; i < BNDefault.GetNodes( ); ++i )
    for( j = 0; j < BNDefault.GetCPT( i ).GetRows( ); ++j )
    for( k = 0; k < BNDefault.GetCPT( i ).GetColumns( ); ++k )
    cout << i << '\t' << j << '\t' << k << '\t' << BNDefault.GetCPT( i ).Get( j, k ) << endl;
    for( m = 0; m < vecpBNs.size( ); ++m )
    for( i = 0; i < vecpBNs[ m ]->GetNodes( ); ++i )
    for( j = 0; j < vecpBNs[ m ]->GetCPT( i ).GetRows( ); ++j )
    for( k = 0; k < vecpBNs[ m ]->GetCPT( i ).GetColumns( ); ++k )
    cout << i << '\t' << j << '\t' << k << '\t' << vecpBNs[ m ]->GetCPT( i ).Get( j, k ) << endl;
    return 0;
    //*/

    if( !sArgs.genome_arg ) {
        cerr << "No genes given" << endl;
        return 1;
    }
    {
        CPCL	PCLGenes( false );

        if( !PCLGenes.Open( sArgs.genome_arg, 1 ) ) {
            cerr << "Could not open: " << sArgs.genome_arg << endl;
            return 1;
        }
        for( i = 0; i < PCLGenes.GetGenes( ); ++i )
            Genome.AddGene( PCLGenes.GetFeature( i, 1 ) );
    }

    vecpGenes.resize( sArgs.inputs_num ? sArgs.inputs_num : 1 );
    vecpYes.resize( vecpGenes.size( ) );
    vecpNo.resize( vecpGenes.size( ) );
    vecstrTmps.resize( vecpNo.size( ) );
    for( i = 0; i < vecpGenes.size( ); ++i ) {
        char*	szTemp;

        vecpGenes[ i ]  = new CGenes( Genome );
        if( sArgs.inputs_num ) {
            ifstream	ifsm;

            ifsm.open( sArgs.inputs[ i ] );
            if( !vecpGenes[ i ]->Open( ifsm, false ) ) {
                cerr << "Couldn't open: " << sArgs.inputs[ i ] << endl;
                return 1;
            }
        }
        else
            vecpGenes[ i ]->Open( Genome.GetGeneNames( ), false );
        vecpYes[ i ] = new CDat( );
        vecpNo[ i ] = new CDat( );
        if( sArgs.memmapout_flag ) {
            vecpYes[ i ]->Open( Genome.GetGeneNames( ), false,
                                ((string) sArgs.output_arg + '/' + (sArgs.inputs_num ? CMeta::Basename(
                                            sArgs.inputs[ i ] ) : "global") + c_acDab).c_str( ) );
            /*FIXME: _tempnam can result in a filename that ends up in use later via race condition */
            if( !(szTemp = _tempnam( sArgs.temporary_arg, NULL )) ) {
                cerr << "Could not generate temporary file name in: " << sArgs.temporary_arg << endl;
                return 1;
            }
            cout << "File: " << szTemp << endl;
            vecstrTmps[ i ] = szTemp;
            free( szTemp );
            vecpNo[ i ]->Open( Genome.GetGeneNames( ), false, vecstrTmps[ i ].c_str( ) );
        } else {
            vecpYes[ i ]->Open( Genome.GetGeneNames( ), false, NULL );
            vecpNo[ i ]->Open( Genome.GetGeneNames( ), false, NULL );
        }
    }

    veciGenes.resize( vecpYes[ 0 ]->GetGenes( ) );
    vecsData.resize( vecpGenes.size( ) );
    vecpthdThreads.resize( vecsData.size( ) );
    for( fFirst = true,iterDataset = mapstriDatasets.begin( ); iterDataset != mapstriDatasets.end( );
            fFirst = false,++iterDataset ) {
        CDataPair	Dat;
        string		strFile;

        if( !Dat.Open( ( strFile = ( (string)sArgs.directory_arg + '/' + iterDataset->first +
                                     c_acDab ) ).c_str( ), false, !!sArgs.memmap_flag ) && !Dat.Open( ( strFile = ( (string)sArgs.directory_arg + '/' + iterDataset->first +
                                             c_acQDab ) ).c_str( ), false, !!sArgs.memmap_flag ) ) {
            cerr << "Couldn't open: " << strFile << endl;
            return 1;
        }
        cerr << "Processing: " << strFile << endl;
        for( i = 0; i < veciGenes.size( ); ++i )
            veciGenes[ i ] = Dat.GetGene( vecpYes[ 0 ]->GetGene( i ) );
        for( iTerm = 0; iTerm < vecpGenes.size( ); iTerm += iThread ) {
            for( iThread = 0; ( ( sArgs.threads_arg == -1 ) || ( iThread < (size_t)sArgs.threads_arg ) ) &&
                    ( ( iTerm + iThread ) < vecpGenes.size( ) ); ++iThread ) {
                i = iTerm + iThread;
                if( sArgs.inputs_num ) {
                    for( j = 0; j < vecpBNs.size( ); ++j )
                        if( vecpBNs[ j ]->GetID( ) == CMeta::Deextension( CMeta::Basename(
                                    sArgs.inputs[ i ] ) ) )
                            break;
                    if( j >= vecpBNs.size( ) ) {
                        cerr << "Could not locate Bayes net for: " << sArgs.inputs[ i ] << endl;
                        return 1;
                    }
                    vecsData[ i ].m_pBN = vecpBNs[ j ];
                }
                else
                    vecsData[ i ].m_pBN = &BNDefault;
                vecsData[ i ].m_pDat = &Dat;
                vecsData[ i ].m_pGenes = vecpGenes[ i ];
                vecsData[ i ].m_pYes = vecpYes[ i ];
                vecsData[ i ].m_pNo = vecpNo[ i ];
                vecsData[ i ].m_iZero = ( ( iterZero = mapstriZeros.find( iterDataset->first ) ) ==
                                          mapstriZeros.end( ) ) ? vecsData[ i ].m_pBN->GetDefault( iterDataset->second + 1 ) :
                                        iterZero->second;
                if( vecsData[ i ].m_iZero == 0xFF )
                    vecsData[ i ].m_iZero = -1;
                vecsData[ i ].m_iNode = iterDataset->second + 1;
                vecsData[ i ].m_pveciGenes = &veciGenes;
                vecsData[ i ].m_fFirst = fFirst;
                vecsData[ i ].m_strName = sArgs.inputs_num ? sArgs.inputs[ i ] : "global";
                if( pthread_create( &vecpthdThreads[ i ], NULL, evaluate, &vecsData[ i ] ) ) {
                    cerr << "Couldn't create evaluation thread: " << sArgs.inputs[ i ] << endl;
                    return 1;
                }
            }
            for( i = 0; i < iThread; ++i )
                pthread_join( vecpthdThreads[ iTerm + i ], NULL );
        }
    }

    for( iTerm = 0; iTerm < vecpGenes.size( ); iTerm += iThread ) {
        for( iThread = 0; ( ( sArgs.threads_arg == -1 ) || ( iThread < (size_t)sArgs.threads_arg ) ) &&
                ( ( iTerm + iThread ) < vecpGenes.size( ) ); ++iThread ) {
            i = iTerm + iThread;
            if( sArgs.inputs_num ) {
                for( j = 0; j < vecpBNs.size( ); ++j )
                    if( vecpBNs[ j ]->GetID( ) == CMeta::Deextension( CMeta::Basename( sArgs.inputs[ i ] ) ) )
                        break;
                if( j >= vecpBNs.size( ) ) {
                    cerr << "Could not locate Bayes net for: " << sArgs.inputs[ i ] << endl;
                    return 1;
                }
                vecsData[ i ].m_pBN = vecpBNs[ j ];
            }
            else
                vecsData[ i ].m_pBN = &BNDefault;
            vecsData[ i ].m_pYes = vecpYes[ i ];
            vecsData[ i ].m_pNo = vecpNo[ i ];
            vecsData[ i ].m_strName = sArgs.inputs_num ? sArgs.inputs[ i ] : "global";
	    vecsData[ i ].m_bLogratio = sArgs.logratio_flag;
            if( pthread_create( &vecpthdThreads[ i ], NULL, finalize, &vecsData[ i ] ) ) {
                cerr << "Couldn't create finalization thread: " << sArgs.inputs[ i ] << endl;
                return 1;
            }
        }
        for( i = 0; i < iThread; ++i )
            pthread_join( vecpthdThreads[ iTerm + i ], NULL );
    }

    if( sArgs.memmapout_flag ) {
        for( i = 0; i < vecstrTmps.size( ); ++i ) {
            _unlink( vecstrTmps[ i ].c_str( ) );
        }
    }
    for( i = 0; i < vecpBNs.size( ); ++i )
        delete vecpBNs[ i ];
    for( i = 0; i < vecpGenes.size( ); ++i ) {
        delete vecpNo[ i ];
        if( !sArgs.memmapout_flag ) {
            vecpYes[ i ]->Save( ((string) sArgs.output_arg + '/' + (sArgs.inputs_num ? CMeta::Basename(
                                     sArgs.inputs[ i ] ) : "global") + c_acDab).c_str( ) );
        }
        delete vecpYes[ i ];
        delete vecpGenes[ i ];
    }

    return 0;
}

void* evaluate( void* pData ) {
    SEvaluate*	psData;
    size_t		i, j, iOne, iTwo, iBin, iIndex;
    float*		adYes;
    float*		adNo;
    float		dNo, dYes;

    psData = (SEvaluate*)pData;
    if( psData->m_fFirst ) {
        float*	adBuffer;

        adBuffer = new float[ psData->m_pYes->GetGenes( ) ];
        for( i = 0; i < psData->m_pNo->GetGenes( ); ++i )
            adBuffer[ i ] = CMeta::GetNaN( );
        for( i = 0; i < psData->m_pNo->GetGenes( ); ++i ) {
            if( !( i % 1000 ) )
                cerr << "IN: " << psData->m_strName << ", " << i << endl;
            psData->m_pNo->Set( i, adBuffer );
        }
        for( i = 0; i < psData->m_pYes->GetGenes( ); ++i ) {
            if( !( i % 1000 ) )
                cerr << "IY: " << psData->m_strName << ", " << i << endl;
            psData->m_pYes->Set( i, adBuffer );
        }
        delete[] adBuffer;
    }

    dNo = log( psData->m_pBN->GetCPT( 0 ).Get( 0, 0 ) );
    dYes = log( psData->m_pBN->GetCPT( 0 ).Get( 1, 0 ) );
    adYes = new float[ psData->m_pYes->GetGenes( ) ];
    adNo = new float[ psData->m_pNo->GetGenes( ) ];
    for( i = 0; i < psData->m_pYes->GetGenes( ); ++i ) {
        if( !( i % 1000 ) )
            cerr << "C: " << psData->m_strName << ", " << i << endl;
        if( ( ( iOne = (*psData->m_pveciGenes)[ i ] ) == -1 ) && ( psData->m_iZero == -1 ) )
            continue;
        memcpy( adYes, psData->m_pYes->Get( i ), ( psData->m_pYes->GetGenes( ) - i - 1 ) * sizeof(*adYes) );
        memcpy( adNo, psData->m_pNo->Get( i ), ( psData->m_pNo->GetGenes( ) - i - 1 ) * sizeof(*adNo) );
        for( j = ( i + 1 ); j < psData->m_pYes->GetGenes( ); ++j ) {
            if( ( ( iTwo = (*psData->m_pveciGenes)[ j ] ) == -1 ) && ( psData->m_iZero == -1 ) )
                continue;
            if( ( iOne == -1 ) || ( iTwo == -1 ) ||
                    ( ( iBin = psData->m_pDat->Quantize( psData->m_pDat->Get( iOne, iTwo ) ) ) == -1 ) )
                iBin = psData->m_iZero;
            if( iBin == -1 )
                continue;
            if( CMeta::IsNaN( adYes[ iIndex = ( j - i - 1 ) ] ) ) {
                adYes[ iIndex ] = dYes;
                adNo[ iIndex ] = dNo;
            }
            adNo[ iIndex ] += log( psData->m_pBN->GetCPT( psData->m_iNode ).Get( iBin, 0 ) );
            adYes[ iIndex ] += log( psData->m_pBN->GetCPT( psData->m_iNode ).Get( iBin, 1 ) );
        }
        psData->m_pNo->Set( i, adNo );
        psData->m_pYes->Set( i, adYes );
    }
    delete[] adYes;
    delete[] adNo;

    return NULL;
}

void* finalize( void* pData ) {
    SEvaluate*	psData;
    size_t		i, j;
    float*		adYes;
    float*		adNo;
    float		dPrior, dYes, dNo;

    psData = (SEvaluate*)pData;
    dPrior = psData->m_pBN->GetCPT( 0 ).Get( 1, 0 );
    adYes = new float[ psData->m_pYes->GetGenes( ) ];
    adNo = new float[ psData->m_pNo->GetGenes( ) ];

    dYes = log( dPrior );
    dNo = log( psData->m_pBN->GetCPT( 0 ).Get( 0, 0 ) );

    for( i = 0; i < psData->m_pYes->GetGenes( ); ++i ) {
        if( !( i % 1000 ) )
            cerr << "F: " << psData->m_strName << ", " << i << endl;
        memcpy( adYes, psData->m_pYes->Get( i ), ( psData->m_pYes->GetGenes( ) - i - 1 ) * sizeof(*adYes) );
        memcpy( adNo, psData->m_pNo->Get( i ), ( psData->m_pNo->GetGenes( ) - i - 1 ) * sizeof(*adNo) );
        for( j = 0; j < ( psData->m_pYes->GetGenes( ) - i - 1 ); ++j ) {
	    if( psData->m_bLogratio ) {
		adYes[ j ] = CMeta::IsNaN( adYes[ j ] ) ? 0.0 :
		             (float)( (double)(adYes[ j ] - dYes) - (double)(adNo[ j ] - dNo) );
	    }
	    else {
		adYes[ j ] = CMeta::IsNaN( adYes[ j ] ) ? dPrior :
		             (float)( 1 / ( 1 + exp( (double)adNo[ j ] - (double)adYes[ j ] ) ) );
	    }
	}
        psData->m_pYes->Set( i, adYes );
    }
    delete[] adNo;
    delete[] adYes;

    return NULL;
}

////////////////////////

int main_inference2( const gengetopt_args_info& sArgs, const map<string, size_t>& mapstriZeros,
                     const map<string, size_t>& mapstriDatasets ) {
    map<string, size_t>::const_iterator	iterDataset;
    vector<size_t>						veciGenes, veciBNs;
    size_t								i, j, iThread;
    vector<CGenes*>						vecpGenes;
    vector<CDat*>						vecpYes, vecpNo;
    vector<string>						vecstrTmps;
    CGenome								Genome;
    vector<SEvaluate2>					vecsData;
    CBayesNetMinimal					BNDefault;
    vector<CBayesNetMinimal*>			vecpBNs;
    ifstream							ifsm;
    uint32_t							iSize;
    map<string, size_t>::const_iterator	iterZero;
    vector<pthread_t>					vecpthdThreads;
    bool								fFirst;

    ifsm.open( sArgs.networks_arg, ios_base::binary );
    if( !BNDefault.Open( ifsm ) ) {
        cerr << "Could not open: " << sArgs.networks_arg << endl;
        return 1;
    }
    ifsm.read( (char*)&iSize, sizeof(iSize) );
    vecpBNs.resize( iSize );
    for( i = 0; i < vecpBNs.size( ); ++i ) {
        vecpBNs[ i ] = new CBayesNetMinimal( );
        if( !vecpBNs[ i ]->Open( ifsm ) ) {
            cerr << "Could not open: " << sArgs.networks_arg << endl;
            return 1;
        }
    }
    ifsm.close( );

    if( !sArgs.genome_arg ) {
        cerr << "No genes given" << endl;
        return 1;
    }
    {
        CPCL	PCLGenes( false );

        if( !PCLGenes.Open( sArgs.genome_arg, 1 ) ) {
            cerr << "Could not open: " << sArgs.genome_arg << endl;
            return 1;
        }
        for( i = 0; i < PCLGenes.GetGenes( ); ++i )
            Genome.AddGene( PCLGenes.GetFeature( i, 1 ) );
    }

    vecpGenes.resize( sArgs.inputs_num ? sArgs.inputs_num : 1 );
    vecpYes.resize( 1 );
    vecpNo.resize( vecpYes.size( ) );
    vecstrTmps.resize( vecpNo.size( ) );
    for( i = 0; i < vecpYes.size( ); ++i ) {
        char*	szTemp;

        vecpYes[ i ] = new CDat( );
        vecpYes[ i ]->Open( Genome.GetGeneNames( ), false, ( (string)sArgs.output_arg + '/' + "global" + c_acDab ).c_str( ) );
        vecpNo[ i ] = new CDat( );
        if( !( szTemp = _tempnam( sArgs.temporary_arg, NULL ) ) ) {
            cerr << "Could not generate temporary file name in: " << sArgs.temporary_arg << endl;
            return 1;
        }
        vecstrTmps[ i ] = szTemp;
        free( szTemp );
        vecpNo[ i ]->Open( Genome.GetGeneNames( ), false, vecstrTmps[ i ].c_str( ) );
    }
    for( i = 0; i < vecpGenes.size( ); ++i ) {
        vecpGenes[ i ]  = new CGenes( Genome );
        if( sArgs.inputs_num ) {
            ifstream	ifsm;

            ifsm.open( sArgs.inputs[ i ] );
            if( !vecpGenes[ i ]->Open( ifsm, false ) ) {
                cerr << "Couldn't open: " << sArgs.inputs[ i ] << endl;
                return 1;
            }
        }
        else
            vecpGenes[ i ]->Open( Genome.GetGeneNames( ), false );
    }

    {
        map<string, size_t>					mapstriBNs;
        map<string, size_t>::const_iterator	iterBN;
        string								strGene;

        for( i = 0; i < vecpBNs.size( ); ++i )
            mapstriBNs[ vecpBNs[ i ]->GetID( ) ] = i;
        veciBNs.resize( vecpYes[ 0 ]->GetGenes( ) );
        for( i = 0; i < veciBNs.size( ); ++i ) {
            strGene = CMeta::Filename( vecpYes[ 0 ]->GetGene( i ) );
            veciBNs[ i ] = ( ( iterBN = mapstriBNs.find( strGene ) ) == mapstriBNs.end( ) ) ? -1 :
                           iterBN->second;
        }
    }

    for( i = 0; i < vecpYes.size( ); ++i ) {
        float*	adBuffer;

        adBuffer = new float[ vecpYes[ i ]->GetGenes( ) ];
        for( j = 0; j < vecpNo[ i ]->GetGenes( ); ++j )
            adBuffer[ j ] = log( BNDefault.GetCPT( 0 ).Get( 0, 0 ) );
        for( j = 0; j < vecpNo[ i ]->GetGenes( ); ++j ) {
            if( !( j % 1000 ) )
                cerr << "IN: global, " << j << endl;
            vecpNo[ i ]->Set( j, adBuffer );
        }
        for( j = 0; j < vecpYes[ i ]->GetGenes( ); ++j )
            adBuffer[ j ] = log( BNDefault.GetCPT( 0 ).Get( 1, 0 ) );
        for( j = 0; j < vecpYes[ i ]->GetGenes( ); ++j ) {
            if( !( j % 1000 ) )
                cerr << "IY: global, " << j << endl;
            vecpYes[ i ]->Set( j, adBuffer );
        }
        delete[] adBuffer;
    }

    veciGenes.resize( vecpYes[ 0 ]->GetGenes( ) );
    vecsData.resize( ( sArgs.threads_arg == -1 ) ? 1 : sArgs.threads_arg );
    vecpthdThreads.resize( vecsData.size( ) );
    for( fFirst = true,iterDataset = mapstriDatasets.begin( ); iterDataset != mapstriDatasets.end( );
            fFirst = false,++iterDataset ) {
        CDataPair		Dat;
        string			strFile;
        vector<size_t>	veciCur;
        size_t			iZero;

        if( !Dat.Open( ( strFile = ( (string)sArgs.directory_arg + '/' + iterDataset->first +
                                     c_acDab ) ).c_str( ), false, !!sArgs.memmap_flag ) ) {
            cerr << "Couldn't open: " << strFile << endl;
            return 1;
        }
        cerr << "Processing: " << strFile << endl;
        for( i = 0; i < veciGenes.size( ); ++i )
            veciGenes[ i ] = Dat.GetGene( vecpYes[ 0 ]->GetGene( i ) );
        iZero = ( ( iterZero = mapstriZeros.find( iterDataset->first ) ) == mapstriZeros.end( ) ) ?
                BNDefault.GetDefault( iterDataset->second + 1 ) : iterZero->second;
        if( iZero == 0xFF )
            iZero = -1;
        iThread = ( vecpYes[ 0 ]->GetGenes( ) + 1 ) / vecsData.size( );
        for( i = j = 0; i < vecsData.size( ); ++i,j += iThread ) {
            vecsData[ i ].m_iBegin = j;
            vecsData[ i ].m_iEnd = min( vecpYes[ 0 ]->GetGenes( ), j + iThread );
            vecsData[ i ].m_pBN = &BNDefault;
            vecsData[ i ].m_pvecpBNs = &vecpBNs;
            vecsData[ i ].m_pDat = &Dat;
            vecsData[ i ].m_pYes = vecpYes[ 0 ];
            vecsData[ i ].m_pNo = vecpNo[ 0 ];
            vecsData[ i ].m_iZero = iZero;
            vecsData[ i ].m_pveciGenes = &veciGenes;
            vecsData[ i ].m_strName = iterDataset->first;
            vecsData[ i ].m_pveciBNs = &veciBNs;
            vecsData[ i ].m_iNode = iterDataset->second + 1;
            if( fFirst )
                pthread_spin_init( &vecsData[ i ].m_sLock, 0 );
            if( pthread_create( &vecpthdThreads[ i ], NULL, evaluate2, &vecsData[ i ] ) ) {
                cerr << "Couldn't create evaluation thread: " << i << endl;
                return 1;
            }
        }
        for( i = 0; i < vecpthdThreads.size( ); ++i )
            pthread_join( vecpthdThreads[ i ], NULL );
    }

    {
        CBayesNetMinimal*	pBN;
        float*				adYes;
        float*				adNo;
        float				dPrior;

        if( sArgs.inputs_num ) {
            for( j = 0; j < vecpBNs.size( ); ++j )
                if( vecpBNs[ j ]->GetID( ) == CMeta::Deextension( CMeta::Basename( sArgs.inputs[ i ] ) ) )
                    break;
            if( j >= vecpBNs.size( ) ) {
                cerr << "Could not locate Bayes net for: " << sArgs.inputs[ i ] << endl;
                return 1;
            }
            pBN = vecpBNs[ j ];
        }
        else
            pBN = &BNDefault;

        dPrior = pBN->GetCPT( 0 ).Get( 1, 0 );
        adYes = new float[ vecpYes[ 0 ]->GetGenes( ) ];
        adNo = new float[ vecpYes[ 0 ]->GetGenes( ) ];
        for( i = 0; i < vecpYes[ 0 ]->GetGenes( ); ++i ) {
            if( !( i % 1000 ) )
                cerr << "F: global, " << i << endl;
            memcpy( adYes, vecpYes[ 0 ]->Get( i ), ( vecpYes[ 0 ]->GetGenes( ) - i - 1 ) * sizeof(*adYes) );
            memcpy( adNo, vecpNo[ 0 ]->Get( i ), ( vecpNo[ 0 ]->GetGenes( ) - i - 1 ) * sizeof(*adNo) );
            for( j = 0; j < ( vecpYes[ 0 ]->GetGenes( ) - i - 1 ); ++j )
                adYes[ j ] = CMeta::IsNaN( adYes[ j ] ) ? dPrior :
                             (float)( 1 / ( 1 + exp( (double)adNo[ j ] - (double)adYes[ j ] ) ) );
            vecpYes[ 0 ]->Set( i, adYes );
        }
    }

    for( i = 0; i < vecsData.size( ); ++i )
        pthread_spin_destroy( &vecsData[ i ].m_sLock );
    for( i = 0; i < vecstrTmps.size( ); ++i )
        _unlink( vecstrTmps[ i ].c_str( ) );
    for( i = 0; i < vecpBNs.size( ); ++i )
        delete vecpBNs[ i ];
    for( i = 0; i < vecpGenes.size( ); ++i )
        delete vecpGenes[ i ];
    for( i = 0; i < vecpYes.size( ); ++i ) {
        delete vecpYes[ i ];
        delete vecpNo[ i ];
    }

    return 0;
}

void* evaluate2( void* pData ) {
    SEvaluate2*		psData;
    size_t			i, j, iBNOne, iBNTwo;
    vector<size_t>	veciCur;
    float			dYes, dNo;

    psData = (SEvaluate2*)pData;
    for( i = psData->m_iBegin; i < psData->m_iEnd; ++i ) {
        size_t	iOne, iTwo, iBin, k;

        if( !( i % 1000 ) )
            cerr << "C: " << psData->m_strName << ", " << i << endl;
        if( ( ( iOne = (*psData->m_pveciGenes)[ i ] ) == -1 ) && ( psData->m_iZero == -1 ) )
            continue;
        veciCur.clear( );
        if( ( iBNOne = (*psData->m_pveciBNs)[ i ] ) != -1 )
            veciCur.push_back( iBNOne );
        for( j = ( i + 1 ); j < psData->m_pYes->GetGenes( ); ++j ) {
            if( ( ( iTwo = (*psData->m_pveciGenes)[ j ] ) == -1 ) && ( psData->m_iZero == -1 ) )
                continue;
            if( ( iOne == -1 ) || ( iTwo == -1 ) || ( ( iBin = psData->m_pDat->Quantize( psData->m_pDat->Get( iOne, iTwo ) ) ) == -1 ) )
                iBin = psData->m_iZero;
            if( iBin == -1 )
                continue;
            if( ( iBNTwo = (*psData->m_pveciBNs)[ j ] ) != -1 )
                veciCur.push_back( iBNTwo );
            for( k = 0; k < max( (size_t)1, veciCur.size( ) ); ++k ) {
                const CBayesNetMinimal*	pBN	= ( k >= veciCur.size( ) ) ? psData->m_pBN: (*psData->m_pvecpBNs)[ veciCur[ k ] ];

                dNo = log( pBN->GetCPT( psData->m_iNode ).Get( iBin, 0 ) );
                dYes = log( pBN->GetCPT( psData->m_iNode ).Get( iBin, 1 ) );
                pthread_spin_lock( &psData->m_sLock );
                psData->m_pNo->Get( i, j ) += dNo;
                psData->m_pYes->Get( i, j ) += dYes;
                pthread_spin_unlock( &psData->m_sLock );
            }
            if( iBNTwo != -1 )
                veciCur.pop_back( );
        }
    }

    return NULL;
}
