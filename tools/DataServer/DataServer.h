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
#ifndef DATASERVER_H
#define DATASERVER_H


struct SDataServerData {
    const CDatabase &m_Database;
    const vector <string> m_vecstrDatasets;
    const CPCL &m_VarPCL;
    const vector <size_t> m_veciPCLGeneIdx;
    const vector <size_t> m_veciPCLDataIdx;
    const string m_strQuant;
    const int m_iThreads;

    SDataServerData(const CDatabase &Database,
                    const vector <string> &vecstrDatasets,
                    const CPCL &VarPCL,
                    vector <size_t> &veciPCLGeneIdx,
                    vector <size_t> &veciPCLDataIdx,
                    const string strQuant,
                    int iThreads) :
            m_Database(Database),
            m_vecstrDatasets(vecstrDatasets),
            m_VarPCL(VarPCL),
            m_veciPCLGeneIdx(veciPCLGeneIdx),
            m_veciPCLDataIdx(veciPCLDataIdx),
            m_strQuant(strQuant),
            m_iThreads(iThreads) {}
};

class CDataServer : public IServerClient {
public:
    CDataServer(SOCKET, const string &, const SDataServerData &);

    ~CDataServer();

    IServerClient *NewInstance(SOCKET, uint32_t, uint16_t);

    void Destroy();

    bool ProcessMessage(const std::vector<unsigned char> &);

    const CDatabase &GetDatabase() const {
        return m_sData.m_Database;
    }

    const size_t GetDatasets() const {
        return GetDatabase().GetDatasets();
    }

    const size_t GetGenes() const {
        return GetDatabase().GetGenes();
    }

    const vector <string> &GetDatasetNames() const {
        return m_sData.m_vecstrDatasets;
    }

    void GetScores(const vector <size_t> &, const vector<float> &,
                   vector<float> &, vector<float> &, vector <size_t> &, CFullMatrix<float> *, CFullMatrix<float> *);

    void GetScores(const vector <size_t> &, const vector <size_t> &,
                   CFullMatrix<float> &, CFullMatrix<float> &);

    void GetDataWeights(size_t, CFullMatrix<float> &, float, float, vector<float> &);

private:

    size_t ProcessDatasetSearch(const vector<unsigned char> &, size_t);

    size_t ProcessDatasetMeasure(const vector<unsigned char> &, size_t);

    size_t ProcessDatasetRetrieve(const vector<unsigned char> &, size_t);

    string GetQuantFile() const { return m_sData.m_strQuant; }

    const CPCL &GetVarPCL() const {
        return m_sData.m_VarPCL;
    }

    const vector <size_t> &GetPCLGeneIdx() const {
        return m_sData.m_veciPCLGeneIdx;
    }

    const vector <size_t> &GetPCLDataIdx() const {
        return m_sData.m_veciPCLDataIdx;
    }

    const float GetGeneVar(size_t iGene, size_t iData) {
        return m_sData.m_VarPCL.Get(m_sData.m_veciPCLGeneIdx[iGene],
                                    m_sData.m_veciPCLDataIdx[iData]);
    }

    const int GetMaxThreads() const {
        return m_sData.m_iThreads;
    }

    typedef size_t (CDataServer::*TPFNProcessor)(const std::vector<unsigned char> &, size_t);

    static const TPFNProcessor c_apfnProcessors[];

    SOCKET m_iSocket;
    const SDataServerData &m_sData;
    string m_strConnection;
};

#endif // DATASERVER_H
