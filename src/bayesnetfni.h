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
#ifndef BAYESNETFNI_H
#define BAYESNETFNI_H

namespace Sleipnir {

    class CBayesNetMinimal;

#ifndef NO_SMILE

    class CBayesNetFNNode {
    protected:
        friend class CBayesNetFN;

        friend class CBayesNetFNImpl;

        static const char c_szType[];

        static CBayesNetFNNode *Open(DSL_node *);

        const std::string &GetName() const;

        unsigned char GetParameters() const;

        void Reverse();

        bool Save(DSL_node *) const;

        bool Learn(const std::vector <size_t> &);

        virtual const char *GetType() const = 0;

        virtual void Randomize() = 0;

        virtual CBayesNetFNNode *New(DSL_node *) const = 0;

        virtual bool Learn(const IDataset *, size_t, size_t) = 0;

        virtual bool Evaluate(float, std::vector<float> &) const = 0;

        virtual bool IsContinuous() const {

            return true;
        }

        std::string m_strName;
        const char *m_szType;
        CFullMatrix<float> m_Params;

        virtual ~CBayesNetFNNode() = default;
    };

    class CBayesNetFNNodeDiscrete : protected CBayesNetFNNode {
    protected:
        friend class CBayesNetFNNode;

        void Randomize();

        bool Learn(const IDataset *, size_t, size_t);

        bool Evaluate(float, std::vector<float> &) const;

        CBayesNetFNNode *New(DSL_node *pNode) const {

            return new CBayesNetFNNodeDiscrete();
        }

        const char *GetType() const {

            return "discrete";
        }

        bool IsContinuous() const {

            return false;
        }
    };

    class CBayesNetFNNodeGaussian : protected CBayesNetFNNode {
    protected:
        friend class CBayesNetFNNode;

        static const size_t c_iMu = 0;
        static const size_t c_iSigma = 1;

        void Randomize();

        bool Learn(const IDataset *, size_t, size_t);

        bool Evaluate(float, std::vector<float> &) const;

        CBayesNetFNNode *New(DSL_node *pNode) const {

            return new CBayesNetFNNodeGaussian();
        }

        const char *GetType() const {

            return "gaussian";
        }
    };

    class CBayesNetFNNodeBeta : protected CBayesNetFNNode {
    protected:
        friend class CBayesNetFNNode;

        static const size_t c_iMin = 0;
        static const size_t c_iMax = 1;
        static const size_t c_iAlpha = 2;
        static const size_t c_iBeta = 3;

        void Randomize();

        bool Learn(const IDataset *, size_t, size_t);

        bool Evaluate(float, std::vector<float> &) const;

        CBayesNetFNNode *New(DSL_node *pNode) const {

            return new CBayesNetFNNodeBeta();
        }

        const char *GetType() const {

            return "beta";
        }
    };

    class CBayesNetFNNodeExponential : protected CBayesNetFNNode {
    protected:
        friend class CBayesNetFNNode;

        static const size_t c_iMin = 0;
        static const size_t c_iBeta = 1;

        void Randomize();

        bool Learn(const IDataset *, size_t, size_t);

        bool Evaluate(float, std::vector<float> &) const;

        CBayesNetFNNode *New(DSL_node *pNode) const {

            return new CBayesNetFNNodeExponential();
        }

        const char *GetType() const {

            return "exponential";
        }
    };

    class CBayesNetFNNodeMOG : protected CBayesNetFNNode {
    protected:
        friend class CBayesNetFNNode;

        static const size_t c_iMu = 0;
        static const size_t c_iSigma = 1;

        void Randomize();

        bool Learn(const IDataset *, size_t, size_t);

        bool Evaluate(float, std::vector<float> &) const;

        CBayesNetFNNode *New(DSL_node *pNode) const {

            return new CBayesNetFNNodeMOG();
        }

        const char *GetType() const {

            return "mog";
        }
    };

    class CBayesNetFNImpl : protected CBayesNetImpl {
    protected:
        CBayesNetFNImpl();

        ~CBayesNetFNImpl();

        void Reset();

        bool Evaluate(const IDataset *, CDat *, std::vector <std::vector<float>> *, bool) const;

        bool Evaluate(const IDataset *, size_t, size_t, bool, std::vector<float> &) const;

        size_t m_iNodes;
        CBayesNetFNNode **m_apNodes;
        bool m_fSmileNet;
        DSL_network m_SmileNet;
    };

#endif // NO_SMILE

    class CBayesNetMinimalNode {
    public:
        CBayesNetMinimalNode() : m_bDefault(0xFF) {}

        unsigned char m_bDefault;
        CDataMatrix m_MatCPT;
    };

    class CBayesNetMinimalImpl : protected CBayesNetImpl, protected CFile {
    protected:
        static bool Counts2Probs(const std::vector <std::string> &, std::vector<float> &, float dAlpha = 1,
                                 float = HUGE_VAL, const CBayesNetMinimal * = NULL, size_t = 0, size_t = 0);

        CBayesNetMinimalImpl() : CBayesNetImpl(true), m_adNY(NULL) {}

        std::string m_strID;
        long double *m_adNY;
        CDataMatrix m_MatRoot;
        std::vector <CBayesNetMinimalNode> m_vecNodes;
    };

}

#endif // BAYESNETFNI_H
