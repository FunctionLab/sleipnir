if WITH_TOOLS
  TOOLS_SUBDIRS = \
	  Answerer \
	  Cliquer \
	  Clinician \
	  Clusterer \
	  Clusters2Dab \
	  COALESCE \
	  Combiner \
	  NetworkCombiner \
	  SparseNetCombiner \
	  DChecker \
	  Dat2Dab \
	  Dat2Graph \
	  Dat2PCL \
	  Data2Bnt \
	  Data2Features \
	  Data2Sql \
	  DataDumper \
	  DataServer \
	  Distancer \
	  Explainer \
	  Filterer \
	  Funcaeologist \
	  Funcifier \
	  Funcographer \
	  Hubber \
	  KNNImputer \
	  Mat2Txt \
	  MCluster \
	  Matcher \
	  MIer \
	  MIed \
	  NetworkRanker \
	  Normalizer \
	  Orthologer \
	  Overlapper \
	  PCLPlotter \
	  PCL2Bin \
	  Randomizer \
	  Seqs2Ngrams \
	  SMRF \
	  SVDer \
	  Synthesizer \
	  NetMiner \
	  Txt2Bin
if WITH_SMILE_TOOLS
    SMILE_TOOLS_SUBDIRS = \
	  BNConverter \
	  BNCreator \
	  SpeciesConnector \
	  BNEvaluator \
	  BNFunc \
	  BNTester \
	  BNTruster \
	  BNs2Txt \
	  BNUnraveler \
	  BNWeaver \
	  Contexter \
	  Counter \
	  Data2DB \
      DBCombiner \
	  DSLConverter \
	  Dab2Dad \
	  Edges2Posteriors \
	  MEFIT
if WITH_GSL_TOOLS
    GSL_TOOLS_SUBDIRS = \
	  SeekTest \
      SeekPrep \
      SeekReader \
      SeekMiner \
	  SeekServer \
	  SeekEvaluator \
	  SeekGeneRecommender \
	  SeekAggregatedDataset \
      SeekPValue \
      SeekIterative \
	  PCLServer
endif

if WITH_BOOST_TOOLS
      BOOST_TOOLS_SUBDIRS = BNServer
endif
endif
if WITH_SVM_TOOLS
    SVM_TOOLS_SUBDIRS = \
	  Data2Svm \
	  SVMer \
	  SVMperfing \
	  SVMperfer \
	SVMfe
endif
if WITH_SVM_MULTICLASS_TOOLS
    SVM_MULTICLASS_SUBDIRS = \
		SVMmulticlass
endif

if WITH_SVM_HIERARCHY_TOOLS
    SVM_HIERARCHY_SUBDIRS = \
		SVMhierarchy
endif

if WITH_LIBSVM_TOOLS
    LIBSVM_TOOLS_SUBDIRS = \
	  LibSVMer
endif
if WITH_VW_TOOLS
    VW_TOOLS_SUBDIRS = VWer
endif
if WITH_READLINE_TOOLS
    READLINE_TOOLS_SUBDIRS = OntoShell
endif
endif

SUBDIRS = $(TOOLS_SUBDIRS) $(SMILE_TOOLS_SUBDIRS) $(BOOST_TOOLS_SUBDIRS) $(SVM_TOOLS_SUBDIRS) $(READLINE_TOOLS_SUBDIRS) $(VW_TOOLS_SUBDIRS) $(GSL_TOOLS_SUBDIRS) $(LIBSVM_TOOLS_SUBDIRS) $(SVM_MULTICLASS_TOOLS_SUBDIRS) $(SVM_HIERARCHY_TOOLS_SUBDIRS)
