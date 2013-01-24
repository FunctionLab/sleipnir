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
	  DChecker \
	  Dat2Dab \
	  Dat2Graph \
	  Data2Bnt \
	  Data2Features \
	  Data2Sql \
	  DataDumper \
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
	  Normalizer \
	  Orthologer \
	  Overlapper \
	  PCLPlotter \
	  PCLServer \
	  PCL2Bin \
	  Randomizer \
	  Seqs2Ngrams \
	  SMRF \
	  SVDer \
	  Synthesizer \
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
	  SeekTest \
      SeekPrep \
      SeekReader \
      SeekMiner \
	  SeekServer \
	  SeekEvaluator \
	  SeekGeneRecommender \
	  SeekAggregatedDataset \
	  DSLConverter \
	  Dab2Dad \
	  Edges2Posteriors \
	  MEFIT
if WITH_BOOST_TOOLS
      BOOST_TOOLS_SUBDIRS = BNServer
endif
endif
if WITH_SVM_TOOLS
    SVM_TOOLS_SUBDIRS = \
          Data2Svm \
	  SVMer \
	  SVMperfer
endif
if WITH_VW_TOOLS
    VW_TOOLS_SUBDIRS = VWer
endif
if WITH_READLINE_TOOLS
    READLINE_TOOLS_SUBDIRS = OntoShell
endif
endif

SUBDIRS = $(TOOLS_SUBDIRS) $(SMILE_TOOLS_SUBDIRS) $(BOOST_TOOLS_SUBDIRS) $(SVM_TOOLS_SUBDIRS) $(READLINE_TOOLS_SUBDIRS) $(VW_TOOLS_SUBDIRS)
