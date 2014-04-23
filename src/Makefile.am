lib_LIBRARIES			= libSleipnir.a

AM_CPPFLAGS = $(BOOST_CFLAGS) $(GSL_CFLAGS) $(LOG4CPP_CFLAGS) $(SMILE_CFLAGS) $(SVM_PERF_CFLAGS) $(SVM_MULTICLASS_CFLAGS) $(SVM_HIERARCHY_CFLAGS) $(PTHREAD_CFLAGS) $(VOWPAL_WABBIT_CFLAGS) $(LIBSVM_CFLAGS)

#LDADD = $(LOG4CPP_LIBS) $(SMILE_LIBS) $(SVM_PERF_LIBS) $(SVM_MULTICLASS_LIBS) $(SVM_HIERARCHY_LIBS) $(PTHREAD_LIBS)

if WITH_SVM_TOOLS
libSleipnir_SVM_SOURCES = svm.cpp svmperf.cpp
libSleipnir_SVM_INCLUDES = svm.h svmi.h svmperf.h
endif

if WITH_SVM_MULTICLASS_TOOLS
libSleipnir_SVM_MULTICLASS_SOURCES = svmstruct.cpp
libSleipnir_SVM_MULTICLASS_INCLUDES = svmstruct.h
endif

if WITH_SVM_HIERARCHY_TOOLS
libSleipnir_SVM_HIERARCHY_SOURCES = svmstructtree.cpp
libSleipnir_SVM_HIERARCHY_INCLUDES = svmstructtree.h
endif

if WITH_LIBSVM_TOOLS
libSleipnir_LIBSVM_SOURCES = libsvm.cpp
libSleipnir_LIBSVM_INCLUDES = libsvm.h
endif

if WITH_GSL_TOOLS
libSleipnir_GSL_SOURCES = \
	seekmap.cpp				\
	seekdataset.cpp			\
	seekweight.cpp			\
	seekevaluate.cpp		\
	seekquery.cpp			\
	seekreader.cpp			\
	seekwriter.cpp			\
	seekplatform.cpp		\
	seekcentral.cpp		\
	seeknetwork.cpp
libSleipnir_GSL_INCLUDES = \
	seekplatform.h			\
	seekcentral.h			\
	seekreader.h			\
	seekwriter.h			\
	seekquery.h				\
	seekevaluate.h			\
	seekweight.h			\
	seekreader.h			\
	seekwriter.h			\
	seekdataset.h			\
	seekmap.h				\
	seeknetwork.h			
endif
libSleipnir_a_SOURCES	=	\
	annotation.cpp			\
	annotationobo.cpp		\
	annotationkegg.cpp		\
	annotationmips.cpp		\
	bayesnet.cpp			\
	bayesnetfn.cpp			\
	bayesnetmixed.cpp		\
	bayesnetpnl.cpp			\
	bayesnetsmile.cpp		\
	bayesnetsmileelr.cpp	\
	clusthierarchical.cpp	\
	clustkmeans.cpp			\
	clustpivot.cpp			\
	clustqtc.cpp			\
	coalesce.cpp			\
	coalescecluster.cpp		\
	coalescemotifs.cpp		\
	coalescestructs.cpp		\
	color.cpp				\
	compactmatrix.cpp		\
	database.cpp			\
	datapair.cpp			\
	datasetcompact.cpp		\
	dataset.cpp				\
	dat.cpp					\
	example.cpp				\
	fasta.cpp				\
	file.cpp				\
	genome.cpp				\
	hmm.cpp					\
	mathb.cpp				\
	measure.cpp				\
	meta.cpp				\
	orthology.cpp			\
	pcl.cpp					\
	pclset.cpp				\
	pst.cpp					\
	serverclient.cpp		\
	server.cpp				\
	statistics.cpp			\
	stdafx.cpp				\
	$(libSleipnir_SVM_SOURCES) \
	$(libSleipnir_GSL_SOURCES) \
	$(libSleipnir_LIBSVM_SOURCES) \
	$(libSleipnir_SVM_MULTICLASS_SOURCES)\
	$(libSleipnir_SVM_HIERARCHY_SOURCES)\
	vwb.cpp					\
	strassen.cpp

include_HEADERS			=	\
	annotation.h			\
	annotationi.h			\
	bayesnetfni.h			\
	bayesnet.h				\
	bayesneti.h				\
	bayesnetint.h			\
	clusthierarchical.h		\
	clusthierarchicali.h	\
	clustkmeans.h			\
	clustpivot.h			\
	clustqtc.h				\
	clustqtci.h				\
	coalesce.h				\
	coalescei.h				\
	coalescebasei.h			\
	coalescecluster.h		\
	coalesceclusteri.h		\
	coalescemotifs.h		\
	coalescemotifsi.h		\
	coalescestructsi.h		\
	color.h					\
	colori.h				\
	compactmatrix.h			\
	compactmatrixi.h		\
	database.h				\
	databasei.h				\
	datapair.h				\
	datapairi.h				\
	dataset.h				\
	dataseti.h				\
	dat.h					\
	dati.h					\
	examplei.h				\
	fasta.h					\
	fastai.h				\
	file.h					\
	filei.h					\
	fullmatrix.h			\
	genome.h				\
	genomei.h				\
	halfmatrix.h			\
	halfmatrixi.h			\
	hmm.h					\
	hmmi.h					\
	mathb.h					\
	mathbi.h				\
	measure.h				\
	measurei.h				\
	meta.h					\
	metai.h					\
	orthology.h				\
	orthologyi.h			\
	pcl.h					\
	pcli.h					\
	pclset.h				\
	pclseti.h				\
	pst.h					\
	psti.h					\
	serverclient.h			\
	serverclienti.h			\
	server.h				\
	serveri.h				\
	sparsematrix.h			\
	statistics.h			\
	statisticsi.h			\
	stdafx.h				\
	strassen.h				\
	$(libSleipnir_SVM_INCLUDES)  \
	$(libSleipnir_GSL_SOURCES) \
	$(libSleipnir_LIBSVM_INCLUDES) \
	$(libSleipnir_SVM_MULTICLASS_INCLUDES)\
	$(libSleipnir_SVM_HIERARCHY_INCLUDES)\
	trie.h					\
	triei.h					\
	typesi.h				\
	vwb.h					\
	vwbi.h
