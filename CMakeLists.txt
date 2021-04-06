cmake_minimum_required(VERSION 3.14.5)
project(sleipnir)

add_link_options(-lpthread)

set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Generating release config")
    set(CMAKE_BUILD_TYPE Release)
    add_compile_options(-O3)
elseif(${CMAKE_BUILD_TYPE} MATCHES Debug)
    message(STATUS "Generating debug config")
    add_compile_options(-O0)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include_directories(src)

option(PACKAGE_TESTS "Build the tests" ON)
if(PACKAGE_TESTS)
    enable_testing()
    include(GoogleTest)
    add_subdirectory(tests)
endif()

set(SLEIPNIR_SOURCES
    src/example_gtest.cpp
    src/example_gtest.h
    src/annotation.cpp
    src/annotation.h
    src/annotationi.h
    src/annotationkegg.cpp
    src/annotationmips.cpp
    src/annotationobo.cpp
    src/bayesnet.cpp
    src/bayesnet.h
    src/bayesnetfn.cpp
    src/bayesnetfni.h
    src/bayesneti.h
    src/bayesnetint.h
    src/bayesnetmixed.cpp
    src/bayesnetpnl.cpp
    src/bayesnetsmile.cpp
    src/bayesnetsmileelr.cpp
    src/clusthierarchical.cpp
    src/clusthierarchical.h
    src/clusthierarchicali.h
    src/clustkmeans.cpp
    src/clustkmeans.h
    src/clustkmeansi.h
    src/clustpivot.cpp
    src/clustpivot.h
    src/clustqtc.cpp
    src/clustqtc.h
    src/clustqtci.h
    src/coalesce.cpp
    src/coalesce.h
    src/coalescebasei.h
    src/coalescecluster.cpp
    src/coalescecluster.h
    src/coalesceclusteri.h
    src/coalescei.h
    src/coalescemotifs.cpp
    src/coalescemotifs.h
    src/coalescemotifsi.h
    src/coalescestructs.cpp
    src/coalescestructsi.h
    src/color.cpp
    src/color.h
    src/colori.h
    src/compactmatrix.cpp
    src/compactmatrix.h
    src/compactmatrixi.h
    src/dat.cpp
    src/dat.h
    src/database.cpp
    src/database.h
    src/databasei.h
    src/datapair.cpp
    src/datapair.h
    src/datapairi.h
    src/dataset.cpp
    src/dataset.h
    src/datasetcompact.cpp
    src/dataseti.h
    src/dati.h
    src/example.cpp
    src/examplei.h
    src/fasta.cpp
    src/fasta.h
    src/fastai.h
    src/file.cpp
    src/file.h
    src/filei.h
    src/fullmatrix.h
    src/genome.cpp
    src/genome.h
    src/genomei.h
    src/halfmatrix.h
    src/halfmatrixi.h
    src/hmm.cpp
    src/hmm.h
    src/hmmi.h
    src/mathb.cpp
    src/mathb.h
    src/mathbi.h
    src/measure.cpp
    src/measure.h
    src/measurei.h
    src/meta.cpp
    src/meta.h
    src/metai.h
    src/orthology.cpp
    src/orthology.h
    src/orthologyi.h
    src/pcl.cpp
    src/pcl.h
    src/pcli.h
    src/pclset.cpp
    src/pclset.h
    src/pclseti.h
    src/pst.cpp
    src/pst.h
    src/psti.h
    src/seekbasic.h
    src/seekcentral.cpp
    src/seekcentral.h
    src/seekdataset.cpp
    src/seekdataset.h
    src/seekerror.cpp
    src/seekerror.h
    src/seekevaluate.cpp
    src/seekevaluate.h
    src/seekhelper.cpp
    src/seekhelper.h
    src/seekmap.cpp
    src/seekmap.h
    src/seeknetwork.cpp
    src/seeknetwork.h
    src/seekplatform.cpp
    src/seekplatform.h
    src/seekquery.cpp
    src/seekquery.h
    src/seekreader.cpp
    src/seekreader.h
    src/seekweight.cpp
    src/seekweight.h
    src/seekwriter.cpp
    src/seekwriter.h
    src/server.cpp
    src/server.h
    src/serverclient.cpp
    src/serverclient.h
    src/serverclienti.h
    src/serveri.h
    src/sparsematrix.h
    src/sparsematrixi.h
    src/statistics.cpp
    src/statistics.h
    src/statisticsi.h
    src/stdafx.cpp
    src/stdafx.h
    src/strassen.cpp
    src/strassen.h
    src/trie.h
    src/triei.h
    src/typesi.h
    src/vwb.cpp
    src/vwb.h
    src/vwbi.h)


# compiler-specific logic for building multithreaded code

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    find_library(LIBOMP libomp.a)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    add_compile_options(-fopenmp)
    add_link_options(-lgslcblas)
    list(APPEND EXTRA_LIBRARIES stdc++fs)
    set(LIBOMP -lgomp)
else()
    message(SEND_INFO "Found some weird compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()


option(USE_GENGETOPT "Use gengetopt" ON)
if(USE_GENGETOPT)
    find_program(GENGETOPT gengetopt)
    set(GENGETOPT_PROGRAM GENGETOPT)
endif()


option(USE_READLINE "Use readline" ON)
if(USE_READLINE)
    find_library(READLINE_LIBRARY readline)
    find_path(READLINE_INCLUDE readline)
    list(APPEND EXTRA_LIBRARIES ${READLINE_LIBRARY})
    list(APPEND EXTRA_INCLUDES ${READLINE_INCLUDE})
endif()


option(USE_LOG4CPP "Use log4cpp" OFF)
if(USE_LOG4CPP)
    find_library(LOG4CPP_LIBRARY log4cpp)
    find_path(LOG4CPP_INCLUDE log4cpp)
    list(APPEND EXTRA_LIBRARIES ${LOG4CPP_LIBRARY})
    list(APPEND EXTRA_INCLUDES ${LOG4CPP_INCLUDE})
else()
    list(APPEND SLEIPNIR_COMPILE_DEFS USE_LOG4CPP_STUB)
endif()

option(USE_SMILE "Use SMILE" OFF)
if(USE_SMILE)
    find_library(SMILE_LIBRARY smile)
    find_path(SMILE_INCLUDE smile.h)
    list(APPEND EXTRA_LIBRARIES ${SMILE_LIBRARY})
    list(APPEND EXTRA_INCLUDES ${SMILE_INCLUDE})
else()
    list(APPEND SLEIPNIR_COMPILE_DEFS NO_SMILE=1)
endif()

find_path(BOOST_INCLUDE boost)
if(BOOST_INCLUDE)
    list(APPEND EXTRA_INCLUDES ${BOOST_INCLUDE})
else()
    message(SEND_ERROR "Missing required Boost headers")
endif()


find_library(BOOST_GRAPH_LIBRARY boost_graph)
if(BOOST_GRAPH_LIBRARY AND BOOST_INCLUDE)
    list(APPEND EXTRA_LIBRARIES ${BOOST_GRAPH_LIBRARY})
else()
    message(SEND_ERROR "Missing required library boost_graph")
endif()


find_library(BOOST_REGEX_LIBRARY boost_regex)
if(BOOST_REGEX_LIBRARY AND BOOST_INCLUDE)
    list(APPEND EXTRA_LIBRARIES ${BOOST_REGEX_LIBRARY})
else()
    message(SEND_ERROR "Missing required library boost_regex")
endif()


include(FindGSL)
if(GSL_FOUND)
    list(APPEND EXTRA_INCLUDES ${GSL_INCLUDE_DIR})
    list(APPEND EXTRA_LIBRARIES ${GSL_LIBRARY})
else()
    message(SEND_ERROR "Missing required library gsl")
endif()


option(USE_OPENBLAS "Use open_blas" OFF)
if(USE_OPENBLAS)
    find_package(OpenBLAS REQUIRED)
    list(APPEND EXTRA_LIBRARIES ${OpenBLAS_LIBRARIES})
    list(APPEND EXTRA_INCLUDES ${OpenBLAS_INCLUDE_DIRS})
endif()


option(USE_SVMLIB "Use svm_lib" ON)
if(USE_SVMLIB)
    list(APPEND SLEIPNIR_SOURCES src/svm_local.cpp
	src/svm_local.h
	src/svm_local.cpp
	src/svmi.h
	src/libsvm.cpp
	src/libsvm.h)
    find_library(SVM_LIBRARY svm)
    find_path(SVM_INCLUDE PATH_SUFFIXES libsvm NAMES svm.h)
    list(APPEND EXTRA_LIBRARIES ${SVM_LIBRARY})
    list(APPEND EXTRA_INCLUDES ${SVM_INCLUDE})
else()
    list(APPEND SLEIPNIR_COMPILE_DEFS NO_LIBSVM=1)
endif()


option(USE_SVM_PERF "Use svm_perf" OFF)
if(USE_SVM_PERF)
    list(APPEND SLEIPNIR_SOURCES
	src/svmperf.cpp
	src/svmperf.h)
    find_library(SVM_PERF_LIBRARY libsvmperf.a)
    find_path(SVM_PERF_INCLUDE
	PATH_SUFFIXES svm_perf
	NAMES svm_common.h svm_struct_api.h)
    list(APPEND EXTRA_LIBRARIES ${SVM_PERF_LIBRARY})
    list(APPEND EXTRA_INCLUDES ${SVM_PERF_INCLUDE})
else()
    list(APPEND SLEIPNIR_COMPILE_DEFS NO_SVM_PERF=1 NO_SVM_STRUCTTREE=1 NO_SVM_STRUCT=1)
endif()

# we explicitly disable the hierarchical and multiclass svm_perfs because they don't work

list(APPEND SLEIPNIR_COMPILE_DEFS NO_SVM_MULTICLASS=1 NO_SVM_STRUCTTREE=1)


option(USE_VOWPAL_WABBIT "Use Vowpal Wabbit" ON)
if(USE_VOWPAL_WABBIT)
    find_library(VOWPAL_WABBIT_LIBRARY libvw.a)
    find_path(VOWPAL_WABBIT_INCLUDE PATH_SUFFIXES vowpalwabbit
            NAMES vw.h)
    if(VOWPAL_WABBIT_LIBRARY AND VOWPAL_WABBIT_INCLUDE)
        list(APPEND EXTRA_LIBRARIES ${VOWPAL_WABBIT_LIBRARY})
        list(APPEND EXTRA_INCLUDES ${VOWPAL_WABBIT_INCLUDE})
    else()
        list(APPEND SLEIPNIR_COMPILE_DEFS NO_VOWPAL_WABBIT)
    endif()
else()
    list(APPEND SLEIPNIR_COMPILE_DEFS NO_VOWPAL_WABBIT=1)
endif()


add_library(sleipnir STATIC ${SLEIPNIR_SOURCES})

target_include_directories(sleipnir PUBLIC ${EXTRA_INCLUDES})
target_compile_definitions(sleipnir PUBLIC ${SLEIPNIR_COMPILE_DEFS})
target_link_libraries(sleipnir)



# command line tools

# these tools depend on having libsmile available to link
if(SMILE_LIBRARY)
    list(APPEND SLEIPNIR_TOOLS
	BNConverter
	BNCreator 
	BNEvaluator
	BNServer 
	BNTester 
	BNTruster 
	BNUnraveler
	BNWeaver 
	BNs2Txt 
	Contexter
	Counter 
	DSLConverter
	Dab2Dad 
	Edges2Posteriors
	SpeciesConnector
	MEFIT
	)
endif()

# tools that depend on having svm_perf

if(SVM_PERF_LIBRARY)
    list(APPEND SLEIPNIR_TOOLS
	Data2SVM
	SVMperfer
	SVMperfing
	SVMfe
	SVMer
	)
endif()


list(APPEND SLEIPNIR_TOOLS
    Answerer
    BNFunc
    COALESCE
    Clinician
    Cliquer
    Clusterer
    Clusters2Dab
    Combiner
    DBCombiner
    DChecker
    Dat2Dab
    Dat2Graph
    Dat2PCL
    Data2Bnt
    Data2DB
    Data2Features
    Data2Sql
    DataDumper
    DataServer
    Distancer
    Explainer
    Filterer
    Funcaeologist
    Funcifier
    Funcographer
    Hubber
    KNNImputer
    LibSVMer
    MCluster
    MIed
    MIer
    Mat2Txt
    Matcher
    NetMiner
    NetworkCombiner
    NetworkRanker
    Normalizer
    OntoShell
    Orthologer
    Overlapper
    PCL2Bin
    PCLPlotter
    PCLServer
    Randomizer
    SMRF
    SVDer
    SparseNetCombiner
    SeekAggregatedDataset
    SeekEvaluator
    SeekGeneRecommender
    SeekIterative
    SeekMiner
    SeekPValue
    SeekPrep
    SeekReader
    SeekServer
    SeekTest
    Seqs2Ngrams
    Synthesizer
    Txt2Bin
    VWer
    )

foreach(TOOL ${SLEIPNIR_TOOLS})
    include_directories(tools/${TOOL})
    if (${TOOL} MATCHES OntoShell)
        add_executable(${TOOL}
            tools/${TOOL}/${TOOL}.cpp
            tools/${TOOL}/cmdline.c
            tools/${TOOL}/stdafx.cpp
            tools/${TOOL}/parser.cpp
            tools/${TOOL}/parserconsole.cpp)
    elseif (${TOOL} MATCHES BNServer)
        add_executable(${TOOL}
            tools/${TOOL}/${TOOL}.cpp
            tools/${TOOL}/cmdline.c
            tools/${TOOL}/stdafx.cpp
            tools/${TOOL}/dot.cpp)
    else ()
        add_executable(${TOOL}
            tools/${TOOL}/${TOOL}.cpp
            tools/${TOOL}/cmdline.c)
	  if (EXISTS tools/${TOOL}/stdafx.cpp)
            add_executable(${TOOL} tools/${TOOL}/stdafx.cpp)
	  endif()
    endif ()
    add_dependencies(${TOOL} sleipnir)
    target_include_directories(${TOOL} PUBLIC ${EXTRA_INCLUDES})
    target_link_libraries(${TOOL} sleipnir ${EXTRA_LIBRARIES} ${LIBOMP})
endforeach()


list(APPEND SINGLE_FILE_TOOLS
    verifyMergedDBFiles
)

foreach(TOOL ${SINGLE_FILE_TOOLS})
    add_executable(${TOOL} tools/${TOOL}.cpp)
    add_dependencies(${TOOL} sleipnir)
    target_include_directories(${TOOL} PUBLIC ${EXTRA_INCLUDES})
    target_link_libraries(${TOOL} sleipnir ${EXTRA_LIBRARIES} ${LIBOMP})
endforeach()
