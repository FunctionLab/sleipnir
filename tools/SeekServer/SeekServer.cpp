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
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>

#define BACKLOG 10   // how many pending connections queue will hold
char *PORT;
int NUM_DSET_MEMORY = 50;
CPCL **pcl;
char *available;
map<string, int> DNAME_MAP;
map<int, string> DNAME_RMAP;

pthread_mutex_t mutexGet;

void sigchld_handler(int s){
    while(waitpid(-1, NULL, WNOHANG) > 0);
}

// get sockaddr, IPv4 or IPv6:
void *get_in_addr(struct sockaddr *sa){
    if (sa->sa_family == AF_INET) {
        return &(((struct sockaddr_in*)sa)->sin_addr);
    }
    return &(((struct sockaddr_in6*)sa)->sin6_addr);
}

#define NUM_THREADS 1
char THREAD_OCCUPIED[NUM_THREADS];

int send_msg(int new_fd, char *c, int size){
	int r = send(new_fd, c, size, 0);
	if(r==-1){
		printf("client exists");
	}
	return r;
}

void cl(char *b, int size){
	int i = 0;
	for(i=0; i<size; i++){
		b[i] = '\0';
	}
}

struct thread_data{
    vector<string> queryName;
	vector<CSeekDataset*> *vc;
	vector<string> *vecstrGenes;
	vector<string> *vecstrDatasets;
	gsl_rng *rnd;
	CDatabase *DB;
    int threadid;
    int new_fd;
};

int cpy(char *d, char *s, int beg, int num){
    int i;
    for(i=0; i<num; i++){
        d[beg+i] = s[i];
    }
    return beg+num;
}

void *do_query(void *th_arg){
	struct thread_data *my = (struct thread_data *) th_arg;
	vector<string> queryName = my->queryName;
	vector<CSeekDataset*> &vc = *(my->vc);
	vector<string> &vecstrGenes = *(my->vecstrGenes);
	vector<string> &vecstrDatasets = *(my->vecstrDatasets);
	gsl_rng *rnd = my->rnd;
	CDatabase &DB = *(my->DB);
	int new_fd = my->new_fd;
	int tid = my->threadid;

	vector< vector<string> > vecstrAllQuery;
	vecstrAllQuery.push_back(queryName);

	CSeekTools::ReadDatabaselets(DB, vecstrAllQuery, vc);

	ushort i;
	ushort j;
	ushort d;
	float RATE = 0.95;
	ushort FOLD = 5;
	enum PartitionMode PART_M = CUSTOM_PARTITION;

	ushort numThreads = omp_get_max_threads();
	bool DEBUG = false;
	size_t iGenes = vecstrGenes.size();
	size_t iDatasets = vecstrDatasets.size();

	for(i=0; i<vecstrAllQuery.size(); i++){

		vector<ushort> queryGenes;
		for(j=0; j<vecstrAllQuery[i].size(); j++){
			size_t m = DB.GetGene(vecstrAllQuery[i][j]);
			if(m==-1) continue;
			queryGenes.push_back(m);
		}
		queryGenes.resize(queryGenes.size());

		vector<char> cQuery;
		CSeekTools::CreatePresenceVector(queryGenes, cQuery, iGenes);

		for(j=0; j<iDatasets; j++){
			vc[j]->InitializeQuery(queryGenes);
		}

		//fprintf(stderr, "Start creating CV partitions\n"); system("date +%s%M 1>&2");
		CSeekQuery query;
		query.InitializeQuery(queryGenes);
		query.CreateCVPartitions(rnd, PART_M, FOLD);
		//fprintf(stderr, "Done creating CV partitions\n"); system("date +%s%M 1>&2");

		ushort iQuery = query.GetQuery().size();
		//if(DEBUG) fprintf(stderr, "Query size: %d\n", iQuery);

		ushort ***rData = new ushort**[numThreads];
		for(j=0; j<numThreads; j++){
			rData[j] = CSeekTools::Init2DArray(iGenes, iQuery, (ushort)0);
		}

		vector<float> master_rank;
		CSeekTools::InitVector(master_rank, iGenes, (float) 0);

		vector<float> sum_weight;
		CSeekTools::InitVector(sum_weight, iGenes, (float) 0);

		vector<ushort> counts;
		CSeekTools::InitVector(counts, iGenes, (ushort) 0);

		vector<float> weight;
		CSeekTools::InitVector(weight, iDatasets, (float) 0);

		float **master_rank_threads = CSeekTools::Init2DArray(numThreads, iGenes, (float) 0);
		float **sum_weight_threads = CSeekTools::Init2DArray(numThreads, iGenes, (float) 0);
		ushort **counts_threads = CSeekTools::Init2DArray(numThreads, iGenes, (ushort) 0);
		vector<ushort> *rank_threads = new vector<ushort>[numThreads];
		vector<ushort> *rank_normal_threads = new vector<ushort>[numThreads];
		for(j=0; j<numThreads; j++){
			rank_threads[j].resize(iGenes);
			rank_normal_threads[j].resize(iGenes);
		}

		fprintf(stderr, "Entering search\n");
		system("date +%s%N 1>&2");

		#pragma omp parallel for \
		shared(weight, query, vc, rData, master_rank_threads, \
				sum_weight_threads, counts_threads, rank_threads) \
		private(d, j) \
		firstprivate(iDatasets, iGenes, iQuery) \
		schedule(dynamic)

		for(d=0; d<iDatasets; d++){
			ushort tid = omp_get_thread_num();
			if(DEBUG) fprintf(stderr, "Dataset %d, %s\n", d, vecstrDatasets[d].c_str());

			CSeekIntIntMap *mapQ = vc[d]->GetQueryMap();
			CSeekIntIntMap *mapG = vc[d]->GetGeneMap();

			vector<ushort> this_q;
			for(j=0; j<mapQ->GetNumSet(); j++){
				this_q.push_back(mapQ->GetReverse(j));
			}

			if(mapQ->GetNumSet()==0){
				if(DEBUG) fprintf(stderr, "This dataset is skipped\n");
				continue;
			}

			if(DEBUG) fprintf(stderr, "Initializing %d\n", this_q.size());

			vc[d]->InitializeDataMatrix(rData[tid], iGenes, iQuery);

			if(DEBUG) fprintf(stderr, "Weighting dataset\n");

			CSeekWeighter::CVWeighting(query, *vc[d], &rank_threads[tid], false);
			float w = vc[d]->GetDatasetSumWeight();

			if(w==-1){
				if(DEBUG) fprintf(stderr, "Bad weight\n");
				continue;
			}

			if(DEBUG) fprintf(stderr, "Doing linear combination\n");

			CSeekWeighter::LinearCombine(rank_normal_threads[tid], this_q, *vc[d], false);

			vc[d]->DeleteQuery();

			if(DEBUG) fprintf(stderr, "Adding contribution of dataset to master ranking: %.5f\n", w);

			ushort iGeneSet = mapG->GetNumSet();
			const vector<ushort> &allRGenes = mapG->GetAllReverse();
			vector<ushort>::const_iterator iterR = allRGenes.begin();
			vector<ushort>::const_iterator endR = allRGenes.begin() + iGeneSet;

			vector<ushort> &Rank_Normal = rank_normal_threads[tid];
			float* Master_Rank = &master_rank_threads[tid][0];
			float* Sum_Weight = &sum_weight_threads[tid][0];
			ushort* Counts = &counts_threads[tid][0];

			for(; iterR!=endR; iterR++){
				if(Rank_Normal[*iterR]==0){
					continue;
				}
				Master_Rank[*iterR] += (float) Rank_Normal[*iterR] * w;
				Sum_Weight[*iterR] += w;
				Counts[*iterR]++;
			}

			weight[d] = w;
		}
		//omp finishes

		for(j=0; j<numThreads; j++){
			ushort k;
			for(k=0; k<iGenes; k++){
				master_rank[k] += master_rank_threads[j][k];
				counts[k] += counts_threads[j][k];
				sum_weight[k]+=sum_weight_threads[j][k];
			}
		}

		CSeekTools::Free2DArray(master_rank_threads);
		CSeekTools::Free2DArray(counts_threads);
		CSeekTools::Free2DArray(sum_weight_threads);

		for(j=0; j<numThreads; j++){
			CSeekTools::Free2DArray(rData[j]);
		}
		delete[] rData;
		delete[] rank_threads;
		delete[] rank_normal_threads;

		if(DEBUG) fprintf(stderr, "Aggregating genes\n");
		for(j=0; j<iGenes; j++){
			if(counts[j]<(int)(0.5*iDatasets)){
				master_rank[j] = -320;
			}else if(sum_weight[j]==0){
				master_rank[j] = -320;
			}else{
				master_rank[j] = (master_rank[j] / sum_weight[j] - 320) / 100.0;
			}
			if(DEBUG) fprintf(stderr, "Gene %d %.5f\n", j, master_rank[j]);
		}

		if(DEBUG) fprintf(stderr, "Sorting genes\n");
		vector<AResultFloat> a;
		a.clear();
		a.resize(iGenes);
		for(j=0; j<iGenes; j++){
			a[j].i = j;
			a[j].f = master_rank[j];
		}
		if(DEBUG) fprintf(stderr, "Begin Sorting genes\n");
		sort(a.begin(), a.end());

		if(DEBUG) fprintf(stderr, "Results:\n");
		ushort jj;
		ushort ii;
		for(ii=0, jj=0; jj<500; ii++){
			if(cQuery[a[ii].i]==1) continue;
			//fprintf(stderr, "%s %.5f\n", DB.GetGene((size_t)a[ii].i).c_str(), a[ii].f);
			jj++;
		}

		if(DEBUG) fprintf(stderr, "Begin sorting dataset weight\n");
		vector<AResultFloat> af;
		af.clear();
		af.resize(iDatasets);
		for(j=0; j<iDatasets; j++){
			af[j].i = j;
			af[j].f = weight[j];
		}
		sort(af.begin(), af.end());
		

		fprintf(stderr, "Done search\n"); system("date +%s%N 1>&2");

		int *si;
		float *uf;

		char *ss = (char*)malloc(8);
		si = (int*)&ss[0]; *si = iDatasets;
		si++; *si = iGenes;
		send_msg(new_fd, ss, 8);
		free(ss);

		char *dw = (char*)malloc(4*iDatasets);
		uf = (float*)&dw[0];
		for(j=0; j<iDatasets; j++, uf++){
			*uf = af[j].f;
		}
		send_msg(new_fd, dw, 4*iDatasets);

		si = (int*)&dw[0];
		for(j=0; j<iDatasets; j++, si++){
			*si = (int) af[j].i;
		}
		send_msg(new_fd, dw, 4*iDatasets);
		free(dw);

		char *gsc = (char*)malloc(4*iGenes);
		uf = (float*)&gsc[0];
		for(j=0; j<iGenes; j++, uf++){
			*uf = a[j].f;
		}
		send_msg(new_fd, gsc, 4*iGenes);

		si = (int*)&gsc[0];
		for(j=0; j<iGenes; j++, si++){
			*si = (int) a[j].i;
		}
		send_msg(new_fd, gsc, 4*iGenes);
		free(gsc);
		
		
		/*sprintf(acBuffer, "results/%d.query", i);
		CSeekTools::WriteArrayText(acBuffer, vecstrAllQuery[i]);

		sprintf(acBuffer, "results/%d.dweight", i);
		CSeekTools::WriteArray(acBuffer, weight);

		sprintf(acBuffer, "results/%d.gscore", i);
		CSeekTools::WriteArray(acBuffer, master_rank);
		*/
	}

	for(j=0; j<iDatasets; j++){
		vc[j]->DeleteQueryBlock();
	}

	THREAD_OCCUPIED[tid] = 0;
	int ret = 0;
	close(new_fd);

	pthread_exit((void*)ret);

}



int main( int iArgs, char** aszArgs ) {
	static const size_t	c_iBuffer	= 1024;
#ifdef WIN32
	pthread_win32_process_attach_np( );
#endif // WIN32
	gengetopt_args_info	sArgs;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }

	ifstream			ifsm;
	istream*			pistm;
	vector<string>		vecstrGenes;
	char				acBuffer[ c_iBuffer ];
	ushort				i;

	/* Random Number Generator Initializations */
	const gsl_rng_type *T;
	gsl_rng *rnd;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rnd = gsl_rng_alloc(T);

	//Reading gene mapping
	if( sArgs.input_arg ) {
		string strGeneInput = sArgs.input_arg;
		vector<string> vecstrGeneID;
		if(!CSeekTools::ReadListTwoColumns(strGeneInput, vecstrGeneID, vecstrGenes)){
			return false;
		}
	}

	bool useNibble = false;
	if(sArgs.is_nibble_flag==1){
		useNibble = true;
	}

	CDatabase DB(useNibble);
	omp_set_num_threads(8);

	if(!sArgs.db_arg){
		fprintf(stderr, "Error!\n"); 
		return false;
	}

	string strDBInput = sArgs.db_arg;
	vector<string> vecstrDatasets, vecstrDP;
	if(!CSeekTools::ReadListTwoColumns(strDBInput, vecstrDatasets, vecstrDP)){
		return false;
	}
	map<string, string> mapstrstrDatasetPlatform;
	for(i=0; i<vecstrDatasets.size(); i++){
		mapstrstrDatasetPlatform[vecstrDatasets[i]] = vecstrDP[i];
	}


	fprintf(stderr, "Start reading platform\n"); system("date +%s%N 1>&2");
	string strPlatformDirectory = sArgs.dir_platform_arg;
	vector<CSeekPlatform> vp;
	map<string, ushort> mapstriPlatform;
	vector<string> vecstrPlatforms;
	CSeekTools::ReadPlatforms(strPlatformDirectory, vp, vecstrPlatforms,
		mapstriPlatform);

	fprintf(stderr, "Done reading platform\n"); system("date +%s%N 1>&2");

	string strInputDirectory = sArgs.dir_in_arg;
	string strPrepInputDirectory = sArgs.dir_prep_in_arg;
	size_t iNumDBs = sArgs.num_db_arg;
	size_t iDatasets = vecstrDatasets.size();
	size_t iGenes = vecstrGenes.size();

	fprintf(stderr, "Start reading CDatabase header\n"); system("date +%s%N 1>&2");
	DB.Open(strInputDirectory, vecstrGenes, iDatasets, iNumDBs);
	fprintf(stderr, "Done reading CDatabase header\n"); system("date +%s%N 1>&2");

	//printf("Done opening"); getchar();
	vector<CSeekDataset*> vc;
	CSeekTools::LoadDatabase(DB, strPrepInputDirectory, 
		vecstrDatasets, mapstrstrDatasetPlatform, mapstriPlatform, vp, vc);


	signal(SIGPIPE, SIG_IGN);
	//ushort i;
	for(i=0; i<NUM_THREADS; i++){
		THREAD_OCCUPIED[i] = 0;
	}

	PORT = sArgs.port_arg;

	int sockfd, new_fd;
	struct addrinfo hints, *servinfo, *p;
	struct sockaddr_storage their_addr;
	socklen_t sin_size;
	struct sigaction sa;
	char s[INET6_ADDRSTRLEN];
	char buf[10];
	int rv;
	int yes = 1;

	memset(&hints, 0, sizeof(hints));
	hints.ai_family=AF_UNSPEC;
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_flags = AI_PASSIVE;

	if((rv=getaddrinfo(NULL, PORT, &hints, &servinfo))!=0){
		fprintf(stderr, "getaddrinfo: %s\n", gai_strerror(rv));
		return 1;
	}

	// loop through all the results and bind to the first we can
	for(p = servinfo; p != NULL; p = p->ai_next) {
		if ((sockfd = socket(p->ai_family, p->ai_socktype, p->ai_protocol)) == -1) {
			perror("server: socket");
			continue;
		}
		if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(int)) == -1) {
			perror("setsockopt");
			exit(1);
		}
		if (bind(sockfd, p->ai_addr, p->ai_addrlen) == -1) {
			close(sockfd);
			perror("server: bind");
			continue;
		}
		break;
	}

	if (p == NULL)  {
		fprintf(stderr, "server: failed to bind\n");
		return 2;
	}

	freeaddrinfo(servinfo); // all done with this structure

	if (listen(sockfd, BACKLOG) == -1) {
		perror("listen");
		exit(1);
	}

	sa.sa_handler = sigchld_handler; // reap all dead processes
	sigemptyset(&sa.sa_mask);
	sa.sa_flags = SA_RESTART;
	if (sigaction(SIGCHLD, &sa, NULL) == -1) {
		perror("sigaction");
		exit(1);
	}

	printf("server: waiting for connections...\n");
	struct thread_data thread_arg[NUM_THREADS];
	pthread_t th[NUM_THREADS];
	int d = 0;

	pthread_mutex_init(&mutexGet, NULL);

	while(1){
		sin_size = sizeof their_addr;
		new_fd = accept(sockfd, (struct sockaddr *) &their_addr, &sin_size);
		if(new_fd==-1){
			perror("accept");
			continue;
		}
		inet_ntop(their_addr.ss_family, get_in_addr((struct sockaddr *)&their_addr), s, sizeof s);
		printf("server, got connection from %s\n", s);
		for(d=0; d<NUM_THREADS; d++){
			if(THREAD_OCCUPIED[d]==0) break;
		}

		if(d==NUM_THREADS){
			close(new_fd); 
			continue;
		}

		THREAD_OCCUPIED[d] = 1;

		char ar[4];
		recv(new_fd, ar, 4, 0);
		int *cn = (int*) ar; 
		int strLen = *cn;
		char *qname = (char*)malloc(strLen);
		recv(new_fd, qname, strLen, 0);

		fprintf(stderr, "%s\n", qname);

		vector<string> queryName;
		CMeta::Tokenize(qname, queryName);
	
		thread_arg[d].threadid = d;
		thread_arg[d].new_fd = new_fd;
		thread_arg[d].queryName = queryName;
		thread_arg[d].vc = &vc;
		thread_arg[d].vecstrGenes = &vecstrGenes;
		thread_arg[d].vecstrDatasets = &vecstrDatasets;
		thread_arg[d].rnd = rnd;
		thread_arg[d].DB = &DB;

		//fprintf(stderr, "Arguments: %d %d %s %s\n", d, new_fd, qname);

		free(qname);

		int ret;
		pthread_create(&th[d], NULL, do_query, (void *) &thread_arg[d]);
		/*pthread_join(th[d], (void **)&ret);
		if(ret==0){
			close(new_fd);
		}*/
	}

#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; }
