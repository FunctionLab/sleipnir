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

#define BACKLOG 10   // how many pending connections queue will hold
char *PORT;
int NUM_DSET_MEMORY = 50;
CPCL **pcl;
CSeekCentral *csfinal;
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

#define NUM_THREADS 8
char THREAD_OCCUPIED[NUM_THREADS];

struct thread_data{
	//CSeekCentral *csf;
	string strQuery;
	string strOutputDir;
	string strSearchDatasets;

	float rbp_p;
	float query_fraction_required;
	float genome_fraction_required;
	string distanceMeasure; //Correlation, Zscore, or ZscoreHubbinessCorrected
	string searchMethod; //RBP, EqualWeighting, or OrderStatistics

	string correlationSign; //positive, negative

    int threadid;
    int new_fd;
};

void *do_query(void *th_arg){
	struct thread_data *my = (struct thread_data *) th_arg;
	int new_fd = my->new_fd;
	int threadid = my->threadid;
	//CSeekCentral *csf = my->csf;
	string strQuery = my->strQuery;
	string strOutputDir = my->strOutputDir;
	string strSearchDatasets = my->strSearchDatasets;
	string distanceMeasure = my->distanceMeasure;
	string searchMethod = my->searchMethod;
	//if "negative", then rank datasets and genes by most negative correlations
	string correlationSign = my->correlationSign;

	float rbp_p = my->rbp_p;
	float query_fraction_required = my->query_fraction_required;
	float genome_fraction_required = my->genome_fraction_required;

	CSeekCentral *csu = new CSeekCentral();
	enum CSeekDataset::DistanceMeasure eDM = CSeekDataset::Z_SCORE;
	bool bSubtractGeneAvg = false;
	bool bNormPlatform = false;

	bool bNegativeCor = false;
	if(correlationSign=="positive"){
		bNegativeCor = false;
	}else{
		bNegativeCor = true;
	}

	if(distanceMeasure=="Correlation"){
		eDM = CSeekDataset::CORRELATION;
	}else if(distanceMeasure=="Zscore"){
		//do nothing
	}else if(distanceMeasure=="ZscoreHubbinessCorrected"){
		bSubtractGeneAvg = true;
		bNormPlatform = true;
		//bNormPlatform = false;
	}

	//fprintf(stderr, "%s\n%s\n%s\n%.2f\n", strOutputDir.c_str(), strQuery.c_str(), strSearchDatasets.c_str(), query_fraction_required);

	bool r = csu->Initialize(strOutputDir, strQuery, strSearchDatasets, csfinal,
		new_fd, query_fraction_required, genome_fraction_required, eDM, bSubtractGeneAvg,
		bNormPlatform, bNegativeCor);

	//if r is false, then one of the query has no datasets 
	//containing any of the query (because of CheckDatasets() in Initialize()),
	//exit in this case
	if(r){
		if(searchMethod=="EqualWeighting"){
			csu->EqualWeightSearch();
		}else if(searchMethod=="OrderStatistics"){
			csu->OrderStatistics();
		}else{
			const gsl_rng_type *T;
			gsl_rng *rnd;
			gsl_rng_env_setup();
			T = gsl_rng_default;
			rnd = gsl_rng_alloc(T);
			gsl_rng_set(rnd, 100);
			utype FOLD = 5;
			//enum PartitionMode PART_M = CUSTOM_PARTITION;
			enum CSeekQuery::PartitionMode PART_M = CSeekQuery::LEAVE_ONE_IN;
			csu->CVSearch(rnd, PART_M, FOLD, rbp_p);
			gsl_rng_free(rnd);
		}
	}

	csu->Destruct();
	delete csu;

	//fprintf(stderr, "Done search\n"); system("date +%s%N 1>&2");

	/*char *sm = (char*)malloc(12);
	sprintf(sm, "Done search");
	sm[11] = '\0';
	send_msg(new_fd, sm, 12);
	free(sm);
	*/
	//Sending back to client, still need to be completed===============
	//pthread_mutex_lock(&mutexGet);
	//THREAD_OCCUPIED[threadid] = 0;
	//close(new_fd);
	//pthread_mutex_lock(&mutexGet);

	pthread_mutex_lock(&mutexGet);
	close(new_fd);
	THREAD_OCCUPIED[threadid] = 0;
	pthread_mutex_unlock(&mutexGet);

	int ret = 0;
	pthread_exit((void*)ret);
}

int main( int iArgs, char** aszArgs ) {
	static const size_t	c_iBuffer	= 1024;
#ifdef WIN32
	pthread_win32_process_attach_np( );
#endif // WIN32
	gengetopt_args_info	sArgs;
	int lineSize = 1024;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1; }

	bool useNibble = false;

	if(sArgs.is_nibble_flag==1){
		fprintf(stderr, "Nibble integration is not supported! Please use a non-nibble CDatabase.\n");
		useNibble = true;
		return 1;
	}

	// Random Number Generator Initializations
	utype i, j;
	bool bOutputWeightComponent = true;
	bool bSimulateWeight = true;
	bool bSubtractAvg = false;
	bool bNormPlatform = false;
	bool bLogit = false;
	bool bVariance = false;

	string strGvar = sArgs.dir_gvar_arg;
	if(strGvar!="NA"){
		bVariance = true;
	}

	csfinal = new CSeekCentral();
	CSeekDBSetting *dbSetting = new CSeekDBSetting(sArgs.dir_gvar_arg,
		sArgs.dir_sinfo_arg, sArgs.dir_platform_arg, sArgs.dir_prep_in_arg,
		sArgs.dir_in_arg, sArgs.input_arg, sArgs.quant_arg, sArgs.dset_arg,
		sArgs.num_db_arg);
	vector<CSeekDBSetting*> cc;
	cc.push_back(dbSetting);

	string add_db = sArgs.additional_db_arg;
	if(add_db!="NA"){
		ifstream ifsm;
		ifsm.open(add_db.c_str());
		if(!ifsm.is_open()){
			fprintf(stderr, "Error opening file %s\n", add_db.c_str());
			return false;
		}
		char acBuffer[lineSize];
		utype c_iBuffer = lineSize;
		vector<map<string,string> > parameters; //an array of CDatabase's
		while(!ifsm.eof()){
			ifsm.getline(acBuffer, c_iBuffer-1);
			if(acBuffer[0]==0) break;
			acBuffer[c_iBuffer-1]=0;
			string strB = acBuffer;
			if(strB=="START"){
				map<string,string> p;
				while(!ifsm.eof()){
					ifsm.getline(acBuffer, c_iBuffer-1);
					if(acBuffer[0]==0){
						fprintf(stderr, "Invalid line (empty)\n");
						return 1;
					}
					strB = acBuffer;
					if(strB=="END") break;
					vector<string> tok;
					CMeta::Tokenize(acBuffer, tok); //separator is tab
					p[tok[0]] = tok[1];
				}
				parameters.push_back(p);
			}
		}
		ifsm.close();
		if(parameters.size()==0){
			fprintf(stderr, "Error, extra_db setting file must begin with START and end with END lines\n");
			return 1;
		}

		for(i=0; i<parameters.size(); i++){
		string sinfo_dir = "NA";
		string gvar_dir = "NA";
		string platform_dir = "NA";
		string prep_dir = "NA";
		string db_dir = "NA";
		string dset_map_file = "NA";
		string gene_map_file = "NA";
		string quant_file = "NA";
		int num_db = -1;

		if(parameters[i].find("SINFO_DIR")->second=="NA"){
			fprintf(stderr, "Please specify an sinfo directory for the extra db\n");
			return false;
		}
		sinfo_dir = parameters[i].find("SINFO_DIR")->second;
		if(parameters[i].find("GVAR_DIR")!=parameters[i].end())
			gvar_dir = parameters[i].find("GVAR_DIR")->second;
		if(parameters[i].find("PREP_DIR")==parameters[i].end() ||
			parameters[i].find("PLATFORM_DIR")==parameters[i].end() ||
			parameters[i].find("DB_DIR")==parameters[i].end() ||
			parameters[i].find("DSET_MAP_FILE")==parameters[i].end() ||
			parameters[i].find("GENE_MAP_FILE")==parameters[i].end() ||
			parameters[i].find("QUANT_FILE")==parameters[i].end() ||
			parameters[i].find("NUMBER_OF_DB")==parameters[i].end()){
			fprintf(stderr, "Some arguments are missing. Please make sure the following are provided:\n");
			fprintf(stderr, "PREP_DIR, DB_DIR, DSET_MAP_FILE, GENE_MAP_FILE, QUANT_FILE, NUMBER_OF_DB\n");
			return false;
		}

		platform_dir = parameters[i].find("PLATFORM_DIR")->second;
		db_dir = parameters[i].find("DB_DIR")->second;
		prep_dir = parameters[i].find("PREP_DIR")->second;
		dset_map_file = parameters[i].find("DSET_MAP_FILE")->second;
		gene_map_file = parameters[i].find("GENE_MAP_FILE")->second;
		quant_file = parameters[i].find("QUANT_FILE")->second;
		num_db = atoi(parameters[i].find("NUMBER_OF_DB")->second.c_str());

		CSeekDBSetting *dbSetting2 = new CSeekDBSetting(gvar_dir, sinfo_dir,
			platform_dir, prep_dir, db_dir, gene_map_file, quant_file, dset_map_file,
			num_db);
		cc.push_back(dbSetting2);
		}
	}

	if(!csfinal->Initialize(cc,
		//"/tmp/ex_query2.txt", 
		sArgs.buffer_arg, !!sArgs.output_text_flag,
		bOutputWeightComponent, bSimulateWeight,
		CSeekDataset::CORRELATION, //to be overwritten by individual search instance's setting
		bVariance, //decide whether or not to load gvar
		bSubtractAvg, bNormPlatform, //to be overwritten by individual search instance's settings
		bLogit, //always false
		sArgs.score_cutoff_arg, 
		0.0, //min query fraction (to be overwrriten)
		0.0, //min genome fraction (to be overwrriten)
		!!sArgs.square_z_flag, //default
		false, 1, false, NULL, useNibble, sArgs.num_threads_arg)) //default
	{
		fprintf(stderr, "Error occurred!\n");
		return -1;
	}

	signal(SIGPIPE, SIG_IGN);
	//utype i;
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
	
		int d = 0;
		pthread_mutex_lock(&mutexGet);
		for(d=0; d<NUM_THREADS; d++){
			if(THREAD_OCCUPIED[d]==0) break;
		}

		if(d==NUM_THREADS){
			close(new_fd); 
			pthread_mutex_unlock(&mutexGet);
			continue;
		}

		THREAD_OCCUPIED[d] = 1;
		pthread_mutex_unlock(&mutexGet);

		//receiving query from client, still need to be completed
		//only needs to receive three strings,
		//queryFile, searchdatasetFile, outputDir

		string strSearchDataset;
		string strQuery;
		string strOutputDir;

		//search parameters
		string strSearchParameter;
		//format: _ delimited		
		//[0]: search_method ("RBP", "OrderStatistics", or "EqualWeighting") 
		//[1]: rbp_p
		//[2]: minimum fraction of query required
		//[3]: minimum fraction of genome required
		//[4]: distance measure ("Correlation", "Zscore", or "ZscoreHubbinessCorrected")
		//[5]: correlation sign ("positive", "negative")

		if(CSeekNetwork::Receive(new_fd, strSearchDataset)==-1){
			fprintf(stderr, "Error receiving from client!\n");
		}

		if(CSeekNetwork::Receive(new_fd, strQuery)==-1){
			fprintf(stderr, "Error receiving from client!\n");
		}

		if(CSeekNetwork::Receive(new_fd, strOutputDir)==-1){
			fprintf(stderr, "Error receiving from client!\n");
		}

		if(CSeekNetwork::Receive(new_fd, strSearchParameter)==-1){
			fprintf(stderr, "Error receiving from client!\n");
		}

		vector<string> searchParameterTokens;
		fprintf(stderr, "%s\n", strSearchParameter.c_str());
		CMeta::Tokenize(strSearchParameter.c_str(), searchParameterTokens, "_");
		thread_arg[d].searchMethod = searchParameterTokens[0];
		thread_arg[d].rbp_p = atof(searchParameterTokens[1].c_str());
		thread_arg[d].query_fraction_required = 
			atof(searchParameterTokens[2].c_str());
		thread_arg[d].genome_fraction_required = 
			atof(searchParameterTokens[3].c_str());
		thread_arg[d].distanceMeasure = searchParameterTokens[4];
		thread_arg[d].correlationSign = searchParameterTokens[5];
		//=========================================================

		thread_arg[d].threadid = d;
		thread_arg[d].new_fd = new_fd;
		thread_arg[d].strQuery = strQuery;
		thread_arg[d].strOutputDir = strOutputDir;
		thread_arg[d].strSearchDatasets = strSearchDataset;
		//thread_arg[d].csfinal = &csfinal;

		int ret;
		pthread_create(&th[d], NULL, do_query, (void*) &thread_arg[d]);

		//int ret;
		//pthread_create(&th[d], NULL, do_query, (void *) &thread_arg[d]);
		//pthread_join(th[d], (void **)&ret);

	}

#ifdef WIN32
	pthread_win32_process_detach_np( );
#endif // WIN32
	return 0; 

}
