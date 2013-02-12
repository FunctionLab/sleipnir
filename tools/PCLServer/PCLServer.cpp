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
#include <list>

#define BACKLOG 10   // how many pending connections queue will hold
char *PORT;
int NUM_DSET_MEMORY = 300;
CPCL **pcl;
list<int> available;
char *loaded;
map<string, int> DNAME_MAP;
map<int, int> return_val;
map<int, string> DNAME_RMAP;

pthread_mutex_t mutexGet;

string strPrepInputDirectory;
string strSinfoInputDirectory;
map<string, int> mapstrintGene;

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
    vector<string> datasetName;
    vector<string> geneName;
	vector<string> queryName;
    int threadid;
    int new_fd;
	bool outputNormalized;
	bool outputCoexpression;
	bool outputQueryCoexpression;
	bool outputExpression;
	bool outputQueryExpression;
};

int cpy(char *d, char *s, int beg, int num){
    int i;
    for(i=0; i<num; i++){
        d[beg+i] = s[i];
    }
    return beg+num;
}

int GetOpenSlot(){
	int i = -1;
	int size = 0;
	while(size<available.size()){
		i = available.front();
		available.pop_front();
		available.push_back(i);
		if(loaded[i]==0) break;
		size++;
	}
	if(size==available.size()){
		return -1;
	}
	return i;
}

void *do_query(void *th_arg){
	struct thread_data *my = (struct thread_data *) th_arg;
	vector<string> datasetName = my->datasetName;
	vector<string> geneName = my->geneName;
	vector<string> queryName = my->queryName;
	bool outputNormalized = my->outputNormalized;
	bool outputCoexpression = my->outputCoexpression;
	bool outputQueryCoexpression = my->outputQueryCoexpression;
	bool outputExpression = my->outputExpression;
	bool outputQueryExpression = my->outputQueryExpression;
	int new_fd = my->new_fd;
	int tid = my->threadid;

	/*
	cl(buf, 5000);
	sprintf(buf, "Begin...\n");
	if(send_msg(new_fd, buf, 5000)==-1){
		THREAD_OCCUPIED[tid] = 0;
		pthread_exit(0);
	}*/
	pthread_mutex_lock(&mutexGet);

	size_t i;

	vector<string>::const_iterator iterS = datasetName.begin();
	vector<CPCL*> vc;

	for(i=0; i<NUM_DSET_MEMORY; i++){
		loaded[i] = 0;
	}

	fprintf(stderr, "start processing...\n");
	for(; iterS!=datasetName.end(); iterS++){
		int n = -1;
		map<string, int>::const_iterator iterM = DNAME_MAP.find(*iterS);
		if(iterM!=DNAME_MAP.end()){
			n = iterM->second;
			fprintf(stderr, "found %d for dataset %s...\n", n, iterS->c_str());
			vc.push_back(pcl[n]);
			loaded[n] = 1;
			continue;
		}
		//pthread_mutex_lock(&mutexGet);
		n = GetOpenSlot();
		map<int, string>::const_iterator iterRM = DNAME_RMAP.find(n);
		if(iterRM!=DNAME_RMAP.end()){
			DNAME_MAP.erase(iterRM->second);
			DNAME_RMAP.erase(n);
		}
		DNAME_MAP[*iterS] = n;
		DNAME_RMAP[n] = *iterS;
		loaded[n] = 1;
		//pthread_mutex_unlock(&mutexGet);

		fprintf(stderr, "acquired %d for dataset %s...\n", n, iterS->c_str());
		//pcl[n]->Reset();
		fprintf(stderr, "dataset reset\n");
		pcl[n]->Open((*iterS).c_str());
		fprintf(stderr, "dataset opened\n");
		vc.push_back(pcl[n]);
	}

	size_t j, k;
	//vector<CFullMatrix<float> *> vC; //gene expression
	//vector<CFullMatrix<float> *> vQ; //query expression (if enabled) (EXTRA)
	//int totNumExperiments = 0;
	int genes = geneName.size();
	int queries = queryName.size(); //(EXTRA)
	int datasets = datasetName.size();

	vector<float> vecG, vecQ, vecCoexpression, vecqCoexpression;
	vector<float> sizeD;

	//Qian added
	//vector<CSeekDataset *> vd; //(EXTRA)
	//vd.resize(vc.size()); //(EXTRA)
	//CFullMatrix<float> *vCoexpression = new CFullMatrix<float>(); //(EXTRA)
	//vCoexpression->Initialize(genes, datasets); //(EXTRA)

	fprintf(stderr, "Reading data...\n");

	float NaN = -9999;	
	//for each dataset
	for(i=0; i<vc.size(); i++){
		CPCL *pp = vc[i];
		int ps = pp->GetExperiments() - 2;
		int gs = pp->GetExperiments();
		CFullMatrix<float> *fq = NULL;
		CFullMatrix<float> *ff = NULL;
		vector<float> vCoexpression;
		vCoexpression.resize(genes);
		vector<float> vqCoexpression;
		vqCoexpression.resize(queryName.size());

		CSeekDataset *vd = NULL;
		sizeD.push_back((float) ps);

		if(outputCoexpression || outputQueryCoexpression){
			vd = new CSeekDataset();
			string strFileStem = datasetName[i].substr(0, datasetName[i].find(".bin"));
			string strAvgPath = strPrepInputDirectory + "/" + strFileStem + ".gavg";
			string strPresencePath = strPrepInputDirectory + "/" + strFileStem + ".gpres";
			string strSinfoPath = strSinfoInputDirectory + "/" + strFileStem + ".sinfo";
			vd->ReadGeneAverage(strAvgPath);
			vd->ReadGenePresence(strPresencePath);
			vd->ReadDatasetAverageStdev(strSinfoPath);
			vd->InitializeGeneMap();

			fq = new CFullMatrix<float>();
			fq->Initialize(queryName.size(), pp->GetExperiments() - 2);
			for(k=0; k<queryName.size(); k++){
				int g = pp->GetGene(queryName[k]);
				if(g==-1){ //does not contain the gene in the dataset
					for(j=2; j<gs; j++){
						fq->Set(k, j-2, NaN);
						vecQ.push_back(fq->Get(k, j-2));
					}
					continue;
				}
				float *vv = pp->Get(g);
				for(j=2; j<gs; j++)
					fq->Set(k, j-2, vv[j]);
				if(!outputNormalized){
					for(j=2; j<gs; j++)
						vecQ.push_back(fq->Get(k, j-2));
				}
				
				//normalize
				float mean = 0;
				for(j=2; j<gs; j++)
					mean+=fq->Get(k, j-2);
				mean/=(float) (gs - 2);
				float stdev = 0;
				for(j=2; j<gs; j++)
					stdev+=(fq->Get(k, j-2) - mean) * (fq->Get(k, j-2) - mean);
				stdev/=(float) (gs - 2);
				stdev = sqrt(stdev);
				for(j=2; j<gs; j++){
					float t1 = fq->Get(k, j-2);
					fq->Set(k, j-2, (t1 - mean) / stdev);
				}

				if(outputNormalized){
					for(j=2; j<gs; j++)
						vecQ.push_back(fq->Get(k, j-2));
				}

			}
		}

		fprintf(stderr, "allocating space %d %d...\n", geneName.size(),
			pp->GetExperiments());
		ff = new CFullMatrix<float>();
		ff->Initialize(genes, ps);
		fprintf(stderr, "done allocating space.\n");

		for(k=0; k<geneName.size(); k++){
			int g = pp->GetGene(geneName[k]);
			if(g==-1){
				for(j=2; j<gs; j++){
					ff->Set(k, j-2, NaN);
					vecG.push_back(ff->Get(k, j-2));
				}
				continue;
			}
			float *vv = pp->Get(g);
			for(j=2; j<gs; j++)
				ff->Set(k, j-2, vv[j]);

			if(!outputNormalized){
				for(j=2; j<gs; j++){
					vecG.push_back(ff->Get(k, j-2));
				}
			}
			
			//normalize	
			float mean = 0;
			for(j=2; j<gs; j++)
				mean+=ff->Get(k, j-2);
			mean/=(float) (gs - 2);
			float stdev = 0;
			for(j=2; j<gs; j++)
				stdev+=(ff->Get(k, j-2) - mean) * (ff->Get(k, j-2) - mean);
			stdev/=(float) (gs - 2);
			stdev = sqrt(stdev);
			for(j=2; j<gs; j++){
				float t1 = ff->Get(k, j-2);
				ff->Set(k, j-2, (t1 - mean) / stdev);
			}

			if(outputNormalized){
				for(j=2; j<gs; j++){
					vecG.push_back(ff->Get(k, j-2));
				}
			}
		}

		if(outputCoexpression){
			int kk;
			for(k=0; k<geneName.size(); k++){
				int g = pp->GetGene(geneName[k]);
				if(g==-1){
					vCoexpression[k] = NaN;
					vecCoexpression.push_back(vCoexpression[k]);
					continue;
				}
				float avgP = 0;
				int totalQueryPresent = queryName.size();
				for(kk=0; kk<queryName.size(); kk++){
					int gg = pp->GetGene(queryName[kk]);
					if(gg==-1){
						totalQueryPresent--;
						continue;
					}
					float p = 0;
					for(j=2; j<gs; j++)
						p+= ff->Get(k, j-2)*fq->Get(kk, j-2);
					p/=(float)(gs-2);
					//fprintf(stderr, "%.2f sinfo: %.2f %.2f\ngene average %d: %.2f\n", p, vd->GetDatasetAverage(), 
					//	vd->GetDatasetStdev(), g, vd->GetGeneAverage(mapstrintGene[queryName[kk]]));
					//p = max((float) min(p - vd->GetGeneAverage(g), (float) 3.2), (float)-3.2);
					//get z-score (dataset wide)
					p = (p - vd->GetDatasetAverage()) / vd->GetDatasetStdev();
					//subtract hubbiness
					p = p - vd->GetGeneAverage(mapstrintGene[queryName[kk]]);
					avgP+=p;
				}
				if(totalQueryPresent==0)
					avgP = NaN;
				else
					avgP/=(float)(totalQueryPresent);
				
				vCoexpression[k] = avgP;
				vecCoexpression.push_back(avgP);
			}
		}

		if(outputQueryCoexpression){
			int kk=0;
			for(k=0; k<queryName.size(); k++){
				int g = pp->GetGene(queryName[k]);
				if(g==-1){
					vqCoexpression[k] = NaN;
					vecqCoexpression.push_back(vqCoexpression[k]);
					continue;
				}
				float avgP = 0;
				int totalQueryPresent = queryName.size() - 1;
				for(kk=0; kk<queryName.size(); kk++){
					if(kk==k) continue;
					int gg = pp->GetGene(queryName[kk]);
					if(gg==-1){
						totalQueryPresent--;
						continue;
					}
					float p = 0;
					for(j=2; j<gs; j++)
						p+= fq->Get(k, j-2)*fq->Get(kk, j-2);
					p/=(float)(gs-2);
					//get z-score (dataset wide)
					p = (p - vd->GetDatasetAverage()) / vd->GetDatasetStdev();
					//subtract hubbiness
					p = p - vd->GetGeneAverage(mapstrintGene[queryName[kk]]);
					avgP+=p;
				}
				if(totalQueryPresent==0)
					avgP = NaN;
				else
					avgP/=(float)(totalQueryPresent);
				
				vqCoexpression[k] = avgP;
				vecqCoexpression.push_back(avgP);
			}
		}
		
		if(outputCoexpression || outputQueryCoexpression){
			delete vd;
			delete fq;
		}

		delete ff;
	}

	if(CSeekNetwork::Send(new_fd, sizeD)!=0){
		fprintf(stderr, "Error sending messages\n");
	}

	if(outputExpression){
		if(CSeekNetwork::Send(new_fd, vecG)!=0){
			fprintf(stderr, "Error sending messages\n");
		}
	}

	if(outputQueryExpression){
		if(CSeekNetwork::Send(new_fd, vecQ)!=0){
			fprintf(stderr, "Error sending messages\n");
		}
	}

	if(outputCoexpression){
		if(CSeekNetwork::Send(new_fd, vecCoexpression)!=0){
			fprintf(stderr, "Error sending messages\n");
		}
	}

	if(outputQueryCoexpression){
		if(CSeekNetwork::Send(new_fd, vecqCoexpression)!=0){
			fprintf(stderr, "Error sending messages\n");
		}
	}

	for(i=0; i<NUM_DSET_MEMORY; i++){
		loaded[i] = 0;
	}

	THREAD_OCCUPIED[tid] = 0;

	pthread_mutex_unlock(&mutexGet);

	int ret = 0;
	close(new_fd);
	pthread_exit((void*)ret);
	//return void;
}

int main( int iArgs, char** aszArgs ) {
	static const size_t	c_iBuffer	= 1024;
#ifdef WIN32
	pthread_win32_process_attach_np( );
#endif // WIN32

	gengetopt_args_info	sArgs;

	if( cmdline_parser( iArgs, aszArgs, &sArgs ) ) {
		cmdline_parser_print_help( );
		return 1;
	}

	signal(SIGPIPE, SIG_IGN);
	size_t i;
	for(i=0; i<NUM_THREADS; i++){
		THREAD_OCCUPIED[i] = 0;
	}

	PORT = sArgs.port_arg;
	
	strPrepInputDirectory = sArgs.prep_arg;
	strSinfoInputDirectory = sArgs.sinfo_arg;

	vector<string> vecstrGeneID;
	vector<string> vecstrGenes;
	if(!CSeekTools::ReadListTwoColumns(sArgs.gene_arg, vecstrGeneID, vecstrGenes))
		return false;
	for(i=0; i<vecstrGenes.size(); i++)
		mapstrintGene[vecstrGenes[i]] = (int) i;

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
		if ((sockfd = socket(p->ai_family, p->ai_socktype,
			p->ai_protocol)) == -1) {
			perror("server: socket");
			continue;
		}
		if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &yes,
			sizeof(int)) == -1) {
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

	fprintf(stderr, "Start init.\n");

	pcl = new CPCL*[NUM_DSET_MEMORY];

	loaded = new char[NUM_DSET_MEMORY];

	available.clear();
	for(i=0; i<NUM_DSET_MEMORY; i++){
		available.push_back(i);
		loaded[i] = 0;
		pcl[i] = new CPCL();
	}

	fprintf(stderr, "Finished initializations.\n");
	while(1){
		sin_size = sizeof their_addr;
		new_fd = accept(sockfd, (struct sockaddr *) &their_addr, &sin_size);
		if(new_fd==-1){
			perror("accept");
			continue;
		}
		inet_ntop(their_addr.ss_family, get_in_addr(
			(struct sockaddr *)&their_addr), s, sizeof s);
		printf("server, got connection from %s\n", s);
		for(d=0; d<NUM_THREADS; d++){
			if(THREAD_OCCUPIED[d]==0) break;
		}
		if(d==NUM_THREADS){
			close(new_fd);
			continue;
		}

		THREAD_OCCUPIED[d] = 1;
		string mode;

		if(CSeekNetwork::Receive(new_fd, mode)!=0){
			fprintf(stderr, "Error: receiving message\n");
			close(new_fd);
			continue;
		}

		//mode, 5 digits
		//1 - output coexpression calculation (with gene hubbiness removed), 
		//require setting two groups of genes, one for query (1 or many), 
		//and one other gene (to calculate coexpression on)
		//2 - output gene-normalized expression instead of original expression
		//3 - output gene expression (can be normalized or unnormalized, depends on 2)
		//4 - output query expression
		//5 - output query coexpression (compared to 1, which is gene-to-query 
		//    coexpression. This is query-to-query coexpression)

		bool outputCoexpression = false;
		bool outputNormalized = false;
		bool outputExpression = false;
		bool outputQueryExpression = false;
		bool outputQueryCoexpression = false;
		if(mode[0]=='1') 
			outputCoexpression = true;
		if(mode[1]=='1') 
			outputNormalized = true;
		if(mode[2]=='1') 
			outputExpression = true;
		if(mode[3]=='1') 
			outputQueryExpression = true;
		if(mode[4]=='1') 
			outputQueryCoexpression = true;

		vector<string> dsetName, geneName, queryName;
		string qname, gname, dname;

		if(outputCoexpression || outputQueryCoexpression){
			if(	CSeekNetwork::Receive(new_fd, qname)!=0 ||
				CSeekNetwork::Receive(new_fd, gname)!=0 ||
				CSeekNetwork::Receive(new_fd, dname)!=0 ){
				fprintf(stderr, "Error: receiving message\n");
				close(new_fd);
				continue;
			}
			CMeta::Tokenize(dname.c_str(), dsetName);
			CMeta::Tokenize(gname.c_str(), geneName);
			CMeta::Tokenize(qname.c_str(), queryName);
		}
		else{
			if(	CSeekNetwork::Receive(new_fd, gname)!=0 ||
				CSeekNetwork::Receive(new_fd, dname)!=0 ){
				fprintf(stderr, "Error: receiving message\n");
				close(new_fd);
				continue;
			}
			CMeta::Tokenize(dname.c_str(), dsetName);
			CMeta::Tokenize(gname.c_str(), geneName);
		}
	
		thread_arg[d].threadid = d;
		thread_arg[d].new_fd = new_fd;
		thread_arg[d].geneName = geneName;
		thread_arg[d].datasetName = dsetName;
		thread_arg[d].queryName = queryName;
		thread_arg[d].outputNormalized = outputNormalized;
		thread_arg[d].outputCoexpression = outputCoexpression;
		thread_arg[d].outputExpression = outputExpression;
		thread_arg[d].outputQueryExpression = outputQueryExpression;
		thread_arg[d].outputQueryCoexpression = outputQueryCoexpression;

		fprintf(stderr, "Arguments: %d %d %s %s\n", d, new_fd, dname.c_str(), gname.c_str());
		if(outputCoexpression){
			fprintf(stderr, "Arguments: %s\n", qname.c_str());
		}

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
