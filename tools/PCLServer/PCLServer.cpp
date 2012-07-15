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
int NUM_DSET_MEMORY = 50;
CPCL **pcl;
list<int> available;
char *loaded;
map<string, int> DNAME_MAP;
map<int, int> return_val;
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
    vector<string> datasetName;
    vector<string> geneName;
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
	int new_fd = my->new_fd;
	int tid = my->threadid;

	char *buf = (char*)malloc(5000);
	char *bb = (char*)malloc(5000);
	/*
	cl(buf, 5000);
	sprintf(buf, "Begin...\n");
	if(send_msg(new_fd, buf, 5000)==-1){
		THREAD_OCCUPIED[tid] = 0;
		pthread_exit(0);
	}*/
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
		pthread_mutex_lock(&mutexGet);
		n = GetOpenSlot();
		map<int, string>::const_iterator iterRM = DNAME_RMAP.find(n);
		if(iterRM!=DNAME_RMAP.end()){
			DNAME_MAP.erase(iterRM->second);
			DNAME_RMAP.erase(n);
		}
		DNAME_MAP[*iterS] = n;
		DNAME_RMAP[n] = *iterS;
		loaded[n] = 1;
		pthread_mutex_unlock(&mutexGet);

		fprintf(stderr, "acquired %d for dataset %s...\n", n, iterS->c_str());
		//pcl[n]->Reset();
		fprintf(stderr, "dataset reset\n");
		pcl[n]->Open((*iterS).c_str());
		fprintf(stderr, "dataset opened\n");
		vc.push_back(pcl[n]);
	}

	size_t j, k;
	vector<CFullMatrix<float> *> vC;
	int totNumExperiments = 0;
	int genes = geneName.size();
	int datasets = datasetName.size();

	fprintf(stderr, "Reading data...\n");
	
	for(i=0; i<vc.size(); i++){
		CPCL *pp = vc[i];
		fprintf(stderr, "allocating space %d %d...\n", geneName.size(),
			pp->GetExperiments());
		CFullMatrix<float> *ff = new CFullMatrix<float>();
		ff->Initialize(geneName.size(), pp->GetExperiments() - 2);
		totNumExperiments += pp->GetExperiments() - 2;
		fprintf(stderr, "done allocating space.\n");

		for(k=0; k<geneName.size(); k++){
			int g = pp->GetGene(geneName[k]);
			int gs = pp->GetExperiments();

			if(g==-1){
				for(j=2; j<gs; j++){
					ff->Set(k, j-2, CMeta::GetNaN());
				}
				continue;
			}
			float *vv = pp->Get(g);
			for(j=2; j<gs; j++){
				ff->Set(k, j-2, vv[j]);
			}
		}
		vC.push_back(ff);
	}

	int est = 4 + 2 + 2 + (genes * totNumExperiments + datasets) * 2;

	int BUFFER = 5000;

	char *c8 = (char*)malloc(8);
	int *p = (int*)c8; *p = est;
	short *u1 = (short*)&c8[4]; *u1 = (short) datasets;
	u1++; *u1 = (short) genes;
	send_msg(new_fd, c8, 8);
	free(c8);
	fprintf(stderr, "%d %d %d\n", est, datasets, genes);

	char *cx = (char*)malloc(2*datasets);
	u1 = (short*) &cx[0];
	for(i=0; i<vc.size(); i++){
		CPCL *pp = vc[i];
		*u1 = (short) pp->GetExperiments()-2;
		fprintf(stderr, "E%d ", *u1);
		u1++;
	}
	send_msg(new_fd, cx, 2*datasets);
	free(cx);
	fprintf(stderr, "\n");

	for(i=0; i<vc.size(); i++){
		CPCL *pp = vc[i];
		CFullMatrix<float> *fm = vC[i];
		int ps = pp->GetExperiments()-2;
		char *cp = (char*)malloc(2*ps*genes);
		short *up = (short*) &cp[0];
		//fprintf(stderr, "Dataset %d\n", i);
		for(k=0; k<genes; k++){
			//fprintf(stderr, "Gene %d: ", k);
			for(j=0; j<ps; j++){
				short sp = -1;
				if(CMeta::IsNaN(fm->Get(k, j))){
					sp = 32767;
				}else{
					sp = fm->Get(k, j) * 100;
					if(sp<-32000){
						sp = -32000;
					}else if(sp>32000){
						sp = 32000;
					}
				}
				*up = sp;
				//fprintf(stderr, "%d ", *up);
				up++;
			}
			//fprintf(stderr, "\n");
		}
		send_msg(new_fd, cp, 2*ps*genes);

		free(cp);
		delete vC[i];
	}
	free(buf);
	free(bb);
	vC.clear();

	for(i=0; i<NUM_DSET_MEMORY; i++){
		loaded[i] = 0;
	}

	THREAD_OCCUPIED[tid] = 0;
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

		char ar[4];
		recv(new_fd, ar, 4, 0);
		int *cn = (int*) ar; 
		int strLen = *cn;
		char *dname = (char*)malloc(strLen);
		recv(new_fd, dname, strLen, 0);

		recv(new_fd, ar, 4, 0);
		strLen = *cn;
		char *gname = (char*)malloc(strLen);
		recv(new_fd, gname, strLen, 0);

		vector<string> dsetName, geneName;
		CMeta::Tokenize(dname, dsetName);
		CMeta::Tokenize(gname, geneName);
	
		thread_arg[d].threadid = d;
		thread_arg[d].new_fd = new_fd;
		thread_arg[d].geneName = geneName;
		thread_arg[d].datasetName = dsetName;

		fprintf(stderr, "Arguments: %d %d %s %s\n", d, new_fd, dname, gname);

		free(dname);
		free(gname);

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
