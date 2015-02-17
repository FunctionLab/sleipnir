/*****************************************************************************
* This file is provided under the Creative Commons Attribution 3.0 license.
*
* You are free to share, copy, distribute, transmit, or adapt this work
* PROVIDED THAT you attribute the work to the authors listed below.
* For more information, please see the following web page:
* http://creativecommons.org/licenses/by/3.0/
*
* This file is a part of SEEK (Search-based exploration of expression compendium)
* which is authored and maintained by: Qian Zhu (qzhu@princeton.edu)
*
* If you use this file, please cite the following publication:
* Qian Zhu, Aaron K Wong, Arjun Krishnan, Miriam R Aure, Alicja Tadych, 
* Ran Zhang, David C Corney, Casey S Greene, Lars A Bongo, 
* Vessela N Kristensen, Moses Charikar, Kai Li & Olga G Troyanskaya
* "Targeted exploration and analysis of large cross-platform human 
* transcriptomic compendia" Nat Methods (2015)
*
* This file is a component of the Sleipnir library for functional genomics,
* authored by:
* Curtis Huttenhower (chuttenh@princeton.edu)
* Mark Schroeder
* Maria D. Chikina
* Olga G. Troyanskaya (ogt@princeton.edu, primary contact)
*
* If you use this library for development, or use any other Sleipnir executable
* tools, please also cite the following publication:
* Curtis Huttenhower, Mark Schroeder, Maria D. Chikina, and
* Olga G. Troyanskaya.
* "The Sleipnir library for computational functional genomics"
*****************************************************************************/

#include "seeknetwork.h"

namespace Sleipnir {

//send message header:
//byte #1-4:
// size of one element of the data structure (in unsigned int)
//byte #5-8:
// (in unsigned int):
// total number of elements to be sent, 1 if a single value, or else
// the size of the array
//byte #9-:
// actual data

int CSeekNetwork::Send(int new_fd, const string &str){
	char *s = (char*)malloc(2*sizeof(unsigned int));
	unsigned int *p = (unsigned int*) s;
	*p = sizeof(char); p++;
	*p = str.length()+1;
	if(CSeekNetwork::Send(new_fd, s, 2*sizeof(unsigned int))==-1){
		fprintf(stderr, "Bad 1\n");
		free(s);
		return -1;
	}
	free(s);
	
	char *c = (char*)malloc(str.length()+1);
	strcpy(c, str.c_str());
	if(CSeekNetwork::Send(new_fd, c, str.length()+1)==-1){
		fprintf(stderr, "Bad 2\n");
		free(c);
		return -1;
	}
	free(c);
	return 0;
}

int CSeekNetwork::Send(int new_fd, const vector<float> &f){
	char *s = (char*)malloc(2*sizeof(unsigned int));
	unsigned int *p = (unsigned int*) s;
	*p = sizeof(float); p++;
	*p = f.size();
	if(CSeekNetwork::Send(new_fd, s, 2*sizeof(unsigned int))==-1){
		fprintf(stderr, "Bad 1\n");
		free(s);
		return -1;
	}
	free(s);
	
	char *c = (char*)malloc(sizeof(float)*f.size());
	float *fp = (float*)c;
	int i;
	for(i=0; i<f.size(); i++){
		*fp = f[i];
		fp++;
	}
	if(CSeekNetwork::Send(new_fd, c, sizeof(float)*f.size())==-1){
		fprintf(stderr, "Bad 2\n");
		free(c);
		return -1;
	}
	free(c);
	return 0;
}

int CSeekNetwork::Send(int new_fd, char *c, int size){
	int tmp_size = size;
	int beg = 0;

	int r  = -1;
	while(1){
		char *p = (char*)malloc(tmp_size);
		Copy(p, c, beg, tmp_size);
		r = send(new_fd, p, tmp_size, 0);
		if(r==-1){
			fprintf(stderr, "client exits");
			break;
		}
		if(r==tmp_size){
			break;
		}
		tmp_size = size - tmp_size;
		beg = beg + r;
		free(p);
	}

	return r;
}

void CSeekNetwork::Clear(char *b, int size){
	int i = 0;
	for(i=0; i<size; i++){
		b[i] = '\0';
	}
}

int CSeekNetwork::Copy(char *d, char *s, int beg, int num){
    int i;
    for(i=0; i<num; i++){
        d[beg+i] = s[i];
    }
    return beg+num;
}

int CSeekNetwork::Receive(int new_fd, vector<float> &f){
	char *ar = (char*)malloc(4);
	char *p = &ar[0];
	int tmp_size = 4;
	int receive_size = -1;
	while(1){
		receive_size = recv(new_fd, p, tmp_size, 0);
		if(receive_size==tmp_size){
			break;
		}
		else if(receive_size==-1){
			fprintf(stderr, "client exits\n");
			break;
		}
		tmp_size = tmp_size - receive_size;
		p+=receive_size;
	}
	if(receive_size==-1)
		return -1;

	int *iP = (int*)ar;
	int length = *iP;
	f.resize(length);

	char *cStr = (char*)malloc(length*sizeof(float));
	tmp_size = length*sizeof(float);
	receive_size = -1;
	p = &cStr[0];

	while(1){
		receive_size = recv(new_fd, p, tmp_size, 0);
		if(receive_size==tmp_size){
			break;
		}
		else if(receive_size==-1){
			fprintf(stderr, "client exits\n");
			break;
		}
		tmp_size = tmp_size - receive_size;
		p+=receive_size;
	}
	if(receive_size==-1)
		return -1;

	int i;
	float *fp = (float*)cStr;
	for(i=0; i<length; i++){
		f[i] = *fp;
		fp++;
	}
	free(ar);
	free(cStr);
	return 0;
}


int CSeekNetwork::Receive(int new_fd, string &s){
	char *ar = (char*)malloc(4);
	char *p = &ar[0];
	int tmp_size = 4;
	int receive_size = -1;
	while(1){
		receive_size = recv(new_fd, p, tmp_size, 0);
		if(receive_size==tmp_size){
			break;
		}
		else if(receive_size==-1){
			fprintf(stderr, "client exits\n");
			break;
		}
		tmp_size = tmp_size - receive_size;
		p+=receive_size;
	}

	if(receive_size==-1){
		return -1;
	}
	int *iP = (int*)ar;
	int length = *iP;
	char *cStr = (char*)malloc(length);
	tmp_size = length;
	receive_size = -1;
	p = &cStr[0];
	
	while(1){
		receive_size = recv(new_fd, p, tmp_size, 0);
		if(receive_size==tmp_size){
			break;
		}
		else if(receive_size==-1){
			fprintf(stderr, "client exits\n");
			break;
		}
		tmp_size = tmp_size - receive_size;
		p+=receive_size;
	}
	if(receive_size==-1){
		return -1;
	}

	s = cStr; //copy result into string
	free(ar);
	free(cStr);
	return 0;
}

}

