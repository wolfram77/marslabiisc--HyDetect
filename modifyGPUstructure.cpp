#include "defs.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "graphHOST.h"
#include "graphGPU.h"
#include "communityGPU.h"
#include <iostream>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <string>
#include "Controller.h"
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <string>
#include <iostream>
#include <limits>
#include <string.h>

#include <unistd.h>
#include <cassert>
#include <inttypes.h>
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <cassert>
#include <inttypes.h>
#include <vector>
#include <omp.h>

Community modifyGPUstructure(Community *dev1_community,unsigned int *statIndices,unsigned int*edges,bool *dirtyg,bool *dirtygpu,long *c,int total,int NV,move1* mo,graph *Gn,unsigned int mid,long *c1,int it,long *cact){
unsigned int newV1=0;
//Community *cc1;
unsigned int edgec1=0;
for(int i=0;i<total-NV;i++)
{
        if(dirtyg[i]==true && !dirtygpu[i]){
        newV1++; 
}

}
cout<<"dirty1"<<" "<<newV1<<endl;
//return *dev1_community;
int k,jj;

 k=0;jj=0;

thrust::host_vector<unsigned int> h_vec((*dev1_community).g.nb_nodes);


unsigned int *stat=new unsigned int[(*dev1_community).g.nb_nodes];
std::copy(((*dev1_community).g.indices).begin(),((*dev1_community).g.indices).end(), h_vec.begin());
	int i=0;
	for(int i=0;i<h_vec.size();i++)
	{
		stat[i]=h_vec[i+1];
	}

unsigned int* edge1=new unsigned int[(*dev1_community).g.nb_links];
thrust::host_vector<unsigned int> h_vec1((*dev1_community).g.nb_links);
std::copy(((*dev1_community).g.links).begin(),((*dev1_community).g.links).end(), h_vec1.begin());

thrust::host_vector<unsigned int> h_vec2((*dev1_community).g.nb_links);

float * weight1=new float[(*dev1_community).g.nb_links];
std::copy(((*dev1_community).g.weights).begin(),((*dev1_community).g.weights).end(), h_vec2.begin());

for(int i=0;i<h_vec1.size();i++)
{
	edge1[i]=h_vec1[i];

}
for(int i=0;i<h_vec1.size();i++)
{
        weight1[i]=h_vec2[i];
}


cout<<"okk done copying"<<endl;
unsigned int edgestoadd;

unsigned    int *stat1=new unsigned int[(*dev1_community).g.nb_nodes+newV1];
	for(int i=0;i<(*dev1_community).g.nb_nodes;i++){
		if(i!=0)
		stat1[i]=stat[i]-stat[i-1];
		else if(i==0)
		stat1[i]=stat[i];

		}

for(int i=(*dev1_community).g.nb_nodes;i<(*dev1_community).g.nb_nodes+newV1;i++)	
		stat1[i]=0;
thrust::host_vector<float> w((*dev1_community).g.nb_links);
int k1=0;int k2=0;
std::copy((*dev1_community).g.weights.begin(),(*dev1_community).g.weights.end(),w.begin());


int s=0;
unsigned int adj1,adj2;

 k=0;
 adj1=0,adj2=0;
#pragma omp parallel for
for(unsigned int i=0;i<total-NV;i++)
{



        if(dirtyg[i]==true && !dirtygpu[i])

        {
             if(i!=0){
       		 adj1=statIndices[i-1];
        	 adj2=statIndices[i];
                	}
        	else if(i==0)
        	{
                	adj1=0;
                	adj2=statIndices[i];

        	}
   

                for(unsigned int j=adj1;j<adj2;j++)
                {
                        unsigned int tail = edges[j];
                        __sync_fetch_and_add(&k, 1);
                }
        }
}

int val=2*k;
unsigned int* edge2=new unsigned int[(*dev1_community).g.nb_links+2*k+mo->edg];
float* weight2=new float[(*dev1_community).g.nb_links+2*k+mo->edg];
for(int i=0;i<(*dev1_community).g.nb_links;i++)
{
        edge2[i]=edge1[i];
        weight2[i]=weight1[i];
}
for(int i=(*dev1_community).g.nb_links;i<(*dev1_community).g.nb_links+edgestoadd;i++)
{
        edge2[i]=0;
        weight2[i]=0;
}

unsigned int cc=(*dev1_community).g.nb_nodes;
unsigned int totalV=(*dev1_community).g.nb_nodes+newV1;
map<unsigned int,unsigned int>** cluPtrIn = (map<unsigned int,unsigned int>**) malloc (totalV*sizeof(map<unsigned int,unsigned int>*));
  assert(cluPtrIn != 0);
map<unsigned int,unsigned int>** count = (map<unsigned int,unsigned int>**) malloc (totalV*sizeof(map<unsigned int,unsigned int>*));
assert(count != 0);
cout<<"*********************"<<endl;
#pragma omp parallel for
  for (long i=0; i<totalV; i++) {
        cluPtrIn[i] = new map<unsigned int,unsigned int>();
        count[i]=new map<unsigned int,unsigned int>();
        if(i<cc){
	if(i!=0){
        adj1=stat[i-1];
        adj2=stat[i];
		}
	else if(i==0)
	{
		adj1=0;
		adj2=stat[i];

	}
       for(unsigned int j=adj1;j<adj2;j++)
                {
                        unsigned int tail=edge1[j];
                        (*(cluPtrIn[i]))[tail] =w[j];

                }
        }
(*(cluPtrIn[i]))[i]=0;

        }
k=0;
cout<<"????????????????????????????"<<endl;
//checked
unsigned int ii=0;
jj=0;
//cout<<stat[cc-1]<<" "<<"......"<<endl;
#pragma omp parallel for
for(unsigned int i=0;i<total-NV;i++)
{



        if(dirtyg[i]==true && !dirtygpu[i])

        {
                int k1=0;
                map<unsigned int, unsigned int>::iterator localIterator;
 		if(i!=0){
        		adj1=statIndices[i-1];
        		adj2=statIndices[i];
               		 }
        	else if(i==0)
        	{
                	adj1=0;
                	adj2=statIndices[i];

        	}
               

       for(unsigned int j=adj1;j<adj2;j++)
                {
                        unsigned int tail = edges[j];

                        localIterator = cluPtrIn[ii+cc]->find(c[edges[j]]);
                        if( localIterator != cluPtrIn[ii+cc]->end() ) {

                                 __sync_fetch_and_add(&(*(cluPtrIn[ii]))[c[edges[j]]] ,1);
                                int pos= (*(count[ii+cc]))[c[edges[j]]];
                                weight2[stat[cc-1]+pos]=(*(cluPtrIn[ii+cc]))[c[edges[j]]];
			        weight2[stat[cc-1]+pos+1]=(*(cluPtrIn[ii+cc]))[c[edges[j]]];

                                }

                        else{
				
                                (*(cluPtrIn[ii+cc]))[c[edges[j]]] = 1;
		
                                edge2[stat[cc-1]+k]=ii+(cc);
				
				weight2[stat[cc-1]+k]=(*(cluPtrIn[ii+cc]))[c[edges[j]]];
				(*(count[ii+cc]))[c[edges[j]]]=k;

				__sync_fetch_and_add(&k,1);
				 __sync_fetch_and_add(&k1,1);
                                edge2[stat[cc-1]+k]=c[edges[j]];

                                weight2[stat[cc-1]+k]=(*(cluPtrIn[ii+cc]))[c[edges[j]]];

                                __sync_fetch_and_add(&k,1);
                                __sync_fetch_and_add(&k1,1);
                                }



                }

                stat1[cc+jj]=k1;
                k1=0;
                __sync_fetch_and_add(&jj, 1);

                __sync_fetch_and_add(&ii, 1);
        }

}

cout<<"value"<< " "<<k<<endl; k2=0;

//return *dev1_community;

 jj=0;

GraphHOST modified_graph;
unsigned int gpunodes = (*dev1_community).g.nb_nodes+newV1;

jj=0;ii=0;


/*#pragma omp parallel for
  for (long i=0; i<totalV; i++){
        delete cluPtrIn[i];
        delete count[i];}
  free(cluPtrIn);
  free(count);*/
//return *dev1_community;
modified_graph.nb_nodes=mo->vertex+gpunodes;
unsigned int*stat2=(unsigned int*)malloc(sizeof(unsigned int)*modified_graph.nb_nodes);
for(unsigned int i=0;i<modified_graph.nb_nodes;i++)
	stat2[i]=0;
modified_graph.degrees.resize(modified_graph.nb_nodes);
cout<<"::::::::::::::::"<<endl;

for(unsigned int i=1;i<gpunodes;i++)
{  
//	stat1[i]+=stat1[i-1];
	stat2[i]=stat1[i];
}
stat2[0]=stat1[0];

 s=0;
for(int i=0;i<5;i++)
	cout<<"stat"<<" "<<stat2[i]<<endl;

modified_graph.nb_links = (*dev1_community).g.nb_links+val+mo->edg;

cout<<"no of previous edges"<<" "<<(*dev1_community).g.nb_links+val<<endl;
modified_graph.links.resize(modified_graph.nb_links);
unsigned int *edge3=(unsigned int *)malloc(sizeof(unsigned int)*modified_graph.nb_links);
unsigned int *weight3=(unsigned int *)malloc(sizeof(unsigned int)*modified_graph.nb_links);

for(unsigned int i=0;i<(*dev1_community).g.nb_links+val;i++){
	edge3[i]=edge2[i];
	weight3[i]=weight2[i];
//	cout<<"edge="<< edge3[i]<<endl;
}
for(int i=(*dev1_community).g.nb_links+val;i<modified_graph.nb_links;i++)    {
                edge3[i]=0;
		weight3[i]=0;
}

int j=0;
int toc=0;
int novt=0;
if(it==1){
for(int i=0;i<mo->edg;i++)
	{	if(mo->edgesa[i]>NV-1){
		edge3[(*dev1_community).g.nb_links+val+j]=c[mo->edgesa[i]-NV];
//		if(c[mo->edgesa[i]-NV]==1) toc++;
		weight3[(*dev1_community).g.nb_links+val+j]=mo->weighta[i];
		j++;
		stat2[c[mo->edgesa[i]-NV]]++;
}
			
		else{
		edge3[(*dev1_community).g.nb_links+val+j]=mo->edgesa[i]+gpunodes;
		                weight3[(*dev1_community).g.nb_links+val+j]=mo->weighta[i];
			j++;
		stat2[mo->edgesa[i]+gpunodes]++;
		
}
	
	}}
else if(it==0){
j=0; novt=0;
cout<<"no of edge"<<" "<<mo->edg<<endl;
for(int i=0;i<Gn->numVertices;i++)
{
edge3[(*dev1_community).g.nb_links+val+j]=gpunodes+mo->edgesa[i];
weight3[(*dev1_community).g.nb_links+val+j]=mo->weighta[i];
 stat2[mo->edgesa[i]+gpunodes]++;
                        j++;

edge3[(*dev1_community).g.nb_links+val+j]=gpunodes+mo->edgesa[i];
weight3[(*dev1_community).g.nb_links+val+j]=mo->weighta[i];
 stat2[mo->edgesa[i]+gpunodes]++;
                        j++;
}
for(int i=Gn->numVertices;i<(mo->edg)/2;i++)
        {      // if(mo->edgesa[i]<=NV-1){
	//	edge3[(*dev1_community).g.nb_links+val+j]=mo->edgesa[i]+gpunodes;j++;
          //      stat2[mo->edgesa[i]+gpunodes]++;
	//	cout<<"check1"<<" "<<mo->edgesa[i]+gpunodes<<endl;
                edge3[(*dev1_community).g.nb_links+val+j]=gpunodes+mo->edgesa[i];
		                weight3[(*dev1_community).g.nb_links+val+j]=mo->weighta[i];
			j++;
	//	cout<<"check2"<<" "<<c[mo->edgesa[i]]<<" "<<mo->edgesa[i]<<endl;
		stat2[mo->edgesa[i]+gpunodes]++; //}
//}
                        
               // else{
               	int firstval=cact[mo->edgesa[i]];
		int act=c[firstval];
                edge3[(*dev1_community).g.nb_links+val+j]=act;
		                weight3[(*dev1_community).g.nb_links+val+j]=mo->weighta[i];
		j++;
                stat2[act]++;
		
//}
        
        }
//modified_graph.nb_nodes=novt+gpunodes;
//cout<<"actual no"<<modified_graph.nb_nodes<<endl;

//modified_graph.degrees.resize(modified_graph.nb_nodes);

}
//modified_graph.nb_nodes=novt+gpunodes;
cout<<"tocheck"<<" "<<toc<<endl;
/*for(int i=0;i<modified_graph.nb_nodes;i++)
	stat2[i]=0;

stat2[i]=0;
for(int i=0;i<modified_graph.nb_links;i++)
{
	stat2[edge3[i]]=stat2[edge3[i]]+1;


}*/
//return *dev1_community;

//for(int i=0;i<10;i++)
  //      cout<<"stat2"<<" "<<stat2[i]<<endl;

cout<<"but working"<<" "<<modified_graph.nb_nodes<<endl;
for(int i=1;i<modified_graph.nb_nodes;i++)
	stat2[i]=stat2[i]+stat2[i-1];

/*for(int i=0;i<modified_graph.nb_nodes;i++)
        cout<<"stat3"<<" "<<stat2[i]<<endl;
for(int i=0;i<modified_graph.nb_links;i++)
	cout<<"edge3"<<" "<<edge3[i]<<endl;*/
cout<<"check final statIndices"<<endl;

/*for(int i=0;i<modified_graph.nb_nodes;i++)
        cout<<stat2[i]<<" ";
cout<<endl;
*/
/*for(int i=0;i<modified_graph.nb_links;i++)
	cout<<"edge"<< " "<<edge3[i]<<endl;*/
std::copy(stat2, stat2 + (modified_graph.nb_nodes), modified_graph.degrees.begin());
//for(int i=0;i<10;i++)
//	cout<<stat2[i]<<" ";

cout<<"::::::::"<<endl;

std::copy(edge3, edge3 + modified_graph.nb_links, modified_graph.links.begin());

for(int i=0;i<10;i++)
	cout<<edge3[i]<<" ";
modified_graph.weights.resize(modified_graph.nb_links);
//std::copy(weight3, weight3 + modified_graph.nb_links, modified_graph.weights.begin());

ofstream fout1("d1.txt");
ofstream fout("e1.txt");

/*statIndices=(unsigned int *)malloc(sizeof(unsigned int)*(modified_graph.nb_nodes));
edges=(unsigned int *)malloc(sizeof(unsigned int)*(modified_graph.nb_links));
for(int i=0;i<modified_graph.nb_nodes;i++)
	statIndices[i]=stat2[i];
for(int i=0;i<modified_graph.nb_links;i++)
	edges[i]=edge3[i];*/
Community *dev11_community;
Community dev12_community(modified_graph, -1, 0.0001);
dev11_community=&dev12_community;
//free(edge3);
//free(stat2);
cout<<"........."<<endl;
double cur_mod = -1.0;
bool TEPS = true;
//return;
	bool islastRound = false;
	int szSmallComm = 100000;
cudaStream_t *streams = NULL;
	int n_streams = 8;


	bool isGauss=true;
bool *dirty1=(bool *)malloc(sizeof(bool)*(modified_graph.nb_nodes));
bool *dirty2=(bool *)malloc(sizeof(bool)*(modified_graph.nb_nodes));
if(it==0){
 c1=(long *)malloc(sizeof(long)*(modified_graph.nb_nodes));
}
for(int i=0;i<modified_graph.nb_nodes;i++)
{	c1[i]=-1;
	dirty1[i]=false;
	dirty2[i]=false;

}
/*for(int i=0;i<modified_graph.nb_nodes;i++)
	cout<<stat2[i]<<" ";
cout<<endl;

for(int i=0;i<modified_graph.nb_links;i++)
	cout<<edge3[i]<<" ";
cout<<endl;
*/
//return *dev1_community;
//return *dev1_community;
int a=gpuonly(modified_graph,c1,stat2,edge3,dev11_community,dirty1,dirty2,1,Gn,mid);
cout<<(*dev11_community).g.nb_nodes<<" "<<(*dev11_community).g.nb_links<<endl;
/*for(int i=0;i<modified_graph.nb_nodes;i++)
	cout<<"communitya"<<" "<<c1[i]<<endl;*/
/*#pragma omp parallel for
  for (long i=0; i<totalV; i++){
        delete cluPtrIn[i];
        delete count[i];}
  free(cluPtrIn);
  free(count);
//return *dev11_community;
free(stat);
free(edge1);
free(weight1);
free(stat2);
free(edge3);
free(edge2);
free(weight2);*/
cout<<"actual one"<<modified_graph.nb_nodes<<endl;
return *dev11_community;
}
