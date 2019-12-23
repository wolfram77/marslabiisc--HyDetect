/*

    Copyright (C) 2016, University of Bergen

    This file is part of Rundemanen - CUDA C++ parallel program for
    community detection

    Rundemanen is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Rundemanen is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Rundemanen.  If not, see <http://www.gnu.org/licenses/>.
    
    */
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










graph* modifyCPUstructure(graph *Gnew,graph *G1,bool *dirty,unsigned int *C_orig,move2 *mo1,unsigned int mid,unsigned int node,bool *dirty3){
//cout<<"inside modification"<<endl;
//cout<<mo1->vertex<<" "<<mo1->edgeno<<endl;
unsigned int newV=0;
unsigned int edgec1=0;
displayGraphCharacteristics(Gnew);

//cout<<G1->numVertices<<" "<<G1->numEdges<<endl;

for(int i=0;i<G1->numVertices;i++)
{
        if(dirty[i]){
        newV++;
     
		}

}
  unsigned int    *vtxPtrIn    = G1->edgeListPtrs;
  edge    *vtxIndIn    = G1->edgeList;
  unsigned int    *vtxPtrOut    = Gnew->edgeListPtrs;
  edge    *vtxIndOut    = Gnew->edgeList;
//  map<unsigned int,unsigned int>** cluPtrIn = (map<unsigned int,unsigned int>**) malloc (Gnew->numVertices*sizeof(map<unsigned int,unsigned int>*));
//  assert(cluPtrIn != 0);
 /* #pragma omp parallel for
  for (long i=0; i<Gnew->numVertices; i++) {
	cluPtrIn[i] = new map<unsigned int,unsigned int>();

	unsigned int adj1=vtxPtrOut[i];
	unsigned int adj2=vtxPtrOut[i+1];
	for(unsigned int j=adj1;j<adj2;j++)
		{
			unsigned int tail=vtxIndOut[j].tail; 
			(*(cluPtrIn[i]))[j] =(unsigned int)vtxIndOut[j].weight;;

		}

	}		*/
	
int k=0;
//cout<<"doubtful vertices ="<<newV<<endl;
#pragma omp parallel for
for(unsigned int i=0;i<G1->numVertices;i++)
{



        if(dirty[i]==true)

        {
		unsigned int adj1=vtxPtrIn[i];
                unsigned int adj2=vtxPtrIn[i+1];
	
		for(unsigned int j=adj1;j<adj2;j++)
                { 	
			unsigned int tail = vtxIndIn[j].tail;
			__sync_fetch_and_add(&k, 1); 
		}
	}
}	

edgec1=G1->numEdges;
unsigned int *edgeP=(unsigned int*)malloc(sizeof(unsigned int)*(newV+Gnew->numVertices+1+mo1->vertex));
edgeP[0]=Gnew->edgeListPtrs[0];
for(unsigned int i=1;i<Gnew->numVertices+1;i++)
{
        edgeP[i]=Gnew->edgeListPtrs[i]-Gnew->edgeListPtrs[i-1];
	
}
for(unsigned int i=Gnew->numVertices+1;i<Gnew->numVertices+1+newV+mo1->vertex;i++)
	edgeP[i]=0;
edge *edge2=(edge*)malloc(sizeof(edge)*((Gnew->numEdges*2)+2*k+(mo1->edgeno)*2));
//unsigned int val=Gnew->numVertices+k+mo1->edgeno;
//edge *edge2=(edge*)malloc(sizeof(edge)*(Gnew->edgeListPtrs[Gnew->numVertices]+edgec1));
/*cout<<"initial situation"<<endl;
for(int i=0;i<Gnew->numVertices+1;i++)
	cout<<Gnew->edgeListPtrs[i]<<" ";
cout<<endl;*/
for(unsigned int i=0;i<(Gnew->edgeListPtrs[Gnew->numVertices]);i++)
        {

        edge2[i].head=Gnew->edgeList[i].head;
        edge2[i].tail=Gnew->edgeList[i].tail;
        edge2[i].weight=Gnew->edgeList[i].weight;
//ut<<Gnew->edgeList[i].head<<" "<<Gnew->edgeList[i].tail<<" "<<Gnew->edgeList[i].weight<<endl;
        }

 k=0;int jj=0;
int *flag;
vector<unsigned int> pos;
pos.resize(Gnew->numVertices+1);
std::fill(pos.begin(),pos.begin()+pos.size(),0);
cout<<"while updating"<<endl;
//ofstream ch("check.txt");
int ii=0;
/*ofstream f1("1.txt");
for(int i=0;i<=Gnew->numVertices;i++)
        f1<<Gnew->edgeListPtrs[i]<<endl;
*/
unsigned int cc=Gnew->numVertices;
Gnew->numVertices=Gnew->numVertices+newV+mo1->vertex;
map<unsigned int,unsigned int>** cluPtrIn = (map<unsigned int,unsigned int>**) malloc (Gnew->numVertices*sizeof(map<unsigned int,unsigned int>*));
  assert(cluPtrIn != 0);
map<unsigned int,unsigned int>** count = (map<unsigned int,unsigned int>**) malloc (Gnew->numVertices*sizeof(map<unsigned int,unsigned int>*));

#pragma omp parallel for
  for (long i=0; i<Gnew->numVertices; i++) {
        cluPtrIn[i] = new map<unsigned int,unsigned int>();
	count[i]=new map<unsigned int,unsigned int>();
	if(i<cc){
        unsigned int adj1=vtxPtrOut[i];
        unsigned int adj2=vtxPtrOut[i+1];
        for(unsigned int j=adj1;j<adj2;j++)
                {
                        unsigned int tail=vtxIndOut[j].tail; 
                        (*(cluPtrIn[i]))[tail] =(unsigned int)vtxIndOut[j].weight;;

                }
	}
(*(cluPtrIn[i]))[i]=0; 

        }
//cout<<"Okk!"<<endl; 
k=0;
#pragma omp parallel for
for(unsigned int i=0;i<G1->numVertices;i++)
{



        if(dirty[i]==true)

        {	
		int k1=0;
		map<unsigned int, unsigned int>::iterator localIterator;	
                unsigned int adj1=vtxPtrIn[i];
                unsigned int adj2=vtxPtrIn[i+1];
		
        for(unsigned int j=adj1;j<adj2;j++)
                {	
			unsigned int tail = vtxIndIn[j].tail; 
                       
			localIterator = cluPtrIn[ii]->find(C_orig[G1->edgeList[j].tail]);
			if( localIterator != cluPtrIn[ii]->end() ) {
				 
				 __sync_fetch_and_add(&(*(cluPtrIn[ii]))[C_orig[G1->edgeList[j].tail]] ,1);
				int pos= (*(count[ii]))[C_orig[G1->edgeList[j].tail]];
				edge2[Gnew->edgeListPtrs[cc]+pos].weight=(*(cluPtrIn[ii]))[C_orig[G1->edgeList[j].tail]];
				}
				 
			else{
				
				(*(cluPtrIn[ii]))[C_orig[G1->edgeList[j].tail]] = 1;
				edge2[Gnew->edgeListPtrs[cc]+k].head=ii+(cc);
                        	edge2[Gnew->edgeListPtrs[cc]+k].tail=C_orig[G1->edgeList[j].tail];
			
	                                edge2[Gnew->edgeListPtrs[cc]+k].weight=(*(cluPtrIn[ii]))[C_orig[G1->edgeList[j].tail]];
				
					(*(count[ii]))[C_orig[G1->edgeList[j].tail]]=k;
				__sync_fetch_and_add(&k,1);
				__sync_fetch_and_add(&k1,1);
				}
			                    

		
                }
	
                edgeP[cc+jj+1]=k1;
		k1=0;
                __sync_fetch_and_add(&jj, 1);

		__sync_fetch_and_add(&ii, 1);
        }

}

cout<<"value"<< " "<<k<<endl;int k2=0;
#pragma omp parallel for
  for (long i=0; i<Gnew->numVertices; i++){
	delete cluPtrIn[i];	
	delete count[i];}
  free(cluPtrIn);				
  free(count);
//cout<<"done1"<<endl;
/*for(int i=Gnew->edgeListPtrs[cc];i<Gnew->edgeListPtrs[cc]+10;i++)
        cout<<edge2[i].head<<" "<<edge2[i].tail<<" "<<edge2[i].weight<<endl;
*/
//cout<<k<<" "<<newV<<endl;
unsigned int prevV=cc;
for(int i=1;i<=newV+prevV;i++)
	edgeP[i]+=edgeP[i-1];

//ut<<edgeP[prevV]<<" "<<edgeP[prevV+newV-1]<<endl;
//cout<<"ok"<<endl;
//cout<<"mo"<<" "<<mo1->edgeno<<endl;
Gnew->numVertices=prevV+newV+mo1->vertex;
//cout<<"check="<<edgeP[prevV+newV];
//cout<<"After transfer"<<endl;
//cout<<"No of vertices="<<" "<<Gnew->numVertices<<" ";
Gnew->numEdges=edgeP[prevV+newV]+mo1->edgeno;

//cout<<"No of edges="<<" "<<Gnew->numEdges<<" "<<edgeP[Gnew->numVertices]<<" "<<Gnew->numVertices<<endl;

//unsigned int val=Gnew->numVertices;
//cout<<"okk?"<<endl;

//free(Gnew->edgeListPtrs);
Gnew->edgeListPtrs=(unsigned int*)malloc(sizeof(unsigned int)*(Gnew->numVertices+1));

//Gnew->edgeListPtrs=(unsigned int*)malloc(sizeof(unsigned int)*(Gnew->numVertices+1));

//cout<<"okk?"<<endl;
//Gnew->edgeListPtrs=(unsigned int*)malloc(sizeof(unsigned int)*(newV+Gnew->numVertices+1));
for(unsigned int i=0;i<prevV+newV+1;i++){
	Gnew->edgeListPtrs[i]=edgeP[i];
//	cout<<edgeP[i]<<endl;
	}
//cout<<"ch"<<" "<<Gnew->edgeListPtrs[Gnew->numVertices]<<endl;
//cout<<"No of edges="<<" "<<Gnew->numEdges<<" "<<k<<endl;;

/*for(unsigned int i=0;i<newV+prevV+1;i++){
        cout<<Gnew->edgeListPtrs[i]<<endl;}*/
//free(edgeP);

//cout<<"ok?"<<endl;*/

//Gnew->edgeList=(edge*)malloc(sizeof(edge)*(k+Gnew->edgeListPtrs[Gnew->numVertices]));
//cout<<"edgeListPtrs"<<" "<<Gnew->edgeListPtrs[Gnew->numVertices]<<endl;
Gnew->edgeList=(edge *)malloc(sizeof(edge)*(Gnew->edgeListPtrs[Gnew->numVertices]));
for(unsigned int i=0;i<Gnew->edgeListPtrs[Gnew->numVertices];i++)
	{

	Gnew->edgeList[i].head=edge2[i].head;
	Gnew->edgeList[i].tail=edge2[i].tail;
	Gnew->edgeList[i].weight=edge2[i].weight;
//	if(i>Gnew->edgeListPtrs[cc] && i<Gnew->edgeListPtrs[cc]+20)
//	cout<<Gnew->edgeList[i].head<<" "<<Gnew->edgeList[i].tail<<" "<<Gnew->edgeList[i].weight<<endl;


	}
//cout<<"donefinal"<<endl;
 double prevMod = 1;
  double currMod = -1;
double tmpTime;
int tmpItr=0;
//int p=Gnew->numVertices;
//Gnew->numVertices=Gnew->numVertices+newV+mo1->vertex;
//Gnew->numEdges=Gnew->edgeListPtrs[Gnew->numVertices]+mo1->edgeno;
//Gnew->edgeListPtrs=(unsigned int*)malloc(sizeof(unsigned int)*(Gnew->numVertices+1));
//for(int i=0;i<=p+newV;i++)
//	Gnew->edgeListPtrs[i]=edgeP[i];
int j=0;
if(mo1->vertex!=0){
//cout<<"ok2"<<mo1->edgeListPtrsM[0]<<endl;
for(int i=cc+newV+1;i<=Gnew->numVertices;i++)
{	Gnew->edgeListPtrs[i]+=mo1->edgeListPtrsM[j];
//	cout<<mo1->edgeListPtrsM[j]<<" ";
	j++;

}
}
/*ofstream f1("1.txt");
for(int i=0;i<=Gnew->numVertices;i++)
	f1<<Gnew->edgeListPtrs[i]<<endl;*/
//cout<<"done2"<<endl;
//cout<<"chance"<<endl;
//Gnew->numEdges=Gnew->edgeListPtrs[Gnew->numVertices]+mo1->edgeno;

j=0;
//cout<<"dirty1="<<" "<<Gnew->edgeListPtrs[prevV+newV]<<" "<<"dirty2"<<Gnew->edgeListPtrs[Gnew->numVertices]<<endl;
if(mo1->vertex>0){
for(unsigned int i=(Gnew->edgeListPtrs[cc+newV]);i<Gnew->numEdges;i++)
        {

        Gnew->edgeList[i].head=mo1->edgesM[j].head;
        Gnew->edgeList[i].tail=mo1->edgesM[j].tail;
        Gnew->edgeList[i].weight=mo1->edgesM[j].weight;
	j++;
        }
}
//cout<<"done3"<<endl;

C_orig=(unsigned int *)malloc(sizeof(unsigned int)*Gnew->numVertices);
unsigned int NV=Gnew->numVertices;
/*for(int i=0;i<NV+1;i++)
	cout<<Gnew->edgeListPtrs[i]<<" ";
cout<<endl;
*/
/*for( int i=0;i<Gnew->edgeListPtrs[NV];i++)
{
	cout<<Gnew->edgeList[i].head<<" "<<Gnew->edgeList[i].tail<<" "<<Gnew->edgeList[i].weight<<endl;


}*/

/*for(int i=0;i<Gnew->numVertices;i++)
{


cout<<Gnew->edgeList[i].head<<" "<<Gnew->edgeList[i].tail<<" "<<Gnew->edgeList[i].weight<< " "<<i<<endl;


}*/
graph *Gnew1=(graph *)malloc(sizeof(graph));
#pragma omp parallel for
	for (long i=0; i<NV; i++) {
  	   C_orig[i] = 0;
	}
bool *dirty1,*dirty2;
//dirty1=(bool *)malloc(sizeof(bool)*NV);
dirty2=(bool*)malloc(sizeof(bool)*NV);
for(int i=0;i<NV;i++)
{
//	dirty1[i]=false;
	dirty2[i]=false;

}
//cout<<"done4"<<endl;
displayGraphCharacteristics(Gnew);
//displayGraphEdgeList(Gnew);
graph *Gnew2=(graph *)malloc(sizeof(graph));
duplicateGivenGraph(Gnew,Gnew2);
	
//Gnew1=runMultiPhaseLouvainAlgorithm(Gnew, C_orig, 0, 1000, 0.0001, 0.0001, 6); 
Gnew1=cpuonly(Gnew,Gnew1,Gnew2,C_orig,dirty3,dirty2,1);
//Gnew2=runMultiPhaseLouvainAlgorithm(Gnew, C_orig,0, 1000, 0.001, 0.001, 6);

unsigned int new1=0;
/*for(int i=0;i<Gnew->numVertices;i++)
	{

		if(dirty3[i])
			new1++;
	}
cout<<"**********"<<endl;
unsigned int NV1=NV+new1;
move1 *mo=new move1[NV1];
//void movefinal(graph *Gnew,graph *Gnew1,bool *dirtycpu,move1 *mo,unsigned int *c,unsigned int mid,unsigned int node)

movefinal(Gnew,Gnew1,dirty1,dirty2,mo,C_orig,mid,node);*/
//void verticesToMoveToGPU(graph *G,bool *dirtycpu,move1 *mo,unsigned int *c,unsigned int mid,unsigned int node)

//final vertices t( move to GPU
//currMod = parallelLouvianMethod(Gnew, C_orig,6, currMod, 0.001, &tmpTime, &tmpItr);
//cout<<"cpu current mod="<<currMod<<endl;
//Gnew->edgeList=edge2;*/
//free(edge2);
//e.clear();
return Gnew1;
//free(flag);
//Gnew->numVertices=Gnew->numVertices+newV-newV11;
//Gnew->numEdges=Gnew->edgeListPtrs[Gnew->numVertices];

pos.clear();
free(edgeP);
free(edge2);




}
