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


unsigned int newV=0;
unsigned int edgec1=0;
cout<<G1->numVertices<<" "<<G1->numEdges<<endl;

for(int i=0;i<G1->numVertices;i++)
{
        if(dirty[i]){
        newV++;
     
		}

}
cout<<"doubtful vertices ="<<newV<<endl;
edgec1=G1->numEdges;
cout<<Gnew->numVertices<<" "<<Gnew->numEdges<<" "<<endl;
unsigned int *edgeP=(unsigned int*)malloc(sizeof(unsigned int)*(newV+Gnew->numVertices+1));
for(unsigned int i=0;i<Gnew->numVertices+1;i++)
{
        edgeP[i]=Gnew->edgeListPtrs[i];
	
}
for(unsigned int i=Gnew->numVertices+1;i<Gnew->numVertices+1+newV;i++)
	edgeP[i]=0;
edge *edge2=(edge*)malloc(sizeof(edge)*(Gnew->edgeListPtrs[Gnew->numVertices]+edgec1));
cout<<"initial situation"<<endl;
for(int i=0;i<Gnew->numVertices+1;i++)
	cout<<Gnew->edgeListPtrs[i]<<" ";
cout<<endl;
for(unsigned int i=0;i<(Gnew->edgeListPtrs[Gnew->numVertices]);i++)
        {

        edge2[i].head=Gnew->edgeList[i].head;
        edge2[i].tail=Gnew->edgeList[i].tail;
        edge2[i].weight=Gnew->edgeList[i].weight;
cout<<Gnew->edgeList[i].head<<" "<<Gnew->edgeList[i].tail<<" "<<Gnew->edgeList[i].weight<<endl;
        }

int k=0;int jj=0;
int *flag;
vector<unsigned int> pos;
pos.resize(Gnew->numVertices+1);
std::fill(pos.begin(),pos.begin()+pos.size(),0);
cout<<"while updating"<<endl;
//ofstream ch("check.txt");
int ii=0;
#pragma omp parallel for
for(unsigned int i=0;i<G1->numVertices;i++)
{



        if(dirty[i]==true)

        {	
		//cout<<"doubtful="<<" "<<i<<endl;
		flag=(int *)malloc(sizeof(int)*Gnew->numVertices);
		for(unsigned int k1=0;k1<Gnew->numVertices;k1++)
			flag[k1]=0;
                unsigned int adj1=G1->edgeListPtrs[i];
                unsigned int adj2=G1->edgeListPtrs[i+1];
		//cout<<adj1<<" "<<adj2<<endl;
	#pragma omp parallel for
         for(unsigned int j=adj1;j<adj2;j++)
                {	if(flag[C_orig[G1->edgeList[j].tail]]==0){
			
                        edge2[Gnew->edgeListPtrs[Gnew->numVertices]+k].head=ii+(Gnew->numVertices);
                        edge2[Gnew->edgeListPtrs[Gnew->numVertices]+k].tail=C_orig[G1->edgeList[j].tail]; //check for the condition where edge connected to same community,weight will be weight+1
			pos.at(C_orig[G1->edgeList[j].tail]	)=k;
			edge2[Gnew->edgeListPtrs[Gnew->numVertices]+k].weight=1;
			                 flag[C_orig[G1->edgeList[j].tail]]+=1;
//cout<<edge2[Gnew->edgeListPtrs[Gnew->numVertices]+k].head<<" "<<edge2[Gnew->edgeListPtrs[Gnew->numVertices]+k].tail<<" "<<edge2[Gnew->edgeListPtrs[Gnew->numVertices]+k].weight<<endl;

					k++;
					}
		else if(flag[C_orig[G1->edgeList[j].tail]]>0)
			{
				edge2[Gnew->edgeListPtrs[Gnew->numVertices]+pos.at(C_orig[G1->edgeList[j].tail])].weight=flag[C_orig[G1->edgeList[j].tail]]+1;
			}
// 	 	cout<<edge2[Gnew->edgeListPtrs[Gnew->numVertices]+k].head<<" "<<edge2[Gnew->edgeListPtrs[Gnew->numVertices]+k].tail<<" "<<edge2[Gnew->edgeListPtrs[Gnew->numVertices]+k].weight<<endl;
			
                }
	//	k++;
		
		unsigned int adj3=Gnew->edgeListPtrs[C_orig[i]];
                unsigned int adj4=Gnew->edgeListPtrs[C_orig[i]+1];
	#pragma omp parallel for
	for(unsigned int j=adj3;j<adj4;j++)
                {
                        if(Gnew->edgeList[j].head!=Gnew->edgeList[j].tail && edge2[j].weight!=0)
                              edge2[j].weight-=1;

                      edgeP[C_orig[i]+1]-=1;
		      edgeP[C_orig[Gnew->edgeList[j].tail]+1]-=1;


                }

		
                edgeP[Gnew->numVertices+jj+1]=edgeP[Gnew->numVertices+jj-1+1]+k;
                jj++;
		ii++;
        }

}
unsigned int prevV=Gnew->numVertices;

//cout<<"ok"<<endl;
//cout<<"mo"<<" "<<mo1->edgeno<<endl;
Gnew->numVertices=prevV+newV+mo1->vertex;
//cout<<"check="<<edgeP[prevV+newV];
cout<<"After transfer"<<endl;
cout<<"No of vertices="<<" "<<Gnew->numVertices<<" ";
Gnew->numEdges=edgeP[prevV+newV]+mo1->edgeno;
//cout<<"No of edges="<<" "<<Gnew->numEdges<<" "<<Gnew->edgeListPtrs[newV+prevV]+k<<endl;;

Gnew->edgeListPtrs=(unsigned int*)malloc(sizeof(unsigned int)*(Gnew->numVertices+1));
//cout<<"okk?"<<endl;
//Gnew->edgeListPtrs=(unsigned int*)malloc(sizeof(unsigned int)*(newV+Gnew->numVertices+1));
for(unsigned int i=0;i<newV+prevV+1;i++){
	Gnew->edgeListPtrs[i]=edgeP[i];
	cout<<edgeP[i]<<endl;
	}
cout<<"No of edges="<<" "<<Gnew->numEdges<<" "<<k<<endl;;

for(unsigned int i=0;i<newV+prevV+1;i++){
        cout<<Gnew->edgeListPtrs[i]<<endl;}
//free(edgeP);

//cout<<"ok?"<<endl;

//Gnew->edgeList=(edge*)malloc(sizeof(edge)*(k+Gnew->edgeListPtrs[Gnew->numVertices]));
Gnew->edgeList=(edge *)malloc(sizeof(edge)*Gnew->edgeListPtrs[Gnew->numVertices]);
for(unsigned int i=0;i<(Gnew->edgeListPtrs[prevV+newV]);i++)
	{

	Gnew->edgeList[i].head=edge2[i].head;
	Gnew->edgeList[i].tail=edge2[i].tail;
	Gnew->edgeList[i].weight=edge2[i].weight;
	}
 double prevMod = -1;
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
cout<<"ok2"<<mo1->edgeListPtrsM[0]<<endl;
for(int i=prevV+newV+1;i<=Gnew->numVertices;i++)
{	//Gnew->edgeListPtrs[i]+=mo1->edgeListPtrsM[j];
	cout<<mo1->edgeListPtrsM[j]<<" ";
	j++;

}
}
//cout<<"chance"<<endl;
//Gnew->numEdges=Gnew->edgeListPtrs[Gnew->numVertices]+mo1->edgeno;

j=0;
cout<<"dirty1="<<" "<<Gnew->edgeListPtrs[prevV+newV]<<" "<<"dirty2"<<Gnew->edgeListPtrs[Gnew->numVertices]<<endl;
for(unsigned int i=(Gnew->edgeListPtrs[prevV+newV]);i<Gnew->numEdges;i++)
        {

        Gnew->edgeList[i].head=mo1->edgesM[j].head;
        Gnew->edgeList[i].tail=mo1->edgesM[j].tail;
        Gnew->edgeList[i].weight=mo1->edgesM[j].weight;
	j++;
        }

C_orig=(unsigned int *)malloc(sizeof(unsigned int)*Gnew->numVertices);
unsigned int NV=Gnew->numVertices;
for(int i=0;i<NV+1;i++)
	cout<<Gnew->edgeListPtrs[i]<<" ";
cout<<endl;

for( int i=0;i<Gnew->edgeListPtrs[NV];i++)
{
	cout<<Gnew->edgeList[i].head<<" "<<Gnew->edgeList[i].tail<<" "<<Gnew->edgeList[i].weight<<endl;


}
graph *Gnew1=(graph *)malloc(sizeof(graph));
#pragma omp parallel for
	for (long i=0; i<NV; i++) {
  	   C_orig[i] = -1;
	}
bool *dirty1,*dirty2;
//dirty1=(bool *)malloc(sizeof(bool)*NV);
dirty2=(bool*)malloc(sizeof(bool)*NV);
for(int i=0;i<NV;i++)
{
//	dirty1[i]=false;
	dirty2[i]=false;

}
displayGraphCharacteristics(Gnew);
graph *Gnew2=(graph *)malloc(sizeof(graph));
duplicateGivenGraph(Gnew,Gnew2);	
//Gnew1=runMultiPhaseLouvainAlgorithm(Gnew, C_orig, 0, 1000, 0.0001, 0.0001, 6); 
Gnew1=cpuonly(Gnew,Gnew1,Gnew2,C_orig,dirty3,dirty2,1);
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
free(flag);
//Gnew->numVertices=Gnew->numVertices+newV-newV11;
//Gnew->numEdges=Gnew->edgeListPtrs[Gnew->numVertices];

pos.clear();
free(edgeP);
free(edge2);




}
