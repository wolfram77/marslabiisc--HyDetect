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







//we need to move subset of doubtful vertices to the other device. This Code find the subset of the vertices that needs to be moved and store those vertices int mo data structure//
void verticesToMoveToGPU(graph *G,bool *dirtycpu,move1 *mo,unsigned int *c,unsigned int mid,unsigned int node)
{
unsigned int NV=G->numVertices;
unsigned int NE=G->numEdges;
//unsigned int *vtxPtr=G->edgeListPtrs;
//edge *vtxInd=G->edgeList;
cout<<NE<<" "<<NV<<endl;
unsigned int edgec=0;
unsigned int newV=0;
for(int i=0;i<NV;i++)
{
        if(dirtycpu[i]){
        newV++;
                       }
}
cout<<"no of doubtful"<<" "<<newV<<endl;
int nn=newV;
bool *flaz=(bool *)malloc(sizeof(bool)*(NV+1));
unsigned int *statIndices1=(unsigned int*)malloc(sizeof(unsigned int)*nn);
for(int i=0;i<nn;i++)
	statIndices1[i]=0;
for(long i=0;i<NV+1;i++)  {
        flaz[i]=false;
                          
}

unsigned int* vtxPtr=(G->edgeListPtrs);
edge *vtxInd=(G->edgeList);
//cout<<"OKK"<<endl;
vector<unsigned int> v;
#pragma omp parallel for
for(long i=0;i<NV;i++)
  {
        if(dirtycpu[i])
        {
        bool bordercheck=false;
        long adj1 = vtxPtr[i];
        long adj2 = vtxPtr[i+1];
    
        for(long j=adj1;j<adj2;j++){
		//if(G->bord[i] && !bordercheck ){
                   

                if(dirtycpu[(vtxInd[j].tail+1)] && !flaz[(vtxInd[j].tail+1)]){
                         v.push_back(vtxInd[j].tail+1);
                         flaz[(vtxInd[j].tail+1)]=true;
                                                                       	}
                                
					}
        }
    }
//cout<<"okk1"<<endl;
std::sort(v.begin(), v.begin()+v.size());
//cout<<"k"<<endl;
cout<<"size="<<v.size()<<endl;
unsigned int totaledgec=0;
//ofstream fout("f1.txt");
int pos1;
int *pos=(int *)malloc(sizeof(int)*nn);
#pragma omp parallel for
for(long i=0;i<=nn;i++)
        pos[i]=0;
vector<unsigned int> *ee=new vector<unsigned int>[nn];
vector<unsigned int> *weight=new vector<unsigned int>[nn];
//cout<<"what"<<endl;

int f;
int ccd=0;
pos1=0;
int qq=0;
vector<unsigned int> pos2;
//cout<<"value"<<" "<<node<<endl;
pos2.resize(node+1);
cout<<"problem size"<<" "<<pos2.size()<<endl;
std::fill(pos2.begin(),pos2.begin()+pos2.size(),0);
unsigned int checkpos=0;
//unsigned int* edgesa1;
//unsigned int *weighta1;
//cout<<"done mys"<<endl;
int *f1=(int*)malloc(sizeof(int)*(node)); 
#pragma omp parallel for
for(long i=0;i<NV;i++)
{
	f=0;
        if(dirtycpu[i])
        {
        ccd++;
//	int *f1=(int*)malloc(sizeof(int)*(node)); //should be total-NV
  //      for(int i1=0;i1<node;i1++)
//	f1[i1]=0;
	bool bordercheck=false;
        long adj1 = vtxPtr[i];
        long adj2 = vtxPtr[i+1];
	
	for(long j=adj1;j<adj2;j++)
        {
                if(G->bord[i] && !bordercheck ){
			//cout<<G->bordno[i]<<"border"<<endl;	
                        statIndices1[qq]=statIndices1[qq]+G->bordno[i];
			//fout<<mo[qq].statIndices<<endl;
                        bordercheck=true;

	
		for(std::vector<unsigned int> ::iterator it=G->bordvalue[i].begin();it<G->bordvalue[i].end();it++){
	//		if(f1[c[*it-mid]]==0){
			ee[qq].push_back(qq);
                        ee[qq].push_back(*it);
			weight[qq].push_back(1);
			weight[qq].push_back(1);
			int b=*it-mid;
		//	cout<<"b="<<c[b]<<endl;
			pos2[c[b]+1]=checkpos;
			checkpos++;
		//	cout<<"done1"<<endl;
			//	f1[c[*it-mid]]+=1;}
	/*		else if(f1[c[*it-mid]]>0)
			{	//unsigned int a=*it-mid;
				unsigned int b=(pos2.at(c[*it-mid]));
				//weighta[qq].push_back(f1[c[*it-mid]]+1);
				weight[qq].at(b+1)+=1;
				weight[qq].at(b)+=1;
		//		cout<<"done2"<<endl;
				f1[c[*it-mid]]+=1;
			}*/
		//	fout<<c[*it]<<endl;
                             pos[qq]++;
                                pos1++;}


                                         }
      if(dirtycpu[(vtxInd[j].tail)+1]){

                statIndices1[qq]=statIndices1[qq]+1;
                std::vector<unsigned int>::iterator itr;
                itr=std::find (v.begin(), v.begin()+v.size(), (vtxInd[j].tail+1));
		ee[qq].push_back(qq);

                ee[qq].push_back((itr - v.begin()) +f);
//		cout<<"done3"<<endl;
		weight[qq].push_back(1);
		weight[qq].push_back(1);
                pos[qq]++;

                pos1++;
                f++;

                                 }
        
	}	
//cout<<"1"<<endl;
//cout<<"size"<<" "<<ee[qq].size()<<endl;
/*edgesa1=(unsigned int*)malloc(sizeof(unsigned int)*(ee[qq].size()));

weighta1=(unsigned int*)malloc(sizeof(unsigned int)*(ee[qq].size()));
//cout<<ee[qq].size()<<endl;
if(ee[qq].size()!=0){
        for(std::vector<unsigned int>::iterator it=ee[qq].begin();it!=ee[qq].end();it++)
        {       int ij=0;
                edgesa1[ij]=*it;
                ij++;
		totaledgec++;
	}
	 for(std::vector<unsigned int>::iterator it=weight[qq].begin();it!=weight[qq].end();it++)
        {       int ij=0;
                weighta1[ij]=*it;
                ij++;
        }

                }*/
        qq++;
        }


}
cout<<"1 done"<<" "<<qq<<endl;
unsigned int k1=0,k2=0;
unsigned int size=0;
for(int i=0;i<qq;i++)
	size+=ee[qq].size();
cout<<"2 done"<<endl;
unsigned int *edgesa1=(unsigned int*)malloc(sizeof(unsigned int)*(size));
unsigned int *weighta1=(unsigned int*)malloc(sizeof(unsigned int)*(2*size));
#pragma omp paraller for
for(int i=0;i<qq;i++){
if(ee[i].size()!=0){
        for(std::vector<unsigned int>::iterator it=ee[i].begin();it!=ee[i].end();it++)
        {       
                edgesa1[k1]=*it;
                k1++;
                totaledgec++;
        }
         for(std::vector<unsigned int>::iterator it=weight[i].begin();it!=weight[i].end();it++)
        {  //     int ij=0;
                weighta1[k2]=*it;
		k2++;
                
        }

                }

}
cout<<"happening"<<" "<<totaledgec<<endl;
mo->vertex=qq;
mo->edg=totaledgec;
mo->statIndices=statIndices1;
mo->edgesa=edgesa1;
mo->weighta=weighta1;
//cout<<"stat"<<mo->edgesa[5]<<" "<<statIndices1[1]<<endl;
//cout<<"1"<<endl;
//return mo;
//free(vtxInd);
//free(vtxPtr);
free(flaz);
free(pos);
free(edgesa1);
free(weighta1);
delete[]ee;
delete[]weight;
pos2.clear();
v.clear();
}
