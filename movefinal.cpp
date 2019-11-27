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




//void verticesToMoveToGPU(dirty3,graph *G1,graph Gnew,bool *dirtycpu,move1 *mo,unsigned int *c,unsigned int mid,unsigned int node)
//


//we need to move subset of doubtful vertices to the other device. This Code find the subset of the vertices that needs to be moved and store those vertices int mo data structure//
void movefinal(graph *Gnew,graph *Gnew1,bool *dirty1,bool *dirtycpu,move1 *mo,unsigned int *c,unsigned int mid,unsigned int node,bool *bord,unsigned int* bordno,vector<unsigned int>* bordvalue)
{
/*unsigned int NV=G1->numVertices;
unsigned int  NE=G1->numEdges;
unsigned int * vtxPtr=(G1->edgeListPtrs);
cout<<vtxPtr<<endl;
edge *vtxInd=(G1->edgeList);
cout<<vtxInd<<endl;
bool *f1=(bool*)malloc(sizeof(bool)*NV);
unsigned int *count1=(unsigned int *)malloc(sizeof(unsigned int )*(NV+1));
vector<long> *borderval=new vector<long>[NV+1];

for(long i=0;i<NV;i++){
        f1[i]=false;
        count1[i]=0;}
for(long i=0;i<NV;i++)
        {
                long adj1=vtxPtr[i];
                long adj2=vtxPtr[i+1];
        for(long j=adj1;j<adj2;j++)
                {
                        if(dirty3[(vtxInd[j].tail+1)]){
                                f1[i]=true;
                                count1[i]=count1[i]+1;
                                borderval[i].push_back((vtxInd[j].tail+1)-1);   }

                }
        }*/
//dirtycpu=(bool *)malloc(sizeof(bool)*
unsigned int NV,NE;
 NV=Gnew->numVertices;
 NE=Gnew->numEdges;
//unsigned int *vtxPtr=G->edgeListPtrs;
//edge *vtxInd=G->edgeList;
cout<<NE<<" "<<NV<<endl;
unsigned int edgec=0;
unsigned int newV=0;
for(int i=0;i<NV;i++)
{
        if(dirty1[i]){
        newV++;
                        }
}
cout<<"no of doubtful vertices="<<newV<<endl;
dirtycpu=(bool *)malloc(sizeof(bool)*(NV+newV));
for(int i=0;i<NV+newV;i++)
	dirtycpu[i]=true;
//cout<<"entry"<<endl;
int nn=newV+NV;
bool *flaz=(bool *)malloc(sizeof(bool)*(NV+newV));
unsigned int *statIndices1=(unsigned int*)malloc(sizeof(unsigned int)*nn);
for(int i=0;i<nn;i++)
	statIndices1[i]=0;
for(long i=0;i<NV+newV;i++)  {
        flaz[i]=false;
                          
}


unsigned int* vtxPtr=(Gnew->edgeListPtrs);
edge *vtxInd=(Gnew->edgeList);
//cout<<"OKK"<<endl;
vector<unsigned int> v;
cout<<"************************************************"<<endl;
#pragma omp parallel for
for(long i=0;i<Gnew->numVertices;i++)
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
cout<<"||************************************************||"<<endl;

//cout<<"okk1"<<endl;
std::sort(v.begin(), v.begin()+v.size());
//cout<<"k"<<endl;
//cout<<"size="<<v.size()<<endl;
for(int i=0;i<v.size();i++)
	cout<<v[i]<<" ";
cout<<endl;
unsigned int totaledgec=0;
//ofstream fout("f1.txt");
int pos1;
int *pos=(int *)malloc(sizeof(int)*nn);
#pragma omp parallel for
for(long i=0;i<nn;i++)
        pos[i]=0;
vector<unsigned int> *ee=new vector<unsigned int>[nn];
vector<unsigned int> *weight=new vector<unsigned int>[nn];
//cout<<"what"<<endl;

int f;
int ccd=0;
pos1=0;
int qq=0;
vector<unsigned int> pos2;
pos2.resize(max(Gnew->numVertices+newV,node));
//ut<<"value"<<" "<<node<<endl;
//s2.resize(node+Gnew->numVertices);
//std::fill(pos2.begin(),pos2.begin()+pos2.size(),0);
unsigned int checkpos=0;
//unsigned int* edgesa1;
//unsigned int *weighta1;
//cout<<"done mys"<<endl;
#pragma omp parallel for
for(long i=0;i<Gnew->numVertices;i++)
{f=0;
        if(dirtycpu[i])
        {
        ccd++;
	int *f1=(int*)malloc(sizeof(int)*(node+newV)); //should be total-NV
        for(int i1=0;i1<node+newV;i1++)
	f1[i1]=0;
	bool bordercheck=false;
        long adj1 = vtxPtr[i];
        long adj2 = vtxPtr[i+1];

        for(long j=adj1;j<adj2;j++)
        {
                if(bord[i] && !bordercheck ){
			//cout<<G->bordno[i]<<"border"<<endl;	
                        statIndices1[qq]=statIndices1[qq]+bordno[i];
			//fout<<mo[qq].statIndices<<endl;
                        bordercheck=true;
                for(std::vector<unsigned int> ::iterator it=bordvalue[i].begin();it!=bordvalue[i].end();it++){
		//	cout<<"edge val"<<" "<<*it<<endl;
			if(f1[c[*it]]==0){
			ee[qq].push_back(qq);
                        ee[qq].push_back(c[*it]);
			weight[qq].push_back(1);
			weight[qq].push_back(1);
			int b=*it;
		//	cout<<"val="<<b<<endl;
		//	cout<<"b="<<c[b]<<endl;
			//cout<<c[b]<<endl;
		//	cout<<"first"<<" "<<c[*it]<<endl;
			int a=c[b];
		//	cout<<"value of a"<<" "<<c[a]<<endl;
		//	pos2.at(a)=checkpos;
			
			checkpos++;
		//	cout<<"done"<<endl;
				f1[c[*it]]+=1;
			//	cout<<"ok1"<<endl;
					}
			else if(f1[c[*it]]>0)
			{	//unsigned int a=*it-mid;
			//	int b=0;
				unsigned int b=(pos2.at(c[*it]));
				//weighta[qq].push_back(f1[c[*it-mid]]+1);
				weight[qq].at(b+1)+=1;
				weight[qq].at(b)+=1;
				f1[c[*it]]+=1;
			}
		//	fout<<c[*it]<<endl;
                               pos[qq]++;
			//	cout<<"ok2"<<endl;
                                pos1++;}


                                         }
       if(dirtycpu[(vtxInd[j].tail)+1]){

                statIndices1[qq]=statIndices1[qq]+1;
                std::vector<unsigned int>::iterator itr;
                itr=std::find (v.begin(), v.begin()+v.size(), (vtxInd[j].tail+1));
		ee[qq].push_back(qq);

                ee[qq].push_back((itr - v.begin()) +f);
		weight[qq].push_back(1);
		weight[qq].push_back(1);
                pos[qq]++;

                pos1++;
                f++;

                                 }
        }
	
cout<<"1"<<endl;
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
unsigned int k1=0,k2=0;
unsigned int size=0;
for(int i=0;i<qq;i++)
	size+=ee[qq].size();
unsigned int *edgesa1=(unsigned int*)malloc(sizeof(unsigned int)*(size));
unsigned int *weighta1=(unsigned int*)malloc(sizeof(unsigned int)*(size));
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
//cout<<"happening"<<endl;
mo->vertex=qq;
mo->edg=totaledgec;
mo->statIndices=statIndices1;
mo->edgesa=edgesa1;
mo->weighta=weighta1;
cout<<"stat"<<mo->statIndices[0]<<" "<<statIndices1[0]<<endl;
//cout<<"1"<<endl;
//return mo;
//free(vtxInd);
//free(vtxPtr);
free(flaz);
free(pos);
delete[]ee;
delete[]weight;
pos2.clear();
v.clear();
}
