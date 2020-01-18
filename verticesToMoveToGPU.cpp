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
void verticesToMoveToGPU(graph *G,bool *dirtycpu,move1 *mo,long *c,unsigned int mid,unsigned int node,long *C_orig)
{
unsigned int NV=G->numVertices;
unsigned int NE=G->numEdges;
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
		
                   

                if(dirtycpu[(vtxInd[j].tail)] && !flaz[(vtxInd[j].tail)]){
                         v.push_back(vtxInd[j].tail);
                         flaz[(vtxInd[j].tail)]=true;
                                                                       	}
                                
					}
        }
    }
cout<<"....."<<endl;
std::sort(v.begin(), v.begin()+v.size());
unsigned int totaledgec=0;
int pos1=0;
int *pos=(int *)malloc(sizeof(int)*nn+1);
#pragma omp parallel for
for(long i=0;i<=nn;i++)
        pos[i]=0;
vector<unsigned int> *ee=new vector<unsigned int>[nn];
/*for(int i=0;i<nn;i++)
	ee[i].resize(node);*/
vector<unsigned int> *weight=new vector<unsigned int>[nn];
/*for(int i=0;i<nn;i++)
        weight[i].resize(node);
*/
cout<<"no no"<<endl;
int f;
int ccd=0;
pos1=0;
int qq=0;
vector<unsigned int> pos2;
pos2.resize(node+1);
std::fill(pos2.begin(),pos2.begin()+pos2.size(),0);
unsigned int checkpos=0;
map<unsigned int,unsigned int>** cluPtrIn = (map<unsigned int,unsigned int>**) malloc (nn*sizeof(map<unsigned int,unsigned int>*));
assert(cluPtrIn != 0);
map<unsigned int,unsigned int>** count = (map<unsigned int,unsigned int>**) malloc (nn*sizeof(map<unsigned int,unsigned int>*));
assert(count!=0);
#pragma omp parallel for
  for (long i=0; i<nn; i++) {
        cluPtrIn[i] = new map<unsigned int,unsigned int>();
        count[i]=new map<unsigned int,unsigned int>();
 /*       if(i<cc){
        unsigned int adj1=vtxPtrOut[i];
        unsigned int adj2=vtxPtrOut[i+1];
        for(unsigned int j=adj1;j<adj2;j++)
                {
                        unsigned int tail=vtxIndOut[j].tail;  
                        (*(cluPtrIn[i]))[tail] =(unsigned int)vtxIndOut[j].weight;;
    
                }
        }*/
(*(cluPtrIn[i]))[i]=0; 
 (*(count[i]))[i]=0;
        }


cout<<"no problem"<<endl;
 int k=0;
#pragma omp parallel for
for(long i=0;i<NV;i++)
{
	f=0;k=0;
        if(dirtycpu[i])
        {
	map<unsigned int, unsigned int>::iterator localIterator;

        ccd++;
//	int *f1=(int*)malloc(sizeof(int)*(node)); //should be total-NV
//      for(int i1=0;i1<node;i1++)
//	f1[i1]=0;
	bool bordercheck=false;
        long adj1 = vtxPtr[i];
        long adj2 = vtxPtr[i+1];
	
//	for(long j=adj1;j<adj2;j++)
  //      {
			
                if(G->bord[i] && !bordercheck ){
			//cout<<G->bordno[i]<<"border"<<endl;	
                       // statIndices1[qq]=statIndices1[qq]+G->bordno[i];
			 __sync_fetch_and_add(&statIndices1[qq],G->bordno[i]);
			//fout<<mo[qq].statIndices<<endl;
                        bordercheck=true;
			
	
	for(std::vector<unsigned int> ::iterator it=G->bordvalue[i].begin();it<G->bordvalue[i].end();it++){
	
                        localIterator = cluPtrIn[qq]->find(c[*it-mid]);
		     //  cout<<"okk"<<endl;	
                       if( localIterator != cluPtrIn[qq]->end() ) {
		//		cout<<"1"<<endl;
                                 __sync_fetch_and_add(&(*(cluPtrIn[qq]))[c[*it-mid]] ,1);
		//		cout<<"2"<<endl;
                                int pos= (*(count[qq]))[c[*it-mid]];
                               // edge2[Gnew->edgeListPtrs[cc]+pos].weight=(*(cluPtrIn[ii]))[*it];
                              // cout<<"pos"<<" "<<pos<<" "<<weight[qq].size()<<endl; 
				if((pos+1)<weight[qq].size()){
                              weight[qq].at(pos)=(*(cluPtrIn[qq]))[c[*it-mid]] ;
		//		cout<<"3"<<endl;
				weight[qq].at(pos+1)=(*(cluPtrIn[qq]))[c[*it-mid]]; }
			//	cout<<".."<<endl;
                                }

                       else{

                                (*(cluPtrIn[qq]))[c[*it-mid]] = 1;
                               // edge2[Gnew->edgeListPtrs[cc]+k].head=ii+(cc);
                               // edge2[Gnew->edgeListPtrs[cc]+k].tail=C_orig[G1->edgeList[j].tail];

                                   //     edge2[Gnew->edgeListPtrs[cc]+k].weight=(*(cluPtrIn[ii]))[C_orig[G1->edgeList[j].tail]];
                        //           	cout<<"1"<<endl;
					 ee[qq].push_back(qq);
                        		 ee[qq].push_back(*it);
                        		 weight[qq].push_back(1);
                        		 weight[qq].push_back(1);

                                        (*(count[qq]))[c[*it-mid]]=k;
                                __sync_fetch_and_add(&k,1);
                              //  __sync_fetch_and_add(&k1,1);
                          //    cout<<"value of k"<<k<<endl;
                                }												
	}
//cout<<"::::::::::"<<endl;
			G->bordvalue[i].clear(); G->bordno[i]=0; G->bord[i]=false;
/*			ee[qq].push_back(qq);
                        ee[qq].push_back(*it);
			weight[qq].push_back(1);
			weight[qq].push_back(1);
			int b=*it-mid;
		//	cout<<"b="<<c[b]<<endl;
	//		pos2[c[b]+1]=checkpos;
		//	checkpos++;
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
               //              pos[qq]++;
                               
			//	 __sync_fetch_and_add(&pos1,1);
					

//}
        
/*      for(long j=adj1;j<adj2;j++)
{ 
      if(dirtycpu[(vtxInd[j].tail)]){

                //statIndices1[qq]=statIndices1[qq]+1;
		 __sync_fetch_and_add(&statIndices1[qq],1);
                std::vector<unsigned int>::iterator itr;
                itr=std::find (v.begin(), v.begin()+v.size(), (vtxInd[j].tail));
		 localIterator = cluPtrIn[qq]->find((itr - v.begin()) +f);
                 if( localIterator != cluPtrIn[qq]->end() ) {
				 __sync_fetch_and_add(&(*(cluPtrIn[qq]))[(itr - v.begin()) +f] ,1);
				int pos= (*(count[qq]))[itr - v.begin() +f];
				weight[qq].at(pos)=(*(cluPtrIn[qq]))[itr - v.begin() +f] ;
                                weight[qq].at(pos+1)=(*(cluPtrIn[qq]))[itr - v.begin() +f];
                                }
	else{
		(*(cluPtrIn[qq]))[(itr - v.begin()) +f] = 1;

		ee[qq].push_back(qq);

                ee[qq].push_back((itr - v.begin()) +f);
		weight[qq].push_back(1);
		weight[qq].push_back(1);
              	(*(count[qq]))[(itr - v.begin()) +f]=k;
                       __sync_fetch_and_add(&k,1);


         //       __sync_fetch_and_add(&pos1,1);
                 __sync_fetch_and_add(&f,1);

                                 }
        
	}}*/
/*	for(long j=adj1;j<adj2;j++)
	{
		G->bord[C_orig[(vtxInd[j].tail)]]=true;
		G->bordno[C_orig[(vtxInd[j].tail)]]+=1;
		G->bordvalue[C_orig[(vtxInd[j].tail)]].push_back(qq+node);
		

	}*/
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
//	if(qq<10)
//	cout<<ee[qq].size()<<" "<<weight[qq].size()<<endl;
         __sync_fetch_and_add(&qq,1);
	k=0;
        }


}
int q=0;
cout<<"1 done"<<" "<<qq<<endl;

for(long i=0;i<NV;i++)
{
        f=0;
       if(dirtycpu[i])
        {
	long adj1 = vtxPtr[i];
        long adj2 = vtxPtr[i+1];

       		for(long j=adj1;j<adj2;j++)
        {
	 if(!dirtycpu[(vtxInd[j].tail)]){
//		cout<<C_orig[(vtxInd[j].tail)]<<endl;
                G->bord[C_orig[(vtxInd[j].tail)]]=true;
                G->bordno[C_orig[(vtxInd[j].tail)]]+=1;
                G->bordvalue[C_orig[(vtxInd[j].tail)]].push_back(q);
			}
                

        }

	}
	q++;
}

cout<<"?"<<endl;
unsigned int k1=0,k2=0;
 int size=0; int size1=0;
/*for(int i=0;i<10;i++)
{
cout<<ee[i].size()<<" "<<weight[i].size()<<endl;
}*/
for(int j=0;j<qq;j++){
	//size+=ee[qq].size();
	__sync_fetch_and_add(&size,ee[j].size());

//	cout<<ee[qq].size()<<" ";
	}
for(int i=0;i<qq;i++)
//        size1+=weight[qq].size();
__sync_fetch_and_add(&size1,weight[i].size());

//cout<<"2 done"<<endl;
//cout<<size<<" "<<size1<<endl;
unsigned int *edgesa1=(unsigned int*)malloc(sizeof(unsigned int)*(size1));
unsigned int *weighta1=(unsigned int*)malloc(sizeof(unsigned int)*(size1));
#pragma omp paraller for
for(int i=0;i<qq;i++){
if(ee[i].size()>0 &&weight[qq].size()>0){
        for(std::vector<unsigned int>::iterator it=ee[i].begin();it!=ee[i].end();it++)
        {       
                edgesa1[k1]=*it;
                k1++;
                totaledgec++;
        }
         for(std::vector<unsigned int>::iterator it=weight[i].begin();it!=weight[i].end();it++)
        {  //     int ij=0;
//	cout<<"val"<< " "<<k2<<" "<<size1<<endl;
                weighta1[k2]=*it;
		k2++;
                
        }

                }

}

cout<<"happening"<<" "<<totaledgec<<" "<<size1<<endl;

mo->vertex=qq;

mo->edg=totaledgec;
//cout<<mo->edg<<endl;
/*for(unsigned int i=0;i<mo->edg;i++)
{
//	mo->statIndices[i]=statIndices1[i];
	
	mo->edgesa[i]=edgesa1[i];
	cout<<"1"<<endl;
	mo->weighta[i]=weighta1[i];
cout<<"2"<<endl;

}*/
mo->statIndices=(unsigned int *)malloc(sizeof(unsigned int)*mo->vertex);
mo->edgesa=(unsigned int* )malloc(sizeof(unsigned int)*mo->edg);
mo->weighta=(unsigned int*)malloc(sizeof(unsigned int)*mo->edg);
for(unsigned int i=0;i<mo->edg;i++)
{
mo->edgesa[i]=edgesa1[i];
  //      cout<<"1"<<endl;
        mo->weighta[i]=weighta1[i];
//cout<<"2"<<endl;
//cout<<mo->edgesa[i]<<" "<<mo->weighta[i]<<endl;

}
for(int i=0;i<mo->vertex;i++)
mo->statIndices[i]=statIndices1[i];
//mo->edgesa=edgesa1;
//mo->weighta=weighta1;
//cout<<"stat"<<mo->edgesa[1]<<" "<<statIndices1[1]<<endl;
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
