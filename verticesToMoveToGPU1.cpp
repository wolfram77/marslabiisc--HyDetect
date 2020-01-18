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
void verticesToMoveToGPU1(graph *G,graph *Gnew,bool *dirtycpu,move1 *mo,long *c,unsigned int mid,unsigned int node,int f1,long * C_orig)
{
unsigned int NV=G->numVertices;
/*for(int i=0;i<NV;i++)
	cout<<G->bord[i]<<" ";*/
cout<<"**************"<<endl;
unsigned int NE=G->numEdges;
//cout<<NE<<" "<<NV<<endl;
unsigned int edgec=0;
unsigned int newV=0;
/*for(int i=0;i<NV;i++)
	dirtycpu[i]=true;
for(int i=0;i<NV;i++)
{
        if(dirtycpu[i]){
        newV++;
                       }
}*/
cout<<"no of doubtful"<<" "<<newV<<endl;
int nn=newV;
bool *flaz=(bool *)malloc(sizeof(bool)*(NV+1));
unsigned int *statIndices1=(unsigned int*)malloc(sizeof(unsigned int)*(nn+G->numVertices));
for(int i=0;i<nn+G->numVertices;i++)
	statIndices1[i]=0;
for(long i=0;i<NV+1;i++)  {
        flaz[i]=false;
                          
}
nn=0;
for(int i=0;i<G->numVertices;i++)
{
	if(G->bord[i])
		nn++;		



}
unsigned int* vtxPtr=(G->edgeListPtrs);
edge *vtxInd=(G->edgeList);
cout<<"OKK"<<endl;
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

cout<<"okk1"<<endl;
std::sort(v.begin(), v.begin()+v.size());
//cout<<"k"<<endl;
//cout<<"size="<<v.size()<<endl;
unsigned int totaledgec=0;
//ofstream fout("f1.txt");
int pos1;
int *pos=(int *)malloc(sizeof(int)*nn);
#pragma omp parallel for
for(long i=0;i<=nn;i++)
        pos[i]=0;
vector<unsigned int> *ee=new vector<unsigned int>[Gnew->numVertices];
/*for(int i=0;i<nn+G->numVertices;i++)
	ee[i].resize(0);*/
vector<unsigned int> *weight=new vector<unsigned int>[Gnew->numVertices];
/*for(int i=0;i<nn+G->numVertices;i++)
        weight[i].resize(0);*/
//cout<<"what"<<endl;

int f;
int ccd=0;
pos1=0;
int qq=0;
vector<unsigned int> pos2;
//cout<<"value"<<" "<<node<<endl;
pos2.resize(node+1);
//cout<<"problem size"<<" "<<pos2.size()<<endl;
std::fill(pos2.begin(),pos2.begin()+pos2.size(),0);
unsigned int checkpos=0;
//unsigned int* edgesa1;
//unsigned int *weighta1;
cout<<"done mys"<<endl;
//int *f1=(int*)malloc(sizeof(int)*(node));
map<unsigned int,unsigned int>** cluPtrIn = (map<unsigned int,unsigned int>**) malloc ((nn)*sizeof(map<unsigned int,unsigned int>*));
assert(cluPtrIn != 0);
map<unsigned int,unsigned int>** count = (map<unsigned int,unsigned int>**) malloc ((nn)*sizeof(map<unsigned int,unsigned int>*));
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
    
        }
/*for(int i=0;i<G->numVertices;i++)
	cout<<G->edgeListPtrs[i]<<endl;

cout<<"no problem"<<endl;*/
cout<<"no prob"<<endl;
 int k=0;
#pragma omp parallel for
for(long i=0;i<Gnew->numVertices;i++)
{
	
//	if(!dirtycpu[i])
//	{

		long adj1 = Gnew->edgeListPtrs[i];
        	long adj2 = Gnew->edgeListPtrs[i+1];
	//	cout<<adj1<<" "<<adj2<<endl;
        	for(long j=adj1;j<adj2;j++)
       		 {

			ee[qq].push_back(i);
			ee[qq].push_back(Gnew->edgeList[j].tail);
			weight[qq].push_back(Gnew->edgeList[j].weight);
                        weight[qq].push_back(Gnew->edgeList[j].weight);			




		}




//	}

 //        __sync_fetch_and_add(&qq,1);

}
cout<<"Ok"<<endl;
qq=0;
//return;
#pragma omp parallel for
for(long i=0;i<NV;i++)
{
//	cout<<"......................"<<endl;
//	f=0;
  //      if(dirtycpu[i])
    //    {
//	cout<<"......................"<<endl;
	map<unsigned int, unsigned int>::iterator localIterator;

        ccd++;
	bool bordercheck=false;
        long adj1 = vtxPtr[i];
        long adj2 = vtxPtr[i+1];
//	cout<<adj1<<" "<<adj2<<" "<<endl;	
//	cout<<"......................"<<endl;
			
          if(G->bord[i] && !bordercheck){
		
	//		cout<<"endtered"<<endl;	
			 __sync_fetch_and_add(&statIndices1[G->numVertices+qq],G->bordno[i]);
                        bordercheck=true;
	//		cout<<"endtered"<<endl;
	
	for(std::vector<unsigned int> ::iterator it=G->bordvalue[i].begin();it<G->bordvalue[i].end();it++){
	//		if(f1[c[*it-mid]]==0){
//		unsigned int tail = vtxIndIn[j].tail;
//			cout<<"......."<<endl;
                        localIterator = cluPtrIn[C_orig[i]]->find(c[*it+f]);
	//	       cout<<"okk"<<endl;	
                        if( localIterator != cluPtrIn[C_orig[i]]->end() ) {

                                 __sync_fetch_and_add(&(*(cluPtrIn[C_orig[i]]))[*it] ,1);
                                int pos= (*(count[C_orig[i]]))[c[*it+f]];
                               // edge2[Gnew->edgeListPtrs[cc]+pos].weight=(*(cluPtrIn[ii]))[*it];
                               weight[C_orig[i]].at(pos)=(*(cluPtrIn[C_orig[i]]))[c[*it+f]] ;
				weight[C_orig[i]].at(pos+1)=(*(cluPtrIn[C_orig[i]]))[c[*it+f]]; 
                                }

                        else{

                                (*(cluPtrIn[C_orig[i]]))[c[*it+f]] = 1;
                               // edge2[Gnew->edgeListPtrs[cc]+k].head=ii+(cc);
                               // edge2[Gnew->edgeListPtrs[cc]+k].tail=C_orig[G1->edgeList[j].tail];

                                   //     edge2[Gnew->edgeListPtrs[cc]+k].weight=(*(cluPtrIn[ii]))[C_orig[G1->edgeList[j].tail]];
                                 //  	cout<<"1"<<endl;
                                // 	cout<<"value of numvertex"<<qq+G->numVertices<<endl;
					 ee[C_orig[i]].push_back(qq);
                        		 ee[C_orig[i]].push_back(*it);
                        		 weight[C_orig[i]].push_back(1);
                        		 weight[C_orig[i]].push_back(1);

                                        (*(count[C_orig[i]]))[c[*it+f]]=k;
                                __sync_fetch_and_add(&k,1);
                              //  __sync_fetch_and_add(&k1,1);
       //                       cout<<"value of k"<<k<<endl;
                                }
		}													}

					

/*else{
     cout<<"..........."<<endl;                                    
    
for(unsigned int j=adj1;j<adj2;j++){
		if(dirtycpu[(vtxInd[j].tail)]){
                //statIndices1[qq]=statIndices1[qq]+1;
		 __sync_fetch_and_add(&statIndices1[qq+G->numVertices],1);
                std::vector<unsigned int>::iterator itr;
                itr=std::find (v.begin(), v.begin()+v.size(), (vtxInd[j].tail));
		 localIterator = cluPtrIn[qq+G->numVertices]->find((itr - v.begin()) +f);
                 if( localIterator != cluPtrIn[qq+G->numVertices]->end() ) {
				 __sync_fetch_and_add(&(*(cluPtrIn[qq+G->numVertices]))[(itr - v.begin()) +f] ,1);
				int pos= (*(count[qq+G->numVertices]))[itr - v.begin() +f];
				weight[qq].at(pos)=(*(cluPtrIn[qq+G->numVertices]))[itr - v.begin() +f] ;
                                weight[qq].at(pos+1)=(*(cluPtrIn[qq+G->numVertices]))[itr - v.begin() +f];
                                }
	else{
		(*(cluPtrIn[qq+G->numVertices]))[(itr - v.begin()) +f] = 1;

		ee[qq+G->numVertices].push_back(qq+G->numVertices);

                ee[qq+G->numVertices].push_back((itr - v.begin()) +f);
		weight[qq+G->numVertices].push_back(1);
		weight[qq+G->numVertices].push_back(1);
              	(*(count[qq+G->numVertices]))[(itr - v.begin()) +f]=k;
                       __sync_fetch_and_add(&k,1);


         //       __sync_fetch_and_add(&pos1,1);
                 __sync_fetch_and_add(&f,1);

                                 }
        
	}
}*/

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
/*	if(qq+G->numVertices<10+G->numVertices)
	cout<<ee[qq+G->numVertices].size()<<" "<<weight[qq+G->numVertices].size()<<endl;*/
         __sync_fetch_and_add(&qq,1);
	k=0;
        


}
cout<<"1 done"<<" "<<qq<<endl;
unsigned int k1=0,k2=0;
 int size=0; int size1=0;
/*for(int i=0;i<10;i++)
{
cout<<ee[i].size()<<" "<<weight[i].size()<<endl;
}*/
for(int j=0;j<Gnew->numVertices;j++){
	//size+=ee[qq].size();
	__sync_fetch_and_add(&size,ee[j].size());

//	cout<<ee[qq].size()<<" ";
	}
for(int i=0;i<Gnew->numVertices;i++)
//        size1+=weight[qq].size();
__sync_fetch_and_add(&size1,weight[i].size());

cout<<"2 done"<<endl;
cout<<size<<" "<<size1<<endl;

unsigned int *edgesa1=(unsigned int*)malloc(sizeof(unsigned int)*(size1));
unsigned int *weighta1=(unsigned int*)malloc(sizeof(unsigned int)*(size1));
#pragma omp paraller for
for(int i=0;i<+Gnew->numVertices;i++){
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

mo->vertex=Gnew->numVertices;

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

cout<<mo->vertex<<" "<<mo->edg<<endl;
mo->statIndices=(unsigned int *)malloc(sizeof(unsigned int)*mo->vertex);
mo->edgesa=(unsigned int* )malloc(sizeof(unsigned int)*mo->edg);
mo->weighta=(unsigned int*)malloc(sizeof(unsigned int)*mo->edg);
for(unsigned int i=0;i<mo->edg;i++)
{
mo->edgesa[i]=edgesa1[i];
  //      cout<<"1"<<endl;
        mo->weighta[i]=weighta1[i];
//cout<<"2"<<endl;
//cout<<"........"<<" "<<mo->edgesa[i]<<" "<<mo->weighta[i]<<endl;

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
