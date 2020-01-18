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


using namespace std;







graph* modifyCPUstructure(graph *Gnew,graph *G1,bool *dirty,long *C_orig,long *C_orig1,move2 *mo1,unsigned int mid,unsigned int node,bool *dirty3,bool *dirtyc){
cout<<"inside modifyCPU"<<endl;
cout<<Gnew->edgeListPtrs[Gnew->numVertices]<<endl;
displayGraphCharacteristics(Gnew);
/*for(int i=0;i<=Gnew->numVertices;i++)
	cout<<Gnew->edgeListPtrs[i]<<" ";
for(int i=0;i<Gnew->edgeListPtrs[Gnew->numVertices];i++)
         cout<<Gnew->edgeList[i].head<<" "<<Gnew->edgeList[i].tail<<" "<<Gnew->edgeList[i].weight<<endl;
*/

unsigned int newV=0;
unsigned int edgec1=0;
for(int i=0;i<G1->numVertices;i++)
{
        if(dirty[i]){
        newV++;
     
		}

}
cout<<"no of doubtful vertices= "<<newV<<endl;
unsigned int    *vtxPtrIn    = G1->edgeListPtrs;
edge    *vtxIndIn    = G1->edgeList;
unsigned int    *vtxPtrOut    = Gnew->edgeListPtrs;
edge    *vtxIndOut    = Gnew->edgeList;
unsigned int* doubtindex=(unsigned int*)malloc(sizeof(unsigned int)*G1->numVertices);
for(int i=0;i<G1->numVertices;i++)
	doubtindex[i]=0;	
int k=0;
bool checkflag=true;
int reduce=0;
int count1=0;int count2=0;
map<unsigned int,unsigned int>** cluPtrIn1 = (map<unsigned int,unsigned int>**) malloc (Gnew->numVertices*sizeof(map<unsigned int,unsigned int>*));
  assert(cluPtrIn1 != 0);
//map<unsigned int,unsigned int>** count1 = (map<unsigned int,unsigned int>**) malloc (Gnew->numVertices*sizeof(map<unsigned int,unsigned int>*));
#pragma omp parallel for
  for (long i=0; i<Gnew->numVertices; i++) {
        cluPtrIn1[i] = new map<unsigned int,unsigned int>();
 //       count1[i]=new map<unsigned int,unsigned int>();
       
        unsigned int adj1=vtxPtrOut[i];
        unsigned int adj2=vtxPtrOut[i+1];
        for(unsigned int j=adj1;j<adj2;j++)
                {
                        unsigned int tail=vtxIndOut[j].tail;
                        (*(cluPtrIn1[i]))[tail] =(unsigned int)vtxIndOut[j].weight;;

                }
        

        }

for(int i=0;i<G1->numVertices;i++)
{

	if(dirty[i])
{
	{
		 /*  if(     (*(cluPtrIn1[C_orig[i]]))[C_orig[i]]>=1)
                                 __sync_fetch_and_sub(&(*(cluPtrIn1[C_orig[i]]))[C_orig[i]] ,1);
		*/
		unsigned int aa=vtxPtrIn[i];
		unsigned int bb=vtxPtrIn[i+1];
		for(unsigned int j=aa;j<bb;j++)
		{
	
               /* for(unsigned int z=adj1;z<adj2;z++)
                {

              //         unsigned int tail = vtxIndOut[z].tail;

			count1++;
                        if(Gnew->edgeList[C_orig[j]].weight>=1){
                        Gnew->edgeList[C_orig[j]].weight-=1;

                        }
                else if(Gnew->edgeList[C_orig[j]].weight==0){
//	cout<<Gnew->edgeList[z].head<<" "<<Gnew->edgeList[z].tail<<" "<<Gnew->edgeList[z].weight<<endl;
                	              count2++;
					 }

                }*/

	//	 localIterator1 = cluPtrIn[C_orig[i]]->find(C_orig[G1->edgeList[j].tail]);
                        
				if(	(*(cluPtrIn1[C_orig[i]]))[C_orig[G1->edgeList[j].tail]]>=1)
                                 __sync_fetch_and_sub(&(*(cluPtrIn1[C_orig[i]]))[C_orig[G1->edgeList[j].tail]] ,1);
				
				else 
					count2++;
                               // int pos= (*(count[ii]))[C_orig[G1->edgeList[j].tail]];
                               // edge2[Gnew->edgeListPtrs[cc]+pos].weight=(*(cluPtrIn[ii]))[C_orig[G1->edgeList[j].tail]];
	//		 if(     (*(cluPtrIn1[C_orig[i]]))[C_orig[i]]>=1)
          //                       __sync_fetch_and_sub(&(*(cluPtrIn1[C_orig[i]]))[C_orig[G1->edgeList[j].tail]] ,1);
                                



		}





}
}


}
cout<<Gnew->edgeListPtrs[Gnew->numVertices]<<endl;
//cout<<"val of count2"<<count_edge<<endl;
unsigned int count_edge=Gnew->edgeListPtrs[Gnew->numVertices]-count2;
//unsigned int count_edge=Gnew->edgeListPtrs[Gnew->numVertices];
//cout<<"val of count2 "<<count2<<" "<<count1<<endl;;

cout<<"val of count2"<<count_edge<<endl;
/*for(int i=0;i<Gnew->numEdges;i++)
	 cout<<Gnew->edgeList[i].head<<" "<<Gnew->edgeList[i].tail<<" "<<Gnew->edgeList[i].weight<<endl;
*/
cout<<endl;
//return Gnew;
edge *edgefinal=(edge *)malloc(sizeof(edge)*count_edge);
for(int i=0;i<count_edge;i++)
{
	edgefinal[i].head=0;
	edgefinal[i].tail=0;
	edgefinal[i].weight=0;


}
//unsigned int *vec=(unsigned int *)malloc(sizeof(unsigned int)*count_edge*2);
//
for(int i=0;i<Gnew->numVertices;i++)
{

        unsigned int adj1=vtxPtrOut[i];
        unsigned int adj2=vtxPtrOut[i+1];
        for(unsigned int j=adj1;j<adj2;j++)
	{
	 unsigned int tail=vtxIndOut[j].tail;
                        vtxIndOut[j].weight= (*(cluPtrIn1[i]))[tail];

	}

}
vector<unsigned int> vec;
int id=0;
int qq=0;
for(int i=0;i<Gnew->numVertices;i++)
{

checkflag=true;
        unsigned int adj1=vtxPtrOut[i];
        unsigned int adj2=vtxPtrOut[i+1];
        for(unsigned int j=adj1;j<adj2;j++)
        {
                if((Gnew->edgeList[j].head !=Gnew->edgeList[j].head)  &&Gnew->edgeList[j].weight>=1){
                        checkflag=false;
	//		cout<<"kkk"<<endl;
                        edgefinal[qq].head=Gnew->edgeList[j].head;
                        edgefinal[qq].tail=Gnew->edgeList[j].tail;
			//vec[id++]=Gnew->edgeList[j].head;
			//vec[id++]=Gnew->edgeList[j].tail;
		//	if(Gnew->edgeList[j].head!=Gnew->edgeList[j].tail){
			vec.push_back(Gnew->edgeList[j].head);
			vec.push_back(Gnew->edgeList[j].tail);
		//	else
			 vec.push_back(Gnew->edgeList[j].head);

                        edgefinal[qq].weight=Gnew->edgeList[j].weight;
                        qq++;
                }
	       if((Gnew->edgeList[j].head ==Gnew->edgeList[j].head)  &&Gnew->edgeList[j].weight>=0){

			
			checkflag=false;
			edgefinal[qq].head=Gnew->edgeList[j].head;
                        edgefinal[qq].tail=Gnew->edgeList[j].tail;

		//	if(Gnew->edgeList[j].head!=Gnew->edgeList[j].tail){
                        vec.push_back(Gnew->edgeList[j].head);
                        vec.push_back(Gnew->edgeList[j].tail);
                  //      else
                         vec.push_back(Gnew->edgeList[j].head);

                        edgefinal[qq].weight=Gnew->edgeList[j].weight;
                        qq++;



			}

        }
        if(checkflag==true)
                reduce++;


}
cout<<Gnew->edgeListPtrs[Gnew->numVertices]<<" "<<qq<<endl;

/*for(int i=0;i<Gnew->edgeListPtrs[Gnew->numVertices];i++)
	cout<<edgefinal[i].head<<" "<<edgefinal[i].tail<<endl;*/
//cout<<"val of reduce ="<<reduce<<endl;
std::sort(vec.begin(),vec.begin()+vec.size());
vector<unsigned int>::iterator ip; 
ip= unique( vec.begin(), vec.end() ) ;
vec.resize(std::distance(vec.begin(), ip)); 
int *index=(int *)malloc(sizeof(int)*Gnew->numVertices);
for(int i=0;i<vec.size();i++)
	index[vec[i]]=i;
cout<<vec.size()<<endl;
/*for(int i=0;i<Gnew->numVertices;i++)
	cout<<index[i]<<"see"<<endl;
return Gnew;*/
Gnew->edgeList=(edge*)malloc(sizeof(edge)*count_edge);
Gnew->numVertices-=reduce;
Gnew->edgeListPtrs=(unsigned int*)malloc(sizeof(unsigned int)*Gnew->numVertices+1);
for(int i=0;i<=Gnew->numVertices;i++)
        Gnew->edgeListPtrs[i]=0;
for(int i=0;i<count_edge;i++)
{

        Gnew->edgeList[i].head=index[edgefinal[i].head];
        Gnew->edgeList[i].tail=index[edgefinal[i].tail];
        Gnew->edgeList[i].weight=edgefinal[i].weight;
        Gnew->edgeListPtrs[index[edgefinal[i].head]]++;
      //  Gnew->edgeListPtrs[index[edgefinal[i].tail]]++;
// cout<<Gnew->edgeList[i].head<<" "<<Gnew->edgeList[i].tail<<" "<<Gnew->edgeList[i].weight<<endl;       

}
//return Gnew;
Gnew->numEdges=count_edge;
/*for(int i=0;i<=Gnew->numVertices;i++)
	cout<<Gnew->edgeListPtrs[i]<<endl;*/
for(int i=1;i<=Gnew->numVertices;i++){
//	cout<<Gnew->edgeListPtrs[i-1]<<"pointer"<<endl;
	Gnew->edgeListPtrs[i]=Gnew->edgeListPtrs[i-1]+Gnew->edgeListPtrs[i];

}
for(int i=Gnew->numVertices;i>=0;i--)
	Gnew->edgeListPtrs[i]=Gnew->edgeListPtrs[i-1];

Gnew->edgeListPtrs[0]=0;
free(edgefinal);
free(index);
displayGraphCharacteristics(Gnew);
//return Gnew;
//cout<<Gnew->numVertices<<" "<<Gnew->numEdges<<endl;
/*for(int i=0;i<=Gnew->numVertices;i++)
	cout<<Gnew->edgeListPtrs[Gnew->numVertices]<<" ";
for(int i=0;i<Gnew->edgeListPtrs[Gnew->numVertices];i++)
         cout<<Gnew->edgeList[i].head<<" "<<Gnew->edgeList[i].tail<<" "<<Gnew->edgeList[i].weight<<endl;
//return Gnew;
//return Gnew;*/
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

for(unsigned int i=0;i<(Gnew->edgeListPtrs[Gnew->numVertices]);i++)
        {

        edge2[i].head=Gnew->edgeList[i].head;
        edge2[i].tail=Gnew->edgeList[i].tail;
        edge2[i].weight=Gnew->edgeList[i].weight;
        }

 k=0;int jj=0;
int *flag;
vector<unsigned int> pos;
pos.resize(Gnew->numVertices+1);
std::fill(pos.begin(),pos.begin()+pos.size(),0);
cout<<"while updating"<<endl;
//ofstream ch("check.txt");
int ii=0;

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
//int count1=0;int count2=0;
#pragma omp parallel for
for(unsigned int i=0;i<G1->numVertices;i++)
{



        if(dirty[i]==true)

        {	doubtindex[i]=ii;
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
		bool checkflag=true;
	
		
		__sync_fetch_and_add(&ii, 1);
        }


}
cout<<"checking"<<endl;
/*
for(int i=0;i<Gnew->edgeListPtrs[Gnew->numVertices];i++)
         cout<<Gnew->edgeList[i].head<<" "<<Gnew->edgeList[i].tail<<" "<<Gnew->edgeList[i].weight<<endl;
return Gnew;
*/
cout<<"value"<< " "<<newV<<endl;int k2=0;
#pragma omp parallel for
  for (long i=0; i<Gnew->numVertices; i++){
	delete cluPtrIn[i];	
	delete count[i];}
  free(cluPtrIn);				
  free(count);

unsigned int prevV=cc;
for(int i=1;i<=newV+prevV;i++)
	edgeP[i]+=edgeP[i-1];


Gnew->numVertices=prevV+newV+mo1->vertex;
cout<<"check"<<" "<<edgeP[prevV+newV]<<" "<<mo1->edgeno<<endl;
Gnew->numEdges=(edgeP[prevV+newV]/2)+mo1->edgeno;

Gnew->edgeListPtrs=(unsigned int*)malloc(sizeof(unsigned int)*(Gnew->numVertices+1));


for(unsigned int i=0;i<prevV+newV+1;i++){
	Gnew->edgeListPtrs[i]=edgeP[i];

	}
for(unsigned int i=prevV+newV+1;i<=Gnew->numVertices;i++)
	Gnew->edgeListPtrs[i]=0;
Gnew->edgeList=(edge *)malloc(sizeof(edge)*(Gnew->edgeListPtrs[Gnew->numVertices]+2*mo1->edgeno));
for(unsigned int i=0;i<edgeP[prevV+newV];i++)
	{

	Gnew->edgeList[i].head=edge2[i].head;
	Gnew->edgeList[i].tail=edge2[i].tail;
	Gnew->edgeList[i].weight=edge2[i].weight;


	}

 double prevMod = 1;
  double currMod = -1;
double tmpTime;
int tmpItr=0;
/*for(int i=0;i<Gnew->edgeListPtrs[Gnew->numVertices];i++)
         cout<<Gnew->edgeList[i].head<<" "<<Gnew->edgeList[i].tail<<" "<<Gnew->edgeList[i].weight<<endl;
*/
int j=0;


j=0;
qq=0;
if(mo1->vertex>0){
for(unsigned int i=(Gnew->edgeListPtrs[cc+newV]);i<Gnew->numEdges;i++)
        {
	qq=Gnew->edgeListPtrs[cc+newV]+qq;
	if(!dirtyc[mo1->edgesM[j].tail] && !dirty[mo1->edgesM[j].tail]){
        Gnew->edgeList[qq].head=mo1->edgesM[j].head+cc+newV;
        Gnew->edgeList[qq].tail=C_orig[mo1->edgesM[j].tail];
        Gnew->edgeList[qq].weight=mo1->edgesM[j].weight;qq++;
	j++;}
	else if(dirty[mo1->edgesM[j].tail])
		{
		Gnew->edgeList[qq].head=mo1->edgesM[j].head+cc+newV;
 		Gnew->edgeList[qq].tail=doubtindex[mo1->edgesM[j].tail]+cc;
        	Gnew->edgeList[qq].weight=mo1->edgesM[j].weight;
		qq++;
		}
        }
}
for(int i=0;i<=Gnew->numVertices;i++)
	Gnew->edgeListPtrs[i]=0;
for(int i=0;i<2*Gnew->numEdges;i++)
{
Gnew->edgeListPtrs[Gnew->edgeList[i].head]++;
//Gnew->edgeListPtrs[Gnew->edgeList[i].tail]++;

}

//Gnew->edgeListPtrs[0]=0;
for(int i=0;i<Gnew->numVertices;i++)
	Gnew->edgeListPtrs[i+1]+=Gnew->edgeListPtrs[i];
for(int i=Gnew->numVertices;i>=0;i--)
	Gnew->edgeListPtrs[i]=Gnew->edgeListPtrs[i-1];
Gnew->edgeListPtrs[0]=0;
cout<<"edgeListPtrs="<<Gnew->edgeListPtrs[cc+newV]<<endl;
Gnew->numEdges=(Gnew->edgeListPtrs[cc+newV]/2)+qq;
//cout<<"done3"<<endl;

C_orig=(long *)malloc(sizeof(long)*Gnew->numVertices);
unsigned int NV=Gnew->numVertices;

graph *Gnew1=(graph *)malloc(sizeof(graph));
#pragma omp parallel for
	for (long i=0; i<NV; i++) {
  	   C_orig[i] = 0;
	}
bool *dirty1,*dirty2;
//cout<<"lll"<<endl;
dirty2=(bool*)malloc(sizeof(bool)*NV);
for(int i=0;i<NV;i++)
{

	dirty2[i]=false;

}

cout<<Gnew->numVertices<<" "<<Gnew->numEdges<<endl;
/*for(int i=0;i<Gnew->numVertices+1;i++)
	cout<<Gnew->edgeListPtrs[i]<<endl;


for(int i=0;i<Gnew->edgeListPtrs[Gnew->numVertices];i++)
cout<<Gnew->edgeList[i].head<<" "<<Gnew->edgeList[i].tail<<" "<<Gnew->edgeList[i].weight<<endl;
*/
//displayGraphCharacteristics(Gnew);

graph *Gnew2=(graph *)malloc(sizeof(graph));
duplicateGivenGraph(Gnew,Gnew2);


//Gnew1=cpuonly(Gnew,Gnew1,Gnew2,C_orig,dirty3,dirty2,1);
cout<<"........."<<endl;
Gnew2=runMultiPhaseLouvainAlgorithm(Gnew, C_orig,0, 1000, 0.001, 0.001, 6);

unsigned int new1=0;

return Gnew2;

pos.clear();
free(edgeP);
free(edge2);




}
