// ***********************************************************************
//
//            Grappolo: A C++ library for graph clustering
//               Mahantesh Halappanavar (hala@pnnl.gov)
//               Pacific Northwest National Laboratory     
//
// ***********************************************************************
//
//       Copyright (2014) Battelle Memorial Institute
//                      All rights reserved.
//
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions 
// are met:
//
// 1. Redistributions of source code must retain the above copyright 
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright 
// notice, this list of conditions and the following disclaimer in the 
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its 
// contributors may be used to endorse or promote products derived from 
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN 
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************



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
double calculateRatioCD(graph graph, unsigned int NUM_VER){
        bool *ghostVertex;
        ghostVertex = (bool*) calloc(graph.nnodes, sizeof(bool));
        //mark other vertices as ghost so that while calculation it doesn't take care of these vertices
        #pragma omp parallel for num_threads(16) schedule(static)
        for(unsigned int i = NUM_VER; i < graph.nnodes; i++)
                ghostVertex[i] = true;

        unsigned int *currCtx, currActiveNodes, currActiveEdges;
        unsigned int *Parents;
        currCtx = (unsigned int*) malloc( NUM_VER * sizeof(unsigned int));
        std::copy(graph.srcsrc,  graph.srcsrc+NUM_VER, currCtx);
        Parents = (unsigned int*) malloc(graph.nnodes * sizeof(unsigned int));
        std::copy(graph.srcsrc, graph.srcsrc + graph.nnodes, Parents);

        currActiveNodes = NUM_VER;
        currActiveEdges = graph.psrc[NUM_VER];

        double cpuTime = 0, gpuTime = 0;
        unsigned int *hostParents, *deviceParents;
        unsigned int *cpuCurrCtx, *gpuCurrCtx;
        omp_set_num_threads(2);
        omp_set_nested(1);
        int part = -1;
	#pragma omp parallel sections
        {
                #pragma omp section
                {
                        hostParents = (unsigned int*) malloc( graph.nnodes * sizeof(unsigned int));
                        std::copy( Parents, Parents + graph.nnodes,  hostParents);                  
                        cpuCurrCtx = (unsigned int*) malloc( currActiveNodes * sizeof(unsigned int));
                        std::copy( currCtx, currCtx + currActiveNodes, cpuCurrCtx);
                        //cpuTime = cdCPU( graph, hostParents, cpuCurrCtx, currActiveNodes, currActiveEdges, ghostVertex, part);
			
                        free( cpuCurrCtx);
                        free( hostParents);
                }
                #pragma omp section
                {
                        deviceParents = (unsigned int*) malloc( graph.nnodes * sizeof(unsigned int));
                        std::copy( Parents, Parents + graph.nnodes,  deviceParents);
                        gpuCurrCtx = (unsigned int*) malloc( currActiveNodes * sizeof(unsigned int));
                        std::copy( currCtx, currCtx + currActiveNodes, gpuCurrCtx);
                        //gpuTime = cdGPU( graph, deviceParents, gpuCurrCtx, currActiveNodes, currActiveEdges, ghostVertex, part);
                        free( gpuCurrCtx);
                        free( deviceParents);

                }
        }
        std::cout<<"CPU Exe Time: "<<cpuTime<<" GPU Exe Time: "<<gpuTime<<"\n";
        free( Parents);
        free( currCtx);
        free( ghostVertex);

return gpuTime/(cpuTime+gpuTime);
}
double partGraphCD( graph graph, int parts, double *sizeofData){
     /*   int NUM_VER = 10000;
        if( graph.nnodes < NUM_VER)
                  NUM_VER =  graph.nnodes / 2;
        double ST = get_time();
        double partitionRatio = calculateRatioCD(graph, NUM_VER);
        std::cout<<"CPU & GPU Exe time: "<<get_time()-ST<<"\n";
        *sizeofData += 3 * graph.nnodes * sizeof(unsigned int);
        *sizeofData += graph.nedges * sizeof(unsigned int);
        *sizeofData/= (1024 * 1024);
        std::cout<<"Size of data: "<<*sizeofData<<" MB\n";
	if( *sizeofData > 3200){
                double otherPart = (double)(1 - partitionRatio);
                if( (otherPart * graph.nedges) < 500000000){
                        while( (otherPart * graph.nedges) < 500000000){
                                partitionRatio -= 0.01;
                                otherPart = (double)(1 - partitionRatio);
                   //     std::cout<<"Ratio: "<<partitionRatio<<"\n";
                        }
                }
        }*/
	double partitionRatio=0.5;
	graph.partGraph( 2, partitionRatio);
        return partitionRatio;
}

/*void indCompCC( Graph graph, unsigned int *&Parents){
        unsigned int *hostParents, *deviceParents;
        omp_set_num_threads(2);
        omp_set_nested(1);
        unsigned int *cpuCurrCtx, *gpuCurrCtx;
        unsigned int cpuCurrNodes, gpuCurrNodes;
        unsigned int cpuCurrEdges, gpuCurrEdges;
        bool *cpuGhostNodes, *gpuGhostNodes;
        int cpuPart = 0, gpuPart = 1;
        #pragma omp parallel sections
        {
                #pragma omp section
                {
                        cpuCurrNodes = graph.PerPartNodes[ cpuPart];
                        cpuCurrEdges = graph.PerPartEdges[ cpuPart];
                        std::cout<<"CPU Part Nodes: "<<cpuCurrNodes<<" Edges: "<<cpuCurrEdges<<"\n";
                        hostParents = (unsigned int *) malloc( graph.nnodes * sizeof(unsigned int));
                        cpuCurrCtx = (unsigned int *) malloc(cpuCurrNodes * sizeof(unsigned int));
                        cpuGhostNodes = (bool*) calloc( graph.nnodes, sizeof(bool));

                        std::copy(Parents, Parents+graph.nnodes, hostParents);
                        std::copy( graph.srcsrc+graph.startVertex[ cpuPart], graph.srcsrc+graph.endVertex[ cpuPart], cpuCurrCtx);
                        #pragma omp parallel for schedule(static) num_threads(8)
                        for(unsigned int i = graph.endVertex[ cpuPart]; i < graph.nnodes; i++)
                                cpuGhostNodes[i] = true;
			ccCPU( graph, hostParents, cpuCurrCtx, cpuCurrNodes, cpuCurrEdges, cpuGhostNodes, cpuPart);
                        std::copy(hostParents, hostParents+cpuCurrNodes, Parents);
                        
			free( hostParents);
                        free( cpuCurrCtx);
                        free( cpuGhostNodes);
                }

		#pragma omp section
                {
                        gpuCurrNodes = graph.PerPartNodes[ gpuPart];
                        gpuCurrEdges = graph.PerPartEdges[ gpuPart];
                        std::cout<<"GPU Part Nodes: "<<gpuCurrNodes<<" Edges: "<<gpuCurrEdges<<"\n";
                        deviceParents = (unsigned int*) malloc( graph.nnodes * sizeof(unsigned int));
                        gpuCurrCtx = (unsigned int *) malloc( gpuCurrNodes * sizeof(unsigned int));
                        gpuGhostNodes =  (bool*) calloc( graph.nnodes, sizeof(bool));

                        std::copy( Parents, Parents+graph.nnodes, deviceParents);
                        std::copy( graph.srcsrc+graph.startVertex[ gpuPart], graph.srcsrc+graph.endVertex[ gpuPart], gpuCurrCtx);
                        #pragma omp parallel for schedule(static) num_threads(8)
                        for(unsigned int i = 0; i < graph.startVertex[ gpuPart]; i++)
                                gpuGhostNodes[i] = true;

//                        deviceFindMst( graph, deviceParents, gpuCurrCtx, gpuCurrNodes, gpuCurrEdges, gpuGhostNodes, gpuPart);
			ccGPU( graph, deviceParents, gpuCurrCtx, gpuCurrNodes, gpuCurrEdges, gpuGhostNodes, gpuPart);
                        std::copy(deviceParents+cpuCurrNodes, deviceParents +graph.nnodes, Parents+cpuCurrNodes);
                        //Free Mem
                        free( deviceParents);
                        free( gpuCurrCtx);
                        free( gpuGhostNodes);
                }
      }
}
*/

int main(int argc, char** argv) {
 


  graph* G = (graph *) malloc (sizeof(graph));
  graph* G1 = (graph *) malloc (sizeof(graph));
  graph* G2 = (graph *) malloc (sizeof(graph));
  graph* Gnew=(graph *) malloc (sizeof(graph));

  graph g1;
  g1.readFromGR(argv[1]);
  double exeTime = get_time();
  double partTime = get_time();
  double sizeofData = 0;
  int parts = 2;
  double partitionRatio=0.92;
  unsigned int mid=g1.partGraph( 2, partitionRatio);

  unsigned int total=g1.nnodes;
  unsigned int NV=mid;
	
  g1.numVertices=mid;
  g1.numEdges=g1.no1/2;
//  cout<<NV<<" "<<g.no1<<endl;
  bool *dirty=new bool[NV];
  bool *dirtycpu=new bool[NV];
  bool *dirtygpu=new bool[total-NV]; //need to change
  bool *dirtyg=new bool[total-NV]; //need to change
  //cout<<"4"<<endl;
  unsigned int *C_orig=(unsigned int *)malloc(sizeof(unsigned int)*NV);
  for(long i=0;i<NV;i++)
	dirtycpu[i]=false;
  for(long i=0;i<total-NV;i++){
	dirtygpu[i]=false;
	dirtyg[i]=false;}
  int coloring = 0;
  G=&g1;
//	cout<<"5"<<endl;
/*	ofstream fout("deg.txt");
	ofstream fout1("edge.txt");
 for(int i=0;i<NV;i++)
	fout<<g.edgeListPtrs[i]<<endl;
	for(int i=0;i<g.no1;i++)
	fout1<<g.edgeList[i].head<<" "<<g.edgeList[i].tail<<" "<<g.edgeList[i].weight<<endl;*/
  displayGraphCharacteristics(&g1);
  duplicateGivenGraph(G,G1);
  duplicateGivenGraph(G,G2);

  for(long i=0;i<NV;i++)
	C_orig[i]=-1;
  /*for(int i=0;i<g.numVertices;i++)
	cout<<g.edgeListPtrs[i]<<" ";*/

   	GraphHOST input_graph;
      	input_graph.nb_nodes = total-mid;
      	input_graph.degrees.resize(input_graph.nb_nodes);

      	for(int i=0;i<input_graph.nb_nodes;i++)
          input_graph.degrees.at(i)=g1.edgeListPtrs1[i];
   
 	unsigned int *statIndices=(unsigned int *)malloc(input_graph.nb_nodes*sizeof(unsigned int));
        unsigned int *edges=(unsigned int*)malloc(g1.no2*sizeof(unsigned int));

 	input_graph.nb_links =g1.no2;
    	input_graph.links.resize(input_graph.nb_links);
    	//input_graph.links=g.links;
	        for(int i=0;i<input_graph.nb_links;i++)
          input_graph.links.at(i)=g1.links.at(i);
	
	
        double threshold = 0.000001;
	Community *dev1_community;
	Community dev_community(input_graph, -1, .000001);
	dev1_community=&dev_community;
	int *c=(int *)malloc((total-mid)*sizeof(int));
	int i=0;

	for(std::vector<unsigned long>::iterator it=input_graph.degrees.begin();it<input_graph.degrees.end();it++)
                {*(statIndices+i)=*it;
                        i++;
                }
	i=0;
//	return 0;
	for(std::vector<unsigned int>::iterator it=input_graph.links.begin();it<input_graph.links.end();it++)
                {*(edges+i)=*it;
                        i++;
                }

	#pragma omp parallel sections
        {
                #pragma omp section
                {
                 
	int f=gpuonly(input_graph,c,statIndices,edges,dev1_community,dirtyg,dirtygpu,0,G1,mid);
//	cout<<(*dev1_community).g.total_weight<<endl;
	print_vector((*dev1_community).g.weights,"weights");
	print_vector((*dev1_community).g.links,"links");
	cout<<"gpu ends"<<endl;
		}

		#pragma omp section
                {
                  Gnew=cpuonly(G, Gnew,G1, C_orig,dirty,dirtycpu);
	cout<<Gnew->numVertices<<endl;
			cout<<Gnew->numEdges<<endl;
		}

	}
return 0;
 // duplicateGivenGraph(&g1,G);

//mo is the structure to move to GPU partition
unsigned int edgec=0;
unsigned int newV=0;
for(int i=0;i<NV;i++)
{
        if(dirtycpu[i]){
        newV++;
			}
}


int nn=newV;
move1 *mo = new move1[nn];

verticesToMoveToGPU(G1,dirtycpu,mo,c,mid);
cout<<"!"<<endl;
//mo1 structure to move from GPU to CPU
//unsigned int edgec=0;
//verticesToMoveToCPU((*dev1_community),dirtygpu,mo1,C_orig,total,NV);
unsigned int newV1=0;
for(int i=0;i<total-NV;i++)
{
        if(dirtygpu[i]){
        newV1++;
	
						}						

}
move2 *mo1 = new move2[newV1];
verticesToMoveToCPU(statIndices,edges,dirtygpu,mo1,C_orig,total,NV,&g1,Gnew,mid);
cout<<"@"<<endl;
cout<<"after"<<G1->edgeListPtrs[0]<<" "<<G1->edgeListPtrs[1]<<endl;
//unsigned int *edgeListPtrsM=(unsigned int *)malloc(newV1*sizeof(unsigned int));
modifyCPUstructure(Gnew,G1,dirty,C_orig);
cout<<"*"<<endl;
//void modifyGPUstructure(Community *dev1_community,unsigned int *statIndices,unsigned int*edges,bool *dirtyg,int *c)
modifyGPUstructure(dev1_community,statIndices,edges,dirtyg,c,total,NV);
cout<<"&"<<endl;
unsigned int newV11=newV;
unsigned int newV12=newV1;
 edgec=0;
 newV1=0;




//edge count 


//Gnew->numVertices=Gnew->numVertices+newV-newV11;
//Gnew->numEdges=Gnew->edgeListPtrs[Gnew->numVertices];
//return 0;
//cpu graph structure changed now change gpu graph structure
 

unsigned int * C_orig1=(unsigned int*)malloc((Gnew->numVertices)*sizeof(unsigned int));

graph* Gnew1 = (graph *) malloc (sizeof(graph));
for(long i=0;i<Gnew->numVertices;i++)
	C_orig1[i]=-1;
//cout<<Gnew->numVertices<<" "<<Gnew->numEdges<<endl;
/*for(int i=0;i<=Gnew->numVertices;i++)
	cout<<Gnew->edgeListPtrs[i]<<" ";*/
//check GPU modified graph and previous cpu weight  update
//according to mo and mo1 change Gnew and statVertices and edges
//Gnew1=runMultiPhaseLouvainAlgorithm(Gnew, C_orig1,0, 1000, 0.0001, .0001, 5);


return 0;


//cout << typeid((*dev_community).g.nb_links).name() << endl;
/*
unsigned int *edgeP1=(unsigned int*)malloc(sizeof(unsigned int)*((*dev_community).g.nb_nodes+newV1));
for(unsigned int i=0;i<(*dev_community).g.nb_nodes;i++)
        edgeP1[i]=statIndices[i];


unsigned int  *edge21=(unsigned int*)malloc(sizeof(unsigned int )*(edgec1+(*dev_community).g.nb_links));
for(unsigned int i=0;i<(statIndices[(*dev_community).g.nb_links]);i++)
        {

     //   edge2[i].head=Gnew->edgeList[i].head;
      //  edge2[i].tail=Gnew->edgeList[i].tail;
      //  edge2[i].weight=Gnew->edgeList[i].weight;
      edge21[i]=edges[i];
        }







//unsigned 


/*Gnew->edgeList1=(edge *)malloc(sizeof(edge)*(Gnew->edgeListPtrs[Gnew->numVertices+1]));
for(int i=0;i<Gnew->numVertices+1;i++)
	cout<<Gnew->edgeListPtrs[i]<<" ";
cout<<endl;
for(int i=0;i<Gnew->edgeListPtrs[Gnew->numVertices];i++)
	cout<<Gnew->edgeList[i].head<<" "<<Gnew->edgeList[i].tail<<" "<<Gnew->edgeList[i].weight<<endl;
int j=0;
unsigned int newV=0;
unsigned int newE=0;
for(unsigned int i=0;i<G->numVertices;i++)
	{

		if(dirty[i])
		{	
			newV++;
			
		}

	}

for(unsigned int i=0;i<(Gnew->edgeListPtrs[Gnew->numVertices+1]);i++)
{
		if(!dirty[i]) //if the vertex is doubtful
		{
		
			//	Gnew->edgeListPtrs[i]=Gnew->edgeListPtrs[i-1];
			Gnew->edgeList1[j].head=Gnew->edgeList[i].head;
			Gnew->edgeList1[j].tail=Gnew->edgeList[i].tail;
			Gnew->edgeList1[j].weight=Gnew->edgeList[i].weight;
			j++;
			
		}
	
unsigned int*edgeListPtrs1=(unsigned int *)malloc(sizeof(unsigned int)*(Gnew->numVertices+1));
for(unsigned int i=0;i<Gnew->numVertices+1;i++)
	edgeListPtrs1[i]=Gnew->edgeListPtrs[i];
Gnew->edgeListPtrs=(unsigned int *)malloc(sizeof(unsigned int)*(Gnew->numVertices+newV+1));
for(int i=0;i<Gnew->numVertices+1;i++)
	Gnew->edgeListPtrs[i]=edgeListPtrs1[i];
for(int i=Gnew->numVertices+2;i<Gnew->numVertices+1+newV;i++)
	Gnew->edgeListPtrs[i]=0;
//add doubtful vertex edges
int cc=1;//to check cumulative degree sequence
for(unsigned int i=0;i<2*G->numEdges;i++)
        {
                        unsigned int*weigh=(unsigned int *)malloc(sizeof(unsigned int)*(Gnew->numVertices));
		int k=0;
                if(dirty[G->edgeList[i].head] && !dirty[G->edgeList[i].tail])
                {
		
			unsigned int adj1=G->edgeListPtrs[i];
			unsigned int adj2=G->edgeListPtrs[i+1];
			for(unsigned int jj=adj1;jj<adj2;jj++)
				{	
                        Gnew->edgeList1[j].head=G->edgeList[i].head;
                        Gnew->edgeList1[j].tail=G->edgeList[c[i]].tail;
                        Gnew->edgeList1[j].weight=1;
                        j++;
			Gnew->edgeListPtrs[c[(Gnew->edgeList1[i].head+1)]+1]-=1;
                        Gnew->edgeListPtrs[c[(Gnew->edgeList1[i].tail+1)]+1]-=1;
			Gnew->edgeListPtrs[Gnew->numVertices+cc]+=1;
			}
			cc++;
                }
		if(dirty[G->edgeList[i].tail] && !dirty[G->edgeList[i].head])
		{
			 Gnew->edgeList1[j].head=G->edgeList[c[i]].head;
                        Gnew->edgeList1[j].tail=G->edgeList[i].tail;
                        Gnew->edgeList1[j].weight=weigh[c[i]];
                        j++;
                        Gnew->edgeListPtrs[(Gnew->edgeList1[i].head+1)-(Gnew->numVertices)]-=1;
                        Gnew->edgeListPtrs[(Gnew->edgeList1[i].tail+1)-(Gnew->numVertices)]-=1;
                        Gnew->edgeListPtrs[cc]+=1;
                        cc++;



		}
		if(dirty[G->edgeList[i].tail] && dirty[G->edgeList[i].head])
                {
                         Gnew->edgeList1[j].head=G->edgeList[i].head;
                        Gnew->edgeList1[j].tail=G->edgeList[i].tail;
                        Gnew->edgeList1[j].weight=1;
                        j++;
                        Gnew->edgeListPtrs[(Gnew->edgeList1[i].head+1)-(Gnew->numVertices)]-=1;
                        Gnew->edgeListPtrs[(Gnew->edgeList1[i].tail+1)-(Gnew->numVertices)]-=1;
                        Gnew->edgeListPtrs[cc]+=1;
                        cc++;



                }



        }

if(dirty[i] && (Gnew->edgeList1[i].head!=Gnew->edgeList1[i].tail) ){
                        Gnew->edgeListPtrs[Gnew->edgeList1[i].head+1]-=1;
                        Gnew->edgeListPtrs[Gnew->edgeList1[i].tail+1]-=1;

}*
else if(dirty[i] && (Gnew->edgeList1[i].head==Gnew->edgeList1[i].tail))
                                        Gnew->edgeListPtrs[Gnew->edgeList1[i].head+1]-=1;

}*/

/*for(int i=0;i<Gnew->numVertices+1;i++)
        cout<<Gnew->edgeListPtrs[i]<<" ";
cout<<endl;
for(int i=0;i<Gnew->edgeListPtrs[Gnew->numVertices];i++)
        cout<<Gnew->edgeList[i].head<<" "<<Gnew->edgeList[i].tail<<endl;*/
//Gnew->edgeList=Gnew->edgeList1;
	/*
	thrust::host_vector<int> h_vec((*dev1_community).g.nb_nodes+1);
	int *stat=new int[(*dev1_community).g.nb_nodes+1];
	
	unsigned int* edge1=new unsigned int[(*dev1_community).g.nb_links];
	
	std::copy(((*dev1_community).g.indices).begin(),((*dev1_community).g.indices).end(), h_vec.begin());
	
	int i=0;
	
	for(int i=0;i<h_vec.size();i++)
	{
		stat[i]=h_vec[i];
	}
	
	thrust::host_vector<unsigned int> h_vec1((*dev1_community).g.nb_links);
//	cout<<"links"<<((*dev1_community).g.nb_links)<<endl;
//	cout<<"no of links"<<((*dev1_community).g.links).size()<<endl;

	std::copy(((*dev1_community).g.links).begin(),((*dev1_community).g.links).end(), h_vec1.begin());
	 i=0;
        
//ofstream d("d.txt");
for(int i=0;i<h_vec1.size();i++)
{
	edge1[i]=h_vec1[i];
//	d<<edge1[i]<<" ";
}



	GraphHOST modified_graph;
	modified_graph.nb_nodes = (*dev1_community).g.nb_nodes;
	modified_graph.degrees.resize(modified_graph.nb_nodes+1);
	std::copy(stat, stat + (modified_graph.nb_nodes), modified_graph.degrees.begin());
	modified_graph.nb_links = (*dev1_community).g.nb_links;
	modified_graph.links.resize(modified_graph.nb_links);
	std::copy(edge1, edge1 + modified_graph.nb_links, modified_graph.links.begin());
	Community *dev11_community;
        Community dev12_community(modified_graph, -1, threshold);
        dev11_community=&dev12_community;

	bool *borcheck=(bool*)malloc(NV1*sizeof(bool));
	for(int i=0;i<NV1;i++)
		borcheck[i]=false;
        cout<<(*dev11_community).g.nb_nodes<<" "<<(*dev11_community).g.nb_links<<endl;
		
	movetogpu(Gnew,G1,dirty,dirtycpu,stat,edge1,dev11_community,c,f);
	
	
	Gnew1= movetocpu(dirtygpu,NV1,statIndices,edges,Gnew,opts,G1,C_orig,borcheck);
	cout<<G2->numVertices<<" "<<G2->numEdges<<endl;
	cout<<Gnew1->numVertices<<" "<<Gnew1->numEdges<<endl;
	
	thrust::host_vector<int> h1_vec((*dev1_community).g.nb_nodes+1);
        int *stat1=new int[(*dev11_community).g.nb_nodes+1];
        unsigned int *edge2=new unsigned int[(*dev11_community).g.nb_links];
        std::copy(((*dev11_community).g.indices).begin(),((*dev11_community).g.indices).end(), h1_vec.begin());
         

        for(int i=0;i<h1_vec.size();i++)
        {
                stat1[i]=h1_vec[i];
        }

        thrust::host_vector<unsigned int> h1_vec1((*dev11_community).g.nb_links);
        std::copy(((*dev11_community).g.links).begin(),((*dev11_community).g.links).end(), h1_vec1.begin());
         i=0;

//ofstream d("d.txt");
	for(int i=0;i<h1_vec1.size();i++)
	{
        edge2[i]=h1_vec1[i];
//	d<<edge2[i]<<" ";
	}
	
	
	GraphHOST modified1_graph;
        modified1_graph.nb_nodes = (*dev11_community).g.nb_nodes;
        modified1_graph.degrees.resize(modified1_graph.nb_nodes);
        std::copy(stat, stat + (modified1_graph.nb_nodes), modified1_graph.degrees.begin());
        modified1_graph.nb_links = (*dev11_community).g.nb_links;
        modified1_graph.links.resize(modified1_graph.nb_links);
        std::copy(edge2, edge2 + modified1_graph.nb_links, modified1_graph.links.begin());
	ofstream b("b.txt");
	for(int i=0;i<(*dev11_community).g.nb_links;i++)
		b<<edge2[i]<<" ";	

        Community *devf_community;
        Community dev1f_community(modified1_graph, -1, threshold);
        devf_community=&dev1f_community;

 	cout<<"new gpu"<<(*devf_community).g.nb_nodes<<" "<<(*devf_community).g.nb_links<<endl;
	
	movefinal(Gnew1,C_orig,G2,statIndices,edge2,devf_community, c,borcheck);
	std::cout<<"Everything is done!!!!"<<std::endl;

*/
//free(borcheck);	
//free(statIndices);
//free(edges); 
//free(dirtycpu);
//free(dirtygpu);
// free(G1); 
  return 0;
}//End of main()
