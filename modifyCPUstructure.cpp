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










void modifyCPUstructure(graph *Gnew,graph *G1,bool *dirty,unsigned int *C_orig){


unsigned int newV=0;
unsigned int edgec1=0;
cout<<G1->numVertices<<" "<<G1->numEdges<<endl;

for(int i=0;i<G1->numVertices;i++)
{
        if(dirty[i]){
        newV++;
     
		}

}
cout<<"newV="<<newV<<endl;
edgec1=G1->numEdges;
cout<<Gnew->numVertices<<" "<<edgec1<<" "<<endl;
unsigned int *edgeP=(unsigned int*)malloc(sizeof(unsigned int)*(newV+Gnew->numVertices+1));
for(unsigned int i=0;i<Gnew->numVertices+1;i++)
        edgeP[i]=Gnew->edgeListPtrs[i];
for(unsigned int i=Gnew->numVertices+1;i<Gnew->numVertices+1+newV;i++)
	edgeP[i]=0;
edge *edge2=(edge*)malloc(sizeof(edge)*2*(Gnew->edgeListPtrs[Gnew->numVertices]+edgec1));
for(unsigned int i=0;i<(Gnew->edgeListPtrs[Gnew->numVertices]);i++)
        {

        edge2[i].head=Gnew->edgeList[i].head;
        edge2[i].tail=Gnew->edgeList[i].tail;
        edge2[i].weight=Gnew->edgeList[i].weight;
        }

cout<<G1->edgeListPtrs[0]<<" "<<G1->edgeListPtrs[1]<<" "<<G1->edgeListPtrs[2]<<endl;
cout<<"done"<<endl;
int k=1;int jj=0;
int *flag;
//ofstream ch("check.txt");
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
			
                        edge2[Gnew->edgeListPtrs[Gnew->numVertices]+k].head=i+(Gnew->numVertices);
                        edge2[Gnew->edgeListPtrs[Gnew->numVertices]+k].tail=C_orig[G1->edgeList[j].tail]; //check for the condition where edge connected to same community,weight will be weight+1
		
			edge2[Gnew->edgeListPtrs[Gnew->numVertices]+k].weight=1;
			                        flag[C_orig[G1->edgeList[j].tail]]+=1;

					}
		else if(flag[C_orig[G1->edgeList[j].tail]]>0)
			{
				edge2[Gnew->edgeListPtrs[Gnew->numVertices]+k].weight=flag[C_orig[G1->edgeList[j].tail]]+1;
			}
 	 k++;	
                }
		
		unsigned int adj3=Gnew->edgeListPtrs[C_orig[i]];
                unsigned int adj4=Gnew->edgeListPtrs[C_orig[i]+1];
	#pragma omp parallel for
	for(unsigned int j=adj3;j<adj4;j++)
                {
                        if(Gnew->edgeList[j].head!=Gnew->edgeList[j].tail)
                              edge2[j].weight-=1;

                      edgeP[C_orig[i]]-=1;
		      edgeP[C_orig[Gnew->edgeList[j].tail]]-=1;


                }

		
                edgeP[Gnew->numVertices+1+jj]=k;
                jj++;

        }

}
cout<<"ok"<<endl;
Gnew->edgeListPtrs=(unsigned int*)malloc(sizeof(unsigned int)*(newV+Gnew->numVertices+1));
for(unsigned int i=0;i<newV+Gnew->numVertices+1;i++)
	Gnew->edgeListPtrs[i]=edgeP[i];
free(edgeP);
/*edge *edge2=(edge*)malloc(sizeof(edge)*(Gnew->edgeListPtrs[Gnew->numVertices]+k));*/
for(unsigned int i=0;i<(Gnew->edgeListPtrs[Gnew->numVertices]);i++)
        {

        edge2[i].head=Gnew->edgeList[i].head;
        edge2[i].tail=Gnew->edgeList[i].tail;
        edge2[i].weight=Gnew->edgeList[i].weight;
        }
/*int cc=0;
for(unsigned int i=Gnew->edgeListPtrs[Gnew->numVertices];i<(Gnew->edgeListPtrs[Gnew->numVertices])+k;i++)
	{
		edge2[i].head=	e.at(cc);
		cc++;
		edge2[i].tail=e.at(cc);
		cc++;
		edge2[i].weight=e.at(cc);
		cc++;	


	}*/
Gnew->edgeList=(edge*)malloc(sizeof(edge)*(k+Gnew->edgeListPtrs[Gnew->numVertices]));
for(unsigned int i=0;i<(k+Gnew->edgeListPtrs[Gnew->numVertices]);i++)
	{

	Gnew->edgeList[i].head=edge2[i].head;
	Gnew->edgeList[i].tail=edge2[i].tail;
	Gnew->edgeList[i].weight=edge2[i].weight;
	}
//Gnew->edgeList=edge2;*/
free(edge2);
//e.clear();
free(flag);
//Gnew->numVertices=Gnew->numVertices+newV-newV11;
//Gnew->numEdges=Gnew->edgeListPtrs[Gnew->numVertices];







}
