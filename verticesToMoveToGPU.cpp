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
void verticesToMoveToGPU(graph *G,bool *dirtycpu,move1 *mo,int *c,unsigned int mid)
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

int nn=newV;
bool *flaz=(bool *)malloc(sizeof(bool)*(NV+1));
for(unsigned int i=0;i<nn;i++){

        mo[i].statIndices=0;

}
for(long i=0;i<NV;i++)  {
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
                if(dirtycpu[(vtxInd[j].tail+1)] && !flaz[(vtxInd[j].tail+1)]){
                         v.push_back(vtxInd[j].tail+1);
                         flaz[(vtxInd[j].tail+1)-1]=true;
                                                                       	}
                                
					}
        }
    }
std::sort(v.begin(), v.end());
cout<<"size="<<v.size();
ofstream fout("f1.txt");
int pos1;
int *pos=(int *)malloc(sizeof(int)*nn);
#pragma omp parallel for
for(long i=0;i<=nn;i++)
        pos[i]=0;
vector<unsigned int> *ee=new vector<unsigned int>[nn];

int f;
int ccd=0;
pos1=0;
int qq=0;
//cout<<"done"<<endl;
#pragma omp parallel for
for(long i=0;i<NV-1;i++)
{f=0;
        if(dirtycpu[i])
        {
        ccd++;

        bool bordercheck=false;
        long adj1 = vtxPtr[i];
        long adj2 = vtxPtr[i+1];

        for(long j=adj1;j<adj2;j++)
        {
                if(G->bord[i] && !bordercheck ){
                        mo[qq].statIndices=mo[qq].statIndices+G->bordno[i];
			//fout<<mo[qq].statIndices<<endl;
                        bordercheck=true;
                for(std::vector<unsigned int> ::iterator it=G->bordvalue[i].begin();it<G->bordvalue[i].end();it++){
                        ee[qq].push_back(c[*it-mid]);
		//	fout<<c[*it]<<endl;
                               pos[qq]++;
                                pos1++;}


                                         }
      if(dirtycpu[(vtxInd[j].tail)+1]){

                mo[qq].statIndices=mo[qq].statIndices+1;
                std::vector<unsigned int>::iterator itr;
                itr=std::find (v.begin(), v.end(), (vtxInd[j].tail+1));
                ee[qq].push_back((itr - v.begin()) +f);

                pos[qq]++;

                pos1++;
                f++;

                                 }
        }
	
//cout<<"1"<<endl;
//cout<<"size"<<" "<<ee[qq].size()<<endl;
mo[qq].edgesa=(unsigned int*)malloc(sizeof(unsigned int)*(ee[qq].size()));
if(ee[qq].size()!=0){
        for(std::vector<unsigned int>::iterator it=ee[qq].begin();it!=ee[qq].end();it++)
        {       int ij=0;
                mo[qq].edgesa[ij]=*it;
                ij++;
	}
                }
        qq++;
        }


}
cout<<"1"<<endl;
//return mo;
//free(vtxInd);
//free(vtxPtr);
free(flaz);
free(pos);
delete[]ee;

}
