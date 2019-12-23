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

void verticesToMoveToCPU(unsigned int *statIndices,unsigned int *edges, bool *dirtygpu,move2* mo1,unsigned int* C_orig,int NV,int total,graph *g1,graph *Gnew,unsigned int mid)
{
int newV1=0;
	for(int i=0;i<total-NV;i++)
	{
        if(dirtygpu[i]){
        newV1++;

                       }

	}
//cout<<"newV1="<<newV1<<endl;
if(newV1==0)
{
mo1->edgeListPtrsM=NULL;
mo1->vertex=0;
mo1->edgesM=NULL;
mo1->edgeno=0;
return ;

}
unsigned int *edgeListPtrsM1=(unsigned int *)malloc(newV1*sizeof(unsigned int));

for(unsigned int i=0;i<newV1;i++)
	edgeListPtrsM1[i]=0;

vector<unsigned int> e1;
int count=0;
int *pos11=(int *)malloc(sizeof(int)*(total-NV));
int j=0;
unsigned int adj1,adj2;
for(int i=0;i<total-NV;i++)
{
        if(dirtygpu[i]==true){
                pos11[i]=j;

                j++;
                }
}

unsigned int *noedge=(unsigned int *)malloc(sizeof(unsigned int)*newV1);
for(int i=0;i<newV1;i++)
        noedge[i]=0;
//cout<<"OKK"<<endl;
for(long i=0;i<total-NV;i++)

{
        int dv=0;
        int pos=0;
        if(dirtygpu[i])
        {
                if(i==0)
                {
                adj1=0;
                adj2=statIndices[(i+1)];
                }
                else{
                adj1 = statIndices[(i)];

                adj2 = statIndices[(i+1)];

                     }
/*	if(g1->bord[i+mid] )
        {       int ch=0;
                for(std::vector<unsigned int> ::iterator it=g1->bordvalue[i+mid].begin();it!=g1->bordvalue[i+mid].end();it++){
		e1.push_back(count);
                e1.push_back(C_orig[*it]);
                e1.push_back(1);
                noedge[dv]++;
                mo1[dv].edgeListPtrsM+=1;
                				}

*/
        
        int ch1=0;
	for(int j=adj1;j<adj2;j++)
        {
		 if(g1->bord[i+mid] && !dirtygpu[edges[j]])
       		 {       int ch=0;
                	for(std::vector<unsigned int> ::iterator it=g1->bordvalue[i+mid].begin();it!=g1->bordvalue[i+mid].end();it++){
               			e1.push_back(pos11[i]+count);
                		e1.push_back(C_orig[*it]);
                		e1.push_back(1);
                		noedge[dv]++;
                		edgeListPtrsM1[dv]+=1;
		//		cout<<"check edge"<<edgeListPtrsM1[dv]<<" ";
                                                }
		}
               	else  if(dirtygpu[edges[j]] && g1->bord[i+mid])
                        {
                       for(std::vector<unsigned int> ::iterator it=g1->bordvalue[i+mid].begin();it!=g1->bordvalue[i+mid].end();it++){

                				edgeListPtrsM1[dv]+=1;
                				noedge[dv]++;
                                                e1.push_back(pos11[i]+count);
                                                e1.push_back(C_orig[*it]);
                                                e1.push_back(1);
                                                }
                                                noedge[dv]++;
                                                e1.push_back(pos11[i]+count);
                                                e1.push_back(pos11[edges[(j)]]+count);
                                                e1.push_back(1);
                                                edgeListPtrsM1[dv]+=1;

                        }
 else if(dirtygpu[edges[(j)]] && !g1->bord[i+mid]){
                        noedge[dv]++;
                        e1.push_back(pos11[i]+count);
                                                e1.push_back(pos11[edges[(j)]]+count);
                                                e1.push_back(1);
                edgeListPtrsM1[dv]+=1;

                                }
        ch1=1;
        }
dv++;
        }
}
//cout<<"OK!"<<endl;
edge *edgesM1=NULL;
for(int i=0;i<newV1;i++)
{
int ij=0;
        for(std::vector<unsigned int>::iterator it=e1.begin();it!=e1.end();it+=3)
{
        edgesM1=(edge*)malloc(sizeof(edge)*noedge[i]);
for(int i=0;i<noedge[i];i++)
{

                edgesM1[ij].head=*it;
                edgesM1[ij].tail=*(it+1);
                edgesM1[ij].weight=*(it+2);
        ij++;
}
}
}
unsigned int edgec=0;
for(unsigned int i=0;i<newV1;i++)
	edgec+=noedge[i];

mo1->edgeListPtrsM=edgeListPtrsM1;
mo1->vertex=newV1;
mo1->edgesM=edgesM1;
mo1->edgeno=edgec;
//cout<<"for checking"<<mo1->edgeListPtrsM[0]<<" "<<newV1<<endl;
//free(edgesM1);

//free(edgeListPtrsM1);
free(pos11);
free(noedge);
e1.clear();























}
