

/*

    Copyright (C) 2016, University of Bergen

    This file is part of Rundemanen - CUDA C++ parallel program for
    community detection

    Rundemanen is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Rundemanen is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Rundemanen.  If not, see <http://www.gnu.org/licenses/>.
    
    */
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

void modifyGPUstructure(Community *dev1_community,unsigned int *statIndices,unsigned int*edges,bool *dirtyg,int *c,int total,int NV)
{
unsigned int newV1=0;
unsigned int edgec1=0;
for(int i=0;i<total-NV;i++)
{
        if(dirtyg[i]){
        newV1++; 
        unsigned int adj1=statIndices[i];

        unsigned int adj2=statIndices[i+1];
//        edgec1+=adj2-adj1;
}

}
cout<<"1"<<endl;
int k,jj;
unsigned int *edgeP=(unsigned int*)malloc(sizeof(unsigned int)*(newV1+(*dev1_community).g.nb_nodes));
for(unsigned int i=0;i<(*dev1_community).g.nb_nodes;i++)
        edgeP[i]=statIndices[i];

 k=0;jj=0;
thrust::host_vector<int> h_vec((*dev1_community).g.nb_nodes);
	int *stat=new int[(*dev1_community).g.nb_nodes];
	std::copy(((*dev1_community).g.indices).begin(),((*dev1_community).g.indices).end(), h_vec.begin());
	int i=0;
	for(int i=0;i<h_vec.size();i++)
	{
		stat[i]=h_vec[i];
	}
cout<<"2"<<endl;

	        unsigned int* edge1=new unsigned int[(*dev1_community).g.nb_links];
cout<<"4"<<endl;
	thrust::host_vector<unsigned int> h_vec1((*dev1_community).g.nb_links);
//	std::copy(((*dev1_community).g.links).begin(),((*dev1_community).g.links).end(), h_vec1.begin());
cout<<"5"<<endl;
	 i=0;
cout<<"3"<<endl;  
//return ;     
for(int i=0;i<h_vec1.size();i++)
{
	edge1[i]=h_vec1[i];
}
cout<<"okk"<<endl;
vector<int> edge2;
        int *stat1=new int[(*dev1_community).g.nb_nodes+1+newV1];
	for(int i=0;i<(*dev1_community).g.nb_nodes;i++)
		stat1[i]=stat[i];
for(int i=(*dev1_community).g.nb_nodes;i<(*dev1_community).g.nb_nodes+newV1;i++)	
		stat[i]=0;
vector<unsigned int> w;	
int k1=0;int k2=0;
std::copy((*dev1_community).g.weights.begin(),(*dev1_community).g.weights.end(),w.begin());
for(unsigned int i=0;i<total-NV;i++)
{

        if(dirtyg[i]==true)
        {
		int *flag=(int*)malloc(sizeof(int)*(*dev1_community).g.nb_nodes);
		for(int i1=0;i1<(*dev1_community).g.nb_nodes;i1++)
			flag[i]=0;
                unsigned int adj1=statIndices[i];
                unsigned int adj2=statIndices[i+1];
                for(unsigned int j=adj1;j<adj2;j++)
                {
			int a=i+(*dev1_community).g.nb_nodes;
			int b=c[edges[j]];
			if(flag[edges[j]]==0){
                        // edgeG[statIndices[(*dev1_community).g.nb_nodes-1]+k]=i;
                        edge2.push_back(i+(*dev1_community).g.nb_nodes);
			k1++;
			w.push_back(1);
                       // edgeG[statIndices[(*dev1_community).g.nb_nodes-1]+k]=c[edges[j]];
                       edge2.push_back(c[edges[j]]);
			flag[edges[j]]+=1;
                        k2++;
			}
			
//
/*		else if(flag[edges[j]]>0)
			w.at(k1);
		stat[c[i]]-=1;
		stat[c[edges[i]]]-=1;*/
		
                }
               // unsigned int adj3=statIndices[c[i]];
               // unsigned int adj4=statIndices[c[i]+1];

  edgeP[(*dev1_community).g.nb_nodes+jj]=k;
                jj++;



        }

}

//statIndices=(unsigned int*)malloc(sizeof(unsigned int)*(newV1+(*dev1_community).g.nb_nodes));
//statIndices=edgeP;
//free(edgeP);
//edges=(unsigned int*)malloc(sizeof(unsigned int)*(edgec1+statIndices[(*dev1_community).g.nb_nodes]));
//edges=edgeG;
//free(edgeG);

}
