
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
#include <bits/stdc++.h>

#include"defs.h"
using namespace std;

void movefinal(graph *Gnew,long *C_orig,graph*G1,unsigned long *statIndices,unsigned int *edges,Community *dev_community, int *cg,bool*borcheck)

{

long NV=G1->numVertices;
long NE=G1->numEdges;
long* vtxPtr=(G1->edgeListPtrs);
cout<<vtxPtr<<endl;
edge *vtxInd=(G1->edgeList);
cout<<vtxInd<<endl;
bool *f1=(bool*)malloc(sizeof(bool)*NV);
long *count1=(long *)malloc(sizeof(long)*(NV+1));
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
			if(borcheck[(vtxInd[j].tail+1)]){
				f1[i]=true;
				count1[i]=count1[i]+1;
				borderval[i].push_back((vtxInd[j].tail+1)-1);	}		

		}
	}

ofstream fout("c.txt");
for(int i=0;i<NV;i++)
{
	for(std::vector<long>::iterator it=borderval[i].begin();it!=borderval[i].end();it++)
	{
		fout<<*it<<endl;
	}
}


cout<<"#################FORMING STRUCTURE FOR MOVEMENT TO GPU#################"<<endl;
int nn=Gnew->numVertices;
cout <<nn<<endl;
move1 *mo = new move1[nn];
unsigned int *edgesa=(unsigned int *)malloc(NV*sizeof(unsigned int));
for(long i=0;i<nn;i++){

mo[i].statIndices=0;
mo[i].edgesa=NULL;
}

bool bordercheck=false;
vtxPtr=(Gnew->edgeListPtrs);
vtxInd=(Gnew->edgeList);


vector<unsigned int> *ee=new vector<unsigned int>[nn];

/*for(long i=0;i<nn;i++)
{

	for(long j=0;j<NV;j++)
	{
		if(C_orig[j]==i &&bord[j])
			{ mo[i].statIndices=mo[i].statIndices+bordno[j];
               		 bordercheck=true;
                	 mo[i].edgesa=addanedge(mo[i].edgesa,i,pos,NV);
                	 pos++;			
			}
		
	
	}
long adj1 = vtxPtr[i];
long adj2 = vtxPtr[i+1];
for(long j=adj1;j<adj2;j++){
	mo[i].statIndices=mo[i].statIndices+1;
        mo[i].edgesa=addanedge(mo[i].edgesa,vtxInd[j].tail+1,pos,NV);}	
		}*/

int qq;
int pos1=0;

int *pos=(int*)malloc(sizeof(int)*nn);
for(long i=0;i<nn;i++)
        pos[i]=0;


vector<long> v;
//cout<<Gnew->numVertices<<" "<<Gnew->numEdges<<endl;

for(long i=0;i<Gnew->numVertices;i++)
	{
	long adj1=vtxPtr[i];
	long adj2=vtxPtr[i+1];
	for(long j=adj1;j<adj2;j++){
//		cout<<(vtxInd[j].tail+1)-1<<" ";
		v.push_back((vtxInd[j].tail+1)-1);}
	}


vtxPtr=(Gnew->edgeListPtrs);
vtxInd=(Gnew->edgeList);
std::sort(v.begin(), v.end());
std::vector<long>::iterator ip = std::unique(v.begin(), v.begin() + v.size()); 
v.resize(std::distance(v.begin(), ip)); 
/*for(std::vector<long>::iterator it=v.begin();it!=v.end();it++)
	cout<<*it<<" ";*/
ofstream a("a.txt");
NV=Gnew->numVertices;
int ccd=0;
for(long i=0;i<NV;i++)
{
	
	qq=0; 
	ccd++;
        long adj1 = vtxPtr[i];
        long adj2 = vtxPtr[i+1];


        for(long j=adj1;j<adj2;j++)
        {
             for(long k=0;k<G1->numVertices;k++){
		        if(C_orig[k]==i && f1[k]){
				//	cout<<"entered"<<endl;
                        	mo[qq].statIndices=mo[qq].statIndices+count1[k];
                        

                       for(std::vector<long> ::iterator it=borderval[i].begin();it<borderval[i].end();it++){
	
                                ee[qq].push_back(cg[*it]);
				a<<cg[*it]<<" ";
                                pos[qq]++;
                                qq++;
                                pos1++;}



                        }                 
      } 
		

               	mo[qq].statIndices=mo[qq].statIndices+1;
                std::vector<long>::iterator itr;
                itr = std::find (v.begin(), v.end(), (vtxInd[j].tail+1)-1);
		a<<((itr - v.begin()) +(*dev_community).g.nb_nodes)<<" ";
		ee[qq].push_back((itr - v.begin()) +(*dev_community).g.nb_nodes);
                pos[qq]++;
                qq++;
                pos1++;

                                 
        }

mo[qq].edgesa=(unsigned int*)malloc(sizeof(unsigned int)*(ee[qq].size()));
cout<<"size="<<ee[qq].size()<<endl;
int ij=0;
if(ee[qq].size()!=0){
        for(std::vector<unsigned int>::iterator it=ee[qq].begin();it<ee[qq].end();it++)
        {
                mo[qq].edgesa[ij]=*it;

                ij++;


        }
		}

qq++;
        


}

//cout<<(*dev_community).g.nb_nodes<<" "<<(*dev_community).g.nb_links<<endl;
//return ;
GraphHOST input1_graph;
cout<<"doubt1"<<endl;
input1_graph.nb_nodes = nn+(*dev_community).g.nb_nodes;
input1_graph.degrees.resize(input1_graph.nb_nodes);
input1_graph.nb_links =(*dev_community).g.nb_links+pos1;
int * statIndices1=(int *)malloc(sizeof(int)*(*dev_community).g.nb_nodes);
int* statIndicesact=(int *)malloc(sizeof(int )*input1_graph.nb_nodes);
cout<<"1"<<endl;
for(int i=0;i<(*dev_community).g.nb_nodes;i++)
        statIndices1[i]=(*dev_community).g.indices[i];

int *statIndicesa=(int*)malloc(sizeof(int)*nn);
for(long i=0;i<nn;i++){
        if(i==0)
        statIndicesa[i]=(*dev_community).g.indices[(*dev_community).g.nb_nodes-1]+mo[i].statIndices;
        else
        statIndicesa[i]=statIndicesa[i-1]+mo[i].statIndices;}

cout<<"2"<<endl;
for(int i=0; i<(*dev_community).g.nb_nodes; i++)
                {
                statIndicesact[i]=statIndices1[i];
                }
int k;
        for(int i=0,  k=(*dev_community).g.nb_nodes; k<input1_graph.nb_nodes && i<nn; i++, k++)
                {
                statIndicesact[k]=statIndicesa[i];
                }
cout<<"3"<<" "<<pos1<<endl;

int cccount=pos1;
cout<<endl;
vector<unsigned int> eeg;
cout<<"5"<<endl;
unsigned int* edges2=(unsigned int*)malloc(sizeof(unsigned int)*cccount);
cout<<"doubt2"<<endl;

int ii=0;
ofstream ff3("ff.txt");
        for(long i=0;i<nn;i++){
                for(std:: vector<unsigned int>:: iterator it=ee[i].begin();it<ee[i].end();it++)
                {       *(edges2+ii)=*it;
                     //   ff3<<*it<<endl;
                        ii++;

                }
                                }


double threshold=0.0001;
unsigned int* edges3=(unsigned int*)malloc(input1_graph.nb_links*sizeof(unsigned int));
for(int i=0; i<(*dev_community).g.nb_links; i++)
                {
                ff3<<edges[i]<<endl;
                edges3[i]=edges[i];
                }

 for(int i=0,  k=(*dev_community).g.nb_links; k<input1_graph.nb_links && i<pos1; i++, k++)
                {
                edges3[k]=edges2[i];
                }

cout<<"doubt3"<<endl;
                std::copy(statIndicesact, statIndicesact + (input1_graph.nb_nodes), input1_graph.degrees.begin());
                input1_graph.links.resize(input1_graph.nb_links);
std::copy(edges3, edges3 + input1_graph.nb_links, input1_graph.links.begin());
cout<<input1_graph.nb_nodes<<endl;
ofstream ff1("pp.txt");
for(long i=0;i<input1_graph.nb_nodes;i++)
        ff1<<edges[i]<<endl;
Community *de;
Community dev_community1(input1_graph, -1, threshold);
de=&dev_community1;
unsigned long *ss;
ss=(unsigned long *)statIndicesact;
int b=1;
bool *dirtycpu=(bool*)malloc(sizeof(bool)*(G1->numVertices));
bool *dirty=(bool*)malloc(sizeof(bool)*(G1->numVertices));

for(int i=0;i<G1->numVertices;i++){
	dirtycpu[i]=false;
	dirty[i]=false;}

gpuonly(input1_graph,cg,ss,edges3,de,dirty,dirtycpu,b);
        cout<<"moved to gpu and computation done"<<endl;


/*GraphHOST input1_graph;

input1_graph.nb_nodes = nn+dev_community.g.nb_nodes;
input1_graph.degrees.resize(input1_graph.nb_nodes);

unsigned long * statIndices1=(unsigned long*)malloc(sizeof(unsigned long)*dev_community.g.nb_nodes);
unsigned long * statIndicesact=(unsigned long*)malloc(sizeof(unsigned long)*input1_graph.nb_nodes);
cout<<"1"<<endl;
for(int i=0;i<dev_community.g.nb_nodes;i++)
        statIndices1[i]=dev_community.g.indices[i];

int *statIndicesa=(int*)malloc(sizeof(int)*nn);
        for(long i=0;i<nn;i++){
        	if(i==0)
        		statIndicesa[i]=dev_community.g.indices[dev_community.g.nb_nodes-1]+mo[i].statIndices;
        	else
        		statIndicesa[i]=statIndicesa[i-1]+mo[i].statIndices;
				}
cout<<"2"<<endl;
for(int i=0; i<dev_community.g.nb_nodes; i++)
                {
                statIndicesact[i]=statIndices1[i];
                }
int k;
for(int i=0,  k=dev_community.g.nb_nodes; k<input1_graph.nb_nodes && i<nn; i++, k++)
{
	statIndicesact[k]=statIndicesa[i];
}
//cout<<"3"<<endl;
int cccount=0;
vector<unsigned int> eeg;
unsigned int* edges1=(unsigned int*)malloc(sizeof(unsigned int)*cccount);
//cout<<"4"<<endl;
unsigned int* edges2=(unsigned int*)malloc(sizeof(unsigned int)*cccount);
cout<<"doubt2"<<endl;
int ii=0;
double threshold=0.0001;
unsigned int* edges3=(unsigned int*)malloc(input1_graph.nb_links*sizeof(unsigned int));
for(int i=0; i<dev_community.g.nb_links; i++)
                {
                edges3[i]=edges[i];
                }
 for(int i=0,  k=dev_community.g.nb_links; k<input1_graph.nb_links && i<pos1; i++, k++)
                {
                edges3[k]=edges2[i];
                }
cout<<"doubt3"<<endl;

std::copy(edges3, edges3 + input1_graph.nb_links, input1_graph.links.begin());
ofstream ff1("pp.txt");
for(long i=0;i<input1_graph.nb_links;i++)
        ff1<<edges3[i]<<endl;
Community dev_community1(input1_graph, -1, threshold);
//int a=gpuonly(input1_graph,cg,statIndicesact,edges3,dev_community1,f1);
*/
 }
