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
vector<string> getNextLineAndSplitIntoTokens(istream& str)
{
    vector<string>   result;
    string                line;
    getline(str,line);

    stringstream          lineStream(line);
    string                cell;

    while(getline(lineStream,cell, ' '))
    {
        result.push_back(cell);
    }
    if (!lineStream && cell.empty())
    {

        result.push_back("");
    }

    
    return result;
}

void movetogpu(graph *Gnew,graph *G1,bool *dirty,bool* dirtycpu,int *statIndices,unsigned int *edges,Community *dev_community,int *c,int f)
{

int x=0;
long NV=G1->numVertices;
long NE=G1->numEdges;
int z1,z2;
//ofstream check1,check2;
//ifstream border("/home/Anwehsa/uk_border1");
int count_doubtful=0;
//bool *bord=(bool *)malloc(sizeof(bool)*(NV+1));
//long *bordno=(long *)malloc(sizeof(long)*(NV+1));
//bool *flaz=(bool *)malloc(sizeof(bool)*(NV+1));
//long *array=(long *)malloc(sizeof(long)*NV*2);
//int qq=0;
vector<long> v;
vector<long> *borderval=new vector<long>[NV+1];
//cout<<"1111"<<endl;
for(long i=0;i<NV;i++)	{
	flaz[i]=false;
	bord[i]=false;
	bordno[i]=0;
			}
for(long i=0;i<(2*NV);i++)
	array[i]=0;
for(long i=0;i<NV;i++)
{
if(dirtycpu[i]==true)
	count_doubtful++;
}
string line;
int in=0;
while(getline(border,line)){
	in++;
}

int bordere=0;
ifstream border22("/home/Anwehsa/uk_border1");
int iij=0;

while(iij<in){
	vector<string> str;
	str = getNextLineAndSplitIntoTokens(border22);
	int n = str.size();
	bordere=0;
	stringstream geek(str[0]);
	geek>>x;
        int border_vertex=x;
        bord[x]=true;
        stringstream geek1(str[1]);
        geek1 >>x;
        bordno[border_vertex]=x;
    	for(int j=2; j<n; j++)
    		{   stringstream geek2(str[j]);
			geek2>>x;
        		borderval[border_vertex].push_back(x);
		}

	iij++;
}

cout<<"#################FORMING STRUCTURE FOR MOVEMENT TO GPU#################"<<endl;
int nn=count_doubtful;
cout<<"no of doubtful"<<count_doubtful<<endl;
move1 *mo = new move1[nn];
unsigned int *edgesa=(unsigned int *)malloc(NV*sizeof(unsigned int));
for(long i=0;i<nn;i++){

mo[i].statIndices=0;
mo[i].edgesa=NULL;
}

long* vtxPtr=(G1->edgeListPtrs);
edge *vtxInd=(G1->edgeList);
for(long i=0;i<NV;i++)
	{
	if(dirtycpu[i])
	{
	bool bordercheck=false;
	long adj1 = vtxPtr[i];
	long adj2 = vtxPtr[i+1];

	for(long j=adj1;j<adj2;j++){
		if(dirtycpu[(vtxInd[j].tail+1)-1] && !flaz[(vtxInd[j].tail+1)-1]){         
          		 v.push_back(vtxInd[j].tail+1);
                   	 flaz[(vtxInd[j].tail+1)-1]=true;                                                       
									}
                    		}
	}
	}
std::sort(v.begin(), v.end());
cout<<"size="<<v.size();  
int pos1;
int *pos=(int *)malloc(sizeof(int)*nn);
ofstream ff2("edge.txt");
for(long i=0;i<=nn;i++)
	pos[i]=0;
vector<unsigned int> *ee=new vector<unsigned int>[nn];
//cout<<"check"<< qq<<" "<<ee[qq].size()<<endl;


int ccd=0;
pos1=0;
for(long i=0;i<NV;i++)
{
	if(dirtycpu[i])
	{
	ccd++;
	
	bool bordercheck=false;
	long adj1 = vtxPtr[i];
	long adj2 = vtxPtr[i+1];

	for(long j=adj1;j<adj2;j++)
        {
        	if(bord[i+1] && !bordercheck ){
                	mo[qq].statIndices=mo[qq].statIndices+bordno[i];
                	bordercheck=true;
		for(std::vector<long> ::iterator it=borderval[i].begin();it<borderval[i].end();it++){
		//	cout<<c[*it]<<endl;
			ee[qq].push_back(c[*it]);
				pos[qq]++;
				ff2<<c[*it]<<endl;					
                		pos1++;}


                               		 }
      if(dirtycpu[(vtxInd[j].tail+1)-1]){
		
             	mo[qq].statIndices=mo[qq].statIndices+1;
		std::vector<long>::iterator itr;
          	itr=std::find (v.begin(), v.end(), (vtxInd[j].tail+1)-1);
	//	cout<<((itr - v.begin()) +(*dev_community).g.nb_nodes)<<endl;
		ee[qq].push_back((itr - v.begin()) +f);
	
		ff2<<((itr - v.begin()) +f)<<endl;
		pos[qq]++;
		
                pos1++;
		
       				 }
	}


mo[qq].edgesa=(unsigned int*)malloc(sizeof(unsigned int)*(ee[qq].size()));
if(ee[qq].size()!=0){
	for(std::vector<unsigned int>::iterator it=ee[qq].begin();it!=ee[qq].end();it++)
	{	int ij=0;
		//cout<<"entered"<<endl;
		mo[qq].edgesa[ij]=*it;
		ij++;

	}
		}
	qq++;
	}

	
}
//cout<<(*dev_community).g.nb_nodes+nn<<endl;

/*for(int i=0;i<nn;i++)
	cout<<mo[i].statIndices<<endl;*/
//cout<<endl;
//cout <<f<<endl;
//return;
//delete[] ee;
//free(ee);

GraphHOST input1_graph;
//cout<<"doubt1"<<endl;
input1_graph.nb_nodes = nn+(*dev_community).g.nb_nodes;
input1_graph.degrees.resize(input1_graph.nb_nodes);
input1_graph.nb_links =(*dev_community).g.nb_links+pos1;
int * statIndices1=(int *)malloc(sizeof(int)*(*dev_community).g.nb_nodes);
int* statIndicesact=(int *)malloc(sizeof(int )*input1_graph.nb_nodes);
//cout<<"1"<<endl;
for(int i=0;i<f;i++)
        statIndices1[i]=(*dev_community).g.indices[i];
int *statIndicesa=(int*)malloc(sizeof(int)*nn);
for(long i=0;i<nn;i++){
        if(i==0)
        statIndicesa[i]=(*dev_community).g.indices[(*dev_community).g.nb_nodes-1]+mo[i].statIndices;
        else
        statIndicesa[i]=statIndicesa[i-1]+mo[i].statIndices;}

//cout<<"2"<<endl;
for(int i=0; i<(*dev_community).g.nb_nodes; i++)
                {
                statIndicesact[i]=statIndices1[i];
                }
int k;
        for(int i=0,  k=(*dev_community).g.nb_nodes; k<input1_graph.nb_nodes && i<nn; i++, k++)
                {
                statIndicesact[k]=statIndicesa[i];
                }
//cout<<"3"<<" "<<pos1<<endl;

int cccount=pos1;
/*for(int i=0;i<input1_graph.nb_nodes;i++)
	cout<<statIndicesact[i]<<" ";*/
cout<<endl;
//vector<unsigned int> eeg;
//unsigned int* edges1=(unsigned int*)malloc(sizeof(unsigned int)*cccount);
//cout<<"4"<<endl;
/*for(int i=0;i<nn;i++)
        {
        for(int j=0;j<pos[i];j++){
                eeg.push_back(mo[i].edgesa[j]);
cccount++;
                }
}*/
//cout<<"5"<<endl;
unsigned int* edges2=(unsigned int*)malloc(sizeof(unsigned int)*cccount);
//cout<<"doubt2"<<endl;

int ii=0;
ofstream ff3("ff.txt");
	for(long i=0;i<nn;i++){
		for(std:: vector<unsigned int>:: iterator it=ee[i].begin();it<ee[i].end();it++)
		{       *(edges2+ii)=*it;
			ff3<<*it<<endl;
			ii++;

		}
				}
//delete[] ee;	
double threshold=0.0001;
unsigned int* edges3=(unsigned int*)malloc(input1_graph.nb_links*sizeof(unsigned int));
for(int i=0; i<(*dev_community).g.nb_links; i++)
                {
		cout<<edges[i]<<endl;
                edges3[i]=edges[i];
                }

 for(int i=0,  k=(*dev_community).g.nb_links; k<input1_graph.nb_links && i<pos1; i++, k++)
                {
                edges3[k]=edges2[i];
                }

//cout<<"doubt3"<<endl;
		std::copy(statIndicesact, statIndicesact + (input1_graph.nb_nodes), input1_graph.degrees.begin());
		input1_graph.links.resize(input1_graph.nb_links);
std::copy(edges3, edges3 + input1_graph.nb_links, input1_graph.links.begin());
//cout<<input1_graph.nb_nodes<<endl;
ofstream ff1("pp.txt");
for(long i=0;i<input1_graph.nb_nodes;i++)
        ff1<<edges[i]<<endl;
Community *de;
Community dev_community1(input1_graph, -1, threshold);
        de=&dev_community1;
unsigned long *ss;
ss=(unsigned long *)statIndicesact;
int b=1;
       gpuonly(input1_graph,c,ss,edges3,de,dirty,dirtycpu,b);
        cout<<"moved to gpu and computation done"<<endl;
 
free(bord);
free(bordno);
free(flaz);
free(array);
free(edgesa);
delete[] borderval;
delete[] mo;
delete[] ee;
}

