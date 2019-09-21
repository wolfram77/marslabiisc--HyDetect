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
vector<string> getNextLineAndSplitIntoTokens3(istream& str)
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

graph* cpuonly(graph *G, graph *Gnew,graph* G1, unsigned int *C_orig,bool *dirty1,bool *dirty2)
{//clock_t beginc,endc;
//beginc=clock();  
//double time1, time2;
 // graph* Gnew = (graph *) malloc (sizeof(graph));
 //raph* G1 = (graph *) malloc (sizeof(graph));
//graph g;
//G=&g;
unsigned int NE        = G->numEdges;
unsigned int *vtxPtr   = G->edgeListPtrs;
duplicateGivenGraph(G,G1);

cout<<"//////////////"<<vtxPtr<<"//////////////"<<endl;
edge * vtxInd   = G->edgeList;
unsigned int NV = G->numVertices;

//cout<<"within RC "<<NV<<endl;
//cout<<G->edgeListPtrs[NV-1]<<endl;
//C_orig = (unsigned int *) malloc (NV * sizeof(unsigned int)); assert(C_orig != 0);
double *deg = (double *) malloc (NV* sizeof(double)); assert(deg != 0);
double *val = (double *) malloc (NV * sizeof(double)); assert(val != 0);
 double *int_deg = (double *) malloc (NV* sizeof(double));
double *relk= (double *) malloc (NV * sizeof(double));
double *rcw= (double *) malloc (NV * sizeof(double));  /*relative commitment */
bool *flag = (bool *) malloc (NV* sizeof(bool));
bool *flag1 = (bool *) malloc (NV* sizeof(bool));
bool *flag2  = (bool *) malloc (NV* sizeof(bool));
//double *max = (double *) malloc (NV * sizeof(double)); assert(max != 0);
/*for(unsigned int i=0;i<NV;i++)
max[i]=0;*/
cout<<"allocation done"<<endl;
//duplicateGivenGraph(G,G1);
/*if(1 ) {
    	printf("Vertex following is enabled.\n");
      //  time1 = omp_get_wtime();
        unsigned int numVtxToFix = 0; //Default zero
        unsigned int *C = (unsigned int *) malloc (G->numVertices * sizeof(unsigned int)); assert(C != 0);
	numVtxToFix = vertexFollowing(G,C); //Find vertices that follow other vertices
	if( numVtxToFix > 0) {  //Need to fix things: build a new graph		      
	        printf("Graph will be modified -- %ld vertices need to be fixed.\n", numVtxToFix);          
		graph *Gnew = (graph *) malloc (sizeof(graph));
		unsigned int numClusters = renumberClustersContiguously(C, G->numVertices);		         
		buildNewGraphVF(G, Gnew, C, numClusters);
  		//Get rid of the old graph and store the new graph
		free(G->edgeListPtrs);
		free(G->edgeList);
		free(G);
		G = Gnew;		
	}
	free(C); //Free up memory
	printf("Graph after modifications:\n");
	displayGraphCharacteristics(G);
   }*/
NV = G1->numVertices;
vtxInd   = G->edgeList;

cout<<"vertex Following"<<endl;

   //Call the clustering algorithm: 
  //No strong scaling -- run once with max threads
	#pragma omp parallel for
	for (unsigned int i=0; i<NV; i++) {
  	   C_orig[i] = -1;
	}	

Gnew=runMultiPhaseLouvainAlgorithm(G, C_orig,0, 1000, 0.001, 0.001, 6); 

/*for(int i=0;i<9;i++)
	cout<<C_orig[i]<<" ";*/
//cout<<"#######%%%%%%"<<Gnew->numVertices<<" "<<Gnew->numEdges<<endl;  
double *max = (double *) malloc (Gnew->numVertices * sizeof(double)); assert(max != 0);
 for(unsigned int i=0;i<Gnew->numVertices;i++)
max[i]=0;
cout<<"############################################compute doubtful###############################################"<<endl;

//cout<<G1->numVertices<<endl;
//cout<<G->numVertices<<endl;

//clock_t begin2 = clock();

vtxPtr   = G1->edgeListPtrs;
vtxInd   = G1->edgeList;
//cout<<G1->numVertices<<endl;
  //displayGraphCharacteristics(G1);
//cout<<vtxPtr[1]<<endl;
//cout<<vtxInd[8].head<<endl;
#pragma omp parallel for
for (unsigned int i = 0; i < NV; i++) {
	
    unsigned int adj1 = vtxPtr[i];
    unsigned int adj2 = vtxPtr[i+1];
//	cout<<adj1<<" "<<adj2<<endl;
	for(unsigned int j=adj1; j<adj2; j++) {
//	cout<<"entered"<<endl;
//	cout<<(vtxInd[j].head+1)-1<<endl;
//	cout<<"enter1"<<endl;
        if(C_orig[i]==C_orig[(vtxInd[j].tail+1)-1]   ){
                deg[i]=deg[i]+1;
//		cout<<"check"<<endl;
                int_deg[(vtxInd[j].tail+1)-1]=int_deg[(vtxInd[j].tail+1)-1]+1;
                if(max[C_orig[(vtxInd[j].tail+1)-1]]<int_deg[(vtxInd[j].tail+1)-1])
                max[C_orig[(vtxInd[j].tail+1)-1]]=int_deg[(vtxInd[j].tail+1)-1];
                }
 }
  }

cout<<"ok1"<<endl;
for(unsigned int i=0;i<NV;i++){
if(deg[i]==0 && int_deg[i]==0){
flag[i]=true;
relk[i]=0;}
else if(int_deg[i]==0){
flag[i]=true;
relk[i]=0;}
else {
relk[i]=log(int_deg[i]+1)/log(max[C_orig[i]]+1);}
//cout<<relk[i]<<endl;
}

//cout<<"ok2"<<endl;
double thres=0.9;
//dirty1 = (bool *) malloc ((NV+1)* sizeof(bool)); 
//dirty2 = (bool *) malloc ((NV+1)* sizeof(bool));
//bool *dirty3 = (bool *) malloc ((NV+1)* sizeof(bool));

for(unsigned int i=0;i<NV;i++)
{
dirty1[i]=false;
dirty2[i]=false;
//dirty3[i]=false;
}
for(unsigned int i=0;i<NV;i++)

{
if(relk[i]<thres)
dirty1[i]=true;
}
int count =0;
for(int i=0;i<NV;i++)
{
	if(dirty1[i]==true)
		count++;


}
cout<<"No of doubtful vertices="<<count<<endl;
//clock_t end2 = clock();
 // double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;
// cout << "time2="<<elapsed_secs2<<endl;

cout<<"&&&&&&&&&&&&&doubtful end&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;

// endc = clock();

//double elapsed_secs1 = double(endc - beginc) / CLOCKS_PER_SEC;
// cout << "time1="<<elapsed_secs1<<endl;
 

cout<<"####################similarity measure#########################"<<endl;

/*unsigned int *cn=(unsigned int*)malloc((NV+1)*sizeof(unsigned int));
for(unsigned int i=0;i<NV+1;i++){
cn[i]=0;
}
bool *d=new bool[NV];
vector<int>* gc;
gc = new vector<int>[NV+1];
for(unsigned int i=0;i<NV;i++)
 {
  if(dirty1[i+1]==true)
    {
    int c=0;
    unsigned int adj1 = vtxPtr[i];
    unsigned int adj2 = vtxPtr[i+1];
        for(unsigned int j=adj1; j<adj2; j++) {
                gc[i].push_back((vtxInd[j].tail+1)-1);

                                        }


    }
 }
int counter=0;
unsigned int maxcn=0;
#pragma omp parallel for num_threads(16) collapse(3)
for(unsigned int i=0;i<NV;i++)
 {

  if(dirty1[i+1]==true)
    {
    int c=0;
    unsigned int adj1 = vtxPtr[i];
    unsigned int adj2 = vtxPtr[i+1];
        for(unsigned int j=adj1; j<adj2; j++) {
        unsigned int adj3=vtxPtr[(vtxInd[j].tail+1)-1];
                unsigned int adj4=vtxPtr[(vtxInd[j].tail+1)];
                for(unsigned int k=adj3;k<adj4;k++){
                        for(std::vector<int>::iterator it=gc[i].begin();it!=gc[i].end();it++)
                                {	counter++;
					cout<<*it<<" "<<"counter="<<counter<<endl;
                                        if(*it==(vtxInd[j].tail+1)-1)
                                        cn[i]=cn[i]+1;
					
                                }
                                             }
                                        }
        }
if(cn[i]>maxcn)
maxcn=cn[i];



}
//cout<<"okk"<<endl;
//maxcn=1;
for(unsigned int i=0;i<NV;i++)
{
double thres=cn[i]/maxcn;
if(thres<0.5)
dirty2[i]=true;
}
*/
/*bool *bord=(bool *)malloc(sizeof(bool)*(NV+1)); 
unsigned int *bordno=(unsigned int *)malloc(sizeof(unsigned int)*(NV+1));

string line;
int in=0;
ifstream border("/home/Anwesha/uk_border1");
while(getline(border,line)){
        in++;
}

int bordere=0;
ifstream border22("/home/Anwesha/uk_border1");
int iij=0;
int x=0;
while(iij<in){
        vector<string> str;
        str = getNextLineAndSplitIntoTokens3(border22);
        int n = str.size();
        bordere=0;
        stringstream geek(str[0]);
        geek>>x;
        int border_vertex=x;
        bord[x]=true;
        stringstream geek1(str[1]);
        geek1 >>x;
        bordno[border_vertex]=x;
       

        iij++;
}
*/

NV=G->numVertices;
NE=G->numEdges;

vtxPtr   = G1->edgeListPtrs;
vtxInd   = G1->edgeList;
/*cout<<"border"<<endl;
for(int i=0;i<G1->numVertices;i++)
{

	if(dirty1[i])
	cout<<i<<" "<<G->bordno[i]<<endl;

}*/
//cout<<G->bordno[1]<<endl;
//cout<<G->bord[1]<<endl;
//return Gnew;
double thres2=0.01;
 count=0;
for(int i=0;i<NV;i++)
	{
		
		if(dirty1[i]==true && G->bord[i])

		{//	cout<<G1->bordno[i]<<endl;
		//	cout<<"entered"<<endl;
			int nonborderedge=0;
			unsigned int borderedge=G->bordno[i];		
		//	cout<<borderedge<<endl;	
			unsigned int adj1 = vtxPtr[i];
    			unsigned int adj2 = vtxPtr[i+1];
		//	cout<<adj1<<" "<<adj2<<endl;
			for(unsigned int j=adj1; j<adj2; j++) {
        			if(C_orig[i]==C_orig[(vtxInd[j].tail+1)])
					nonborderedge++;
			}
	if(i==6)	cout<<"nonborderedge"<<nonborderedge<<endl;
//		cout<<(borderedge/nonborderedge)<<endl;	
		if(nonborderedge!=0 &&(borderedge/nonborderedge)>thres2){
			dirty2[i]=true;

			}

		}
	//	else if(dirty1[i]==true && !G1->bord[i])
	//		dirty2[i]=true;
	

	}
//cout<<count<<endl;
cout<<"#####################similarity ends##########################"<<endl;

free(deg);
free(val);
free(int_deg);
free(relk);
free(rcw);
free(flag1);
free(flag2);
free(max);
free(dirty1);
return Gnew;
//free(dirty2);
//free(cn);


}

