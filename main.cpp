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

#include "thrust/device_vector.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <algorithm>


#include"fstream"
#include "iostream"
#include "graphHOST.h"
#include "graphGPU.h"
#include "communityGPU.h"
#include"list"
#include<math.h>
using namespace std;

int main(int argc, char** argv) {


	char* file_w = NULL;
	int type = UNWEIGHTED;

	std::cout << "#Args: " << argc << std::endl;
	for (int i = 0; i < argc; i++) {
		std::cout << i << " : " << argv[i] << std::endl;
	}

	if (argc == 3) {
		file_w = argv[2];
		type = WEIGHTED;
		if (type == WEIGHTED)
			std::cout << "Weighted Graph \n";
	}

	if (file_w)
		std::cout << "inputGraph: " << argv[1] << " Corresponding Weight: " << file_w << std::endl;
	else if (argc==2)
		std::cout << "inputGraph: " << argv[1] << std::endl;
	else 
		std::cout<<"No input graph provided, creating a sample graph"<<std::endl;
double partratio=0.1;
long NV=39459925 ;// vertices and 100805268 links
//long noedge=83427536;
long noedge=1842690156;
  unsigned long *statIndices=(unsigned long *)malloc((NV)*sizeof(unsigned long));
  unsigned int *edges=(unsigned int*)malloc(noedge*sizeof(unsigned int));

for(long i=0;i<(NV);i++)
        statIndices[i]=0;
 cout<<"check 1"<<endl;
for(long i=0;i<noedge;i++)
        edges[i]=0;
        cout<<"check2"<<endl;
  GraphHOST input_graph(argv[1], file_w, type,statIndices,edges);
        cout<<"graphhost okk"<<endl;
	// Read Graph in  host memory
//	GraphHOST input_graph(argv[1], file_w, type);

	//Create a graph in host memory
/*	GraphHOST input_graph; // Sample graph
	if (1) {

		input_graph.nb_nodes = 7;
		input_graph.degrees.resize(input_graph.nb_nodes);
		int statIndices[] = {3, 4, 5, 8, 10, 13, 14};
		std::copy(statIndices, statIndices + (input_graph.nb_nodes), input_graph.degrees.begin());

		input_graph.nb_links = 14;
		input_graph.links.resize(input_graph.nb_links);
		unsigned int edges[] = {1, 2, 3, 0, 0, 0, 4, 5, 3, 5, 3, 4, 6, 5};
		std::copy(edges, edges + input_graph.nb_links, input_graph.links.begin());

	}*/

	//input_graph.display();

	double threshold = 0.000001;
	if(argc==4) threshold = atof(argv[2]);
	double binThreshold = 0.01;
	if(argc==4) binThreshold=atof(argv[3]);
	//binThreshold=threshold;
	//Copy Graph to Device
	Community dev_community(input_graph, -1, threshold);
	double cur_mod = -1.0, prev_mod = 1.0;
	bool improvement = false;

	std::cout << "threshold: " << threshold << " binThreshold: " << binThreshold << std::endl;

	//Read Prime numbers
	dev_community.readPrimes("fewprimes.txt");

	cudaStream_t *streams = NULL;
	int n_streams = 8;

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
int node=dev_community.g.nb_nodes;

std::vector<int>::iterator it;
vector<int> *pp=new vector<int>[node];

	//clock_t clk_decision, clk_contraction;
	std::vector<clock_t> clkList_decision;
	std::vector<clock_t> clkList_contration;
 int *nc=(int*)malloc((dev_community.g.nb_nodes)*sizeof(int));
int *c=(int *)malloc(node*sizeof(int));
        int *nc1=(int*)malloc(node*sizeof(int));
                int *index=(int *)malloc(node*sizeof(int));
                int counter=0;
                int node1;

	clock_t t1, t2, t3;
	t1 = clock();

	/*			
				dev_community.preProcess();
				dev_community.gatherStatistics(true); // MUST be "true" to filter out isolated vertices at the beginning
				dev_community.compute_next_graph();
				dev_community.set_new_graph_as_current();
	 */
	bool TEPS = true;
	bool islastRound = false;
	int szSmallComm = 100000;
	bool isGauss =true;// false;

	if(isGauss)
		std::cout<<"\n Update method:  Gauss–Seidel (in batch) \n";
	else
		std::cout<<"\n Update method: Jacobi\n";

	int stepID = 1;

	do {

		std::cout << "---------------Calling method for modularity optimization------------- \n";
		t2 = clock();
		prev_mod = cur_mod;

		cur_mod = dev_community.one_levelGaussSeidel(cur_mod, islastRound,
				szSmallComm, binThreshold, isGauss &&(dev_community.community_size > szSmallComm),
				streams, n_streams, start, stop);

		t2 = clock() - t2;

		clkList_decision.push_back(t2); // push the clock for the decision

		std::cout<< "step: " <<stepID <<", Time for modularity optimization: " << ((float) t2) / CLOCKS_PER_SEC << std::endl;
		stepID++;
		if (TEPS == true) {
			std::cout<<binThreshold<<"_"<<threshold<< " #E:" <<  dev_community.g.nb_links << "  TEPS: " << dev_community.g.nb_links / (((float) t2) / CLOCKS_PER_SEC) << std::endl;
			TEPS = false;
		}

		//break;
		std::cout << "Computed modularity: " << cur_mod << " ( init_mod = " << prev_mod << " ) " << std::endl;

		if ((cur_mod - prev_mod) > threshold) {

			t2 = clock();   
			t3 = t2;
			dev_community.gatherStatistics();
			t2 = clock() - t2;
			//std::cout << "T_gatherStatistics: " << ((float) t2) / CLOCKS_PER_SEC << std::endl;

			t2 = clock();
		nc=dev_community.compute_next_graph(streams, n_streams, start, stop);
			t2 = clock() - t2;
if(counter==0)
                          {
                for(int i=0;i<node;i++)
                                 c[i]=nc[i];

                        for(int i=0;i<node;i++){
                                     nc1[i]=nc[i];

                                                }
                                }

                        else if(counter>0)
                         {
                        cout<<"no_of nodes in this iteration"<<dev_community.g.nb_nodes<<endl;
 for(int i=0;i<dev_community.g.nb_nodes;i++)
                                        {
                                for (std::vector<int>::iterator itr = pp[i].begin(); itr != pp[i].end(); ++itr)
                                        {  c[*itr]=nc[i];

                                        }
                                                }
 int *b=(int*)malloc(sizeof(int)*node);int c1=0;
                                for(int i=0;i<node;i++)
                                        {
                                                b[i]=nc1[i];
                                                c1=nc[b[i]];
                                                b[i]=c1;
                                                nc1[i]=b[i];
                                        }
                        }

                        for(int i=0;i<node;i++)

                        {
                        pp[nc1[i]].push_back(i);
                        }
node1=dev_community.g.nb_nodes;
counter=counter+1;
t2 = clock() - t2;

			std::cout << "Time to compute next graph: " << ((float) t2) / CLOCKS_PER_SEC << std::endl;

			//clkList_contration.push_back(t2); // push the clock for the contraction

			t2 = clock();
			dev_community.set_new_graph_as_current();
			t2 = clock() - t2;
			t3 = clock() -t3;

			clkList_contration.push_back(t3); // push the clock for the contraction

			//std::cout << "T_new_graph_as_current: " << ((float) t2) / CLOCKS_PER_SEC << std::endl;
			// break;

			//std::cout << "\n Back to Main \n";
			//int sc;
			//std::cin>>sc;


		} else {
			if (islastRound == false) {
				islastRound = true;
			} else {
				break;
			}
		}
		//improvement = false;
	} while (true);

	std::cout<< "#phase: "<<stepID<<std::endl; 

	t2 = clock();
	float diff = ((float) t2 - (float) t1);
	float seconds = diff / CLOCKS_PER_SEC;

	if( argc ==1){
		std::cout <<  binThreshold<<"_"<<threshold<<" Running Time: " << seconds << " ;  Final Modularity: "
			<< prev_mod  << std::endl;
	}else{

		std::cout <<  binThreshold<<"_"<<threshold<<" Running Time: " << seconds << " ;  Final Modularity: "
			<< prev_mod << " inputGraph: " << argv[1] << std::endl;
	}


	std::cout << "#Record(clk_optimization): " << clkList_decision.size()
		<< " #Record(clk_contraction):" << clkList_contration.size() << std::endl;

	int nrPhase = std::min(clkList_decision.size(), clkList_contration.size());

	std::ofstream ofs ("time.txt", std::ofstream::out);

	//------------------------------------To Plot------------------------------------//
	/*  
	    int nrphaseToPlot = 40;
	    for (int i = 0; i < std::min(nrphaseToPlot, nrPhase); i++) {
	    std::cout << (float) clkList_decision[i] / CLOCKS_PER_SEC << " ";
	    }
	    std::cout << std::endl;

	    for (int i = 0; i < std::min(nrphaseToPlot, nrPhase); i++) {
	    std::cout << (float) clkList_contration[i] / CLOCKS_PER_SEC << " ";
	    }
	    std::cout << std::endl;
	 */


	//----------------------------------------------------//

	float t_decision = 0, t_contraction = 0;

	for (int i = 0; i < clkList_decision.size(); i++) {

		t_decision += (float) clkList_decision[i] / CLOCKS_PER_SEC;
		if(i<nrPhase) ofs<< (float) clkList_decision[i] / CLOCKS_PER_SEC<<" ";
		else std::cout<<  (float) clkList_decision[i] / CLOCKS_PER_SEC<<" -> "<<std::endl;

	}

	ofs<<"\n";

	for (int i = 0; i < clkList_contration.size(); i++) {
		t_contraction += (float) clkList_contration[i] / CLOCKS_PER_SEC;
		if(i<nrPhase) ofs<< (float) clkList_contration[i] / CLOCKS_PER_SEC<<" ";

	}

	ofs<<"\n";
	ofs.close();

	std::cout<< " Optimization and contraction time  ratio:"
		<< (100 * t_decision)/(t_decision + t_contraction) << " " << (100 * t_contraction)/(t_decision+t_contraction) << std::endl;


	/*
	   for (int i = nrPhase; i < clkList_decision.size(); i++) {

	   std::cout << clkList_decision[i] << " "
	   << (float) clkList_decision[i] / CLOCKS_PER_SEC << std::endl;
	   }
	 */

	std::cout << "(graph):      #V  " << dev_community.g.nb_nodes << " #E   " << dev_community.g.nb_links << std::endl;
	std::cout << "(new graph)  #V  " << dev_community.g_next.nb_nodes << " #E  " << dev_community.g_next.nb_links << std::endl;
	/* 
	   for (int i = 0; i < n_streams; i++) {
	//CHECK(cudaStreamDestroy(streams[i]));
	}
	 */

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	//free(streams);
bool *flag=(bool*)malloc(node*sizeof(bool));
for(long i=0;i<node;i++)
        flag[i]=false;
long* int_deg=(long *)malloc(node*sizeof(long));
for(long i=0;i<node;i++)
        int_deg[i]=0;
long *max=(long *)malloc(dev_community.g.nb_nodes*sizeof(long));
for(long i=0;i<dev_community.g.nb_nodes;i++)
        max[i]=0;
long p,q;
for(long i=0;i<node;i++)
{
if(i==0)
{p=0;q=statIndices[0];}
else{
 p=statIndices[i-1];
 q=statIndices[i];}
for(long j=p;j<q;j++)
{
//cout<<"is it ok?"<<endl;
unsigned int f=edges[j];
//cout<<"yes"<<endl;
if(c[f]==c[i])
{
//cout<<"thanks"<<endl;
flag[f]=flag[i]=true;
int_deg[i]=int_deg[i]+1;
int_deg[f]=int_deg[f]+1;
//cout<<"god"<<endl;
if(max[c[f]]<int_deg[f])
                max[c[f]]=int_deg[f];
//cout<<"god"<<endl;

if(max[c[i]]<int_deg[i])
                max[c[i]]=int_deg[i];
}


}


}
cout<<"loop checked"<<endl;
bool *dirty=(bool *)malloc(sizeof(bool)*node);
for(long i=0;i<node;i++)
        dirty[i]=false;
double thresd=0.59;
double *rel=(double *)malloc(node*sizeof(double));
for(long i=0;i<node;i++)
{
rel[i]=log(int_deg[i]+1)/log(max[c[i]]+1);
if(rel[i]<thresd)
dirty[i]=true;

}

//double elapsed_secs1 = double(endg - beging) / CLOCKS_PER_SEC;
 //cout << "time1="<<elapsed_secs1<<endl;

cout<<"###########doubtful computation done###############"<<endl;
for(long i=0;i<dev_community.g.nb_links;i++)
	cout<<"edge"<<edges[i]<<endl;
/*cout<<"############similarity measure###################"<<endl;

long *cn=(long*)malloc((node+1)*sizeof(long));
for(long i=0;i<node;i++){
cn[i]=0;
}

vector<int>* gc;
cout <<"!"<<endl;
gc = new vector<int>[node];
for(long i=0;i<node;i++)
 {
  if(dirty[i+1]==true)
    {
    int c=0;
    int adj1 = statIndices[i];
    int adj2 = statIndices[i+1];
        for(int j=adj1; j<adj2; j++) {
 gc[i].push_back(edges[j]);

                                        }


    }
 }
long adj1,adj2;
//cout <<statIndices[0]<<" "<<statIndices[1]<<endl;
for(long i=0;i<node;i++)
 {
if(dirty[i+1]==true)
    {int c=0;
        if(i==0)
        {adj1=0;        adj2=statIndices[0];}


else{
     adj1 = statIndices[i-1];
     adj2 = statIndices[i];}
        for(long j=adj1; j<adj2; j++) {
long adj3=statIndices[edges[j]];
                long adj4=statIndices[edges[j+1]];
                for(long k=adj3;k<adj4;k++){
 for(std::vector<int>::iterator it=gc[i].begin();it<gc[i].end();it++)
                                {
                                        if(*it==edges[i])
                                        cn[i]=cn[i]+1;
                                }
                                             }
                                        }
        }
}
cout<<"#####################similarity ends############################"<<endl;
*/
free(dirty);
free(flag);
free(rel);
//free(cn);
free(max);
free(int_deg);

	return 0;
}


