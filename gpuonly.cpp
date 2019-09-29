#include <bits/stdc++.h>
#include"defs.h"
#include <thrust/for_each.h>
#include <thrust/iterator/counting_iterator.h>
#include<thrust/host_vector.h>
#include<thrust/device_vector.h>
#include <thrust/transform.h>
#include<thrust/functional.h>
using namespace std;
vector<string> getNextLineAndSplitIntoTokens4(istream& str)
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
/*struct doubtful_functor
{

  __host__ __device__
  void operator()(unsigned long  x)
  {
    
              printf("%d\n");
   }
};*/
struct sum_Functor {
    unsigned long *sum;
    sum_Functor(){}
    __host__ __device__
    void operator()(unsigned long i)
    {
        //*sum+=i;
        printf("In functor: i %d",i);

    }
//printf("inside first functor");

};
int  gpuonly(GraphHOST input_graph,unsigned int*c,unsigned int *statIndices,unsigned int *edges,Community *dev_community,bool *dirty1,bool*dirty2,int b,graph* G, unsigned int mid)
{

cout<<"performing independent gpu computation"<<endl;
clock_t beging,endg;
beging=clock();
double threshold = 0.000001;
double binThreshold = 0.01;
double cur_mod,prev_mod;

cur_mod = -1.0; prev_mod = 0.0;

//cur_mod=0;prev_mod=0.88;
bool improvement = false;
bool TEPS = true;
bool islastRound = false;
int szSmallComm = 10000;
bool isGauss =true;
std::cout << "threshold: " << threshold << " binThreshold: " << binThreshold << std::endl;
//cout<<"Check"<<endl;
(*dev_community).readPrimes("fewprimes.txt");
//cout<<"1"<<endl;
cudaStream_t *streams = NULL;
int n_streams = 8;
cudaEvent_t start, stop;
//cout<<"2"<<endl;
cudaEventCreate(&start);
//cout<<"3"<<endl;
cudaEventCreate(&stop);
std::vector<clock_t> clkList_decision;
std::vector<clock_t> clkList_contration;
clock_t t1, t2, t3;
t1 = clock();
//cout<<"what"<<endl;
int node=(*dev_community).g.nb_nodes;
std::vector<int>::iterator it;
vector<int> *pp=new vector<int>[node];
int che;
if(isGauss)
	std::cout<<"\n Update method:  Gauss–Seidel (in batch) \n";
else
	std::cout<<"\n Update method: Jacobi\n";

//cout<<"1"<<endl;
int stepID = 1;
int *nc=(int*)malloc(((*dev_community).g.nb_nodes)*sizeof(int));
int *nc1=(int*)malloc(node*sizeof(int));
int *index=(int *)malloc(node*sizeof(int));
int counter=0;
int node1;
cout<<"###################initialization done##############################"<<endl;
 do {
	std::cout << "---------------Calling method for modularity optimization------------- \n";
	t2 = clock();
	prev_mod = cur_mod;
	cur_mod = (*dev_community).one_levelGaussSeidel(cur_mod, islastRound,szSmallComm, binThreshold, isGauss &&((*dev_community).community_size > szSmallComm),streams, n_streams, start, stop);
	t2 = clock() - t2;
	cout<<cur_mod<<endl;
	clkList_decision.push_back(t2); 
	std::cout<< "step: " <<stepID <<", Time for modularity optimization: " << ((float) t2) / CLOCKS_PER_SEC << std::endl;
	stepID++;
	if (TEPS == true) {
	std::cout<<binThreshold<<"_"<<threshold<< " #E:" << (* dev_community).g.nb_links << "  TEPS: " <<( *dev_community).g.nb_links 		/(((float) t2) / CLOCKS_PER_SEC) << std::endl;
 	TEPS = false;
                	  }
std::cout << "Computed modularity: " << cur_mod << " ( init_mod = " << prev_mod << " ) " << std::endl;
cout <<"*******************************************************"<<(*dev_community).g.nb_nodes<<endl;
	if ((cur_mod - prev_mod) > threshold) {
		t2 = clock();   
		t3 = t2;
		(*dev_community).gatherStatistics();                                                                                                 t2 = clock() - t2;
		t2 = clock();
                nc=(*dev_community).compute_next_graph(streams, 8, start, stop);
		if(counter==0)
                	{
                        for(int i=0;i<node;i++)
                                { c[i]=nc[i];
			//		cout<<"c"<<c[i]<<endl;
				}
			for(int i=0;i<node;i++){
                                     nc1[i]=nc[i];
			//		cout<<"nc1"<<nc1[i]<<endl;
                                                }
			}

              	else if(counter>0)
                         {
			for(int i=0;i<(*dev_community).g.nb_nodes;i++)
                        	{
					for (std::vector<int>::iterator itr = pp[i].begin(); itr != pp[i].end(); ++itr)
                                        	{      
							c[*itr]=nc[i];

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
//cout<<nc1[3]<<endl;
             	for(int i=0;i<node;i++)
		{
			pp[nc1[i]].push_back(i);
              	}
//	cout<<"okk"<<endl;
	node1=(*dev_community).g.nb_nodes;
	counter=counter+1;
	t2 = clock() - t2;
	std::cout << "Time to compute next graph: " << ((float) t2) / CLOCKS_PER_SEC << std::endl;
	t2 = clock();
	(*dev_community).set_new_graph_as_current();
	t2 = clock() - t2;
	t3 = clock() -t3;
	} else 
		{	
		if (islastRound == false) {
                	islastRound = true;
                        		} 
			else {
                                break;
                              }
                }
		che++;
        } while (true);
/*cout<<"community"<<endl;
for(int i=0;i<node;i++)
	cout<<c[i]<<" ";*/
        std::cout<< "#phase: "<<stepID<<std::endl;

        t2 = clock();
        float diff = ((float) t2 - (float) t1);
        float seconds = diff / CLOCKS_PER_SEC;

     /*   if( argc ==1){
                std::cout <<  binThreshold <<"_"<<threshold<<" Running Time: " << seconds << " ;  Final Modularity: "
                        << prev_mod  << std::endl;
        }else{

                std::cout <<  binThreshold<<"_"<<threshold<<" Running Time: " << seconds << " ;  Final Modularity: "
                        << prev_mod << " inputGraph: " << argv[1] << std::endl;
        }*/


        std::cout << "#Record(clk_optimization): " << clkList_decision.size()
                << " #Record(clk_contraction):" << clkList_contration.size() << std::endl;

        int nrPhase = std::min(clkList_decision.size(), clkList_contration.size());

        std::ofstream ofs ("time.txt", std::ofstream::out);



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

        std::cout << "(graph):      #V  " <<(* dev_community).g.nb_nodes << " #E   " << (*dev_community).g.nb_links << std::endl;
        std::cout << "(new graph)  #V  " <<(* dev_community).g_next.nb_nodes << " #E  " <<(* dev_community).g_next.nb_links << std::endl;


cout<<"#########################GPU computation ends############################################"<<endl;

//return (* dev_community).g.nb_nodes;
/*for(int i=0;i<node;i++)
	cout<<c[i]<<" ";*/
cout <<"########################calculate doubtful vertices  in the gpu##########################"<<endl;

bool *flag=(bool*)malloc(node*sizeof(bool));
for(long i=0;i<node;i++)
	flag[i]=false;
long* int_deg=(long *)malloc(node*sizeof(long));
for(long i=0;i<node;i++)
	int_deg[i]=0;
long *max=(long *)malloc((*dev_community).g.nb_nodes*sizeof(long));
for(long i=0;i< (*dev_community).g.nb_nodes;i++)
	max[i]=0;
long p,q;
for(long i=0;i<node;i++)
{
	if(i==0)
	{
		p=0;
		q=statIndices[0];
	}
	else{
 		p=statIndices[i-1];
 		q=statIndices[i];
    	}
	for(long j=p;j<q;j++)
	{
	unsigned int f=edges[j];
		if(c[f]==c[i])
		{

		flag[f]=flag[i]=true;
		int_deg[i]=int_deg[i]+1;
		int_deg[f]=int_deg[f]+1;
		if(max[c[f]]<int_deg[f])
			max[c[f]]=int_deg[f];

		if(max[c[i]]<int_deg[i])
                	max[c[i]]=int_deg[i];
		}		


	}


}


//copy statIndices in a host_vector

//thrust::device_vector<unsigned long> vec(statIndices, statIndices+ node);
//thrust::device_vector<unsigned long> dirtygpu;
//thrust::device_vector<unsigned long> dd;
//for(unsigned long i=0;i<5;i++)
//thrust::transform(vec.begin(),vec.end(),dd,dirtygpu,sum_Functor());
//unsigned long i=0;
//sum_Functor(&i);

int sum=0;
//dirty1=(bool *)malloc(sizeof(bool)*node);
for(long i=0;i<node;i++)
	dirty1[i]=false;
double thresd=0.5;
double *rel=(double *)malloc(node*sizeof(double));
for(long i=0;i<node;i++)
{
rel[i]=log(int_deg[i]+1)/log(max[c[i]]+1);
if(rel[i]<thresd)
dirty1[i]=true;	

}

double elapsed_secs1 = double(endg - beging) / CLOCKS_PER_SEC;
 cout << "time1="<<elapsed_secs1<<endl;
if(b==1)
return 0;
cout<<"###########doubtful computation done###############"<<endl;
cout<<"############similarity measure###################"<<endl;
//return (* dev_community).g.nb_nodes;

/*long *cn=(long*)malloc((node+1)*sizeof(long));
for(long i=0;i<node;i++){
cn[i]=0;
}

vector<int>* gc;
//cout <<"!"<<endl;
gc = new vector<int>[node];
for(long i=0;i<node;i++)
 {
  if(dirty[i+1]==true)
    {
    int c=0;
    int adj1 = statIndices[i];
    int adj2 = statIndices[i+1];
        for(int j=adj1; j<adj2; j++) {
           //     long adj3=statIndices[edges[j]];
             //   long adj4=statIndices[edges[j+1]];
                gc[i].push_back(edges[j]);

                                        }


    }
 }
//cout<<"@"<<endl;
long adj1,adj2;
cout <<statIndices[0]<<" "<<statIndices[1]<<endl;
for(long i=0;i<node;i++)
 {
//cout<<"#"<<endl;
  if(dirty[i+1]==true)
    {int c=0;
	if(i==0)
	{adj1=0;	adj2=statIndices[0];}
	
   
else{
     adj1 = statIndices[i-1];
     adj2 = statIndices[i];}
        for(long j=adj1; j<adj2; j++) {
		//cout<<"$"<<endl;
        long adj3=statIndices[edges[j]];
                long adj4=statIndices[edges[j+1]];
                for(long k=adj3;k<adj4;k++){
				//cout<<"%"<<endl;
                          for(std::vector<int>::iterator it=gc[i].begin();it<gc[i].end();it++)
                                {
                                        if(*it==edges[i])
                                        cn[i]=cn[i]+1;
                                }
                                             }
                                        }
        }

  }
*/
/*bool *bord=(bool *)malloc(sizeof(bool)*(node+1));
long *bordno=(long *)malloc(sizeof(long)*(node+1));
dirty2=new bool[(*dev_community).g.nb_nodes];
string line;
int in=0;
ifstream border("/home/Anwesha/arabic_border2");
while(getline(border,line)){
        in++;
}

int bordere=0;
ifstream border22("/home/Anwesha/arabic_border2");
int iij=0;
int x=0;
while(iij<in){
        vector<string> str;
        str = getNextLineAndSplitIntoTokens4(border22);
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
}*/
cout<<"check1"<<endl;
/*for(long i=0;i<(*dev_community).g.nb_nodes;i++)
{

	if(dirty1[i]==true)
		dirty2[i]=true;
}*/
/*
int NV=(*dev_community).g.nb_nodes;
//NE=Gnew->numEdges;
double thres2=0.5;
for(int i=0;i<NV;i++)
        {
                
                if(dirty1[i+1]==true && bord[i+1])

                {
                        int nonborderedge=0;
                        int borderedge=bordno[i+1];                     
                        long adj1 = statIndices[i];
                        long adj2 = statIndices[i+1];
                        for(long j=adj1; j<adj2; j++) {
                                if(c[i]==c[edges[j]])
                                        nonborderedge++;
		                        }

		
		
                if(nonborderedge!=0 && (borderedge/nonborderedge)<thres2)
                        dirty2[i]=true;

                }
                else if(dirty1[i+1]==true && !bord[i+1])
                        dirty2[i]=true;
	}*/
unsigned int borderedge=0;
unsigned int nonborderedge=0;
double thres2=0.6;
for(long i=0;i<node;i++)
{
	if(dirty1[i] && G->bord[i+mid]){
        if(i==0)
        {
                p=0;
                q=statIndices[0];
        }
        else{
                p=statIndices[i-1];
                q=statIndices[i];
        }
	borderedge=G->bordno[i+mid];
        for(long j=p;j<q;j++)
	{
		if(c[i]==c[edges[j]])
			nonborderedge++;

	}
      if(nonborderedge!=0 &&(borderedge/nonborderedge)>thres2){
                        dirty2[i]=true;

                        }
 
	
}

}

cout<<"#####################similarity ends############################"<<endl;
//return dirty;

//return (* dev_community).g.nb_nodes ;
free(flag);
free(rel);
free(max);
free(int_deg);
return (* dev_community).g.nb_nodes ;

}
