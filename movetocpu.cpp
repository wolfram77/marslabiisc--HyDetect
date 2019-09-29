#include <bits/stdc++.h>

#include"defs.h"
using namespace std;
vector<string> getNextLineAndSplitIntoTokens1(istream& str)
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

void addedge(graph *Gnew,long c,long d,int ch)
{
long    *vtxPtr1   = Gnew->edgeListPtrs;
edge    *vtxInd1   = Gnew->edgeList;

for (long i = 0; i < Gnew->numVertices; i++) {
	if(i==c){
    long adj1 = vtxPtr1[i];
    long adj2 = vtxPtr1[i+1];
	
long f=(vtxInd1[adj2-adj1+1].tail+1);
    f=d;

		}		
						}
if(ch==0)
Gnew->numVertices=Gnew->numVertices+1;
Gnew->numEdges=Gnew->numEdges+1;
}

graph* movetocpu(bool *dirtygpu,int actualnode,unsigned long *statIndices,unsigned int *edges,graph *Gnew,clustering_parameters opts,			graph *G1,long *C_cpu,bool *borcheck){
cout<<"inside movement to cpu"<<endl;
int countgc=0;
int x=0;             
long NV=G1->numVertices;
long* vtxPtr=G1->edgeListPtrs;
edge* vtxInd=G1->edgeList;
int *pos1=(int *)malloc(sizeof(int)*actualnode);
int j=1;
for(int i=0;i<actualnode;i++)
{
	if(dirtygpu[i]==true){
		pos1[i]=j;
		countgc++;
		j++;
		}
}
string line;
int in=0;
bool *bord=(bool *)malloc(sizeof(bool)*NV);
long *bordno=(long *)malloc(sizeof(long)*NV);
bool *flaz=(bool *)malloc(sizeof(bool)*NV);
long *array=(long *)malloc(sizeof(long)*NV*2);
int qq=0;
 vector<long> v;
vector<long> *borderval=new vector<long>[NV+1];

for(long i=0;i<NV;i++)  {
        flaz[i]=false;
        bord[i]=false;
        bordno[i]=0;
                        }
 for(long i=0;i<(2*NV);i++)
        array[i]=0;

ifstream border("/home/Anwehsa/uk_border2");
while(getline(border,line)){
        in++;
}

//cout<<"hello"<<endl;
int bordere=0;
ifstream border22("/home/Anwehsa/uk_border2");
int iij=0;
while(iij<in){

        vector<string> str;
        str = getNextLineAndSplitIntoTokens1(border22);

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


//cout<<"hello"<<endl;

int adj1,adj2;
int eecount=0;
/*for(long i=0;i<actualnode;i++)
{
	int pos=0;
		if(dirtygpu[i]==true)
		{
			if(i==0)
			{	
			adj1=0;
			adj2=statIndices[i+1];
			}	
			else{
			adj1 = statIndices[i];
			adj2 = statIndices[i+1];
			}
	
			for(int j=adj1;j<adj2;j++)
			{
//vector<unsigned int> mm;
//mm.push_back(edges[j]);
				eecount++;
				addedge(Gnew,cg[edges[j]],(Gnew->numVertices+pos));
			}
		}

	pos=pos+1;
}
*/

cout<<Gnew->numVertices<<" "<<Gnew->numEdges<<endl;

for(long i=0;i<actualnode;i++)
{
	int pos=0;
	if(dirtygpu[i])
	{
		if(i==0)
		{
		adj1=0;
		adj2=statIndices[i+1];
		}
		else{
		adj1 = statIndices[i];
                adj2 = statIndices[i+1];
			}
	//case where the doubtful vertex has connection with the border
	if(bord[i+1])
	{	int ch=0;
		for(std::vector<long> ::iterator it=borderval[i+1].begin();it!=borderval[i+1].end();it++){	
			addedge(Gnew,C_cpu[*it],Gnew->numVertices+pos,ch);	
			ch=1;						}


	}
	int ch1=0;
	//case where there is connection between doubtful vertices
	for(int j=adj1;j<adj2;j++)
	{	
		if(dirtygpu[edges[j]] && bord[i+1])
			addedge(Gnew,pos1[edges[j]],Gnew->numVertices+pos,1);
		else if(dirtygpu[edges[j]] && !bord[i+1])
			 addedge(Gnew,pos1[edges[j]],Gnew->numVertices+pos,ch1);
		ch1=1;
	}	


}
pos=pos+1;

}


graph *gf=(graph *)malloc(sizeof(graph));
//ng newV=Gnew->numVertices+countgc;
long * C_orig1=(long*)malloc((Gnew->numVertices)*sizeof(long));
//gf->numVertices=newV;
//long newE=Gnew->numEdges+eecount;
//->numEdges=newE;
cout<<Gnew->numVertices<<" "<<Gnew->numEdges<<endl;
graph* Gnew1 = (graph *) malloc (sizeof(graph));
for(long i=0;i<Gnew->numVertices;i++)
	C_orig1[i]=-1;
  displayGraphCharacteristics(Gnew);

Gnew1=runMultiPhaseLouvainAlgorithm(Gnew, C_orig1,0, opts.minGraphSize, opts.threshold, opts.C_thresh, 5); 
ifstream border3("/home/Anwehsa/uk_border2");
 in=0;
while(getline(border3,line)){
        in++;
}
//forming new border set 
borcheck=(bool*)malloc(sizeof(bool)*(actualnode+1));
ifstream bord32("/home/Anwehsa/uk_border2");
for(int i=0;i<actualnode;i++)
        borcheck[i]=false;
while(iij<in){
        vector<string> str;
        str = getNextLineAndSplitIntoTokens1(bord32);

        int n = str.size();


        bordere=0;

        stringstream geek(str[0]);

        geek>>x;
        borcheck[x]=true;
        }

for(int i=0;i<actualnode;i++)
{
	if(borcheck[i+1] && dirtygpu[i])
	{
		if(i==0)
                {
                adj1=0;
                adj2=statIndices[i+1];
                }
                else{
                adj1 = statIndices[i];
                adj2 = statIndices[i+1];
                        }
		for(int j=adj1;j<adj2;j++)
		{
			borcheck[edges[j]+1]=true;
		}
	borcheck[i+1]=false;
	}

}
return Gnew1;

}

