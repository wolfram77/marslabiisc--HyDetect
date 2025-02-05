
#include <bits/stdc++.h>
//#include"getNextLineAndSplitIntoTokens.cpp"
#ifndef _DEFS_H
#define _DEFS_H
#include"graphHOST.h"
#include"communityGPU.h"
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <unistd.h> //For getopts()
#include"vector"
#include"assert.h"
#include"commonconstants.h"
#include "graphGPU.h"
#include "cuda.h"
#include"cuda_runtime_api.h"
#include "graphHOST.h"
#include"commonconstants.h"
#include"hostconstants.h"
#include"string"
#define MilanRealMax HUGE_VAL       // +INFINITY
#define MilanRealMin -MilanRealMax  // -INFINITY
#define PRINT_DETAILED_STATS_
 using namespace std; 
double get_time();



typedef struct move_gpu
{
unsigned int id;
int statindices;
int edge;
}move_gpu;

typedef struct
{
unsigned int vertex;
unsigned int edg;
unsigned int *statIndices;
unsigned int *edgesa;
unsigned int*weighta; 
} move1;

typedef struct comm
{
  unsigned int size;
  unsigned int degree;
}Comm;

typedef struct /* the edge data structure */
{
  unsigned int head;
  unsigned int tail;
  //unsigned int weight;
  unsigned weight;
} edge;

typedef struct
{
unsigned int vertex;
unsigned int edgeno;
unsigned int *edgeListPtrsM;
edge *edgesM;
} move2;

typedef struct /* the graph data structure */
{
  unsigned int numVertices;        /* Number of columns                                */
  unsigned int sVertices;          /* Number of rows: Bipartite graph: number of S vertices; T = N - S */
  unsigned int numEdges;           /* Each edge stored twice, but counted once        */
  unsigned int * edgeListPtrs;     /* start vertex of edge, sorted, primary key        */
  unsigned int *edgeListPtrs1;
  edge * edgeList;edge *edgeList1;         /* end   vertex of edge, sorted, secondary key      */
  unsigned read(char file[]);
  unsigned allocOnHost();
  unsigned deallocOnHost();
  unsigned readFromGR(char file[]);
  unsigned partGraph(int, double);
  void getPartitionedData( int, int *);
  void findSubGraph0(unsigned,unsigned,unsigned);
  void findSubGraph1(unsigned,unsigned,unsigned);
  //for gpugraph
  vector<unsigned int> degree;
  vector<unsigned int >  links;
	vector<unsigned int> weights;
  //
  unsigned nnodes, nedges, maxDegNode, secondMaxDegNode;
  unsigned *noutgoing, *nincoming, *srcsrc, *psrc, *edgessrcdst;
  unsigned *edgessrcwt;
  unsigned *maxDeg;
  double dataStructureSpace;

  //for part Graph
  unsigned middleVertex;
	unsigned no1,no2;
  unsigned *PerPartEdges, *PerPartNodes, *startVertex, *endVertex, *startEdge, *endEdge;
  bool *bord;
  unsigned int *bordno;
  vector<unsigned int>* bordvalue;
//  unsigned findSubGraph();

  //for subgarph
  unsigned int *subGraphPsrc, *subGraphOutgoing, *subGraphEdgesDst, *subGraphEdgesWt;
  unsigned int subGraphNodes, subGraphEdges;
  unsigned findSubGraph();
  unsigned *verPart, *borderVer;
  unsigned **partVer, *partNumVer;

  unsigned part0Edges, part1Edges;
  unsigned *part0Psrc, *part0Outgoing, *part0EdgesDst;
  unsigned *part1Psrc, *part1Outgoing, *part1EdgesDst;

  unsigned *partNumEdges;
  unsigned **partPsrc, **partOutgoing, **partEdgesSrcDst;
} graph;

struct clustering_parameters 
{
  const char *inFile; //Input file
  int ftype;  //File type

  bool strongScaling; //Enable strong scaling
  bool output; //Printout the clustering data
  bool VF; //Vertex following turned on
  bool coloring; //If coloring is turned on 

  double C_thresh; //Threshold with coloring on
  unsigned int minGraphSize; //Min |V| to enable coloring
  double threshold; //Value of threshold
       
  clustering_parameters();
  void usage();    
  bool parse(int argc, char *argv[]);
};

//vector<string> getNextLineAndSplitIntoTokens(istream& str);

/////////////////// FUNCTION CALLS ////////////////////
//ctor<string> getNextLineAndSplitIntoTokens(istream& str);

void displayGraphCharacteristics(graph *G);
void displayGraph(graph *G);
void displayGraphEdgeList(graph *G);
void displayGraphEdgeList(graph *G, FILE* out);
//Graph Clustering (Community detection 
int  gpuonly(GraphHOST input_graph,long *c,unsigned int *statIndices,unsigned int *edges,Community *dev_community,bool *dirty,bool *dirty1,int b,graph *G,unsigned int mid);
graph* cpuonly(graph *G, graph *Gnew,graph* G1,long *c,bool *d,bool *d1,int b);
void verticesToMoveToGPU(graph *G,bool *dirtycpu,move1 *mo,long *c,unsigned int mid,unsigned int mid1,long *c1);
void verticesToMoveToCPU(unsigned int*statIndices,unsigned int* edges,bool * dirtygpu,move2* mo1,long *C_orig,int total,int NV,graph *G1,graph *G2,unsigned int mid);
graph * modifyCPUstructure(graph *Gnew,graph *G1,bool *dirty,long *C_orig,long *C11,move2 *mo1,unsigned int mid,unsigned int node,bool *d,bool *dc);
Community modifyGPUstructure(Community *dev1_community,unsigned int *statIndices,unsigned int*edges,bool *dirtyg,bool *dirtygpu,long  *c,int total,int NV,move1 *mo,graph *Gnew,unsigned int mid, long *c1,int it,long *cact);
void verticesToMoveToGPU1(graph *G,graph * gn,bool *dirtycpu,move1 *mo,long *c,unsigned int mid,unsigned int node,int n,long *cc);

void movetogpu(graph *Gnew,graph *G1,bool *d,bool* dirtycpu,int *statIndices,unsigned int *edges,Community* dev_community,int *c,int f);
//void movefinal(graph *Gnew,unsigned int *C_orig,graph*G1, unsigned int *statIndices,unsigned int *edges,Community *dev_community, int *cg,bool *c);
//void movefinal(graph *Gnew,graph *Gnew1,bool *dirty,bool *dirtycpu,move1 *mo,unsigned int *c,unsigned int mid,unsigned int node);
void movefinal(graph *Gnew,graph *Gnew1,bool *dirty1,bool *dirtycpu,move1 *mo,unsigned int *c,unsigned int mid,unsigned int node,bool *bord,unsigned int* bordno,vector<unsigned int>* bordervalue);
graph* movetocpu(bool *dirtygpu,int actualnode,unsigned int *statIndices,unsigned int *edges,graph *Gnew,clustering_parameters opts,graph *G1,unsigned int *c,bool *d);

double parallelLouvianMethod(graph *G, unsigned int *C, int nThreads, double Lower, 
				double thresh, double *totTime, int *numItr);
double algoLouvainWithDistOneColoring(graph* G, unsigned int *C, int nThreads, int* color, 
			int numColor, double Lower, double thresh, double *totTime, int *numItr);
graph* runMultiPhaseLouvainAlgorithm(graph *G, long *C_orig, int coloring, unsigned int minGraphSize, 
			double threshold, double C_threshold, int numThreads);


//***  Clustering Utility Functions ***//
//Distance-1 Coloring
int algoDistanceOneVertexColoring(graph *G, int *vtxColor, int nThreads, double *totTime);
int algoDistanceOneVertexColoringOpt(graph *G, int *vtxColor, int nThreads, double *totTime);

//Other 
inline void Visit(unsigned int v, unsigned int myCommunity, short *Visited, unsigned int *Volts, 
				  unsigned int* vtxPtr, edge* vtxInd, unsigned int *C);
unsigned int buildCommunityBasedOnVoltages(graph *G, unsigned int *Volts, unsigned int *C, unsigned int *Cvolts);
void buildNextLevelGraph(graph *Gin, graph *Gout, unsigned int *C, unsigned int numUniqueClusters);
unsigned int renumberClustersContiguously(unsigned int *C, unsigned int size);
double buildNextLevelGraphOpt(graph *Gin, graph *Gout, unsigned int *C, unsigned int numUniqueClusters, int nThreads);
//Vertex following functions:
unsigned int vertexFollowing(graph *G, unsigned int *C);
double buildNewGraphVF(graph *Gin, graph *Gout, unsigned int *C, unsigned int numUniqueClusters);

//***  Utility Functions ***//
void duplicateGivenGraph(graph *Gin, graph *Gout);
void writeEdgeListToFile(graph *G, FILE* out);

//Random Number Generation:
void generateRandomNumbers(double *RandVec, unsigned int size);

void displayGraph(graph *G);
void displayGraphCharacteristics(graph *G);
graph * convertDirected2Undirected(graph *G);

void segregateEdgesBasedOnVoltages(graph *G, unsigned int *Volts);
void writeGraphPajekFormat(graph* G, char * filename);
void writeGraphPajekFormatWithNodeVolts(graph* G, unsigned int *Cvolts, char * filename);
void writeGraphBinaryFormat(graph* G, char * filename); //Binary (each edge once)
void writeGraphMetisSimpleFormat(graph* G, char *filename); //Metis format; no weights

//File parsers:
void parse_Dimacs9FormatDirectedNewD(graph * G, char *fileName);
unsigned int removeEdges(unsigned int NV, unsigned int NE, edge *edgeList);
void SortNodeEdgesByIndex2(unsigned int NV, edge *list1, edge *list2, unsigned int *ptrs);
void SortEdgesUndirected2(unsigned int NV, unsigned int NE, edge *list1, edge *list2, unsigned int *ptrs);

void loadMetisFileFormat(graph *G, const char* filename); //Metis (DIMACS#10)
void parse_MatrixMarket(graph * G, char *fileName);       //Matrix-Market
void parse_MatrixMarket_Sym_AsGraph(graph * G, char *fileName);

void parse_Dimacs1Format(graph * G, char *fileName);      //DIMACS#1 Challenge format
void parse_Dimacs9FormatDirectedNewD(graph * G, char *fileName); //DIMACS#9 Challenge format
void parse_PajekFormat(graph * G, char *fileName);        //Pajek format (each edge stored only once
void parse_PajekFormatUndirected(graph * G, char *fileName);
void parse_DoulbedEdgeList(graph * G, char *fileName);

void parse_EdgeListBinary(graph * G, char *fileName); //Binary: Each edge stored only once 
void parse_SNAP(graph * G, char *fileName);

//For reading power grid data
unsigned int* parse_MultiKvPowerGridGraph(graph * G, char *fileName); //Four-column format


//Graph partitioning with Metis:
void MetisGraphPartitioner( graph *G, unsigned int *VertexPartitioning, int numParts );

#endif

