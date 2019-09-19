

#include "defs.h"

using namespace std;
//WARNING: This will overwrite the original graph data structure to 
//         minimize memory footprint
// Return: C_orig will hold the cluster ids for vertices in the original graph
//         Assume C_orig is initialized appropriately
//WARNING: Graph G will be destroyed at the end of this routine
graph* runMultiPhaseLouvainAlgorithm(graph *G, unsigned int *C_orig, int coloring, unsigned int minGraphSize,
			double threshold, double C_threshold, int numThreads) 
{
  double totTimeClustering=0, totTimeBuildingPhase=0, totTimeColoring=0, tmpTime;
  int tmpItr=0, totItr = 0;  
  unsigned int NV = G->numVertices;

  int *colors;
  int numColors = 0;
  if(coloring==1) {
	colors = (int *) malloc (G->numVertices * sizeof(int)); assert (colors != 0);
#pragma omp parallel for
	for (unsigned int i=0; i<G->numVertices; i++) {
		colors[i] = -1;
	}
	numColors = algoDistanceOneVertexColoringOpt(G, colors, numThreads, &tmpTime)+1;
	totTimeColoring += tmpTime;
	//printf("Number of colors used: %d\n", numColors);
  }
cout<<"Inside multiphase run"<<endl;
//cout<<G->edgeListPtrs[NV]<<endl;	
  /* Step 3: Find communities */
  double prevMod = -1;
  double currMod = -1;
  unsigned int phase = 1;

  graph *Gnew,*GFinal; //To build new hierarchical graphs
  unsigned int numClusters;
  unsigned int *C = (unsigned int *) malloc (NV * sizeof(unsigned int));
  assert(C != 0);
#pragma omp parallel for
  for (unsigned int i=0; i<NV; i++) {
  	C[i] = -1;
  }
cout<<"inside 2"<<endl;	
  bool nonColor = false; //Make sure that at least one phase with lower threshold runs
  while(1){
        printf("===============================\n");
	printf("Phase %ld\n", phase);
        printf("===============================\n");
    	prevMod = currMod;
	//Compute clusters
	if((coloring == 1)&&(G->numVertices > minGraphSize)&&(nonColor == false)) {
		//Use higher modularity for the first few iterations when graph is big enough
		currMod = algoLouvainWithDistOneColoring(G, C, numThreads, colors, numColors, currMod, C_threshold, &tmpTime, &tmpItr);
		totTimeClustering += tmpTime;
                totItr += tmpItr;
	} else {
		currMod = parallelLouvianMethod(G, C, 6, currMod, threshold, &tmpTime, &tmpItr);
		totTimeClustering += tmpTime;
                totItr += tmpItr;
		nonColor = true; 
	}
  	//Renumber the clusters contiguiously
	numClusters = renumberClustersContiguously(C, G->numVertices);
	printf("Number of unique clusters: %ld\n", numClusters);
 	//printf("About to update C_orig\n");
	//Keep track of clusters in C_orig
	if(phase == 1) {
#pragma omp parallel for
	  for (unsigned int i=0; i<NV; i++) {
	  	C_orig[i] = C[i]; //After the first phase
	  }	
	} else {
#pragma omp parallel for
	  for (unsigned int i=0; i<NV; i++) {
                assert(C_orig[i] < G->numVertices);
                if (C_orig[i] >=0)
	  		C_orig[i] = C[C_orig[i]]; //Each cluster in a previous phase becomes a vertex
	  }	
	}
        printf("Done updating C_orig\n");
	//Break if too many phases or iterations
	if((phase > 200)||(totItr > 10000)) {
		break;
	}
        //Check for modularity gain and build the graph for next phase
	//In case coloring is used, make sure the non-coloring routine is run at least once
	if( (currMod - prevMod) > threshold ) {
		Gnew = (graph *) malloc (sizeof(graph)); assert(Gnew != 0);
				//cout <<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
		tmpTime = buildNextLevelGraphOpt(G, Gnew, C, numClusters, numThreads);
	//	cout<<"done"<<endl;
		totTimeBuildingPhase += tmpTime;
		//Free up the previous graph		
	//	cout<<"done2"<<endl;
		free(G->edgeListPtrs);	
		free(G->edgeList);
	//	free(G);
//		cout<<"done2"<<endl;

		G = Gnew; //Swap the pointers
		G->edgeListPtrs = Gnew->edgeListPtrs;
                G->edgeList = Gnew->edgeList;
		//Free up the previous cluster & create new one of a different size
		cout<<"done1"<<endl;
		unsigned int    *vtxPtr    = Gnew->edgeListPtrs;
		edge    *vtxInd    = Gnew->edgeList;
	/*	for (unsigned int i = 0; i <(Gnew->numVertices) ; i++) {
    			unsigned int adj1 = vtxPtr[i];
    			unsigned int adj2 = vtxPtr[i+1];
			//cout <<i<< " ";
   			// printf("\nVtx: %ld [%ld]: ",i+1,adj2-adj1);
    			for(unsigned int j=adj1; j<adj2; j++) { 
				cout<<(vtxInd[j].tail+1)-1<<" ";
							}
				cout <<endl;
						}

	*/	
		

		free(C);
		C = (unsigned int *) malloc (numClusters * sizeof(unsigned int)); assert(C != 0);
#pragma omp parallel for
		for (unsigned int i=0; i<numClusters; i++) {
			C[i] = -1;
		}
		phase++; //Increment phase number
		//If coloring is enabled & graph is of minimum size, recolor the new graph
		if((coloring == 1)&&(G->numVertices > minGraphSize)&&(nonColor = false)){
#pragma omp parallel for
			for (unsigned int i=0; i<G->numVertices; i++){
				colors[i] = -1;
			}
			numColors = algoDistanceOneVertexColoringOpt(G, colors, numThreads, &tmpTime)+1;
			totTimeColoring += tmpTime;
		}
	} else {
		if ( (coloring == 1)&&(nonColor == false) ) {
			nonColor = true; //Run at least one loop of non-coloring routine 
		}
		else {
			break; //Modularity gain is not enough. Exit.
		}
	}
		//GFinal = (graph *) malloc (sizeof(graph)); assert(GFinal != 0);
		//	GFinal=Gnew;

  } //End of while(1)
 cout <<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<Gnew->numVertices<<" "<<Gnew->numEdges<<endl;	
//  printf("********************************************\n"); 
 // printf("*********    Compact Summary   *************\n");
 // printf("********************************************\n");
 // printf("Number of threads              : %ld\n", numThreads);
 // printf("Total number of phases         : %ld\n", phase);
 // printf("Total number of iterations     : %ld\n", totItr);
  printf("Final number of clusters       : %ld\n", numClusters);
  printf("Final modularity               : %lf\n", prevMod);
 // printf("Total time for clustering      : %lf\n", totTimeClustering);
 // printf("Total time for building phases : %lf\n", totTimeBuildingPhase);
  if(coloring == 1) {
     printf("Total time for coloring        : %lf\n", totTimeColoring);
  }
//  printf("********************************************\n");
  printf("TOTAL TIME                     : %lf\n", (totTimeClustering+totTimeBuildingPhase+totTimeColoring) );
//  printf("********************************************\n");

  //Clean up:
  free(C);
/*  if(G != 0) {
    free(G->edgeListPtrs);
    free(G->edgeList);
    free(G);}*/
return Gnew;
  

  if(coloring==1) {
    if(colors != 0) free(colors);
  }
//return Gnew;
}//End of runMultiPhaseLouvainAlgorithm()
