

#include "defs.h"
#include "utilityClusteringFunctions.h"

using namespace std;

double parallelLouvianMethod(graph *G, unsigned int *C, int nThreads, double Lower, 
				double thresh, double *totTime, int *numItr) {
#ifdef PRINT_DETAILED_STATS_  
  printf("Within parallelLouvianMethod()\n");
#endif
  
    omp_set_num_threads(6);
  int nT;
#pragma omp parallel
  {
    nT = 6;
  }
#ifdef PRINT_DETAILED_STATS_
  printf("Actual number of threads: %d (requested: %d)\n", nT, nThreads);
#endif
  double time1, time2, time3, time4; //For timing purposes  
  double total = 0, totItr = 0;
  
  unsigned int    NV        = G->numVertices;
  unsigned int    NS        = G->sVertices;      
  unsigned int    NE        = G->numEdges;
  unsigned int    *vtxPtr   = G->edgeListPtrs;
  edge    *vtxInd   = G->edgeList;
 
  /* Variables for computing modularity */
  unsigned int totalEdgeWeightTwice;
  double constantForSecondTerm;
  double prevMod=-1;
  double currMod=-1;
  //double thresMod = 0.000001;
  double thresMod = thresh; //Input parameter
  int numItrs = 0;
  
  /********************** Initialization **************************/
  time1 = omp_get_wtime();
  //Store the degree of all vertices
  unsigned int* vDegree = (unsigned int *) malloc (NV * sizeof(unsigned int)); assert(vDegree != 0);
  //Community info. (ai and size)
  Comm *cInfo = (Comm *) malloc (NV * sizeof(Comm)); assert(cInfo != 0);
  //use for updating Community
  Comm *cUpdate = (Comm*)malloc(NV*sizeof(Comm)); assert(cUpdate != 0);
  //use for Modularity calculation (eii)
  unsigned int* clusterWeightInternal = (unsigned int*) malloc (NV*sizeof(unsigned int)); assert(clusterWeightInternal != 0);

  sumVertexDegree(vtxInd, vtxPtr, vDegree, NV , cInfo);	// Sum up the vertex degree
  
  /*** Compute the total edge weight (2m) and 1/2m ***/
  constantForSecondTerm = calConstantForSecondTerm(vDegree, NV); // 1 over sum of the degree
  
  //Community assignments:
  //Store previous iteration's community assignment
  unsigned int* pastCommAss = (unsigned int *) malloc (NV * sizeof(unsigned int)); assert(pastCommAss != 0);
  //Store current community assignment
  unsigned int* currCommAss = (unsigned int *) malloc (NV * sizeof(unsigned int)); assert(currCommAss != 0);  
  //Store the target of community assignment  
  unsigned int* targetCommAss = (unsigned int *) malloc (NV * sizeof(unsigned int)); assert(targetCommAss != 0);
// cout<<"ok1"<<endl; 
  //Initialize each vertex to its own cluster
  initCommAss(pastCommAss, currCommAss, NV); 
  
  time2 = omp_get_wtime();
  printf("Time to initialize: %3.3lf\n", time2-time1);
	
#ifdef PRINT_DETAILED_STATS_
  printf("========================================================================================================\n");
  printf("Itr      E_xx            A_x2           Curr-Mod         Time-1(s)       Time-2(s)        T/Itr(s)\n");
  printf("========================================================================================================\n");
#endif
#ifdef PRINT_TERSE_STATS_
  printf("=====================================================\n");
  printf("Itr      Curr-Mod         T/Itr(s)      T-Cumulative\n");
  printf("=====================================================\n");
#endif
  //Start maximizing modularity
  while(true) {
    numItrs++;    
    time1 = omp_get_wtime();
    /* Re-initialize datastructures */
#pragma omp parallel for
    for (unsigned int i=0; i<NV; i++) {
      clusterWeightInternal[i] = 0; 
      cUpdate[i].degree =0;
      cUpdate[i].size =0;
    }
//cout<<G->edgeListPtrs[12]<<endl;
//ut<<vtxPtr[NV-1]<<endl;
//cout<<"doubt clarified"<<endl;
// cout<<NV<<endl;   
#pragma omp parallel for
    for (unsigned int i=0; i<NV; i++) {
      unsigned int adj1 = vtxPtr[i];
      unsigned int adj2 = vtxPtr[i+1];
//	cout<<adj1<<" "<<adj2<<endl;
	  unsigned int selfLoop = 0;
	  //Build a datastructure to hold the cluster structure of its neighbors      	
	  map<unsigned int, unsigned int> clusterLocalMap; //Map each neighbor's cluster to a local number
	  map<unsigned int, unsigned int>::iterator storedAlready;
	  vector<double> Counter; //Number of edges in each unique cluster
	  //Add v's current cluster:
	  if(adj1 != adj2){	
              clusterLocalMap[currCommAss[i]] = 0;
	      Counter.push_back(0); //Initialize the counter to ZERO (no edges incident yet)
	      //Find unique cluster ids and #of edges incident (eicj) to them
	      selfLoop = buildLocalMapCounter(adj1, adj2, clusterLocalMap, Counter, vtxInd, currCommAss, i);
	      // Update delta Q calculation
	      clusterWeightInternal[i] += (unsigned int)Counter[0]; //(e_ix)
	      //Calculate the max
	      targetCommAss[i] = max(clusterLocalMap, Counter, selfLoop, cInfo, vDegree[i], currCommAss[i], constantForSecondTerm);
              //assert((targetCommAss[i] >= 0)&&(targetCommAss[i] < NV));
      	  } else {
		  targetCommAss[i] = -1;	
          }          

         //Update
         if(targetCommAss[i] != currCommAss[i]  && targetCommAss[i] != -1) {
	  __sync_fetch_and_add(&cUpdate[targetCommAss[i]].degree, vDegree[i]);
	  __sync_fetch_and_add(&cUpdate[targetCommAss[i]].size, 1);
	  __sync_fetch_and_sub(&cUpdate[currCommAss[i]].degree, vDegree[i]);
	  __sync_fetch_and_sub(&cUpdate[currCommAss[i]].size, 1);
        }//End of If()      
        clusterLocalMap.clear();      
        Counter.clear();
    }//End of for(i)
    time2 = omp_get_wtime();
 
    time3 = omp_get_wtime();    
    double e_xx = 0;
    double a2_x = 0;	

#pragma omp parallel for \
  reduction(+:e_xx) reduction(+:a2_x)
    for (unsigned int i=0; i<NV; i++) {
      e_xx += clusterWeightInternal[i];
      a2_x += (cInfo[i].degree)*(cInfo[i].degree);
    }
    time4 = omp_get_wtime();

    currMod = (e_xx*(double)constantForSecondTerm) - (a2_x*(double)constantForSecondTerm*(double)constantForSecondTerm);
    totItr = (time2-time1) + (time4-time3);
    total += totItr;
#ifdef PRINT_DETAILED_STATS_
    printf("%d \t %g \t %g \t %lf \t %3.3lf \t %3.3lf  \t %3.3lf\n",numItrs, e_xx, a2_x, currMod, (time2-time1), (time4-time3), totItr );
#endif
#ifdef PRINT_TERSE_STATS_
   printf("%d \t %lf \t %3.3lf  \t %3.3lf\n",numItrs, currMod, totItr, total);
#endif
 
    //Break if modularity gain is not sufficient
    if((currMod - prevMod) < thresMod) {
      break;
    }
    
    //Else update information for the next iteration
    prevMod = currMod;
    if(prevMod < Lower)
	prevMod = Lower;
#pragma omp parallel for 
    for (unsigned int i=0; i<NV; i++) {
      cInfo[i].size += cUpdate[i].size;
      cInfo[i].degree += cUpdate[i].degree;
    }
    
    //Do pointer swaps to reuse memory:
    unsigned int* tmp;
    tmp = pastCommAss;
    pastCommAss = currCommAss; //Previous holds the current
    currCommAss = targetCommAss; //Current holds the chosen assignment
    targetCommAss = tmp;      //Reuse the vector
    
  }//End of while(true)
  *totTime = total; //Return back the total time for clustering
  *numItr  = numItrs;

#ifdef PRINT_DETAILED_STATS_
  printf("========================================================================================================\n");
  printf("Total time for %d iterations is: %lf\n",numItrs, total);  
  printf("========================================================================================================\n");
#endif  
#ifdef PRINT_TERSE_STATS_
  printf("========================================================================================================\n");
  printf("Total time for %d iterations is: %lf\n",numItrs, total);  
  printf("========================================================================================================\n");
#endif

  //Store back the community assignments in the input variable:
  //Note: No matter when the while loop exits, we are interested in the previous assignment
#pragma omp parallel for 
  for (unsigned int i=0; i<NV; i++) {
    C[i] = pastCommAss[i];
  }
  //Cleanup
  free(pastCommAss);
  free(currCommAss);
  free(targetCommAss);
  free(vDegree);
  free(cInfo);
  free(cUpdate);
  free(clusterWeightInternal);

  return prevMod;
}
