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

#include "defs.h"

using namespace std;

unsigned int vertexFollowing(graph *G, unsigned int *C)
{
	unsigned int    NV        = G->numVertices;
	unsigned int    *vtxPtr   = G->edgeListPtrs;
	edge    *vtxInd   = G->edgeList;
	unsigned int numNode = 0;
	double time1 = omp_get_wtime();
// Initialize the Communities
#pragma omp parallel for  //Parallelize on the outer most loop
  	for (unsigned int i=0; i<NV; i++) {
		C[i] = i; //Initialize each vertex to its own cluster
	}

// Remove Isolated and degree-one vertices
#pragma omp parallel for
    	for (unsigned int i=0; i<NV; i++) {
		unsigned int adj1 = vtxPtr[i];
		unsigned int adj2 = vtxPtr[i+1];
		if(adj1 == adj2) {	// Isolated vertex
			__sync_fetch_and_add(&numNode, 1);
			C[i] = -1;
		} else {
			if( (adj2-adj1) == 1 ) { //Degree one
			    //Check if the tail has degree greater than one:
			    unsigned int tail = vtxInd[adj1].tail;
			    unsigned int adj11 = vtxPtr[tail];
			    unsigned int adj12 = vtxPtr[tail+1];
                            if( ((adj12-adj11) > 1)||(i > tail) ) { //Degree of tail greater than one
				__sync_fetch_and_add(&numNode, 1);
				C[i] = tail;
			    } //else don't do anything
			}//End of if(degree one)
		}//End of else
	}//End of for(i)

        time1 = omp_get_wtime() - time1;
#ifdef PRINT_DETAILED_STATS_
        printf("Time to determine number of vertices (numNode) to fix: %lf\n", time1);	
#endif
	return numNode; //These are nodes that need to be removed
}//End of vertexFollowing()

//WARNING: Will assume that the cluster id have been renumbered contiguously
//Return the total time for building the next level of graph
//This will not add any self-loops
double buildNewGraphVF(graph *Gin, graph *Gout, unsigned int *C, unsigned int numUniqueClusters) {
  int nT;
#pragma omp parallel
  {
    nT = omp_get_num_threads();
  }
#ifdef PRINT_DETAILED_STATS_
  printf("Within buildNewGraphVF(): # of unique clusters= %ld\n",numUniqueClusters);
  printf("Actual number of threads: %d \n", nT);
#endif

  double time1, time2, TotTime=0; //For timing purposes  
  double total = 0, totItr = 0;  
  //Pointers into the input graph structure:
  unsigned int    NV_in        = Gin->numVertices;
  unsigned int    NE_in        = Gin->numEdges;
  unsigned int    *vtxPtrIn    = Gin->edgeListPtrs;
  edge    *vtxIndIn    = Gin->edgeList;
  
  time1 = omp_get_wtime();
  // Pointers into the output graph structure
  unsigned int NV_out = numUniqueClusters;
  unsigned int NE_self = 0; //Not all vertices get self-loops
  unsigned int NE_out = 0;  //Cross edges
  unsigned int *vtxPtrOut = (unsigned int *) malloc ((NV_out+1)*sizeof(unsigned int));
  assert(vtxPtrOut != 0);
  vtxPtrOut[0] = 0; //First location is always a zero
  /* Step 1 : Regroup the node into cluster node */
  map<unsigned int,unsigned int>** cluPtrIn = (map<unsigned int,unsigned int>**) malloc(numUniqueClusters*sizeof(map<unsigned int,unsigned int>*));
  assert(cluPtrIn != 0);	

#pragma omp parallel for
  for (unsigned int i=0; i<numUniqueClusters; i++) {
	cluPtrIn[i] = new map<unsigned int,unsigned int>();
	//Do not add self-loops
        //(*(cluPtrIn[i]))[i] = 0; //Add for a self loop with zero weight
  }
#pragma omp parallel for
  for (unsigned int i=1; i<=NV_out; i++)
	vtxPtrOut[i] = 0; 

  //Create an array of locks for each cluster
  omp_lock_t *nlocks = (omp_lock_t *) malloc (numUniqueClusters * sizeof(omp_lock_t));
  assert(nlocks != 0);
#pragma omp parallel for
  for (unsigned int i=0; i<numUniqueClusters; i++) {
    omp_init_lock(&nlocks[i]); //Initialize locks
  }
  time2 = omp_get_wtime();
  TotTime += (time2-time1);
#ifdef PRINT_DETAILED_STATS_
  printf("Time to initialize: %3.3lf\n", time2-time1);
#endif
  time1 = omp_get_wtime();
#pragma omp parallel for
  for (unsigned int i=0; i<NV_in; i++) {
	if((C[i] < 0)||(C[i]>numUniqueClusters))
		continue; //Not a valid cluster id
	unsigned int adj1 = vtxPtrIn[i];
	unsigned int adj2 = vtxPtrIn[i+1];
	map<unsigned int, unsigned int>::iterator localIterator;
        assert(C[i] < numUniqueClusters);
  	//Now look for all the neighbors of this cluster
	for(unsigned int j=adj1; j<adj2; j++) {
		unsigned int tail = vtxIndIn[j].tail; 
		assert(C[tail] < numUniqueClusters);			
		//Add the edge from one endpoint	
		if(C[i] >= C[tail]) {
                        omp_set_lock(&nlocks[C[i]]);  // Locking the cluster
	
			localIterator = cluPtrIn[C[i]]->find(C[tail]); //Check if it exists			
			if( localIterator != cluPtrIn[C[i]]->end() ) {	//Already exists
                                localIterator->second += (unsigned int)vtxIndIn[j].weight;
			} else {
				(*(cluPtrIn[C[i]]))[C[tail]] = (unsigned int)vtxIndIn[j].weight; //Self-edge
				__sync_fetch_and_add(&vtxPtrOut[C[i]+1], 1);
				if(C[i] == C[tail])
                                	__sync_fetch_and_add(&NE_self, 1); //Keep track of self #edges 
				if(C[i] > C[tail]) {
					__sync_fetch_and_add(&NE_out, 1); //Keep track of non-self #edges
					__sync_fetch_and_add(&vtxPtrOut[C[tail]+1], 1); //Count edge j-->i
				}
			}
                        
                        omp_unset_lock(&nlocks[C[i]]); // Unlocking the cluster
		} //End of if
	}//End of for(j)
  }//End of for(i)  
  //Prefix sum:
  for(unsigned int i=0; i<NV_out; i++) {
	vtxPtrOut[i+1] += vtxPtrOut[i];
  }
  
  time2 = omp_get_wtime();
  TotTime += (time2-time1);
  printf("NE_out= %ld   NE_self= %ld\n", NE_out, NE_self);
  printf("These should match: %ld == %ld\n",(2*NE_out + NE_self), vtxPtrOut[NV_out]);
#ifdef PRINT_DETAILED_STATS_
  printf("Time to count edges: %3.3lf\n", time2-time1);
#endif
  assert(vtxPtrOut[NV_out] == (NE_out*2+NE_self)); //Sanity check
  
  time1 = omp_get_wtime();
  // Step 3 : build the edge list:
  unsigned int numEdges   = vtxPtrOut[NV_out];
  unsigned int realEdges  = NE_out + NE_self; //Self-loops appear once, others appear twice
  edge *vtxIndOut = (edge *) malloc (numEdges * sizeof(edge));
  assert (vtxIndOut != 0);
  unsigned int *Added = (unsigned int *) malloc (NV_out * sizeof(unsigned int)); //Keep track of what got added
  assert (Added != 0);

#pragma omp parallel for
  for (unsigned int i=0; i<NV_out; i++) {
	Added[i] = 0;
  }  
  //Now add the edges in no particular order
#pragma omp parallel for
  for (unsigned int i=0; i<NV_out; i++) {
	unsigned int Where;
	map<unsigned int, unsigned int>::iterator localIterator = cluPtrIn[i]->begin();
	//Now go through the other edges:
	while ( localIterator != cluPtrIn[i]->end()) {
		Where = vtxPtrOut[i] + __sync_fetch_and_add(&Added[i], 1);
		vtxIndOut[Where].head = i; //Head
		vtxIndOut[Where].tail = localIterator->first; //Tail
		vtxIndOut[Where].weight = localIterator->second; //Weight
		if(i != localIterator->first) {		
			Where = vtxPtrOut[localIterator->first] + __sync_fetch_and_add(&Added[localIterator->first], 1);
			vtxIndOut[Where].head = localIterator->first;
			vtxIndOut[Where].tail = i; //Tail
			vtxIndOut[Where].weight = localIterator->second; //Weight
			//printf("%d\n",localIterator->first);
		}
		localIterator++;
	}	
  }//End of for(i)  
  time2 = omp_get_wtime();
  TotTime += (time2-time1);
#ifdef PRINT_DETAILED_STATS_
  printf("Time to build the graph: %3.3lf\n", time2-time1);
  printf("Total time: %3.3lf\n", TotTime);
#endif
#ifdef PRINT_TERSE_STATS_
  printf("Total time to build next phase: %3.3lf\n", TotTime);
#endif
  // Set the pointers
  Gout->numVertices  = NV_out;
  Gout->sVertices    = NV_out;
  //Note: Self-loops are represented ONCE, but others appear TWICE
  Gout->numEdges     = realEdges; //Add self loops to the #edges
  Gout->edgeListPtrs = vtxPtrOut;
  Gout->edgeList     = vtxIndOut;
	
  //Clean up
  free(Added);
#pragma omp parallel for
  for (unsigned int i=0; i<numUniqueClusters; i++)
	delete cluPtrIn[i];	
  free(cluPtrIn);

#pragma omp parallel for
  for (unsigned int i=0; i<numUniqueClusters; i++) {
    omp_destroy_lock(&nlocks[i]); 
  }
  free(nlocks);
  
  return TotTime;
}//End of buildNextLevelGraph2()


