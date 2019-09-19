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
// ************************************************************************//#include "graphReaderFromGr.h"
#include"defs.h"
#include <iostream>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <string>
#include "Controller.h"
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <string>
#include <iostream>
#include <limits>
#include <string.h>

#include <unistd.h>
#include <cassert>
#include <inttypes.h>
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <cassert>
#include <inttypes.h>
#include <vector>
#include <omp.h>
using namespace std;
double get_time()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}


unsigned graph::allocOnHost() {
        edgessrcdst = (unsigned int *)malloc((nedges) * sizeof(unsigned int));        // first entry acts as null.
        edgessrcwt = (unsigned int *)malloc((nedges) * sizeof(unsigned int)); // first entry acts as null.
        psrc = (unsigned int *)calloc(nnodes+1, sizeof(unsigned int));  // init to null.
        psrc[nnodes] = nedges;  // last entry points to end of edges, to avoid thread divergence in drelax.
        noutgoing = (unsigned int *)calloc(nnodes, sizeof(unsigned int));       // init to 0.
  //      nincoming = (unsigned int *)calloc(nnodes, sizeof(unsigned int));       // init to 0.
        srcsrc = (unsigned int *)malloc(nnodes * sizeof(unsigned int));
      edgeListPtrs = (unsigned int *)  malloc((nnodes+1) * sizeof(unsigned int)); assert(edgeListPtrs != NULL);
	      edgeListPtrs1 = (unsigned int *)  malloc((nnodes+1) * sizeof(unsigned int)); assert(edgeListPtrs != NULL);

        edgeList = (edge *) malloc(2*nedges*sizeof(edge));
	       edgeList1 = (edge *) malloc(2*nedges*sizeof(edge));

        maxDeg = (unsigned *)malloc(sizeof(unsigned));
//	maxDegNode = (unsigned *)malloc(sizeof(unsigned));
  //      maxInDegree = (unsigned *)malloc(sizeof(unsigned));
        *maxDeg = 0;
//        *maxInDegree = 0;
	
	verPart  = (unsigned int*) malloc( (nnodes) * sizeof(int));
	borderVer = (unsigned int*) calloc( (nnodes), sizeof(int));
	dataStructureSpace = (2 * nedges + 3 * nnodes) * sizeof(unsigned int);
	dataStructureSpace = dataStructureSpace / (1024 * 1024);
//        memory = AllocatedOnHost;
	PerPartEdges = (unsigned int*) malloc( 2 * sizeof(unsigned int));
        PerPartNodes = (unsigned int*) malloc( 2 * sizeof(unsigned int));
        startVertex  = (unsigned int*) malloc( 2 * sizeof(unsigned int));
        endVertex    = (unsigned int*) malloc( 2 * sizeof(unsigned int));       
	startEdge    = (unsigned int*) malloc( 2 * sizeof(unsigned int));
	endEdge      = (unsigned int*) malloc( 2 * sizeof(unsigned int));
	bord=(bool *)malloc(sizeof(bool)*nnodes);
		bordno=(unsigned int *)malloc(sizeof(unsigned int)*nnodes);
		bordvalue=new vector<unsigned int> [nnodes];
        return 0;
		}


unsigned graph::deallocOnHost() {
        free(noutgoing);
//        free(nincoming);
        free(srcsrc);
        free(psrc);
        free(edgessrcdst);
        free(edgessrcwt);

        free(maxDeg);
        //free(maxInDegree);
        free( PerPartEdges);
	free( PerPartNodes);
	free( startVertex);
	free( endVertex);
	free( startEdge);
	free(endEdge);
        return 0;
}


unsigned graph::readFromGR(char file[]) {
        std::ifstream cfile;
        cfile.open(file);

        int masterFD = open(file, O_RDONLY);
        if (masterFD == -1) {
                printf("FileGraph::structureFromFile: unable to open %s.\n", file);
                return 1;
       			 }
        struct stat buf;
        int f = fstat(masterFD, &buf);
        if (f == -1) {
             printf("FileGraph::structureFromFile: unable to stat %s.\n", file);
             abort();
        }
        size_t masterLength = buf.st_size;
        int _MAP_BASE = MAP_PRIVATE;

	void* m = mmap(0, masterLength, PROT_READ, _MAP_BASE, masterFD, 0);
        if (m == MAP_FAILED) {
                m = 0;
                printf("FileGraph::structureFromFile: mmap failed.\n");
                abort();
        }
	
	//parse file
	uint64_t* fptr = (uint64_t*)m;
    	__attribute__((unused)) uint64_t version = le64toh(*fptr++);
    	assert(version == 1);
        uint64_t sizeEdgeTy = le64toh(*fptr++);
        uint64_t numNodes = le64toh(*fptr++);
        uint64_t numedges = le64toh(*fptr++);
        uint64_t *outIdx = fptr;
        fptr += numNodes;
        uint32_t *fptr32 = (uint32_t*)fptr;
        uint32_t *outs = fptr32;
        fptr32 += numedges;
        if (numedges % 2) fptr32 += 1;
        unsigned  *edgeData = (unsigned *)fptr32;


	// cuda.
	nnodes = numNodes;
    	nedges = numedges;
    	numVertices=nnodes;
    	numEdges=nedges;
	//	unsigned int secondMax = 0;
        printf("nnodes=%ld, nedges=%ld.\n", numVertices, numEdges);
        allocOnHost();
        int j=0;

//	bool flag[nnodes][nnodes];
        for (unsigned ii = 0; ii < nnodes; ++ii) {
		// fill unsigned *noutgoing, *nincoming, *srcsrc, *psrc, *edgessrcdst; foru *edgessrcwt
		srcsrc[ii] = ii;
		if (ii > 0) {
			
                        psrc[ii] = le64toh(outIdx[ii - 1]) ;
                        noutgoing[ii] = le64toh(outIdx[ii]) - le64toh(outIdx[ii - 1]);
			/*if( *maxDeg < noutgoing[ii] ){
			//	secondMax = *maxDeg;
				secondMaxDegNode = maxDegNode;
				*maxDeg = noutgoing[ii];
				maxDegNode = ii;
			}*/
                } else {
                        psrc[0] = 0;
                        noutgoing[0] = le64toh(outIdx[0]);
						*maxDeg = noutgoing[0];
                }
                for (unsigned jj = 0; jj < noutgoing[ii]; ++jj) {
                        unsigned edgeindex = psrc[ii] + jj + 1;
                        unsigned dst = le32toh(outs[edgeindex - 1]);
			 if (dst >= nnodes) printf("\tinvalid edge from %d to %d at index %d(%d).\n", ii, dst, jj, edgeindex);
                        edgessrcdst[edgeindex-1] = dst;
                        edgessrcwt[edgeindex-1] = edgeData[edgeindex - 1];
		//	edgeData[edgeindex-1]=1;
		//	cout<<ii<<" "<<dst<<" "<< edgessrcwt[edgeindex-1]<<endl;
		 }

	
 

        }
	cfile.close();  
	//endtime = rtclock();
/*for(unsigned int i=0;i<nnodes+1;i++)
			 {
			 	cout<<"i"<<i<<" "<<edgeListPtrs[i]<<endl;
			 }

	for(unsigned int i=0;i<nedges;i++)
	{
		cout<<edgeList[i].head<< " "<<edgeList[i].tail<<" "<<edgeList[j].weight<<endl;
        //printf("read %lld bytes in %0.2f ms (%0.2f MB/s)\n", masterLength, 1000 * (endtime - starttime), (masterLength / 1048576) / (endtime - starttime));
	}*/
//	for(unsigned int i=0;i<2*nedges;i++)
//		cout<<"edge="<<links.at(i)<<" ";
	printf("Space required to store data structure ( %lf MB )\n",dataStructureSpace);
//	printf("Max Dedgree: %ld\n",*maxDeg);
        return 0;
}

unsigned graph::partGraph(int parts, double partRatio){
//	if( dataStructureSpace < 2548 )
	unsigned int middleEdge = nedges * partRatio;
	unsigned int i = 0;
	for( i = 0; i < nnodes; i++)
		if( psrc[ i] > middleEdge)
			break;
	middleVertex = i;
//	std::cout<<"MiddleVertex: "<<middleVertex<<"\n";
	//middleVertex = nnodes * partRatio;
//	middleVertex = nnodes * partRatio;

//	else
//		middleVertex = nnodes * (1 - (948 / dataStructureSpace) );

	

	PerPartEdges[0] = psrc[ middleVertex];
	PerPartEdges[1] = nedges - psrc[ middleVertex];

	PerPartNodes[0] = middleVertex;
	PerPartNodes[1] = nnodes - middleVertex; 

	startVertex[0] = 0;
	endVertex[0]   = middleVertex;

	startVertex[1] = middleVertex;
	endVertex[1]   = nnodes;

	startEdge[0] = 0;
	endEdge[0]   = PerPartEdges[0];

	startEdge[1] = PerPartEdges[0];
	endEdge[1]   = nedges;
unsigned int borderNodes = 0;
      edgeListPtrs1 = (unsigned int *)  malloc((nnodes-middleVertex+1) * sizeof(unsigned int)); assert(edgeListPtrs1 != NULL);
    //  edgeListPtrs = (unsigned int *)  malloc((middleVertex+1) * sizeof(unsigned int)); assert(edgeListPtrs != NULL);
    for(unsigned int i=0;i<middleVertex+1;i++)
		edgeListPtrs[i]=0;
edgeListPtrs[0]=0;
	no1=no2=0;
	int jj=0;
	unsigned int *noutgoingp=(unsigned int *)malloc(nnodes*sizeof(unsigned int));
	for(unsigned int i=0;i<nnodes;i++)
		noutgoingp[i]=0;
//cout<<"2"<<endl;
bool flag=false;
	#pragma omp parallel num_threads(16)
	{
	unsigned int localBorderNodes = 0;
	#pragma omp for schedule(dynamic,100) 
	for( unsigned int i = 0; i < middleVertex; i++){	
		unsigned int start = psrc[i];
                unsigned int end   = psrc[i] + noutgoing[i];
		while( start != end ){
			unsigned int dst = edgessrcdst[ start];
			if( dst >= middleVertex ){
				localBorderNodes++;
				//partoneborderno++;
				//flag1[i]=true;
				//partoneborderedges.push_back(dst);
				bord[i]=true;
				bordno[i]+=1;
				bordvalue[i].push_back(dst);
			}
			else
			{
				edgeList[jj].head=i;
                        	edgeList[jj].tail=dst;
                        	edgeList[jj].weight=1;
				noutgoingp[i]+=1;	
				flag=true;		
			//	cout<<edgeList[jj].head<<" "<<edgeList[jj].tail<<" "<<edgeList[jj].weight<<endl;
				no1++;
				jj++;
				edgeListPtrs[i+1]=noutgoingp[i]+edgeListPtrs[i];
			}
		start++;
		}
if(flag==false){
	edgeListPtrs[i+1]=edgeListPtrs[i];}
flag=false;
	}

//cout<<"4"<<endl;
/*	 for(unsigned int i=0;i<middleVertex+1;i++)
                         {
                                edgeListPtrs[i+1]=noutgoingp[i];
                         }
                         for(unsigned int i=0;i<middleVertex+1;i++)
                         {
                                edgeListPtrs[i+1]+=edgeListPtrs[i];
                                edgeListPtrs[0]=0;
                                

                         }*/
	#pragma omp for schedule(dynamic,100)
        for( unsigned int i = middleVertex; i < nnodes; i++){ 
                unsigned int start = psrc[i];
                unsigned int end   = psrc[i] + noutgoing[i];
                while( start != end ){
                        unsigned int dst = edgessrcdst[ start];
                        if( dst < middleVertex ){
                                localBorderNodes++;
         			bord[i]=true;
				bordno[i]+=1;
				bordvalue[i].push_back(dst);	               	
                        	//parttwoborderedges.push_back(dst);
                        }
			else

			{
				links.push_back(i-middleVertex);
				links.push_back(dst-middleVertex);
				noutgoingp[i]+=1;
				no2++;
			}
		start++;
                }
        }
		__sync_fetch_and_add(&borderNodes, localBorderNodes);
	}
//cout<<"3"<<endl;
/*for(int i=0;i<middleVertex;i++)
	cout<<noutgoingp[i]<<" ";	*/
/*	for(unsigned int i=0;i<middleVertex+1;i++)
                         {
                                edgeListPtrs[i+1]=noutgoingp[i];
                         }
                         for(unsigned int i=0;i<middleVertex+1;i++)
                         {
                                edgeListPtrs[i+1]=edgeListPtrs[i+1]+edgeListPtrs[i];
                                edgeListPtrs[0]=0;


                         }*/
	int k=0;	for(unsigned int i=middleVertex;i<=nnodes;i++)
		{
			edgeListPtrs1[k]=noutgoingp[i];
			k++;
		}
	for(unsigned int i=0;i<(nnodes-middleVertex);i++)
                         {
                                edgeListPtrs1[i+1]=edgeListPtrs1[i+1]+edgeListPtrs1[i];
                                edgeListPtrs1[0]=0;


                    }




	unsigned int part0Edges = 0, part1Edges = 0;
 	#pragma omp parallel for num_threads(16) schedule(dynamic,100)
        for( unsigned int i = 0; i < middleVertex; i++){
                __sync_fetch_and_add(&part0Edges, noutgoing[ i]);

        }
	#pragma omp parallel for num_threads(16) schedule(dynamic,100)
	for( unsigned int i = middleVertex; i < nnodes; i++){
		__sync_fetch_and_add(&part1Edges, noutgoing[ i]);

	}
/*	for(int j=0;j<no1;j++)
		cout<<edgeList[j].head<<" "<<edgeList[j].tail<<endl;*/
	printf("Part 0 Ver: %ld and Part 1 Ver: %ld\n",PerPartNodes[0], PerPartNodes[1]);
//	printf("partoneborder no= %ld\n ",partoneborderno );
//	printf("parttwoborder no= %ld\n ",parttwoborderno );
	printf("Total border Nodes: %ld\n",borderNodes);
	printf("Total part0 Edges: %ld\n", part0Edges);
	printf("Total part1 Edges: %ld\n", part1Edges);
free(noutgoingp);
//free(noutgoingp);
free(edgessrcdst);
free(edgessrcwt);
free(psrc);
free(noutgoing);
free(srcsrc);
	return middleVertex;
}

unsigned graph::findSubGraph(){
	subGraphNodes = nnodes/20;
	unsigned int index = psrc[ subGraphNodes+1];
//	unsigned int *subGraphPsrc, *subGraphOutgoing, *subGraphEdgesDst, *subGraphEdgesWt;
	subGraphPsrc = (unsigned int*) malloc( (subGraphNodes+1) * sizeof(unsigned int));
	subGraphOutgoing = (unsigned int*) malloc( subGraphNodes * sizeof(unsigned int));
	subGraphEdgesDst = (unsigned int*) malloc( index * sizeof(unsigned int));
	subGraphEdgesWt  = (unsigned int*) malloc( index * sizeof(unsigned int));
	unsigned int jj = 0;
	subGraphPsrc[ 0] = 0;
	for( unsigned int i = 0; i < subGraphNodes; i++){
		unsigned int start = psrc[i];
		unsigned int end   = psrc[i] + noutgoing[i];
		unsigned int outgoing = 0;
		while ( start != end ){
			unsigned int dst = edgessrcdst[ start];
			unsigned int  wt = edgessrcwt[ start];
			if( dst  < subGraphNodes ){
				subGraphEdgesDst[jj] = dst;
				subGraphEdgesWt [jj] = wt;
				jj++;
				outgoing++;
			}		
		subGraphPsrc[ i+1] = jj;
		subGraphOutgoing[i] = outgoing;
		start++;
		}
	}
	subGraphEdges = jj;
	return 0;
}

void graph::getPartitionedData( int parts, int *partvec){
	partNumVer = (unsigned int*) calloc( parts, sizeof(unsigned int));
	partVer    = (unsigned int**) malloc( parts * sizeof(unsigned int*));
//	#pragma omp parallel num_threads(16)
//	{
		std::vector< std::vector<unsigned int> > localPartVer(parts);
		unsigned int *startLoc;
		startLoc = (unsigned int*) calloc( parts, sizeof(unsigned int));
//		#pragma omp for schedule(static)
		for( unsigned int i = 0; i < nnodes; i++){
			int partNum = partvec[ i];
			if( partNum != -1)
				localPartVer[ partNum].push_back(i);
		}

		for( int j = 0; j < parts; j++)
			startLoc[ j] = __sync_fetch_and_add(&partNumVer[j], localPartVer[j].size());
		
//		#pragma omp barrier
//		#pragma omp single
//		{
			for( int j = 0; j < parts; j++)
				partVer[ j] = (unsigned int*) malloc( partNumVer[j] * sizeof(unsigned int));
//		}

		for( int j = 0; j < parts; j++)
			std::copy(localPartVer[j].begin(), localPartVer[j].end(), partVer[j]+startLoc[ j]);
/*		for( int j = 0; j < parts; j++)
			localPartVer[j].clear();
		LOcalPartVer.clear();
		free(startLoc);*/
//	}

/*	for( int i = 0; i < parts; i++){
		for(unsigned int j = 0; j < partNumVer[ i]; j++){
			std::cout<<partVer[i][j]<<"\t";
		}
		std::cout<<"\n";
	}*/
	
	//#pragma omp parallel num_threads(16)
        //{
//	unsigned int *partEdges;
	partNumEdges = (unsigned int*) calloc( parts, sizeof(unsigned int));
	//#pragma omp for schedule(static)
	for( int i = 0; i < parts; i++){
               	for(unsigned int j = 0; j < partNumVer[ i]; j++){
			partNumEdges[ i] += noutgoing[partVer[i][j]];
		}
	}
	for( int i = 0; i < parts; i++){	
		std::cout<<"Vertices: "<<partNumVer[i]<<" Edges: "<<partNumEdges[i]<<"\n";	
	}

//	unsigned int **partPsrc, **partOutgoing, **partEdgesSrcDst;
	partPsrc = (unsigned int**) malloc( parts * sizeof(unsigned int));
        partOutgoing = (unsigned int**) malloc( parts * sizeof(unsigned int));
        partEdgesSrcDst = (unsigned int**) malloc( parts * sizeof(unsigned int));

	for( int i = 0; i < parts; i++){
		partPsrc[i] = (unsigned int*) calloc( nnodes, sizeof(unsigned int));	
		partOutgoing[i] = (unsigned int*) calloc( nnodes, sizeof(unsigned int));
		partEdgesSrcDst[i] = (unsigned int*) malloc( partNumEdges[i] * sizeof(unsigned int));
	}
	omp_set_num_threads(2);
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			unsigned int newPsrc = 0;
			for(unsigned int j = 0; j < partNumVer[ 0]; j++){
				unsigned int vertexId = partVer[0][ j];
				partPsrc[ 0][ vertexId] = newPsrc;
				partOutgoing[ 0][ vertexId] = noutgoing[ vertexId];
				std::copy(edgessrcdst + psrc[ vertexId], edgessrcdst + psrc[ vertexId] + noutgoing[ vertexId], partEdgesSrcDst[0] + newPsrc);
				newPsrc += noutgoing[ vertexId];
			}
		}
		#pragma omp section
                {
			unsigned int newPsrc = 0;
			for(unsigned int j = 0; j < partNumVer[ 1]; j++){
				unsigned int vertexId = partVer[1][ j];
			        partPsrc[ 1][ vertexId] = newPsrc;
                                partOutgoing[ 1][ vertexId] = noutgoing[ vertexId];
				std::copy(edgessrcdst + psrc[ vertexId], edgessrcdst + psrc[ vertexId] + noutgoing[ vertexId], partEdgesSrcDst[1] + newPsrc);
                                newPsrc += noutgoing[ vertexId];
			}

                }
	}
	
	/*for( int i = 0; i < parts; i++){

		std::cout<<"Part: "<<i<<"------------\n";
		for(unsigned int j = 0; j < partNumVer[ i]; j++){
			unsigned int vertexId = partVer[i][j];
	//		unsigned int start = partPsrc[ i][ vertexId];
	//		unsigned int end = partPsrc[ i][ vertexId] + partOutgoing[ i][ vertexId];
	//		while( start != end ){
				std::cout<<"src: "<<vertexId<<"\n"; //" dst: "<<partEdgesSrcDst[ i][start]<<"\n";
	//		start++;
	//		}
		}
	}*/
}	



unsigned graph::read(char file[]) {
        return readFromGR(file);
}

void graph::findSubGraph0(unsigned int currNodes, unsigned int startNode, unsigned int endNode){

        part0Outgoing = (unsigned int*) calloc(nnodes, sizeof(unsigned int));
        std::copy(noutgoing, noutgoing+currNodes, part0Outgoing);
        //traverse the adjacency of nodes
        #pragma omp parallel for schedule(dynamic,700) num_threads(7)
        for(unsigned i = currNodes; i < nnodes; i++){
                unsigned start = psrc[i];
                unsigned end   = start + noutgoing[i];
                while( start != end){
                        unsigned dst = edgessrcdst[ start];
                        if( dst < endNode)
                                part0Outgoing[i]++;
                start++;
                }
        }
	//prefix sum
	part0Psrc = (unsigned int*) calloc( (nnodes+1),sizeof(unsigned int));
        std::copy(psrc, psrc+currNodes+1, part0Psrc);
        for(unsigned i = currNodes; i < nnodes; i++)
                part0Psrc[i+1] = part0Psrc[i] + part0Outgoing[i];

        part0Edges = part0Psrc[ nnodes];
        part0EdgesDst = (unsigned int*) malloc(part0Edges * sizeof(unsigned int));
        std::copy(edgessrcdst, edgessrcdst + part0Psrc[ currNodes], part0EdgesDst);
        #pragma omp parallel for schedule(dynamic,700) num_threads(7)
        for(unsigned i = currNodes; i < nnodes; i++){
                unsigned start = psrc[i];
                unsigned end   = start + noutgoing[i];
                unsigned storePos = part0Psrc[i];
                if( part0Outgoing[i] > 0){
                        while( start != end){
                                unsigned dst = edgessrcdst[ start];
                                if( dst < endNode)
                                        part0EdgesDst[ storePos++] = dst;
                        start++;
                        }
                }
        }
        //std::cout<<"subGraph0 Update Time: "<<get_time()-traveTime<<"\n";
        printf("Number of Part 0 Edges: %ld\n", part0Edges);
}

void graph::findSubGraph1(unsigned int currNodes, unsigned int startNode, unsigned int endNode){
        part1Outgoing = (unsigned int*) calloc(nnodes, sizeof(unsigned int));
        std::copy(noutgoing+startNode, noutgoing+endNode, part1Outgoing+startNode);
        //traverse the adjacency of nodes
        #pragma omp parallel for schedule(dynamic,700) num_threads(7)
        for(unsigned i = 0; i < startNode; i++){
                     unsigned start = psrc[i];
                     unsigned end   = start + noutgoing[i];
                     while( start != end){
                            unsigned dst = edgessrcdst[ start];
                            if( dst >= startNode)
                                    part1Outgoing[i]++;
                     start++;
                     }
        }
	part1Psrc = (unsigned int*) calloc( (nnodes+1),sizeof(unsigned int));
        for(unsigned i = 0; i < startNode; i++)
                part1Psrc[i+1] = part1Psrc[i] + part1Outgoing[i];
        unsigned addedEdges = part1Psrc[ startNode];
        unsigned minusEdges = psrc[startNode];
        for(unsigned i = startNode+1; i <= endNode; i++)
                part1Psrc[i]   = psrc[i] -minusEdges + addedEdges;
        //for ghostVertex store the incident edges
        part1Edges = part1Psrc[ nnodes];
        part1EdgesDst = (unsigned int*) malloc(part1Edges * sizeof(unsigned int));
        std::copy(edgessrcdst+PerPartEdges[0], edgessrcdst + nedges, part1EdgesDst+addedEdges);
        #pragma omp parallel for schedule(dynamic,700) num_threads(7)
        for(unsigned i = 0; i < startNode; i++){
                unsigned start = psrc[i];
                unsigned end   = start + noutgoing[i];
                unsigned storePos = part1Psrc[i];
                if( part1Outgoing[i] > 0){
                        while( start != end){
                                unsigned dst = edgessrcdst[ start];
                                if( dst >= startNode)
                                        part1EdgesDst[ storePos++] = dst;
                        start++;
                        }
                }
        }
//        std::cout<<"subGraph1 Update Time: "<<get_time()-traveTime<<"\n";
	printf("Number of Part 1 Edges: %ld\n", part1Edges);
}


options parse_cmdline( int argc, char *argv[]){
	int c;
	optarg = NULL;
	options opt;
	int apps, mode;
	while( (c = getopt( argc, argv, "a:i:m:h:d:")) != EOF ){
		switch(c)
		{
			case 'a':
				apps = atoi( optarg);
				if( (apps > 4) || (apps < 0)){
					std::cout<<"Please check the manual\n";
					exit(0);
				}
				opt.apps = (application)apps;
				break;

			case 'i':
				if( optarg == NULL){
					std::cout<<"Please specify an input file\n";
					exit(0);
				}
				opt.fileName = optarg;
				break;

			case 'm':
				mode = atoi( optarg);
				if( (mode >= 3) || ( mode < 0)){
					std::cout<<"Please check the manual\n";
                                        exit(0);
				}
				opt.mode = (arch)mode;
				break;

			default:
				std::cout<<"wrong option\n";
				exit(0);
		}
	}
return opt;
}


/*int main(int argc, char *argv[]){
	options opt = parse_cmdline( argc, argv);

	std::cout<<"Input file: "<<opt.fileName<<"\n";
	std::cout<<"Application: "<<opt.apps<<"\n";
	std::cout<<"Architecture: "<<opt.mode<<"\n";
//unsigned Graph::partGraph(int parts, double partRatio)
//unsigned Graph::readFromGR(char file[]) {
graph g;
g.readFromGR("/home/anwesha/graphread/hypar/input/mst2_G.gr");
g.partGraph(2,0.3);	
//void Graph::getPartitionedData( int parts, int *partvec){
g.numVertices=g.nnodes;
g.numEdges=g.nedges;
return 0;	

}*/
