computeInternals(int *indices, unsigned int *links, float *weights, int *n2c, float *in, unsigned int nrComms, int graphType) {

    unsigned int vid = threadIdx.x / PHY_WRP_SZ;
    unsigned int laneId = threadIdx.x % PHY_WRP_SZ;
    vid = blockIdx.x * (blockDim.x / PHY_WRP_SZ) + vid;
    while (vid < nrComms) {
	if(dirty[vid]){
        unsigned int startNbr = indices[vid];
        unsigned int endNbr = indices[vid + 1];
        for (unsigned int i = startNbr + laneId; i < endNbr; i = i + PHY_WRP_SZ) {
            unsigned int nbr = links[i];
             {

                float toAdd = 0.0;

                if (graphType == UNWEIGHTED) {
                    toAdd = 1.0;
                } else {
                    toAdd = weights[i];
                }
                atomicAdd(&in[vid], toAdd);
            }
        }
        vid = vid + (blockDim.x * gridDim.x) / PHY_WRP_SZ;
    }
}


void editEdgeList(int* indices, unsigned int* links, float* weights, int gType,
        int* dirty, unsigned int dirtyvertex, unsigned int mark,
        int* vtsForPostProcessing) {

    unsigned int wid = threadIdx.x / PHY_WRP_SZ;
    unsigned int laneId = threadIdx.x % PHY_WRP_SZ;
  while (wid < nrUniDegVrts) {

        int vid = vtsForPostProcessing[wid];
        unsigned int dirtyvertex = dirty[wid];
	int startOfNbrs = indices[vid];
        int endOfNbrs = indices[vid + 1];

        for (int j = startOfNbrs + laneId; j <= endOfNbrs; j = j + PHY_WRP_SZ) {
            if (links[j] == uniDegVtx) {
                links[j] = vid;
}
        }

        wid = wid + (blockDim.x * gridDim.x) / PHY_WRP_SZ;
    }

}
