
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
#ifndef GRAPHHOST_H
#define	GRAPHHOST_H

//#include"thrust/host_vector.h"
#include"vector"
#include"assert.h"
#include"commonconstants.h"
class GraphHOST {
public:
    unsigned int nb_nodes;
    unsigned int nb_links;
    double total_weight;
    /*
    thrust::host_vector<unsigned long> degrees;
    thrust::host_vector<unsigned int> links;
    thrust::host_vector<float> weights;
     */

    std::vector<unsigned long> degrees;
    std::vector<unsigned int> links;
    std::vector<unsigned int> weights;

    GraphHOST();

    GraphHOST(char *filename, char *filename_w, int type,unsigned int* statIndices,unsigned int*edges);

    // return the weighted degree of the node
    inline double weighted_degree(unsigned int node);

    // return the number of neighbors (degree) of the node
    inline unsigned int nb_neighbors(unsigned int node);

    //inline thrust::pair<thrust::host_vector<unsigned int>::iterator, thrust::host_vector<float>::iterator >neighbors(unsigned int node);
    inline std::pair<std::vector<unsigned int>::iterator, std::vector<unsigned int>::iterator >neighbors(unsigned int node);

    void display();

};

inline unsigned int
GraphHOST::nb_neighbors(unsigned int node) {
    assert(node >= 0 && node < nb_nodes);

    if (node == 0)
        return degrees[0];
    else
        return degrees[node] - degrees[node - 1];
}

//inline thrust::pair<thrust::host_vector<unsigned int>::iterator, thrust::host_vector<float>::iterator >

inline std::pair<std::vector<unsigned int>::iterator, std::vector<unsigned int>::iterator >
GraphHOST::neighbors(unsigned int node) {
    assert(node >= 0 && node < nb_nodes);

    if (node == 0)
        return std::make_pair(links.begin(), weights.begin());
    else if (weights.size() != 0)
        return std::make_pair(links.begin() + degrees[node - 1], weights.begin() + degrees[node - 1]);
    else
        return std::make_pair(links.begin() + degrees[node - 1], weights.begin());
}

inline double
GraphHOST::weighted_degree(unsigned int node) {
    assert(node >= 0 && node < nb_nodes);

    if (weights.size() == 0)
        return (double) nb_neighbors(node);
    else {
        //thrust::pair<thrust::host_vector<unsigned int>::iterator, thrust::host_vector<float>::iterator > p = neighbors(node);
        std::pair<std::vector<unsigned int>::iterator, std::vector<unsigned int>::iterator> p = neighbors(node);
        double res = 0;
        for (unsigned int i = 0; i < nb_neighbors(node); i++) {
            res += (double) *(p.second + i);
        }
        return res;
    }
}

#endif	/* GRAPHHOST_H */

