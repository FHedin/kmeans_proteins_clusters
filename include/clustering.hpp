/*
 * kmeans : c++ implementation of the kmeans algorithm (http://en.wikipedia.org/wiki/K-means_clustering)
 * for finding clusters of stable microstates when studying ligand migration in proteins
 * 
 * Copyright (c) 2015, Florent Hédin, Pierre-André Cazade, and the University of Basel.
 * All rights reserved.
 * 
 * The 3-clause BSD license is applied to this software.
 * See LICENSE.txt
 */

#ifndef CLUSTERING_HPP_INCLUDED
#define CLUSTERING_HPP_INCLUDED

#include<vector>

namespace kmeans
{
    
void bubble_sort(std::vector<int>& a, std::vector<int>& b,int n);

void findCl(std::vector< std::vector<double> >& cl,double& x, double& y,
            double& z, double& d, int& p,int skip);

void lumpCenters(std::vector< std::vector<double> >& cl,std::vector< std::vector<double> >& ave,
                    std::vector<int>& nStates,std::vector<int>& nAve,double& rExclude);

void zeroArrays(std::vector< std::vector<double> >& ave,std::vector<int>& nStates,std::vector<int>& nAve);
    
}

#endif // CLUSTERING_HPP_INCLUDED