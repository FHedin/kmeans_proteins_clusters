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
#include <cmath>

#include "clustering.hpp"

namespace KMEANS
{
    
using namespace std;
    
// for sorting 2 vectors
void bubble_sort(vector<int>& a, vector<int>& b,int n)
{
    int t(0),s(1);
    while (s)
    {
        s = 0;
        for (int i(1); i < n; i++)
        {
            if (a[i] > a[i - 1])
            {
                t = a[i];
                a[i] = a[i - 1];
                a[i - 1] = t;
                
                t = b[i];
                b[i] = b[i - 1];
                b[i - 1] = t;
                
                s = 1;
            }
        }
    }
}

/*
 * For given x,y,z microstate coordinates, find closest cluster
 * store in p index of the cluster
 * and in d distance to the cluster
 * If skip =-1 work done for whole cl vector
 * otherwise if skip>=0 skip the corresponding entry from cl
 */
void findCl(vector< vector<double> >& cl,double& x, double& y,
            double& z, double& d, int& p,int skip)
{
    double r(0.);
    double rMin(1.e20);
    
    for(size_t i=0; i<cl.size(); i++)
    {
        //for skipping a given component of the vector cl
        if((int)i==skip)
            continue;
        
        r  = (x-cl.at(i).at(0))*(x-cl.at(i).at(0));
        r += (y-cl.at(i).at(1))*(y-cl.at(i).at(1));
        r += (z-cl.at(i).at(2))*(z-cl.at(i).at(2));
        
        r=sqrt(r);
        
        if(r<rMin)
        {
            rMin=r;
            
            p=i;
            d=r;
            
        }
    }
    
}

// merging some clusters if necessary during the cycling in main()
void lumpCenters(vector< vector<double> >& cl,vector< vector<double> >& ave,
                 vector<int>& nStates,vector<int>& nAve,double& rExclude)
{
    int p(0);
    double d(0.);
    
    //cout<<"Lumping clusters."<<endl;
    
    //starting from the end of the cluster vector
    for(int i(cl.size()-1); i>0; i--)
    {
        findCl(cl,cl.at(i).at(0),cl.at(i).at(1),cl.at(i).at(2),d,p,i);
        //cout<<"cluster and shortest distance: "<<i+1<<" "<<d<<" "<<p<<endl;
        if(d<=rExclude)
        {
            // merges the i and p clusters
            cl.at(p).at(0) = ( (cl.at(p).at(0)+cl.at(i).at(0)) ) / 2.0;
            cl.at(p).at(1) = ( (cl.at(p).at(1)+cl.at(i).at(1)) ) / 2.0;
            cl.at(p).at(2) = ( (cl.at(p).at(2)+cl.at(i).at(2)) ) / 2.0;
            
            cl.erase(cl.begin()+i);
            ave.erase(ave.begin()+i);
            nStates.erase(nStates.begin()+i);
            nAve.erase(nAve.begin()+i);
        }
        else if(nStates.at(i)==0 || nAve.at(i)==0)
        {
            //if cluster was empty remove it
            cl.erase(cl.begin()+i);
            ave.erase(ave.begin()+i);
            nStates.erase(nStates.begin()+i);
            nAve.erase(nAve.begin()+i);
        }
    }
}

// initializing 1d and 2d vectors with zeroes
void zeroArrays(vector< vector<double> >& ave,vector<int>& nStates,vector<int>& nAve)
{
    for(size_t i=0; i<ave.size(); i++)
    {
        ave.at(i).at(0)=0.0;
        ave.at(i).at(1)=0.0;
        ave.at(i).at(2)=0.0;
        
        nStates.at(i)=0;
        nAve.at(i)=0;
    }
}

}//end of namespace