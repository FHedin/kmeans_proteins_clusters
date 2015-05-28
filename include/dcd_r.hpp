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


#include "dcd.hpp"

#ifndef DCD_R_HPP
#define	DCD_R_HPP

class DCD_R : public DCD
{

private:
    //no private attributes
    //private methods
    virtual void alloc();
    void read_header();
    bool is_allocated;
    
public:
    
    // no public attributes
    // public methods
    DCD_R(const char filename[]); //constructor
    DCD_R(const DCD_R& d);
    

    void read_oneFrame();
    void printHeader() const;
        
    ~DCD_R();

};

#endif	/* DCD_R_HPP */

