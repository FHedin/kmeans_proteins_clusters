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

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "dcd_r.hpp"
#include "align.hpp"
#include "clustering.hpp"

using namespace std;
using namespace KMEANS;


int main(int argc, char* argv[])
{
    if(argc<10)
    {
        cout << "Error, not enough arguments, usage is : " << argv[0] << " -idx {index of atom to study (taken from a PSF for example)} -dcd {number of dcd files} {paths to DCDs} -out {path to outputFile} -xyz {path to outputXYZ}"<<endl;
        cout << "inputDCDs and outputFile and outputXYZ are necessary fileNames" << endl << endl;
        cout << "optional arguments : " << endl;
        cout << "-align {first} {last} \t if present, align all frames from all dcds relative to first frame of first dcd before starting clustering" << endl;
        cout << "-interactive \t if present the user will have to provide parameters interactively" << endl;
        cout << "-cycles \t provide number of cycles" << endl;
        cout << "-cutoff \t provide cutoff value for cluster determination" << endl;
        cout << "-mult \t provide the Multiplicator factor of the cutoff for exclusion of lonely microstates" << endl;
        cout << "-thresh \t provide the Threshold to consider a microstate is in water" << endl;
        cout << "-tolerance \t provide the Tolerance for convergence" << endl;
        return EXIT_FAILURE;
    }


    FILE *outfile=nullptr, *xyzf=nullptr;

    int dcdIndex(-1);
    
    //maximum number of iterations
    int maxCycle(250);

    //Cutoff for cluster determination
    double rCutoff(1.7);

    //Multiplicator factor of the cutoff for exclusion of lonely microstates
    double mult(5.);

    //Threshold to consider a microstate is in water
    double rThrs(25.);

    //Tolerance for convergence
    double Tol(1e-4);

    bool useDefault(true);
    
    vector<string> dcds_list;
    DCD_R *dcdf = nullptr;
    
    bool needAlign(false);
    int firstAlign(-1);
    int lastAlign(-1);
    
    // arguments parsing
    for (int i=1; i<argc; i++)
    {
        // index of atom to study from dcd
        if (!strcasecmp(argv[i],"-idx"))
        {
            dcdIndex = atoi(argv[++i]);
        }
        // get name of dcd files
        else if (!strcasecmp(argv[i],"-dcd"))
        {
            int ndcd = atoi(argv[++i]);
            cout << "User provided " << ndcd << " dcd files : " << endl;
            
            for(int it=0;it<ndcd;it++)
            {
                dcds_list.push_back(string(argv[++i]));
                cout << dcds_list.at(it) << endl;
            }    
        }
        // output file
        else if (!strcasecmp(argv[i],"-out"))
        {
            outfile = fopen(argv[++i],"wt");
        }
        // output xyz
        else if (!strcasecmp(argv[i],"-xyz"))
        {
            xyzf = fopen(argv[++i],"wt");
        }
        // get number of cycles
        else if (!strcasecmp(argv[i],"-cycles"))
        {
            maxCycle = atoi(argv[++i]);
        }
        else if (!strcasecmp(argv[i],"-cutoff"))
        {
            rCutoff = strtod(argv[++i],nullptr);
        }
        else if (!strcasecmp(argv[i],"-mult"))
        {
            mult = strtod(argv[++i],nullptr);
        }
        else if (!strcasecmp(argv[i],"-thresh"))
        {
            rThrs = strtod(argv[++i],nullptr);
        }
        else if (!strcasecmp(argv[i],"-tolerance"))
        {
            Tol = strtod(argv[++i],nullptr);
        }
        else if (!strcasecmp(argv[i],"-interactive"))
        {
            useDefault = false;
        }
        // optionnall align all frames
        else if (!strcasecmp(argv[i],"-align"))
        {
            needAlign = true;
            firstAlign = atoi(argv[++i]);
            lastAlign = atoi(argv[++i]);
        }
    }
    
    if(!useDefault)
    {
        cout<<"Cutoff for cluster determination?"<<endl;
        cin>>rCutoff;
        cout<<"Multiplicator factor of the cutoff for exclusion of lonely microstates?"<<endl;
        cin>>mult;
        cout<<"Threshold to consider a microstate is in water?"<<endl;
        cin>>rThrs;
        cout<<"Tolerance for convergence?"<<endl;
        cin>>Tol;
        cout<<"Max cycles?"<<endl;
        cin>>maxCycle;
    }

    cout << "Values used for parameters :" << endl;
    cout << "\t rCutoff : " << rCutoff << endl;
    cout << "\t mult : " << mult << endl;
    cout << "\t rThrs : " << rThrs << endl;
    cout << "\t Tol : " << Tol << endl;
    cout << "\t maxCycle : " << maxCycle << endl;

    double rExclude(mult*rCutoff);

    cout<<"Global exclusion (Cutoff*MultiplicatorFactor) is rExclude = "<<rExclude<<endl;
    
    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<float> lx, lxr, ly, lyr, lz, lzr;
    vector<bool> selection;
    bool refRead=false;
  
    // iterate over all provided dcds
    for (string st : dcds_list)
    {
        cout << endl << "Reading coordinates from dcd : " << st << endl;
        dcdf = new DCD_R(st.c_str());
        dcdf->printHeader();
        
        if(needAlign)
        {
            selection.assign(dcdf->getNATOM(),false);
//             cout << "Dump of aligning selection : " << endl;
            for(int it=firstAlign-1;it<=lastAlign-1;it++)
            {
                selection.at(it) = true;
            }
//             for(bool bl : selection)
//                 cout << bl << '\t';
//             cout << endl;
        }
        
        // in this loop the coordinates are read frame by frame
        for(int i=0; i<dcdf->getNFILE(); i++)
        {
            const float *dcdx=nullptr,*dcdy=nullptr,*dcdz=nullptr;
            
            dcdf->read_oneFrame();
            
            dcdx=dcdf->getX();
            dcdy=dcdf->getY();
            dcdz=dcdf->getZ();
            
            if(needAlign)
            {
                // first frame of first dcd used as reference for aligning all frames
                if(refRead==false)
                {

                    for (int it = 0; it < dcdf->getNATOM(); it++)
                    {
                        lxr.push_back(dcdx[it]);
                        lyr.push_back(dcdy[it]);
                        lzr.push_back(dcdz[it]);
                    }
                    
                    refRead=true;
                    
                    x.push_back( lxr.at(dcdIndex) );
                    y.push_back( lyr.at(dcdIndex) );
                    z.push_back( lzr.at(dcdIndex) );
                    
                }//first read
                else
                {
                    for (int it = 0; it < dcdf->getNATOM(); it++)
                    {
                        lx.push_back(dcdx[it]);
                        ly.push_back(dcdy[it]);
                        lz.push_back(dcdz[it]);
                    }
                    
                    // now align coordinates from dcd to lxr lyr lzr
                    fprintf(stdout,"Reading dcd %s and reading+aligning frame %d of %d\r",st.c_str(),i+1,dcdf->getNFILE());
                    fflush(stdout);
                    ALIGN::align_to_ref(lx,ly,lz,lxr,lyr,lzr,selection);
                    
                    //now aligned coordinates are ready for being added to real x y z vector used later
                    x.push_back( lx.at(dcdIndex) );
                    y.push_back( ly.at(dcdIndex) );
                    z.push_back( lz.at(dcdIndex) );

                    
                    lx.clear();
                    ly.clear();
                    lz.clear();
                }//others
            }//need align
            else
            {
                x.push_back( dcdx[dcdIndex] );
                y.push_back( dcdy[dcdIndex] );
                z.push_back( dcdz[dcdIndex] );
            }//not needed align
            
        }//loop on reading frames
        
        cout << "Done for dcd : " << st << endl;
        
        delete dcdf;
        dcdf=nullptr;
        
    }
    
    cout << endl << "Total number of frames read from the " << dcds_list.size() << " dcds is : " << x.size() << endl ;
    
//     cout << "Dump of vectors x y z of size : " << x.size() << '\t' << y.size() << '\t' << z.size() << endl;
//     for (uint i=0; i<x.size(); i++)
//     {
//         printf("%lf\t%lf\t%lf\n",x[i],y[i],z[i]);
//     }
    
    selection.clear();
    lx.clear();
    ly.clear();
    lz.clear();
    lxr.clear();
    lyr.clear();
    lzr.clear();
    
//     exit(0);
 
    /*
     * vectors for storing data read from input file
     */
    vector<int> idTraj;
    vector<double> t;
    vector<double> r;
    
    idTraj.resize(x.size(),0);
    t.resize(x.size(),0.0);
    r.resize(x.size(),0.0);

    vector<int> p(x.size(),0);
    vector<int> nStates,nAve;
    vector< vector <double> > cl;
    vector< vector <double> > ave;

    int nCycle(0);
    bool isConverged(false);
    double drMin(1.e20);
    double drMax(0.0);
    
    cout << endl << "Now performing the clustering task ..." << endl << endl ;

    while(!isConverged && nCycle<maxCycle )
    {
        //cout<<"Cycle: "<<nCycle+1<<" nClusters: "<<cl.size()<<endl;
        // remove some clusters if required
        if(nCycle>0)
        {
            //cout<<"Cycle: "<<nCycle+1<<" nClusters: "<<cl.size()<<endl;
            lumpCenters(cl,ave,nStates,nAve,rExclude);
            zeroArrays(ave,nStates,nAve);
        }

        //iteration over number of microstates from x or y or z vectors
        for(size_t i(0); i<x.size(); i++)
        {
            double d(0.);
            int t_p(0);

            // if cl not empty
            if(cl.size()>0)
            {
                findCl(cl,x.at(i),y.at(i),z.at(i),d,t_p,-1);
                if(d>rCutoff)
                {
                    if(d<=rExclude)
                    {
                        p.at(i)=t_p+1;
                        nStates.at(t_p)++;
                    }
                    else
                    {
                        cl.push_back(vector<double>(3));
                        cl.at(cl.size()-1).at(0)=x.at(i);
                        cl.at(cl.size()-1).at(1)=y.at(i);
                        cl.at(cl.size()-1).at(2)=z.at(i);

                        nStates.push_back(1);

                        p.at(i)=cl.size();

                        ave.push_back(vector<double>(3));
                        ave.at(ave.size()-1).at(0)=x.at(i);
                        ave.at(ave.size()-1).at(1)=y.at(i);
                        ave.at(ave.size()-1).at(2)=z.at(i);

                        nAve.push_back(1);
                    }
                }
                else
                {
                    p.at(i)=t_p+1;
                    nStates.at(t_p)++;

                    ave.at(t_p).at(0)+=x.at(i);
                    ave.at(t_p).at(1)+=y.at(i);
                    ave.at(t_p).at(2)+=z.at(i);

                    nAve.at(t_p)++;
                }

            }
            // if empty at beginning of programm
            else
            {
                cl.push_back(vector<double>(3,0.));
                cl.at(cl.size()-1).at(0)=x.at(i);
                cl.at(cl.size()-1).at(1)=y.at(i);
                cl.at(cl.size()-1).at(2)=z.at(i);

                nStates.push_back(1);

                p.at(i)=cl.size();

                ave.push_back(vector<double>(3,0.));
                ave.at(ave.size()-1).at(0)=x.at(i);
                ave.at(ave.size()-1).at(1)=y.at(i);
                ave.at(ave.size()-1).at(2)=z.at(i);

                nAve.push_back(1);
            }
        }//for iteration over coordinates

        double dr(0.);
        drMin=1.e20;
        drMax=0.;
        isConverged=true;
        for(size_t i(0); i<cl.size(); i++)
        {
            if(nAve.at(i)>0)
            {
                ave.at(i).at(0)/=double(nAve.at(i));
                ave.at(i).at(1)/=double(nAve.at(i));
                ave.at(i).at(2)/=double(nAve.at(i));

                dr=(cl.at(i).at(0)-ave.at(i).at(0))*(cl.at(i).at(0)-ave.at(i).at(0));
                dr+=(cl.at(i).at(1)-ave.at(i).at(1))*(cl.at(i).at(1)-ave.at(i).at(1));
                dr+=(cl.at(i).at(2)-ave.at(i).at(2))*(cl.at(i).at(2)-ave.at(i).at(2));

                dr=sqrt(dr);

                //cout << dr << "\t" << rThrs << endl;

                if(dr>rThrs)
                    continue;

                if(dr>drMax)
                    drMax=dr;

                if(dr<drMin)
                    drMin=dr;

                if(dr>Tol)
                {
                    isConverged=false;
                    cl.at(i).at(0)=ave.at(i).at(0);
                    cl.at(i).at(1)=ave.at(i).at(1);
                    cl.at(i).at(2)=ave.at(i).at(2);
                }
            }
        }

        nCycle++;

    }// end while loop

    for(size_t j(0); j<x.size(); j++)
    {
        int i(p.at(j)-1);

        double dr(cl.at(i).at(0)*cl.at(i).at(0));
        dr+=(cl.at(i).at(1)*cl.at(i).at(1));
        dr+=(cl.at(i).at(2)*cl.at(i).at(2));

        dr=sqrt(dr);

        if(dr>rThrs)
            p.at(j)=-1;
    }

    for(int i(cl.size()-1); i>=0; i--)
    {
        double dr(cl.at(i).at(0)*cl.at(i).at(0));
        dr+=(cl.at(i).at(1)*cl.at(i).at(1));
        dr+=(cl.at(i).at(2)*cl.at(i).at(2));

        dr=sqrt(dr);

        if(dr>rThrs)
        {
            cl.erase(cl.begin()+i);
            nStates.erase(nStates.begin()+i);
        }
    }

    vector<int> mask(nStates.size());
    for(size_t i(0); i<mask.size(); i++)
        mask.at(i)=i;

    bubble_sort(nStates,mask,nStates.size());

    if(!isConverged)
    {
        cout<<"Not Converged! Number of clusters: "<<cl.size()<<" drMin: "<<drMin<<" drMax "<<drMax<<endl;

        fprintf(xyzf,"%d\n",(int)cl.size()+1);
        fprintf(xyzf,"clusters\n");
        fprintf(xyzf,"Fe  0.0000 0.00000 0.00000\n");
        for(size_t i(0); i<nStates.size(); i++)
        {
            cout<<mask.at(i)<<" "<<nStates.at(i)<<" "<<double(nStates.at(i))/double(x.size())*100.<<endl;
            fprintf(xyzf,"%s %lf %lf %lf\n","Ar  ",cl.at(mask.at(i)).at(0),cl.at(mask.at(i)).at(1),cl.at(mask.at(i)).at(2));
        }

        for(size_t i(0); i<x.size(); i++)
            fprintf(outfile,"%d %lf %lf %lf %lf %lf %d\n",idTraj.at(i),t.at(i),x.at(i),y.at(i),z.at(i),r.at(i),p.at(i));
    }
    else
    {
        fprintf(xyzf,"%d\n",(int)cl.size()+1);
        fprintf(xyzf,"clusters\n");
        fprintf(xyzf,"Fe  0.0000 0.00000 0.00000\n");
        for(size_t i(0); i<nStates.size(); i++)
        {
            cout<<mask.at(i)<<" "<<nStates.at(i)<<" "<<double(nStates.at(i))/double(x.size())*100.<<endl;
            fprintf(xyzf,"%s %lf %lf %lf\n","Ar  ",cl.at(mask.at(i)).at(0),cl.at(mask.at(i)).at(1),cl.at(mask.at(i)).at(2));
        }

        for(size_t i(0); i<x.size(); i++)
            fprintf(outfile,"%d %lf %lf %lf %lf %lf %d\n",idTraj.at(i),t.at(i),x.at(i),y.at(i),z.at(i),r.at(i),p.at(i));
    }

    fclose(outfile);
    fclose(xyzf);

    return(0);

}

