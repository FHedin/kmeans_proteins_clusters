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
#include <cmath>

#include "dcd_r.hpp"

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

int main(int argc, char* argv[])
{
    
    if(argc<9)
    {
        cout << "Error, not enough arguments, usage is : " << argv[0] << " -idx {index of atom to study (taken from a PSF for example)} -dcd {path to DCD} -out {path to outputFile} -xyz {path to outputXYZ}"<<endl;
        cout << "inputFile and outputFile and outputXYZ are necessary fileNames" << endl << endl;
        cout << "optional arguments : " << endl;
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
    
    DCD_R *dcdf = nullptr;
    
    // arguments parsing
    for (int i=1; i<argc; i++)
    {
        // index of atom to study from dcd
        if (!strcasecmp(argv[i],"-idx"))
        {
            dcdIndex = atoi(argv[++i]);
        }
        // get name of dcd file
        else if (!strcasecmp(argv[i],"-dcd"))
        {
            dcdf = new DCD_R(argv[++i]);
            dcdf->read_header();
            dcdf->printHeader();
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
    }

//     if(argc>=4)
//     {
//         inp=fopen(argv[1],"r");
//         outfile=fopen(argv[2],"w");
//         // coordinates file for storing clusters location
//         xyzf=fopen(argv[3],"w");
//         if(argc>=5)
//             useDefault=(strcmp("--no-default",argv[4])==0)?false:true;
//         //cout << useDefault << endl;
//     }
//     else
//     {
//         cout << "Error, not enough arguments, usage is : " << argv[0] << " {path to inputFile} {path to outputFile} {path to outputXYZ} [--no-default]"<<endl;
//         cout << "inputFile and outputFile and outputXYZ are necessary, --no-default is optional if not given default values for some internal variables are used,"
//              "otherwise the user will have to provide those values from command line." << endl;
//         return EXIT_FAILURE;
//     }


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
    
//     exit(0);

    /*
     * vectors sor storing data read from input file
     */
    vector<int> idTraj;
    vector<double> t;
    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<double> r;

//     // How was generated the input file ??
// //     while(fscanf(inp,"%d %lf %lf %lf %lf %lf %d",&t_idTraj,&t_t,&t_x,&t_y,&t_z,&t_r,&t_p)!=EOF)
//     FILE *inp;
//     while(fscanf(inp,"%lf %lf %lf",&t_x,&t_y,&t_z)!=EOF)
//     {
// //         idTraj.push_back(t_idTraj);
// //         t.push_back(t_t);
//         x.push_back(t_x);
//         y.push_back(t_y);
//         z.push_back(t_z);
// //         r.push_back(t_r);
//     }

//     fclose(inp);
   
    
    // in this loop the coordinates are read frame by frame
    for(int i=0;i<dcdf->getNFILE();i++)
    {
//         dcdf->
        
        const float *lx,*ly,*lz;
        
        dcdf->read_oneFrame();

        lx=dcdf->getX();
        ly=dcdf->getY();
        lz=dcdf->getZ();
        
        x.push_back(lx[dcdIndex-1]);
        y.push_back(ly[dcdIndex-1]);
        z.push_back(lz[dcdIndex-1]);
    }

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

// int main2(int argc, char* argv[])
// {
//     // instance of a new object DCD_R attached to a dcd file 
//     DCD_R dcdf("dyna.dcd");
//     
//     // read the header and print it
//     dcdf.read_header();
//     dcdf.printHeader();
//     
//     const float *x,*y,*z;
//     
//     // in this loop the coordinates are read frame by frame
//     for(int i=0;i<dcdf.getNFILE();i++)
//     {
//         dcdf.read_oneFrame();
//         
//         /* your code goes here */
//         
//         x=dcdf.getX();
//         y=dcdf.getY();
//         z=dcdf.getZ();
//         
//         /* ... */
//         
//     }
//     
//     return EXIT_SUCCESS;
// }
