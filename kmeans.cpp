/**
 * The k-means clustering algorithm
 *
 * http://en.wikipedia.org/wiki/K-means_clustering
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include <cmath>

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
 * For given x,y,z coordinates, find closest cluster
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
#ifdef _OPENMP
    #pragma omp parallel default(none) shared(x,y,z,cl,skip,p,d) private(r) firstprivate(rMin)
    {
        #pragma omp for nowait schedule(dynamic)
#endif //_OPENMP
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
#ifdef _OPENMP
                #pragma omp critical
                {
#endif //_OPENMP
                    p=i;
                    d=r;
#ifdef _OPENMP
                }
#endif //_OPENMP
            }
        }
#ifdef _OPENMP
    }//end parallel section
#endif //_OPENMP
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
#ifdef _OPENMP
    #pragma omp parallel default(none) shared(ave,nStates,nAve)
    {
        #pragma omp for nowait
#endif //_OPENMP
        for(size_t i=0; i<ave.size(); i++)
        {
            ave.at(i).at(0)=0.0;
            ave.at(i).at(1)=0.0;
            ave.at(i).at(2)=0.0;

            nStates.at(i)=0;
            nAve.at(i)=0;
        }
#ifdef _OPENMP
    }
#endif //_OPENMP

}

int main(int argc, char* argv[])
{

    FILE *inp,*out,*outc;

    //maximum number of iterations
    int maxCycle(250);

    //Cutoff for cluster determination
    double rCutoff(2.);

    //Multiplicator factor of the cutoff for exclusion of lonely microstates
    double mult(1.);

    //Threshold to consider a microstate is in water
    double rThrs(25.);

    //Tolerance for convergence
    double Tol(1e-4);

    bool useDefault(true);

    if(argc>=4)
    {
        inp=fopen(argv[1],"r");
        out=fopen(argv[2],"w");
        // coordinates file for storing clusters location
        outc=fopen(argv[3],"w");
        if(argc>=5)
            useDefault=(strcmp("--no-default",argv[4])==0)?false:true;
        //cout << useDefault << endl;
    }
    else
    {
        cout << "Error, not enough arguments, usage is : " << argv[0] << " {path to inputFile} {path to outputFile} {path to outputXYZ} [--no-default]"<<endl;
        cout << "inputFile and outputFile and outputXYZ are necessary, --no-default is optional if not given default values for some internal variables are used,"
             "otherwise the user will have to provide those values from command line." << endl;
        return EXIT_FAILURE;
    }


    if(useDefault)
    {
        cout << "Using default values for :" << endl;
        cout << "\t rCutoff : " << rCutoff << endl;
        cout << "\t mult : " << mult << endl;
        cout << "\t rThrs : " << rThrs << endl;
        cout << "\t Tol : " << Tol << endl;
        cout << "\t maxCycle : " << maxCycle << endl;
    }
    else
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

    double rExclude(mult*rCutoff);

    cout<<"Global exclusion (Cutoff*MultiplicatorFactor) is rExclude = "<<rExclude<<endl;

    /*
     * vectors sor storing data read from input file
     */
    vector<int> idTraj;
    vector<double> t;
    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<double> r;

    //temporary variables for reading from file and then storing in vectors
    int t_idTraj(0),t_p(0);
    double t_t(0.),t_x(0.),t_y(0.),t_z(0.),t_r(0.);

    // How was generated the input file ??
//     while(fscanf(inp,"%d %lf %lf %lf %lf %lf %d",&t_idTraj,&t_t,&t_x,&t_y,&t_z,&t_r,&t_p)!=EOF)
    while(fscanf(inp,"%lf %lf %lf",&t_x,&t_y,&t_z)!=EOF)
    {
//         idTraj.push_back(t_idTraj);
//         t.push_back(t_t);
        x.push_back(t_x);
        y.push_back(t_y);
        z.push_back(t_z);
//         r.push_back(t_r);
    }

    fclose(inp);

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
        // remove some clusters if required
        if(nCycle>0)
        {
            //cout<<"Cycle: "<<nCycle+1<<" nClusters: "<<cl.size()<<endl;
            lumpCenters(cl,ave,nStates,nAve,rExclude);
            zeroArrays(ave,nStates,nAve);
        }

        //iteration over number of atoms from x or y or z vectors
        for(size_t i(0); i<x.size(); i++)
        {
            double d(0.);

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

        fprintf(outc,"%d\n",(int)cl.size()+1);
        fprintf(outc,"clusters\n");
        fprintf(outc,"Fe  0.0000 0.00000 0.00000\n");
        for(size_t i(0); i<nStates.size(); i++)
        {
            cout<<mask.at(i)<<" "<<nStates.at(i)<<" "<<double(nStates.at(i))/double(x.size())*100.<<endl;
            fprintf(outc,"%s %lf %lf %lf\n","Ar  ",cl.at(mask.at(i)).at(0),cl.at(mask.at(i)).at(1),cl.at(mask.at(i)).at(2));
        }

        for(size_t i(0); i<x.size(); i++)
            fprintf(out,"%d %lf %lf %lf %lf %lf %d\n",idTraj.at(i),t.at(i),x.at(i),y.at(i),z.at(i),r.at(i),p.at(i));
    }
    else
    {
        fprintf(outc,"%d\n",(int)cl.size()+1);
        fprintf(outc,"clusters\n");
        fprintf(outc,"Fe  0.0000 0.00000 0.00000\n");
        for(size_t i(0); i<nStates.size(); i++)
        {
            cout<<mask.at(i)<<" "<<nStates.at(i)<<" "<<double(nStates.at(i))/double(x.size())*100.<<endl;
            fprintf(outc,"%s %lf %lf %lf\n","Ar  ",cl.at(mask.at(i)).at(0),cl.at(mask.at(i)).at(1),cl.at(mask.at(i)).at(2));
        }

        for(size_t i(0); i<x.size(); i++)
            fprintf(out,"%d %lf %lf %lf %lf %lf %d\n",idTraj.at(i),t.at(i),x.at(i),y.at(i),z.at(i),r.at(i),p.at(i));
    }

    fclose(out);
    fclose(outc);

    return(0);

}
