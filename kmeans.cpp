#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


// for sorting an array
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


void findCl(vector< vector<double> >& cl,double& x, double& y,
	    double& z, double& d, int& p,int skip)
{
  double r(0.);
  double rMin(1.e20);
  
  p=0;
  d=0.0;
  for(int i(0);i<cl.size();i++)
  {
    if(i==skip)
      continue;
    
    r=(x-cl.at(i).at(0))*(x-cl.at(i).at(0));
    r+=(y-cl.at(i).at(1))*(y-cl.at(i).at(1));
    r+=(z-cl.at(i).at(2))*(z-cl.at(i).at(2));
    
    r=sqrt(r);
    
    if(r<rMin)
    {
      rMin=r;
      p=i;
      d=r;
    }
  }
}

void lumpCenters(vector< vector<double> >& cl,vector< vector<double> >& ave,
		 vector<int>& nStates,vector<int>& nAve,double& rExclude)
{
  int p(0);
  double d(0.);
  
  //cout<<"Lumping clusters."<<endl;
  for(int i(cl.size()-1);i>0;i--)
  {
    findCl(cl,cl.at(i).at(0),cl.at(i).at(1),cl.at(i).at(2),d,p,i);
    //cout<<"cluster and shortest distance: "<<i+1<<" "<<d<<" "<<p<<endl;
    if(d<=rExclude)
    {
      //cl.at(p).at(0)=((nStates.at(p)*cl.at(p).at(0))+(nStates.at(i)*cl.at(i).at(0)))/double(nStates.at(p)+nStates.at(i));
      //cl.at(p).at(1)=((nStates.at(p)*cl.at(p).at(1))+(nStates.at(i)*cl.at(i).at(1)))/double(nStates.at(p)+nStates.at(i));
      //cl.at(p).at(2)=((nStates.at(p)*cl.at(p).at(2))+(nStates.at(i)*cl.at(i).at(2)))/double(nStates.at(p)+nStates.at(i));
      
      cl.at(p).at(0)=((cl.at(p).at(0)+cl.at(i).at(0)))/2.0;
      cl.at(p).at(1)=((cl.at(p).at(1)+cl.at(i).at(1)))/2.0;
      cl.at(p).at(2)=((cl.at(p).at(2)+cl.at(i).at(2)))/2.0;
      
      cl.erase(cl.begin()+i);
      ave.erase(ave.begin()+i);
      nStates.erase(nStates.begin()+i);
      nAve.erase(nAve.begin()+i);
    }
    else if(nStates.at(i)==0 || nAve.at(i)==0)
    {
      cl.erase(cl.begin()+i);
      ave.erase(ave.begin()+i);
      nStates.erase(nStates.begin()+i);
      nAve.erase(nAve.begin()+i);
    }
  }
}

void zeroArrays(vector< vector<double> >& ave,vector<int>& nStates,vector<int>& nAve)
{
  for(int i(0);i<ave.size();i++)
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
  
  FILE *inp,*out,*outc;
  
  int maxCycle(250);
  
  double rCutoff(2.);
  double mult(3.);
  double rThrs(25.);
  double Tol(1e-4);
  
  outc=fopen("clusters.xyz","w");
  
  if(argc>=3)
  {
    inp=fopen(argv[1],"r");
    out=fopen(argv[2],"w");
  }
  else
    cout<<argv[0]<<" input output cutoff mult"<<endl;
  
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
  
  double rExclude(mult*rCutoff);
  
  cout<<"rExclude: "<<rExclude<<endl;
  
  vector<int> idTraj;
  vector<double> t;
  vector<double> x;
  vector<double> y;
  vector<double> z;
  vector<double> r;
  
  int t_idTraj(0),t_p(0);
  double t_t(0.),t_x(0.),t_y(0.),t_z(0.),t_r(0.);
  
  while(fscanf(inp,"%d %lf %lf %lf %lf %lf %d",&t_idTraj,&t_t,&t_x,&t_y,&t_z,&t_r,&t_p)!=EOF)
  {
    idTraj.push_back(t_idTraj);
    t.push_back(t_t);
    x.push_back(t_x);
    y.push_back(t_y);
    z.push_back(t_z);
    r.push_back(t_r);
  }
  
  fclose(inp);
  
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
    if(nCycle>0)
    {
      //cout<<"Cycle: "<<nCycle+1<<" nClusters: "<<cl.size()<<endl;
      lumpCenters(cl,ave,nStates,nAve,rExclude);
      zeroArrays(ave,nStates,nAve);
    }
    
    for(int i(0);i<x.size();i++)
    {
      double d(0.);
      
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
    }
    
    double dr(0.);
    drMin=1.e20;
    drMax=0.;
    isConverged=true;
    for(int i(0);i<cl.size();i++)
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
  }
  
  for(int j(0);j<x.size();j++)
  {
    int i(p.at(j)-1);
    
    double dr(cl.at(i).at(0)*cl.at(i).at(0));
    dr+=(cl.at(i).at(1)*cl.at(i).at(1));
    dr+=(cl.at(i).at(2)*cl.at(i).at(2));
    
    dr=sqrt(dr);
    
    if(dr>rThrs)
      p.at(j)=-1; 
  }
  
  for(int i(cl.size()-1);i>=0;i--)
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
  for(int i(0);i<mask.size();i++)
    mask.at(i)=i;
  
  bubble_sort(nStates,mask,nStates.size());
  
  if(!isConverged)
  {
    cout<<"Not Converged! Number of clusters: "<<cl.size()<<" drMin: "<<drMin<<" drMax "<<drMax<<endl;
    
    fprintf(outc,"%d\n",cl.size()+1);
    fprintf(outc,"clusters\n");
    fprintf(outc,"Fe  0.0000 0.00000 0.00000\n");
    for(int i(0);i<nStates.size();i++)
    {
      cout<<mask.at(i)<<" "<<nStates.at(i)<<" "<<double(nStates.at(i))/double(x.size())*100.<<endl;
      fprintf(outc,"%s %lf %lf %lf\n","Ar  ",cl.at(mask.at(i)).at(0),cl.at(mask.at(i)).at(1),cl.at(mask.at(i)).at(2));
    }
    
    for(int i(0);i<x.size();i++)
      fprintf(out,"%d %lf %lf %lf %lf %lf %d\n",idTraj.at(i),t.at(i),x.at(i),y.at(i),z.at(i),r.at(i),p.at(i));
  }
  else
  {
    fprintf(outc,"%d\n",cl.size()+1);
    fprintf(outc,"clusters\n");
    fprintf(outc,"Fe  0.0000 0.00000 0.00000\n");
    for(int i(0);i<nStates.size();i++)
    {
      cout<<mask.at(i)<<" "<<nStates.at(i)<<" "<<double(nStates.at(i))/double(x.size())*100.<<endl;
      fprintf(outc,"%s %lf %lf %lf\n","Ar  ",cl.at(mask.at(i)).at(0),cl.at(mask.at(i)).at(1),cl.at(mask.at(i)).at(2));
    }
    
    for(int i(0);i<x.size();i++)
      fprintf(out,"%d %lf %lf %lf %lf %lf %d\n",idTraj.at(i),t.at(i),x.at(i),y.at(i),z.at(i),r.at(i),p.at(i));
  }
  
  fclose(out);
  fclose(outc);
  
  return(0);
  
}
