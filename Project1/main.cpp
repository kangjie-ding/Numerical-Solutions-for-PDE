#include"BVP.h"
#include"function.h"
#include<jsoncpp/json/json.h>
#include<fstream>
#include<ctime>

int main(int argc,char* argv[])
{
  clock_t start,end;
  int n=64;
  string domain;
  string boundary;
  double center[2];
  double radius;
  double alpha;
  double beta;
  Json::Reader reader;
  Json::Value root;

  ifstream in("input.json",ios::binary);

  if(!in.is_open())
    {
      cout<<"Error opening file"<<endl;
      return 1;
    }

  if(reader.parse(in,root))
    {
      for(int i=0;i<root["problem_domain"].size();i++)
	{
	  domain=root["problem_domain"][i].asString();
	  for(int j=0;j<root["boundary_condition"].size();j++)
	    {
	      boundary=root["boundary_condition"][j].asString();
	      if(domain=="regular"&&boundary=="Dirichlet")
		{
		  start=clock();
		Square_Dir(n,u2,f2);
		end=clock();
		cout<<(end-start)/1000000<<endl;
		}
	      else if(domain=="regular"&&boundary=="Neumann")
		Square_Neu(n,u1,f1,u1x,u1y);
	      else if(domain=="regular"&&boundary=="mixed")
		{
		  alpha=root["parameters_mixed"]["alpha"].asDouble();
		  beta=root["parameters_mixed"]["beta"].asDouble();
		  Square_mixed(n,alpha,beta,u1,f1,u1x,u1y);
		}
	      else
		{
		  center[0]=root["disk"]["center"][0].asDouble();
		  center[1]=root["disk"]["center"][1].asDouble();
		  radius=root["disk"]["radius"].asDouble();
		  if(Valid_Check(n,center,radius))
		    {
		      if(boundary=="Dirichlet")
			Disk_Dir(n,center,radius,u1,f1);
		      else if(boundary=="Neumann")
			Disk_Neu(n,center,radius,u3,f3,u3x,u3y);
		      else
			Disk_mixed(n,alpha,beta,center,radius,u3,f3,u3x,u3y);
		    }
		  else
		    cerr<<"Invalid problem doamin!Please retype the parameters of the disk."<<endl;
		}
	    }
	}
    }
  else
    {
      cout<<"parse error!"<<endl;
    }
  return 0;
}
