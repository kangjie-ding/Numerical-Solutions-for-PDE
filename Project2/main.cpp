#include<ctime>
#include"Multigrid.h"
#include<jsoncpp/json/json.h>
#include<fstream>

int main(int argc,char* argv[])
{
  clock_t start,end;
  int dim;
  int boundary;
  bool restriction;
  bool interpolation;
  bool cycle;
  int iteration;
  double epsilon;
  Json::Reader reader;
  Json::Value root;

  ifstream in("input.json",ios::binary);

  if(!in.is_open())
    {
      cout<<"Error opening file!"<<endl;
      return 1;
    }

  if(reader.parse(in,root))
    {
      if(root["dimension"].asInt()==1)
	 dim=1;
      else if(root["dimension"].asInt()==2)
	dim=2;
      else
	cout<<"Invalid dimension!"<<endl;
      
      if(root["boundary_condition"].asString()=="Dirichlet")
	boundary=1;
      else if(root["boundary_condition"].asString()=="Neumann")
	boundary=2;
      else if(root["boundary_condition"].asString()=="mixed")
	boundary=3;
      else
	cout<<"Invalid boundary condition!"<<endl;
      
      if(root["restriction_operator"].asString()=="full weighting")
	restriction=true;
      else if(root["restriction_operator"].asString()=="injection")
	restriction=false;
      else
	cout<<"Invalid restriction operator!"<<endl;

      if(root["interpolation_operator"].asString()=="linear")
	interpolation=true;
      else if(root["interpolation_operator"].asString()=="quadratic")
	interpolation=false;
      else
	cout<<"Invalid interpolation operator!"<<endl;

      if(root["cycles"].asString()=="V-cycle")
	cycle=true;
      else if(root["cycles"].asString()=="FMG")
	cycle=false;
      else
	cout<<"Invalid cycle!"<<endl;

      iteration=root["stopping_criteria"]["max_iteration"].asInt();
      epsilon=root["stopping_criteria"]["accuracy"].asDouble();
    }
  else
    cout<<"parse error!"<<endl;
  Multigrid_Solver test{dim,(Boundary_Condition)boundary,restriction,interpolation,cycle,iteration,epsilon};
  start=clock();
  test.Solver();
  end=clock();
  cout<<"the CPU running time is:"<<(end-start)/1000000<<" seconds"<<endl;
  return 0;
}


