#include<ctime>
#include<jsoncpp/json/json.h>
#include<fstream>
#include"generator.h"

int main(int argc,char* argv[])
{
  //some constant parameters
  const double T1=17.06521656015796;
  const double T2=19.14045706162071;
  const double u01[4]={0.994,0,0,-2.0015851063790825224};
  const double u02[4]={0.87978,0,0,-0.3797};

  clock_t start,end;
  int method;
  bool initial;
  int steps;
  double period;
  int accuracy;
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
      if(root["ODE_method"].asString()=="Adams-Bashforth")
	method=1;
      else if(root["ODE_method"].asString()=="Adams-Moulton")
	method=2;
      else if(root["ODE_method"].asString()=="BDFs")
	method=3;
      else if(root["ODE_method"].asString()=="ClassicalRK")
	method=4;
      else if(root["ODE_method"].asString()=="ESDIRK")
	method=5;
      else if(root["ODE_method"].asString()=="Gauss-Legendre")
	method=6;
      else if(root["ODE_method"].asString()=="Fehlberg")
	method=7;
      else if(root["ODE_method"].asString()=="Dormand-Prince")
	method=8;
      else
	cerr<<"Invalid input method!!"<<endl;
      
      if(root["initial_value"].asInt()==1)
	initial=true;
      else if(root["initial_value"].asInt()==2)
	initial=false;
      
      if(root["period"].asInt()==1)
	period=T1;
      else if(root["period"].asInt()==2)
	period=T2;
      else
	cerr<<"Invalid input period!!"<<endl;

      steps=root["steps"].asInt();
      accuracy=root["order_of_accuracy_or_s-stage"].asInt();
    }
  else
    cout<<"parse error!"<<endl;

  //use factory pattern to create a object
  TimeIntegratorFactory* test=new TimeIntegratorFactory();
  TimeIntegrator* solve=NULL;
  if(initial)
    solve=test->CreateMethod((METHOD)method,u01,steps,period,accuracy);
  else
    solve=test->CreateMethod((METHOD)method,u02,steps,period,accuracy);
  start=clock();
  solve->Solver();
  end=clock();
  cout<<endl;
  cout<<"the CPU time for the test is:"<<double(end-start)/double(CLOCKS_PER_SEC)<<" seconds"<<endl;
  delete test;
  delete solve;
  return 0;
}
