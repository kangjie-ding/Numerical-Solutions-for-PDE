#ifndef _GENERATOR_H_
#define _GENERATOR_H_

#include<iostream>
#include<cmath>

using namespace std;

static double mu=1/81.45;
typedef enum ODEMethodType{
  AdBa=1,AdMo,BDFs,ClassicalRK,ESDIRK,GLRK,FehlbergRK,DoPrRK
}METHOD;

class TimeIntegrator
{
 public:
  virtual void Solver()=0;
  virtual ~TimeIntegrator(){};
};

//to declare classes for different methods
class AdBa_Method:public TimeIntegrator
{
 private:
  int p;
  const double* initial;
  int steps;
  double period;

 public:
 AdBa_Method(int p_,const double* initial_,int steps_,double period_):
  p{p_},initial{initial_},steps{steps_},period{period_}{}
  void Solver();
  ~AdBa_Method(){};
};

class AdMo_Method:public TimeIntegrator
{
  private:
  int p;
  const double* initial;
  int steps;
  double period;

 public:
  AdMo_Method(int p_,const double* initial_,int steps_,double period_):
  p{p_},initial{initial_},steps{steps_},period{period_}{}
  void Solver();
  ~AdMo_Method(){};
};

class BDFs_Method:public TimeIntegrator
{
  private:
  int p;
  const double* initial;
  int steps;
  double period;

 public:
  BDFs_Method(int p_,const double* initial_,int steps_,double period_):
  p{p_},initial{initial_},steps{steps_},period{period_}{}
  void Solver();
  ~BDFs_Method(){};
};

class ClassicalRK_Method:public TimeIntegrator
{
 private:
  const double* initial;
  int steps;
  double period;

 public:
 ClassicalRK_Method(const double* initial_,int steps_,double period_):
  initial{initial_},steps{steps_},period{period_}{}
  void Solver();
  ~ClassicalRK_Method(){};
};

class ESDIRK_Method:public TimeIntegrator
{
  private:
  const double* initial;
  int steps;
  double period;

 public:
  ESDIRK_Method(const double* initial_,int steps_,double period_):
  initial{initial_},steps{steps_},period{period_}{}
  void Solver();
  ~ESDIRK_Method(){};
};

class GLRK_Method:public TimeIntegrator
{
 private:
  int s;
  const double* initial;
  int steps;
  double period;

 public:
  GLRK_Method(int s_,const double* initial_,int steps_,double period_):
  s{s_},initial{initial_},steps{steps_},period{period_}{}
  void Solver();
  ~GLRK_Method(){};
};

class FehlbergRK_Method:public TimeIntegrator
{
  private:
  const double* initial;
  int steps;
  double period;
  int order=5;

 public:
  FehlbergRK_Method(const double* initial_,int steps_,double period_):
  initial{initial_},steps{steps_},period{period_}{}
  void Solver();
  ~FehlbergRK_Method(){};
};

class DoPrRK_Method:public TimeIntegrator
{
  private:
  const double* initial;
  int steps;
  double period;
  int order=4;

 public:
  DoPrRK_Method(const double* initial_,int steps_,double period_):
  initial{initial_},steps{steps_},period{period_}{}
  void Solver();
  ~DoPrRK_Method(){};
};

class TimeIntegratorFactory
{
 public:
  TimeIntegrator* CreateMethod(METHOD type,const double* initial_,int steps_,double period_,int p_or_s=1)
  {
    switch(type)
      {
      case AdBa:
	return new AdBa_Method(p_or_s,initial_,steps_,period_);
      case AdMo:
	return new AdMo_Method(p_or_s,initial_,steps_,period_);
      case BDFs:
        return new BDFs_Method(p_or_s,initial_,steps_,period_);
      case ClassicalRK:
	return new ClassicalRK_Method(initial_,steps_,period_);
      case ESDIRK:
	return new ESDIRK_Method(initial_,steps_,period_);
      case GLRK:
        return new GLRK_Method(p_or_s,initial_,steps_,period_);
      case FehlbergRK:
	return new FehlbergRK_Method(initial_,steps_,period_);
      case DoPrRK:
	return new DoPrRK_Method(initial_,steps_,period_);
      default:
	return NULL;
      }
  }
};

#endif
