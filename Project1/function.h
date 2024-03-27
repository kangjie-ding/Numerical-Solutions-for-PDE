#ifndef _INPUT_H_
#define _INPUT_H_

#include<cmath>

double u1(double x,double y)
{
  return exp(x*x+y*y);
}

double f1(double x,double y)
{
  return -4*exp(x*x+y*y)*(1+x*x+y*y);
}

double u1x(double x,double y)
{
  return 2*x*exp(x*x+y*y);
}

double u1y(double x,double y)
{
  return 2*y*exp(x*x+y*y);
}
//
//
double u2(double x,double y)
{
  return exp(y+sin(x));
}

double u2x(double x,double y)
{
  return exp(y+sin(x))*cos(x);
}

double u2y(double x,double y)
{
  return exp(y+sin(x));
}

double f2(double x,double y)
{
  return exp(sin(x))*exp(y)*(sin(x)-cos(x)*cos(x)-1);
}
//
//
double u3(double x,double y)
{
  return exp(x+y);
}

double u3x(double x,double y)
{
  return exp(x+y);
}

double u3y(double x,double y)
{
  return exp(x+y);
}

double f3(double x,double y)
{
  return -2*exp(x+y);
}
#endif
