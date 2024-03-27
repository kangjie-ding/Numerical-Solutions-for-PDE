#ifndef _INPUT_H_
#define _INPUT_H_

#include<cmath>
//the first function

double u(double x)
{
  return exp(sin(x));
}

double f(double x)
{
  return exp(sin(x))*(sin(x)-cos(x)*cos(x));
}

double ux(double x)
{
  return exp(sin(x))*cos(x);
}

double u(double x,double y)
{
  return exp(y+sin(x));
}

double ux(double x,double y)
{
  return exp(y+sin(x))*cos(x);
}

double uy(double x,double y)
{
  return exp(y+sin(x));
}

double f(double x,double y)
{
  return exp(sin(x))*exp(y)*(sin(x)-cos(x)*cos(x)-1);
}

//the second function
/*
double u(double x)
{
  return exp(x*x);
}

double ux(double x)
{
  return 2*x*exp(x*x);
}

double f(double x)
{
  return -exp(x*x)*(4*x*x+2);
}

double u(double x,double y)
{
  return exp(x*x+y*y);
}

double f(double x,double y)
{
  return -4*exp(x*x+y*y)*(1+x*x+y*y);
}

double ux(double x,double y)
{
  return 2*x*exp(x*x+y*y);
}

double uy(double x,double y)
{
  return 2*y*exp(x*x+y*y);
}

//the third function

double u(double x)
{
  return exp(x);
}

double ux(double x)
{
  return exp(x);
}

double f(double x)
{
  return -exp(x);
}

double u(double x,double y)
{
  return exp(x+y);
}

double ux(double x,double y)
{
  return exp(x+y);
}

double uy(double x,double y)
{
  return exp(x+y);
}

double f(double x,double y)
{
  return -2*exp(x+y);
}
*/
#endif
