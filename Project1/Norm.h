#ifndef _NORM_H_
#define _NORM_H_

#include<cmath>

double norm1(double* p,int n)
{
  double norm=0;
  for(int i=0;i<n;i++)
    norm+=abs(p[i]);
  return norm;
}

double norm2(double* p,int n)
{
  double norm=0;
  for(int i=0;i<n;i++)
    norm+=p[i]*p[i];
  return sqrt(norm);
}

double normi(double* p,int n)
{
  double norm=0;
  for(int i=0;i<n;i++)
    if(abs(p[i])>norm)
      norm=abs(p[i]);
  return norm;
}
#endif
