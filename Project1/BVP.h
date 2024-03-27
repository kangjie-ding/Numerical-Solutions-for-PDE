#ifndef _BVP_H_
#define _BVP_H_

#include"Matrix.h"
#include"lapacke.h"
#include"MatrixO.h"

#include<iostream>

using namespace std;
using namespace Numeric_lib;

void Square_Dir(int n,double (*u)(double x,double y),double (*f)(double x,double y));//n is the grid number

void Square_Neu(int n,double (*u)(double x,double y),double (*f)(double x,double y),double (*ux)(double x,double y),double (*uy)(double x,double y));//n is the grid number

void Square_mixed(int n,double alpha,double beta,double (*u)(double x,double y),double (*f)(double x,double y),double (*ux)(double x,double y),double (*uy)(double x,double y));//n is the grid number

void Disk_Dir(int n,double* center,double radius,double (*u)(double x,double y),double (*f)(double x,double y));//n is the grid number

void Disk_Neu(int n,double* center,double radius,double (*u)(double x,double y),double (*f)(double x,double y),double (*ux)(double x,double y),double (*uy)(double x,double y));//n is the grid number

void Disk_mixed(int n,double alpha,double beta,double* center,double radius,double (*u)(double x,double y),double (*f)(double x,double y),double (*ux)(double x,double y),double (*uy)(double x,double y));//n is the grid number

bool Valid_Check(int n,double* center,double radius);

#endif
