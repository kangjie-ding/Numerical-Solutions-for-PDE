#ifndef _MULTIGRID_H_
#define _MULTIGRID_H_

#include<iostream>
#include<cmath>
#include"lapacke.h"
#include"Matrix.h"
#include"MatrixO.h"

using namespace std;
using namespace Numeric_lib;
using Numeric_lib::Index;

const int v1=2;
const int v2=1;
const int N=32; //the interval number of the coarsest grid
static int n=N;
const double w=0.8;
const int alpha=1;
const int beta=1;
static bool flag1=false;
static bool flag2=false;
enum Boundary_Condition{
  Dirichlet=1,Neumann,mixed
};


class Multigrid_Solver
{
 private:
  int dim;
  Boundary_Condition boundary_con;
  bool res_op;//true for full weighting,false for injection
  bool inter_op;//true for linear,false for quadratic
  bool cycle;//true for V-cycle,false for FMG
  int max_iteration;
  double epsilon;
  double guess[N+1]={0};

 public:
  friend Matrix<double,1> VC1(Multigrid_Solver& obj,Matrix<double,1>& v_h,Matrix<double,1>& f_h,int v1,int v2);
  friend Matrix<double,1> FMG1(Multigrid_Solver& obj,Matrix<double,1>& f_h,int v1,int v2);
  friend Matrix<double,1> VC2(Multigrid_Solver& obj,Matrix<double,1>& v_h,Matrix<double,1>& f_h,int v1,int v2);
  friend Matrix<double,1> FMG2(Multigrid_Solver& obj,Matrix<double,1>& f_h,int v1,int v2);
 Multigrid_Solver(int dim_,Boundary_Condition con,bool res,bool inter,bool cycle_,int max_iter,double epsilon_):dim{dim_},boundary_con{con},res_op{res},inter_op{inter},cycle{cycle_},max_iteration{max_iter},epsilon{epsilon_}{}
  void Solver();
 ~Multigrid_Solver(){};
};


#endif
