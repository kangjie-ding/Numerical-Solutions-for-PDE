#include"Multigrid.h"
#include"function.h"

double norm(Matrix<double,1> a);
Matrix<double,1> operator*(const Matrix<double,2>& A,const Matrix<double,1>& b);
Matrix<double,1> operator+(const Matrix<double,1>& a,const Matrix<double,1>& b);
Matrix<double,1> operator-(const Matrix<double,1>& a,const Matrix<double,1>& b);
Matrix<double,1> operator/(const Matrix<double,1>& a,const Matrix<double,1>& b);


Matrix<double,1> VC1(Multigrid_Solver& obj,Matrix<double,1>& v_h,Matrix<double,1>& f_h,int v1,int v2)
{
  double h=1.0/n;
  Matrix<double,1> f_2h(n/2+1);
  Matrix<double,1> v_2h(n/2+1);
  Matrix<double,2> a(n+1,n+1);
  Matrix<double,2> T(n+1,n+1);
  Matrix<double,1> c(n+1);
  Matrix<double,2> I_h_2h(n/2+1,n+1);
  Matrix<double,2> I_2h_h(n+1,n/2+1);
  //define coefficient matrix
  if(obj.boundary_con==Dirichlet)
    {
       for(int i=0;i<=n;i++)
	{
	  if(i==0||i==n)
	      a(i,i)=1;
	  else
	    {
	      a(i,i)=2/(h*h);
	      a(i,i-1)=-1/(h*h);
	      a(i,i+1)=-1/(h*h);
	    }
	 }
    }
  else if(obj.boundary_con==Neumann)
    {
    for(int i=0;i<=n;i++)
      {
	a(i,i)=2/(h*h);
	if(i==0)
	  a(i,i)=1;
	else if(i==n)
	    a(i,i-1)=-2/(h*h);
	else
	  {
	    a(i,i-1)=-1/(h*h);
	    a(i,i+1)=-1/(h*h);
	  }
      }
    }
  else
    {
      for(int i=0;i<=n;i++)
	{
	if(i==0)
	  {
	    a(i,i)=2/(h*h);
	    a(i,i+1)=-2/(h*h);
	  }
	else if(i==n)
	    a(i,i)=1;
	else
	  {
	    a(i,i)=2/(h*h);
	    a(i,i+1)=-1/(h*h);
	    a(i,i-1)=-1/(h*h);
	  }
	}
    }
   //define iteration matrix T and corresponding vector c
   for(int i=0;i<=n;i++)
    {
      T(i,i)=1-w;
      if(i==0)
	T(i,i+1)=-w/a(i,i)*a(i,i+1);
      else if(i==n)
	T(i,i-1)=-w/a(i,i)*a(i,i-1);
      else
	{
	  T(i,i+1)=-w/a(i,i)*a(i,i+1);
	  T(i,i-1)=-w/a(i,i)*a(i,i-1);
	}
    }
  for(int i=0;i<=n;i++)
    c(i)=w/a(i,i)*f_h(i);
  //define restriction operator
    if(obj.res_op)
    {
      for(int i=0;i<=n/2;i++)
	{
	if(i==0)
	  I_h_2h(i,i)=1;
	else if(i==n/2)
	  I_h_2h(i,n)=1;
	else
	  {
	    I_h_2h(i,2*i-1)=0.25;
	    I_h_2h(i,2*i)=0.5;
	    I_h_2h(i,2*i+1)=0.25;
	  }
	}	
    }
  else
    {
      for(int i=0;i<=n/2;i++)
	I_h_2h(i,2*i)=1;
    }
  //define interpolation operator
    if(obj.inter_op)
    {
      for(int i=0;i<=n/2;i++)
	{
	if(i==0)
	  {
	  I_2h_h(i,i)=1;
	  I_2h_h(i+1,i)=0.5;
	  }
	else if(i==n/2)
	  {
	  I_2h_h(n,i)=1;
	  I_2h_h(n-1,i)=0.5;
	  }
	else
	  {
	    I_2h_h(2*i-1,i)=0.5;
	    I_2h_h(2*i,i)=1;
	    I_2h_h(2*i+1,i)=0.5;
	  }
	}
    }
  else
    {
      ;
    }
  for(int i=1;i<=v1;i++)
    v_h=T*v_h+c;
  if(n==1)
    {
      flag1=true;
      for(int i=1;i<=v2;i++)
	v_h=T*v_h+c;
    }
  else
    {
      n=n/2;
      f_2h=I_h_2h*(f_h-a*v_h);
      v_2h={0};
      v_2h=VC1(obj,v_2h,f_2h,v1,v2);
    }
  if(!flag1)
    {
    v_h=v_h+I_2h_h*v_2h;
    n=2*n;
    }
  if(!flag1)
    {
    for(int i=1;i<=v2;i++)
      v_h=T*v_h+c;
    }
  else
    flag1=false;
    return v_h;
}

Matrix<double,1> FMG1(Multigrid_Solver& obj,Matrix<double,1>& f_h,int v1,int v2)
{
  Matrix<double,1> v_h(n+1);
  Matrix<double,1> f_2h(n/2+1);
  Matrix<double,1> v_2h(n/2+1);
  Matrix<double,2> I_h_2h(n/2+1,n+1);
  Matrix<double,2> I_2h_h(n+1,n/2+1);
  
  //define restriction operator
    if(obj.res_op)
    {
      for(int i=0;i<=n/2;i++)
	{
	if(i==0)
	  I_h_2h(i,i)=1;
	else if(i==n/2)
	  I_h_2h(i,n)=1;
	else
	  {
	    I_h_2h(i,2*i-1)=0.25;
	    I_h_2h(i,2*i)=0.5;
	    I_h_2h(i,2*i+1)=0.25;
	  }
	}	
    }
  else
    {
      for(int i=0;i<=n/2;i++)
	I_h_2h(i,2*i)=1;
    }
  //define interpolation operator
    if(obj.inter_op)
    {
      for(int i=0;i<=n/2;i++)
	{
	if(i==0)
	  {
	  I_2h_h(i,i)=1;
	  I_2h_h(i+1,i)=0.5;
	  }
	else if(i==n/2)
	  {
	  I_2h_h(n,i)=1;
	  I_2h_h(n-1,i)=0.5;
	  }
	else
	  {
	    I_2h_h(2*i-1,i)=0.5;
	    I_2h_h(2*i,i)=1;
	    I_2h_h(2*i+1,i)=0.5;
	  }
	}
    }
  else
    {
      ;
    }
  if(n==1)
    {
      flag2=true;
      v_h={0};
      v_h=VC1(obj,v_h,f_h,v1,v2);
    }
  else
    {
      n=n/2;
      f_2h=I_h_2h*f_h;
      v_2h=FMG1(obj,f_2h,v1,v2);
    }
  if(!flag2)
    {
      n=n*2;
      v_h=I_2h_h*v_2h;
    }
  if(!flag2)
    v_h=VC1(obj,v_h,f_h,v1,v2);
  else
    flag2=false;
  return v_h;
}

Matrix<double,1> VC2(Multigrid_Solver& obj,Matrix<double,1>& v_h,Matrix<double,1>& f_h,int v1,int v2)
{
  double h=1.0/n;
  Matrix<double,1> f_2h((n/2+1)*(n/2+1));
  Matrix<double,1> v_2h((n/2+1)*(n/2+1));
  Matrix<double,2> a((n+1)*(n+1),(n+1)*(n+1));
  Matrix<double,2> T((n+1)*(n+1),(n+1)*(n+1));
  Matrix<double,1> c((n+1)*(n+1));
  Matrix<double,2> I_h_2h((n/2+1)*(n/2+1),(n+1)*(n+1));
  Matrix<double,2> I_2h_h((n+1)*(n+1),(n/2+1)*(n/2+1));
  //define coefficient matrix
  if(obj.boundary_con==Dirichlet)
    {
      int k=0;
      for(int j=0;j<=n;j++)
	for(int i=0;i<=n;i++)
	  {
	    if(i==0||j==0||i==n||j==n)
	      {
		a(k,k)=1;
		k++;
	      }
	    else
	      {
		a(k,k)=4/(h*h);
		a(k,k+1)=-1/(h*h);
		a(k,k-1)=-1/(h*h);
		a(k,k-n-1)=-1/(h*h);
		a(k,k+n+1)=-1/(h*h);
		k++;
	      }
	  }
    }
  else if(obj.boundary_con==Neumann)
    {
      int k=0;
    for(int j=0;j<=n;j++)
    for(int i=0;i<=n;i++)
      {
	a(k,k)=4/(h*h);
	if(j==0&&i!=0&&i!=n)
	  {
	    a(k,k)=4/(h*h);
	    a(k,k+n+1)=-2/(h*h);
	    a(k,k-1)=-1/(h*h);
	    a(k,k+1)=-1/(h*h);
	  }
	else if(i==0&&j!=0&&j!=n)
	  {
	    a(k,k)=4/(h*h);
	    a(k,k+1)=-2/(h*h);
	    a(k,k+n+1)=-1/(h*h);
	    a(k,k-n-1)=-1/(h*h);
	  }
	else if(i==n&&j!=0&&j!=n)
	  {
	    a(k,k-1)=-2/(h*h);
	    a(k,k)=4/(h*h);
	    a(k,k+n+1)=-1/(h*h);
	    a(k,k-n-1)=-1/(h*h);
	  }
	else if(j==n&&i!=0&&i!=n)
	  {
	    a(k,k)=4/(h*h);
	    a(k,k-n-1)=-2/(h*h);
	    a(k,k+1)=-1/(h*h);
	    a(k,k-1)=-1/(h*h);
	  }
	else if((i==0&&j==0)||(i==0)&&(j==n)||(i==n&&j==n)||(i==n&&j==0))
	  {
	    a(k,k)=1;
	  }
	else
	  {
	    a(k,k-1)=-1/(h*h);
	    a(k,k+1)=-1/(h*h);
	    a(k,k+n+1)=-1/(h*h);
	    a(k,k-n-1)=-1/(h*h);
	  }
	k++;
    }
    }
  else
    {
      int k=0;
    for(int j=0;j<=n;j++)
    for(int i=0;i<=n;i++)
      {
	if(j==0&&i!=0&&i!=n)
	  {
	    a(k,k)=(4-2*h*alpha/beta)/(h*h);
	    a(k,k+n+1)=-2/(h*h);
	    a(k,k-1)=-1/(h*h);
	    a(k,k+1)=-1/(h*h);
	  }
	else if(i==0&&j!=0&&j!=n)
	  {
	    a(k,k)=(4-2*h*alpha/beta)/(h*h);
	    a(k,k+1)=-2/(h*h);
	    a(k,k+n+1)=-1/(h*h);
	    a(k,k-n-1)=-1/(h*h);
	  }
	else if(i==n&&j!=0&&j!=n)
	  {
	    a(k,k-1)=-2/(h*h);
	    a(k,k)=(4+2*h*alpha/beta)/(h*h);
	    a(k,k+n+1)=-1/(h*h);
	    a(k,k-n-1)=-1/(h*h);
	  }
	else if(j==n&&i!=0&&i!=n)
	  {
	    a(k,k)=(4+2*h*alpha/beta)/(h*h);
	    a(k,k-n-1)=-2/(h*h);
	    a(k,k+1)=-1/(h*h);
	    a(k,k-1)=-1/(h*h);
	  }
	else if((i==0&&j==0)||(i==0)&&(j==n)||(i==n&&j==n)||(i==n&&j==0))
	  {
	    a(k,k)=1;
	  }
	else
	  {
	    a(k,k)=4/(h*h);
	    a(k,k-1)=-1/(h*h);
	    a(k,k+1)=-1/(h*h);
	    a(k,k+n+1)=-1/(h*h);
	    a(k,k-n-1)=-1/(h*h);
	  }
	k++;
    }
    }
  //cout<<a<<endl;
  //define iteration matrix T and corresponding vector c
  for(int i=0;i<(n+1)*(n+1);i++)
    {
      T(i,i)=1-w;
      if(i==0)
	T(i,i+1)=-w/a(i,i)*a(i,i+1);
      else if(i==(n+1)*(n+1)-1)
	T(i,i-1)=-w/a(i,i)*a(i,i-1);
      else
	{
	  T(i,i+1)=-w/a(i,i)*a(i,i+1);
	  T(i,i-1)=-w/a(i,i)*a(i,i-1);
	}
    }
  for(int i=0;i<(n+1)*(n+1);i++)
     c(i)=w/a(i,i)*f_h(i);
  //define restriction operator
 if(obj.res_op)
    {
      int k=0;
      for(int i=0;i<=n/2;i++)
	for(int j=0;j<=n/2;j++)
	  {
	    if(i==0||j==0||i==n/2||j==n/2)
	      {
	      I_h_2h(k,2*i*(n+1)+2*j)=1;
	      k++;
	      }
	    else
	      {
		I_h_2h(k,2*i*(n+1)+2*j)=0.25;
		I_h_2h(k,2*i*(n+1)+2*j+1)=0.125;
		I_h_2h(k,2*i*(n+1)+2*j-1)=0.125;
		I_h_2h(k,(2*i-1)*(n+1)+2*j)=0.125;
		I_h_2h(k,(2*i+1)*(n+1)+2*j)=0.125;
		I_h_2h(k,(2*i+1)*(n+1)+2*j+1)=0.0625;
		I_h_2h(k,(2*i+1)*(n+1)+2*j-1)=0.0625;
		I_h_2h(k,(2*i-1)*(n+1)+2*j+1)=0.0625;
		I_h_2h(k,(2*i-1)*(n+1)+2*j-1)=0.0625;
		k++;
	      }
	  }
    }
  else
    {
      int k=0;
      for(int i=0;i<=n/2;i++)
	for(int j=0;j<=n/2;j++)
	  {
	    I_h_2h(k,2*i*(n+1)+2*j)=1;
	    k++;
	  }
	  }
  //define interpolation operator
  if(obj.inter_op&&n!=1)
    {
      int k=0;
      for(int i=0;i<=n;i++)
	for(int j=0;j<=n;j++)
	  {
	    if(i%2==1)
	      {
		if(j%2==1)
		  {
		    I_2h_h(k,(i+1)/2*(n/2+1)+(j+1)/2)=0.25;
		    I_2h_h(k,(i+1)/2*(n/2+1)+(j-1)/2)=0.25;
		    I_2h_h(k,(i-1)/2*(n/2+1)+(j+1)/2)=0.25;
		    I_2h_h(k,(i-1)/2*(n/2+1)+(j-1)/2)=0.25;
		    k++;
		  }
		else
		  {
		    I_2h_h(k,(i+1)/2*(n/2+1)+j/2)=0.5;
		    I_2h_h(k,(i-1)/2*(n/2+1)+j/2)=0.5;
		    k++;
		  }
	      }
	    else
	      {
		if(j%2==0)
		  {
		    I_2h_h(k,i/2*(n/2+1)+j/2)=1;
		    k++;
		  }
		else
		  {
		    I_2h_h(k,i/2*(n/2+1)+(j+1)/2)=0.5;
		    I_2h_h(k,i/2*(n/2+1)+(j-1)/2)=0.5;
		    k++;
		  }
	      }
	  }
    }
  else
    {
      ;
      }
  for(int i=1;i<=v1;i++)
    v_h=T*v_h+c;
  if(n==1)
    {
      flag1=true;
      for(int i=1;i<=v2;i++)
	v_h=T*v_h+c;
    }
  else
    {
      n=n/2;
      f_2h=I_h_2h*(f_h-a*v_h);
      v_2h={0};
      v_2h=VC2(obj,v_2h,f_2h,v1,v2);
    }
  if(!flag1)
    {
    v_h=v_h+I_2h_h*v_2h;
    n=2*n;
    }
  if(!flag1)
    {
    for(int i=1;i<=v2;i++)
      v_h=T*v_h+c;
    }
  else
    flag1=false;
  return v_h;
 }


Matrix<double,1> FMG2(Multigrid_Solver& obj,Matrix<double,1>& f_h,int v1,int v2)
{
  Matrix<double,1> v_h((n+1)*(n+1));
  Matrix<double,1> f_2h((n/2+1)*(n/2+1));
  Matrix<double,1> v_2h((n/2+1)*(n/2+1));
  Matrix<double,2> I_h_2h((n/2+1)*(n/2+1),(n+1)*(n+1));
  Matrix<double,2> I_2h_h((n+1)*(n+1),(n/2+1)*(n/2+1));
  //define restriction operator
 if(obj.res_op)
    {
      int k=0;
      for(int i=0;i<=n/2;i++)
	for(int j=0;j<=n/2;j++)
	  {
	    if(i==0||j==0||i==n/2||j==n/2)
	      {
	      I_h_2h(k,2*i*(n+1)+2*j)=1;
	      k++;
	      }
	    else
	      {
		I_h_2h(k,2*i*(n+1)+2*j)=0.25;
		I_h_2h(k,2*i*(n+1)+2*j+1)=0.125;
		I_h_2h(k,2*i*(n+1)+2*j-1)=0.125;
		I_h_2h(k,(2*i-1)*(n+1)+2*j)=0.125;
		I_h_2h(k,(2*i+1)*(n+1)+2*j)=0.125;
		I_h_2h(k,(2*i+1)*(n+1)+2*j+1)=0.0625;
		I_h_2h(k,(2*i+1)*(n+1)+2*j-1)=0.0625;
		I_h_2h(k,(2*i-1)*(n+1)+2*j+1)=0.0625;
		I_h_2h(k,(2*i-1)*(n+1)+2*j-1)=0.0625;
		k++;
	      }
	  }
    }
  else
    {
      int k=0;
      for(int i=0;i<=n/2;i++)
	for(int j=0;j<=n/2;j++)
	  {
	    I_h_2h(k,2*i*(n+1)+2*j)=1;
	    k++;
	  }
	  }
  //define interpolation operator
  if(obj.inter_op&&n!=1)
    {
      int k=0;
      for(int i=0;i<=n;i++)
	for(int j=0;j<=n;j++)
	  {
	    if(i%2==1)
	      {
		if(j%2==1)
		  {
		    I_2h_h(k,(i+1)/2*(n/2+1)+(j+1)/2)=0.25;
		    I_2h_h(k,(i+1)/2*(n/2+1)+(j-1)/2)=0.25;
		    I_2h_h(k,(i-1)/2*(n/2+1)+(j+1)/2)=0.25;
		    I_2h_h(k,(i-1)/2*(n/2+1)+(j-1)/2)=0.25;
		    k++;
		  }
		else
		  {
		    I_2h_h(k,(i+1)/2*(n/2+1)+j/2)=0.5;
		    I_2h_h(k,(i-1)/2*(n/2+1)+j/2)=0.5;
		    k++;
		  }
	      }
	    else
	      {
		if(j%2==0)
		  {
		    I_2h_h(k,i/2*(n/2+1)+j/2)=1;
		    k++;
		  }
		else
		  {
		    I_2h_h(k,i/2*(n/2+1)+(j+1)/2)=0.5;
		    I_2h_h(k,i/2*(n/2+1)+(j-1)/2)=0.5;
		    k++;
		  }
	      }
	  }
    }
  else
    {
      ;
      }
  if(n==1)
    {
      flag2=true;
      v_h={0};
      v_h=VC2(obj,v_h,f_h,v1,v2);
    }
  else
    {
      n=n/2;
      f_2h=I_h_2h*f_h;
      v_2h=FMG2(obj,f_2h,v1,v2);
    }
  if(!flag2)
    {
      n=n*2;
      v_h=I_2h_h*v_2h;
    }
  if(!flag2)
    v_h=VC2(obj,v_h,f_h,v1,v2);
  else
    flag2=false;
  return v_h;
}

void Multigrid_Solver::Solver()
{
  Matrix<double,1> F1(N+1);
  Matrix<double,1> u1_(N+1); //to distinguish from the function u(x),here we use u_ to represent the computed solution from iteration
  Matrix<double,2> a1(N+1,N+1);
  Matrix<double,1> F2((N+1)*(N+1));
  Matrix<double,1> u2_((N+1)*(N+1));
  Matrix<double,2> a2((N+1)*(N+1),(N+1)*(N+1));
  double h=1.0/N;
  int ipiv1[N+1];
  int ipiv2[(N+1)*(N+1)];
  if(dim==1)
    {
   if(boundary_con==Dirichlet)
      for(int i=0;i<=N;i++)
	{
	  if(i==0||i==N)
	    {
	      F1(i)=u(i*h);
	      a1(i,i)=1;
	    }
	  else
	    {
	      F1(i)=f(i*h);
	      a1(i,i)=2/(h*h);
	      a1(i,i-1)=-1/(h*h);
	      a1(i,i+1)=-1/(h*h);
	    }
	}
  else if(boundary_con==Neumann)
    for(int i=0;i<=N;i++)
      {
	a1(i,i)=2/(h*h);
	if(i==0)
	  {
	    F1(i)=u(i*h);
	    a1(i,i)=1;
	  }
	else if(i==N)
	  {
	    a1(i,i-1)=-2/(h*h);
	    F1(i)=f(i*h)+2*ux(i*h)/h;
	  }
	else
	  {
	    a1(i,i+1)=-1/(h*h);
	    a1(i,i-1)=-1/(h*h);
	    F1(i)=f(i*h);
	  }
      }
  else
    {
      for(int i=0;i<=N;i++)
	{
	  a1(i,i)=2/(h*h);
	if(i==0)
	  {
	    a1(i,i+1)=-2/(h*h);
	    F1(i)=f(i*h)-2*ux(i*h)/h;
	  }
	else if(i==N)
	  {
	    a1(i,i)=1;
	    F1(i)=u(i*h);
	  }
	else
	  {
	    a1(i,i+1)=-1/(h*h);
	    a1(i,i-1)=-1/(h*h);
	    F1(i)=f(i*h);
	  }
	}
    }
   double* A1=a1.data();
   double* f1=F1.data();
   if(cycle)
     for(int i=1;i<=max_iteration;i++)
       {
	 u1_=VC1(*this,u1_,F1,v1,v2);
	 cout<<"the residual norm after"<<i<<"th V-cycle is:"<<norm(a1*u1_-F1)<<endl;
       }
   else
     u1_=FMG1(*this,F1,v1,v2);
   LAPACKE_dgesv(LAPACK_ROW_MAJOR,N+1,1,A1,N+1,ipiv1,f1,1);
   cout<<endl;
   cout<<"the maximum norm of error vector is:"<<norm(F1-u1_)<<endl;
    }
  else
    {
      if(boundary_con==Dirichlet)
	{
	  int k=0;
    for(int j=0;j<=n;j++)
    for(int i=0;i<=n;i++)
      {
	double x=i*h;
	double y=j*h;
	if(i==0||j==0||i==n||j==n)
	  {
	    a2(k,k)=1;
	    F2(k)=u(x,y);
	  }
	else
	  {
	    a2(k,k)=4/(h*h);
	    a2(k,k+1)=-1/(h*h);
	    a2(k,k-1)=-1/(h*h);
	    a2(k,k-n-1)=-1/(h*h);
	    a2(k,k+n+1)=-1/(h*h);
	    F2(k)=f(x,y);
	  }
	k++;
      }
	}
      else if(boundary_con==Neumann)
	{
	  int k=0;
    for(int j=0;j<=N;j++)
    for(int i=0;i<=N;i++)
      {
	double x=i*h;
	double y=h*j;
	a2(k,k)=4/(h*h);
	if(j==0&&i!=0&&i!=n)
	  {
	    a2(k,k)=4/(h*h);
	    a2(k,k+n+1)=-2/(h*h);
	    a2(k,k-1)=-1/(h*h);
	    a2(k,k+1)=-1/(h*h);
	    F2(k)=f(x,y)-2*uy(x,y)/h;
	  }
	else if(i==0&&j!=0&&j!=n)
	  {
	    a2(k,k)=4/(h*h);
	    a2(k,k+1)=-2/(h*h);
	    a2(k,k+n+1)=-1/(h*h);
	    a2(k,k-n-1)=-1/(h*h);
	    F2(k)=f(x,y)-2*ux(x,y)/h;
	  }
	else if(i==n&&j!=0&&j!=n)
	  {
	    a2(k,k-1)=-2/(h*h);
	    a2(k,k)=4/(h*h);
	    a2(k,k+n+1)=-1/(h*h);
	    a2(k,k-n-1)=-1/(h*h);
	    F2(k)=f(x,y)+2*ux(x,y)/h;
	  }
	else if(j==n&&i!=0&&i!=n)
	  {
	    a2(k,k)=4/(h*h);
	    a2(k,k-n-1)=-2/(h*h);
	    a2(k,k+1)=-1/(h*h);
	    a2(k,k-1)=-1/(h*h);
	    F2(k)=f(x,y)+2*uy(x,y)/h;
	  }
	else if((i==0&&j==0)||(i==0)&&(j==n)||(i==n&&j==n)||(i==n&&j==0))
	  {
	    F2(k)=u(x,y);
	    a2(k,k)=1;
	  }
	else
	  {
	    a2(k,k-1)=-1/(h*h);
	    a2(k,k+1)=-1/(h*h);
	    a2(k,k+n+1)=-1/(h*h);
	    a2(k,k-n-1)=-1/(h*h);
	    F2(k)=f(x,y);
	  }
	k++;
      }
	}
          else
	{
      	  int k=0;
    for(int j=0;j<=n;j++)
    for(int i=0;i<=n;i++)
      {
	double x=i*h;
	double y=h*j;
	F2(k)=f(x,y);
	if(j==0&&i!=0&&i!=n)
	  {
	    a2(k,k)=(4-2*h*alpha/beta)/(h*h);
	    a2(k,k+n+1)=-2/(h*h);
	    a2(k,k-1)=-1/(h*h);
	    a2(k,k+1)=-1/(h*h);
	    F2(k)-=2*(beta*uy(x,y)+alpha*u(x,y))/(beta*h);
	  }
	else if(i==0&&j!=0&&j!=n)
	  {
	    a2(k,k)=(4-2*h*alpha/beta)/(h*h);
	    a2(k,k+1)=-2/(h*h);
	    a2(k,k+n+1)=-1/(h*h);
	    a2(k,k-n-1)=-1/(h*h);
	    F2(k)-=2*(beta*ux(x,y)+alpha*u(x,y))/(h*beta);
	  }
	else if(i==n&&j!=0&&j!=n)
	  {
	    a2(k,k-1)=-2/(h*h);
	    a2(k,k)=(4+2*h*alpha/beta)/(h*h);
	    a2(k,k+n+1)=-1/(h*h);
	    a2(k,k-n-1)=-1/(h*h);
	    F2(k)+=2*(beta*ux(x,y)+alpha*u(x,y))/(h*beta);
	  }
	else if(j==n&&i!=0&&i!=n)
	  {
	    a2(k,k)=(4+2*h*alpha/beta)/(h*h);
	    a2(k,k-n-1)=-2/(h*h);
	    a2(k,k+1)=-1/(h*h);
	    a2(k,k-1)=-1/(h*h);
	    F2(k)+=2*(beta*uy(x,y)+alpha*u(x,y))/(h*beta);
	  }
	else if((i==0&&j==0)||(i==0)&&(j==n)||(i==n&&j==n)||(i==n&&j==0))
	  {
	    F2(k)=u(x,y);
	    a2(k,k)=1;
	  }
	else
	  {
	    a2(k,k)=4/(h*h);
	    a2(k,k-1)=-1/(h*h);
	    a2(k,k+1)=-1/(h*h);
	    a2(k,k+n+1)=-1/(h*h);
	    a2(k,k-n-1)=-1/(h*h);
	  }
	k++;
      }
	}
      if(cycle)
	{
	  for(int i=1;i<=max_iteration;i++)
	    {
	      u2_=VC2(*this,u2_,F2,v1,v2);
	      cout<<"the maximum norm for the residual after "<<i<<"th V-cycle is:"<<norm(a2*u2_-F2)<<endl;
	    }
	}
      else
      	u2_=FMG2(*this,F2,v1,v2);
      double* A2=a2.data();
      double* f2=F2.data();
      LAPACKE_dgesv(LAPACK_ROW_MAJOR,(N+1)*(N+1),1,A2,(N+1)*(N+1),ipiv2,f2,1);
      cout<<"the maximum norm for the error is:"<<norm(u2_-F2)<<endl;
    }
}

//some auxiliary functions
double norm(const Matrix<double,1> a)
{
  const Index n=a.dim1();
  double norm=abs(a(0));
  for(Index i=1;i<n;i++)
    {
      if(abs(a(i))>norm)
	norm=abs(a(i));
    }
  return norm;
}
Matrix<double,1> operator*(const Matrix<double,2>& A,const Matrix<double,1>& b)
{
  const Index n=A.dim1();
  Matrix<double> c(n);
  for(Index i=0;i<n;i++)
    c(i)=dot_product(A[i],b);
  return c;
}
Matrix<double,1> operator+(const Matrix<double,1>& a,const Matrix<double,1>& b)
{
  const Index n=a.dim1();
  Matrix<double> c(n);
  for(int i=0;i<n;i++)
    c(i)=a(i)+b(i);
  return c;
}
Matrix<double,1> operator-(const Matrix<double,1>& a,const Matrix<double,1>& b)
{
  const Index n=a.dim1();
  Matrix<double> c(n);
  for(int i=0;i<n;i++)
    c(i)=a(i)-b(i);
  return c;
}
Matrix<double,1> operator/(const Matrix<double,1>& a,const Matrix<double,1>& b)
{
  const Index n=a.dim1();
  Matrix<double> c(n);
  for(int i=0;i<n;i++)
    c(i)=abs(a(i)/b(i));
  return c;
}
