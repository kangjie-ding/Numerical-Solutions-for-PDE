#include"BVP.h"
#include"Norm.h"
#include<cmath>

void Square_Dir(const  int n,double (*u)(double x,double y),double (*f)(double x,double y))
{
  const int dim=(n+1)*(n+1);
  const double h=1.0/n;
  double x;
  double y;
  Matrix<double,2> a(dim,dim);
  double* b=new double[dim];
  double* U=new double[dim];
  double* error=new double[dim];
  int ipiv[dim];
  int k=0;
  for(int j=0;j<=n;j++)
    for(int i=0;i<=n;i++)
      {
	x=i*h;
	y=j*h;
	if(i==0||j==0||i==n||j==n)
	  {
	    a(k,k)=h*h;
	    b[k]=u(x,y);
	    U[k]=u(x,y);
	    k++;
	  }
	else
	  {
	    a(k,k)=4;
	    a(k,k+1)=-1;
	    a(k,k-1)=-1;
	    a(k,k-n-1)=-1;
	    a(k,k+n+1)=-1;
	    b[k]=f(x,y);
	    U[k]=u(x,y);
	    k++;
	  }
      }
  
  a *= 1.0/(h*h);
  double* A=a.data();
  LAPACKE_dgesv(LAPACK_ROW_MAJOR,dim,1,A,dim,ipiv,b,1);
  for(int i=0;i<dim;i++)
    cout<<b[i]<<endl;
  for(int i=0;i<dim;i++)
    error[i]=b[i]-U[i];
  cout<<"The norm of error with regular-Dirichlet with n="<<n<<endl;
  cout<<"1-norm:"<<norm1(error,dim)<<endl;
  cout<<"2-norm:"<<norm2(error,dim)<<endl;
  cout<<"infinity-norm:"<<normi(error,dim)<<endl;
  delete []b;
  delete []U;
  delete []error;
}

void Square_Neu(int n,double (*u)(double x,double y),double (*f)(double x,double y),double (*ux)(double x,double y),double (*uy)(double x,double y))
{
  const int dim=(n+1)*(n+1);
  const double h=1.0/n;
  double x;
  double y;
  Matrix<double,2> a(dim,dim);
  double* b=new double[dim];
  double* U=new double[dim];
  double* error=new double[dim];
  int ipiv[dim];
  int k=0;
  for(int j=0;j<=n;j++)
    for(int i=0;i<=n;i++)
      {
	x=i*h;
	y=h*j;
	a(k,k)=4;
	U[k]=u(x,y);
	if(j==0&&i!=0&&i!=n)
	  {
	    a(k,k)=4;
	    a(k,k+n+1)=-2;
	    a(k,k-1)=-1;
	    a(k,k+1)=-1;
	    b[k]=f(x,y)-2*uy(x,y)/h;
	  }
	else if(i==0&&j!=0&&j!=n)
	  {
	    a(k,k)=4;
	    a(k,k+1)=-2;
	    a(k,k+n+1)=-1;
	    a(k,k-n-1)=-1;
	    b[k]=f(x,y)-2*ux(x,y)/h;
	  }
	else if(i==n&&j!=0&&j!=n)
	  {
	    a(k,k-1)=-2;
	    a(k,k)=4;
	    a(k,k+n+1)=-1;
	    a(k,k-n-1)=-1;
	    b[k]=f(x,y)+2*ux(x,y)/h;
	  }
	else if(j==n&&i!=0&&i!=n)
	  {
	    a(k,k)=4;
	    a(k,k-n-1)=-2;
	    a(k,k+1)=-1;
	    a(k,k-1)=-1;
	    b[k]=f(x,y)+2*uy(x,y)/h;
	  }
	else if((i==0&&j==0)||(i==0)&&(j==n)||(i==n&&j==n)||(i==n&&j==0))
	  {
	    b[k]=u(x,y);
	    a(k,k)=h*h;
	  }
	else
	  {
	    a(k,k-1)=-1;
	    a(k,k+1)=-1;
	    a(k,k+n+1)=-1;
	    a(k,k-n-1)=-1;
	    b[k]=f(x,y);
	  }
	k++;
      }
  a *= 1.0/(h*h);
  double* A=a.data();
  LAPACKE_dgesv(LAPACK_ROW_MAJOR,dim,1,A,dim,ipiv,b,1);
  for(int i=0;i<dim;i++)
    error[i]=b[i]-U[i];
  cout<<"The norm of error with regular-Neumann with n="<<n<<endl;
  cout<<"1-norm:"<<norm1(error,dim)<<endl;
  cout<<"2-norm:"<<norm2(error,dim)<<endl;
  cout<<"infinity-norm:"<<normi(error,dim)<<endl;
  delete []b;
  delete []U;
  delete []error;
}

void Square_mixed(int n,double alpha,double beta,double (*u)(double x,double y),double (*f)(double x,double y),double (*ux)(double x,double y),double (*uy)(double x,double y))
{
  const int dim=(n+1)*(n+1);
  const double h=1.0/n;
  double x;
  double y;
  Matrix<double,2> a(dim,dim);
  double* b=new double[dim];
  double* U=new double[dim];
  double* error=new double[dim];
  int ipiv[dim];
  int k=0;
  for(int j=0;j<=n;j++)
    for(int i=0;i<=n;i++)
      {
	x=i*h;
	y=h*j;
	U[k]=u(x,y);
	b[k]=f(x,y);
	if(j==0&&i!=0&&i!=n)
	  {
	    a(k,k)=4-2*h*alpha/beta;
	    a(k,k+n+1)=-2;
	    a(k,k-1)=-1;
	    a(k,k+1)=-1;
	    b[k]-=2*(beta*uy(x,y)+alpha*u(x,y))/(beta*h);
	  }
	else if(i==0&&j!=0&&j!=n)
	  {
	    a(k,k)=4-2*h*alpha/beta;
	    a(k,k+1)=-2;
	    a(k,k+n+1)=-1;
	    a(k,k-n-1)=-1;
	    b[k]-=2*(beta*ux(x,y)+alpha*u(x,y))/(h*beta);
	  }
	else if(i==n&&j!=0&&j!=n)
	  {
	    a(k,k-1)=-2;
	    a(k,k)=4+2*h*alpha/beta;
	    a(k,k+n+1)=-1;
	    a(k,k-n-1)=-1;
	    b[k]+=2*(beta*ux(x,y)+alpha*u(x,y))/(h*beta);
	  }
	else if(j==n&&i!=0&&i!=n)
	  {
	    a(k,k)=4+2*h*alpha/beta;
	    a(k,k-n-1)=-2;
	    a(k,k+1)=-1;
	    a(k,k-1)=-1;
	    b[k]+=2*(beta*uy(x,y)+alpha*u(x,y))/(h*beta);
	  }
	else if((i==0&&j==0)||(i==0)&&(j==n)||(i==n&&j==n)||(i==n&&j==0))
	  {
	    b[k]=u(x,y);
	    a(k,k)=h*h;
	  }
	else
	  {
	    a(k,k)=4;
	    a(k,k-1)=-1;
	    a(k,k+1)=-1;
	    a(k,k+n+1)=-1;
	    a(k,k-n-1)=-1;
	  }
	k++;
      }
  a *= 1.0/(h*h);
  double* A=a.data();
  LAPACKE_dgesv(LAPACK_ROW_MAJOR,dim,1,A,dim,ipiv,b,1);
  for(int i=0;i<dim;i++)
    error[i]=b[i]-U[i];
  cout<<"The norm of error with regular-mixed with n="<<n<<endl;
  cout<<"1-norm:"<<norm1(error,dim)<<endl;
  cout<<"2-norm:"<<norm2(error,dim)<<endl;
  cout<<"infinity-norm:"<<normi(error,dim)<<endl;
  delete []b;
  delete []U;
  delete []error;
}

void Disk_Dir(int n,double* center,double radius,double (*u)(double x,double y),double (*f)(double x,double y))
{
  const int dim=(n+1)*(n+1);
  const double h=1.0/n;
  double x;
  double y;
  Matrix<double,2> a(dim,dim);
  double* b=new double[dim];
  double* U=new double[dim];
  double* error=new double[dim];
  int ipiv[dim];
  int k=0;
  for(int j=0;j<=n;j++)
    for(int i=0;i<=n;i++)
      {
	x=i*h;
	y=j*h;
	if(i==0||j==0||i==n||j==n)
	  {
	    a(k,k)=h*h;
	    b[k]=u(x,y);
	    U[k]=u(x,y);
	    k++;
	  }
	else
	  {
	    if(pow(x-center[0],2)+pow(y-center[1],2)<=pow(radius,2))
	      {
		a(k,k)=h*h;
		b[k]=u(x,y);
		U[k]=u(x,y);
		k++;
	      }
	    else
	      {
		bool right=false,left=false,up=false,below=false;
		b[k]=f(x,y);
		U[k]=u(x,y);
		if(pow(h*(i+1)-center[0],2)+pow(y-center[1],2)<=pow(radius,2))
		  {
		    right=true;
		    double c=(abs(center[0]-x)-sqrt(pow(radius,2)-pow(center[1]-y,2)))/h;
		    a(k,k)+=2/c;
		    a(k,k-1)=-2/(1+c);
		    b[k]+=2.0*u(x+c*h,y)/(c*(1+c)*pow(h,2));
		  }
		 if(pow(h*(i-1)-center[0],2)+pow(y-center[1],2)<=pow(radius,2))
		  {
		    left=true;
		    double c=(abs(center[0]-x)-sqrt(pow(radius,2)-pow(center[1]-y,2)))/h;
		    a(k,k)+=2/c;
		    a(k,k+1)=-2/(1+c);
		    b[k]+=2.0*u(x-c*h,y)/(c*(1+c)*pow(h,2));
		  }
		if(pow(x-center[0],2)+pow(h*(j+1)-center[1],2)<=pow(radius,2))
		  {
		    up=true;
		    double c=(abs(center[1]-y)-sqrt(pow(radius,2)-pow(center[0]-x,2)))/h;
		    a(k,k)+=2/c;
		    a(k,k-n-1)=-2/(1+c);
		    b[k]+=2.0*u(x,y+c*h)/(c*(1+c)*pow(h,2));
		  }
		if(pow(x-center[0],2)+pow(h*(j-1)-center[1],2)<=pow(radius,2))
		  {
		    below=true;
		    double c=(abs(center[1]-y)-sqrt(pow(radius,2)-pow(center[0]-x,2)))/h;
		    a(k,k)+=2/c;
		    a(k,k+n+1)=-2/(1+c);
		    b[k]+=2.0*u(x,y-c*h)/(c*(1+c)*pow(h,2));
		  }
		if(!left&&!right&&!up&&!below)
		  {
		    a(k,k)=4;
		    a(k,k+1)=-1;
		    a(k,k-1)=-1;
		    a(k,k+n+1)=-1;
		    a(k,k-n-1)=-1;
		  }
		if(left||right)
		  if(!up&&!below)
		    {
		      a(k,k)+=2;
		      a(k,k-n-1)=-1;
		      a(k,k+n+1)=-1;
		    }
		if(up||below)
		  if(!left&&!right)
		    {
		      a(k,k)+=2;
		      a(k,k-1)=-1;
		      a(k,k+1)=-1;
		    }
		k++;
	      }
	  }
      }
  a *= 1.0/(h*h);
  double* A=a.data();
  LAPACKE_dgesv(LAPACK_ROW_MAJOR,dim,1,A,dim,ipiv,b,1);
  for(int i=0;i<dim;i++)
    error[i]=b[i]-U[i];
  cout<<"The norm of error with irregular_Dirichlet with n="<<n<<endl;
  cout<<"1-norm:"<<norm1(error,dim)<<endl;
  cout<<"2-norm:"<<norm2(error,dim)<<endl;
  cout<<"infinity-norm:"<<normi(error,dim)<<endl;
  delete []b;
  delete []U;
  delete []error;
}

void Disk_Neu(int n,double* center,double radius,double (*u)(double x,double y),double (*f)(double x,double y),double (*ux)(double x,double y),double (*uy)(double x,double y))
{
  const int dim=(n+1)*(n+1);
  const double h=1.0/n;
  double x;
  double y;
  Matrix<double,2> a(dim,dim);
  double* b=new double[dim];
  double* U=new double[dim];
  double* error=new double[dim];
  int ipiv[dim];
  int k=0;
  for(int j=0;j<=n;j++)
    for(int i=0;i<=n;i++)
      {
	x=i*h;
	y=j*h;
	U[k]=u(x,y);
	if(i==0&&j!=0&&j!=n)
	  {
	    a(k,k)=4;
	    a(k,k+1)=-2;
	    a(k,k+n+1)=-1;
	    a(k,k-n-1)=-1;
	    b[k]=f(x,y)-2*ux(x,y)/h;
	    k++;
	  }
	else if(i==n&&j!=0&&j!=n)
	  {
	    a(k,k-1)=-2;
	    a(k,k)=4;
	    a(k,k+n+1)=-1;
	    a(k,k-n-1)=-1;
	    b[k]=f(x,y)+2*ux(x,y)/h;
	    k++;
	  }
	else if(j==0&&i!=0&&i!=n)
	  {
	    a(k,k)=4;
	    a(k,k+n+1)=-2;
	    a(k,k-1)=-1;
	    a(k,k+1)=-1;
	    b[k]=f(x,y)-2*uy(x,y)/h;
	    k++;
	  }
	else if(j==n&&i!=0&&i!=n)
	  {
	    a(k,k)=4;
	    a(k,k-n-1)=-2;
	    a(k,k+1)=-1;
	    a(k,k-1)=-1;
	    b[k]=f(x,y)+2*uy(x,y)/h;
	    k++;
	  }
	else if((i==0&&j==0)||(i==0&&j==n)||(j==0&&i==n)||(i==n&&j==n))
	  {
	    a(k,k)=h*h;
	    b[k]=u(x,y);
	    k++;
	  }
	else
	  {
	    if(pow(x-center[0],2)+pow(y-center[1],2)<=pow(radius,2))
	      {
		a(k,k)=h*h;
		b[k]=u(x,y);
		k++;
	      }
	    else
	      {
		bool right=false,left=false,up=false,below=false;
		b[k]=f(x,y);
		if(pow(h*(i+1)-center[0],2)+pow(y-center[1],2)<=pow(radius,2))
		  {
		    right=true;
		    double c=(abs(center[0]-x)-sqrt(pow(radius,2)-pow(center[1]-y,2)))/h;
		    a(k,k)+=2/c;
		    a(k,k-1)=-2/(1+c);
		    b[k]+=2.0*u(x+c*h,y)/(c*(1+c)*pow(h,2));
		  }
		 if(pow(h*(i-1)-center[0],2)+pow(y-center[1],2)<=pow(radius,2))
		  {
		    left=true;
		    double c=(abs(center[0]-x)-sqrt(pow(radius,2)-pow(center[1]-y,2)))/h;
		    a(k,k)+=2/c;
		    a(k,k+1)=-2/(1+c);
		    b[k]+=2.0*u(x-c*h,y)/(c*(1+c)*pow(h,2));
		  }
		if(pow(x-center[0],2)+pow(h*(j+1)-center[1],2)<=pow(radius,2))
		  {
		    up=true;
		    double c=(abs(center[1]-y)-sqrt(pow(radius,2)-pow(center[0]-x,2)))/h;
		    a(k,k)+=2/c;
		    a(k,k-n-1)=-2/(1+c);
		    b[k]+=2.0*u(x,y+c*h)/(c*(1+c)*pow(h,2));
		  }
		if(pow(x-center[0],2)+pow(h*(j-1)-center[1],2)<=pow(radius,2))
		  {
		    below=true;
		    double c=(abs(center[1]-y)-sqrt(pow(radius,2)-pow(center[0]-x,2)))/h;
		    a(k,k)+=2/c;
		    a(k,k+n+1)=-2/(1+c);
		    b[k]+=2.0*u(x,y-c*h)/(c*(1+c)*pow(h,2));
		  }
		if(!left&&!right&&!up&&!below)
		  {
		    a(k,k)=4;
		    a(k,k+1)=-1;
		    a(k,k-1)=-1;
		    a(k,k+n+1)=-1;
		    a(k,k-n-1)=-1;
		  }
		if(left||right)
		  if(!up&&!below)
		    {
		      a(k,k)+=2;
		      a(k,k-n-1)=-1;
		      a(k,k+n+1)=-1;
		    }
		if(up||below)
		  if(!left&&!right)
		    {
		      a(k,k)+=2;
		      a(k,k-1)=-1;
		      a(k,k+1)=-1;
		    }
		k++;
	      }
	  }
      }
  a *= 1.0/(h*h);
  double* A=a.data();
  LAPACKE_dgesv(LAPACK_ROW_MAJOR,dim,1,A,dim,ipiv,b,1);
  for(int i=0;i<dim;i++)
    error[i]=b[i]-U[i];
  cout<<"The norm of error with irregular_Neumann with n="<<n<<endl;
  cout<<"1-norm:"<<norm1(error,dim)<<endl;
  cout<<"2-norm:"<<norm2(error,dim)<<endl;
  cout<<"infinity-norm:"<<normi(error,dim)<<endl;
  delete []b;
  delete []U;
  delete []error;
}

void Disk_mixed(int n,double alpha,double beta,double* center,double radius,double (*u)(double x,double y),double (*f)(double x,double y),double (*ux)(double x,double y),double (*uy)(double x,double y))
{
  const int dim=(n+1)*(n+1);
  const double h=1.0/n;
  double x;
  double y;
  Matrix<double,2> a(dim,dim);
  double* b=new double[dim];
  double* U=new double[dim];
  double* error=new double[dim];
  int ipiv[dim];
  int k=0;
  for(int j=0;j<=n;j++)
    for(int i=0;i<=n;i++)
      {
	x=i*h;
	y=j*h;
	U[k]=u(x,y);
	if(i==0&&j!=0&&j!=n)
	  {
	    a(k,k)=4-2*h*alpha/beta;
	    a(k,k+1)=-2;
	    a(k,k+n+1)=-1;
	    a(k,k-n-1)=-1;
	    b[k]=f(x,y)-2*(beta*ux(x,y)+alpha*u(x,y))/(h*beta);
	    k++;
	  }
	else if(i==n&&j!=0&&j!=n)
	  {
	    a(k,k-1)=-2;
	    a(k,k)=4+2*h*alpha/beta;
	    a(k,k+n+1)=-1;
	    a(k,k-n-1)=-1;
	    b[k]=f(x,y)+2*(beta*ux(x,y)+alpha*u(x,y))/(h*beta);
	    k++;
	  }
	else if(j==0&&i!=0&&i!=n)
	  {
	    a(k,k)=4-2*h*alpha/beta;
	    a(k,k+n+1)=-2;
	    a(k,k-1)=-1;
	    a(k,k+1)=-1;
	    b[k]=f(x,y)-2*(beta*uy(x,y)+alpha*u(x,y))/(beta*h);
	    k++;
	  }
	else if(j==n&&i!=0&&i!=n)
	  {
	    a(k,k)=4+2*h*alpha/beta;
	    a(k,k-n-1)=-2;
	    a(k,k+1)=-1;
	    a(k,k-1)=-1;
	    b[k]=f(x,y)+2*(beta*uy(x,y)+alpha*u(x,y))/(h*beta);
	    k++;
	  }
	else if((i==0&&j==0)||(i==0&&j==n)||(j==0&&i==n)||(i==n&&j==n))
	  {
	    a(k,k)=h*h;
	    b[k]=u(x,y);
	    k++;
	  }
	else
	  {
	    if(pow(x-center[0],2)+pow(y-center[1],2)<=pow(radius,2))
	      {
		a(k,k)=h*h;
		b[k]=u(x,y);
		k++;
	      }
	    else
	      {
		bool right=false,left=false,up=false,below=false;
		b[k]=f(x,y);
		if(pow(h*(i+1)-center[0],2)+pow(y-center[1],2)<=pow(radius,2))
		  {
		    right=true;
		    double c=(abs(center[0]-x)-sqrt(pow(radius,2)-pow(center[1]-y,2)))/h;
		    a(k,k)+=2/c;
		    a(k,k-1)=-2/(1+c);
		    b[k]+=2.0*u(x+c*h,y)/(c*(1+c)*pow(h,2));
		  }
		 if(pow(h*(i-1)-center[0],2)+pow(y-center[1],2)<=pow(radius,2))
		  {
		    left=true;
		    double c=(abs(center[0]-x)-sqrt(pow(radius,2)-pow(center[1]-y,2)))/h;
		    a(k,k)+=2/c;
		    a(k,k+1)=-2/(1+c);
		    b[k]+=2.0*u(x-c*h,y)/(c*(1+c)*pow(h,2));
		  }
		if(pow(x-center[0],2)+pow(h*(j+1)-center[1],2)<=pow(radius,2))
		  {
		    up=true;
		    double c=(abs(center[1]-y)-sqrt(pow(radius,2)-pow(center[0]-x,2)))/h;
		    a(k,k)+=2/c;
		    a(k,k-n-1)=-2/(1+c);
		    b[k]+=2.0*u(x,y+c*h)/(c*(1+c)*pow(h,2));
		  }
		if(pow(x-center[0],2)+pow(h*(j-1)-center[1],2)<=pow(radius,2))
		  {
		    below=true;
		    double c=(abs(center[1]-y)-sqrt(pow(radius,2)-pow(center[0]-x,2)))/h;
		    a(k,k)+=2/c;
		    a(k,k+n+1)=-2/(1+c);
		    b[k]+=2.0*u(x,y-c*h)/(c*(1+c)*pow(h,2));
		  }
		if(!left&&!right&&!up&&!below)
		  {
		    a(k,k)=4;
		    a(k,k+1)=-1;
		    a(k,k-1)=-1;
		    a(k,k+n+1)=-1;
		    a(k,k-n-1)=-1;
		  }
		if(left||right)
		  if(!up&&!below)
		    {
		      a(k,k)+=2;
		      a(k,k-n-1)=-1;
		      a(k,k+n+1)=-1;
		    }
		if(up||below)
		  if(!left&&!right)
		    {
		      a(k,k)+=2;
		      a(k,k-1)=-1;
		      a(k,k+1)=-1;
		    }
		k++;
	      }
	  }
      }
  a *= 1.0/(h*h);
  double* A=a.data();
  LAPACKE_dgesv(LAPACK_ROW_MAJOR,dim,1,A,dim,ipiv,b,1);
  for(int i=0;i<dim;i++)
    error[i]=b[i]-U[i];
  cout<<"The norm of error with irregular_mixed with n="<<n<<endl;
  cout<<"1-norm:"<<norm1(error,dim)<<endl;
  cout<<"2-norm:"<<norm2(error,dim)<<endl;
  cout<<"infinity-norm:"<<normi(error,dim)<<endl;
  delete []b;
  delete []U;
  delete []error;
}
bool Valid_Check(int n,double* center,double radius)
{
  const double h=1.0/n;
  double x;
  double y;
  int valid_grid=0;
  for(int i=1;i<n;i++)
    for(int j=1;j<n;j++)
      {
	x=i*h;
	y=i*h;
	if((pow(x-center[0],2)+pow(y-center[1],2))>pow(radius,2))
	  valid_grid++;
      }
  if(valid_grid<4)
     return false;
  return true;
}

