#include"generator.h"

double result[4];

double* f(const double x1,const double x2,const double x4,const double x5)
{
  result[0]=x4;
  result[1]=x5;
  result[2]=2*x5+x1-mu*(x1+mu-1)/pow(pow(x2,2)+pow(x1+mu-1,2),1.5)-(1-mu)*(x1+mu)/pow(pow(x2,2)+pow(x1+mu,2),1.5);
  result[3]=-2*x4+x2-mu*x2/pow(pow(x2,2)+pow(x1+mu-1,2),1.5)-(1-mu)*x2/pow(pow(x2,2)+pow(x1+mu,2),1.5);
  return result;
}

void AdBa_Method::Solver()
{
  double u1[steps+1];
  double u2[steps+1];
  double u4[steps+1];
  double u5[steps+1];
  const double step_size=period/steps;
  
  u1[0]=initial[0];u2[0]=initial[1];u4[0]=initial[2];u5[0]=initial[3];
  for(int i=p;i<=steps;i++)
    {
      if(p==1)
	{
	  double* y=f(u1[i-1],u2[i-1],u4[i-1],u5[i-1]);
	  u1[i]=u1[i-1]+step_size*y[0];
	  u2[i]=u2[i-1]+step_size*y[1];
	  u4[i]=u4[i-1]+step_size*y[2];
	  u5[i]=u5[i-1]+step_size*y[3];
	}
      else if(p==2)
	{
	  double* y1=f(u1[i-1],u2[i-1],u4[i-1],u5[i-1]);
	  double* y2=f(u1[i-2],u2[i-2],u4[i-2],u5[i-2]);
	  u1[i]=u1[i-1]+1.5*step_size*y1[0]-0.5*step_size*y2[0];
	  u2[i]=u2[i-1]+1.5*step_size*y1[1]-0.5*step_size*y2[1];
	  u4[i]=u4[i-1]+1.5*step_size*y1[2]-0.5*step_size*y2[2];
	  u5[i]=u5[i-1]+1.5*step_size*y1[3]-0.5*step_size*y2[3];
	}
      else if(p==3)
	{
	  double* y1=f(u1[i-1],u2[i-1],u4[i-1],u5[i-1]);
	  double* y2=f(u1[i-2],u2[i-2],u4[i-2],u5[i-2]);
	  double* y3=f(u1[i-3],u2[i-3],u4[i-3],u5[i-3]);
	  u1[i]=u1[i-1]+23.0/12*step_size*y1[0]-16.0/12*step_size*y2[0]+5.0/12*step_size*y3[0];
	  u2[i]=u2[i-1]+23.0/12*step_size*y1[1]-16.0/12*step_size*y2[1]+5.0/12*step_size*y3[1];
	  u4[i]=u4[i-1]+23.0/12*step_size*y1[2]-16.0/12*step_size*y2[2]+5.0/12*step_size*y3[2];
	  u5[i]=u5[i-1]+23.0/12*step_size*y1[3]-16.0/12*step_size*y2[3]+5.0/12*step_size*y3[3];
	}
      else
	{
	  double* y1=f(u1[i-1],u2[i-1],u4[i-1],u5[i-1]);
	  double* y2=f(u1[i-2],u2[i-2],u4[i-2],u5[i-2]);
	  double* y3=f(u1[i-3],u2[i-3],u4[i-3],u5[i-3]);
	  double* y4=f(u1[i-4],u2[i-4],u4[i-4],u5[i-4]);
	  u1[i]=u1[i-1]+55.0/24*step_size*y1[0]-59.0/24*step_size*y2[0]+37.0/24*step_size*y3[0]-9.0/24*step_size*y4[0];
	  u2[i]=u2[i-1]+55.0/24*step_size*y1[1]-59.0/24*step_size*y2[1]+37.0/24*step_size*y3[1]-9.0/24*step_size*y4[1];
	  u4[i]=u4[i-1]+55.0/24*step_size*y1[2]-59.0/24*step_size*y2[2]+37.0/24*step_size*y3[2]-9.0/24*step_size*y4[2];
	}
    }
  /*
  for(int i=0;i<=steps;i++)
    cout<<u1[i]<<' '<<u2[i]<<endl;
  */
  double temp1=abs(u1[steps]-u1[0]),temp2=abs(u2[steps]-u2[0]);
  cout<<"the max-norm of solution error after "<<steps<<" steps is:";
  if(temp1>temp2)
    cout<<temp1<<endl;
  else
    cout<<temp2<<endl;
  
}
