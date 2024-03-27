#include"generator.h"

double result0[4];

double* f0(const double x1,const double x2,const double x4,const double x5)
{
  result0[0]=x4;
  result0[1]=x5;
  result0[2]=2*x5+x1-mu*(x1+mu-1)/pow(pow(x2,2)+pow(x1+mu-1,2),1.5)-(1-mu)*(x1+mu)/pow(pow(x2,2)+pow(x1+mu,2),1.5);
  result0[3]=-2*x4+x2-mu*x2/pow(pow(x2,2)+pow(x1+mu-1,2),1.5)-(1-mu)*x2/pow(pow(x2,2)+pow(x1+mu,2),1.5);
  return result0;
}

void ClassicalRK_Method::Solver()
{
  double u1[steps+1];
  double u2[steps+1];
  double u4[steps+1];
  double u5[steps+1];
  const double step_size=period/steps;
  
  u1[0]=initial[0];u2[0]=initial[1];u4[0]=initial[2];u5[0]=initial[3];
  for(int i=1;i<=steps;i++)
    {
      double* y=f0(u1[i-1],u2[i-1],u4[i-1],u5[i-1]);
      u1[i]=u1[i-1]+step_size*y[0]/6;u2[i]=u2[i-1]+step_size*y[1]/6;
      u4[i]=u4[i-1]+step_size*y[2]/6;u5[i]=u5[i-1]+step_size*y[3]/6;
      
      y=f0(u1[i-1]+0.5*step_size*y[0],u2[i-1]+0.5*step_size*y[1],u4[i-1]+0.5*step_size*y[2],u5[i-1]+0.5*step_size*y[3]);
      u1[i]=u1[i]+step_size*y[0]/3;u2[i]=u2[i]+step_size*y[1]/3;
      u4[i]=u4[i]+step_size*y[2]/3;u5[i]=u5[i]+step_size*y[3]/3;
      
      y=f0(u1[i-1]+0.5*step_size*y[0],u2[i-1]+0.5*step_size*y[1],u4[i-1]+0.5*step_size*y[2],u5[i-1]+0.5*step_size*y[3]);
      u1[i]=u1[i]+step_size*y[0]/3;u2[i]=u2[i]+step_size*y[1]/3;
      u4[i]=u4[i]+step_size*y[2]/3;u5[i]=u5[i]+step_size*y[3]/3;

      y=f0(u1[i-1]+step_size*y[0],u2[i-1]+step_size*y[1],u4[i-1]+step_size*y[2],u5[i-1]+step_size*y[3]);
      u1[i]=u1[i]+step_size*y[0]/6;u2[i]=u2[i]+step_size*y[1]/6;
      u4[i]=u4[i]+step_size*y[2]/6;u5[i]=u5[i]+step_size*y[3]/6;
    }
  /*
  for(int i=0;i<=steps;i++)
    cout<<u1[i]<<' '<<u2[i]<<endl;
  */
  cout<<"the max-norm of solution error after "<<steps<<" steps is:";
  double temp1=abs(u1[steps]-u1[0]),temp2=abs(u2[steps]-u2[0]);
  if(temp1>temp2)
    cout<<temp1<<endl;
  else
    cout<<temp2<<endl;
  
}
