#include"generator.h"

double result2[4];

double* f2(const double x1,const double x2,const double x4,const double x5)
{
  result2[0]=x4;
  result2[1]=x5;
  result2[2]=2*x5+x1-mu*(x1+mu-1)/pow(pow(x2,2)+pow(x1+mu-1,2),1.5)-(1-mu)*(x1+mu)/pow(pow(x2,2)+pow(x1+mu,2),1.5);
  result2[3]=-2*x4+x2-mu*x2/pow(pow(x2,2)+pow(x1+mu-1,2),1.5)-(1-mu)*x2/pow(pow(x2,2)+pow(x1+mu,2),1.5);
  return result2;
}

void DoPrRK_Method::Solver()
{
  double u1[steps+1];
  double u2[steps+1];
  double u4[steps+1];
  double u5[steps+1];
  const double step_size=period/steps;
  
  u1[0]=initial[0];u2[0]=initial[1];u4[0]=initial[2];u5[0]=initial[3];
  for(int i=1;i<=steps;i++)
    {
      double* y=f2(u1[i-1],u2[i-1],u4[i-1],u5[i-1]);
      double y1[4]={y[0],y[1],y[2],y[3]};

      y=f2(u1[i-1]+0.2*step_size*y1[0],u2[i-1]+0.2*step_size*y1[1],u4[i-1]+0.2*step_size*y1[2],u5[i-1]+0.2*step_size*y1[3]);
      double y2[4]={y[0],y[1],y[2],y[3]};

      y=f2(u1[i-1]+3.0/40*step_size*y1[0]+9.0/40*step_size*y2[0],u2[i-1]+3.0/40*step_size*y1[1]+9.0/40*step_size*y2[1],u4[i-1]+3.0/40*step_size*y1[2]+9.0/40*step_size*y2[2],u5[i-1]+3.0/40*step_size*y1[3]+9.0/40*step_size*y2[3]);
      double y3[4]={y[0],y[1],y[2],y[3]};

      y=f2(u1[i-1]+44.0/45*step_size*y1[0]-56.0/15*step_size*y2[0]+32.0/9*step_size*y3[0],u2[i-1]+44.0/45*step_size*y1[1]-56.0/15*step_size*y2[1]+32.0/9*step_size*y3[1],u4[i-1]+44.0/45*step_size*y1[2]-56.0/15*step_size*y2[2]+32.0/9*step_size*y3[2],u5[i-1]+44.0/45*step_size*y1[3]-56.0/15*step_size*y2[3]+32.0/9*step_size*y3[3]);
      double y4[4]={y[0],y[1],y[2],y[3]};

      y=f2(u1[i-1]+19372.0/6561*step_size*y1[0]-25360.0/2187*step_size*y2[0]+64448.0/6561*step_size*y3[0]-212.0/729*step_size*y4[0],u2[i-1]+19372.0/6561*step_size*y1[1]-25360.0/2187*step_size*y2[1]+64448.0/6561*step_size*y3[1]-212.0/729*step_size*y4[1],u4[i-1]+19372.0/6561*step_size*y1[2]-25360.0/2187*step_size*y2[2]+64448.0/6561*step_size*y3[2]-212.0/729*step_size*y4[2],u5[i-1]+19372.0/6561*step_size*y1[3]-25360.0/2187*step_size*y2[3]+64448.0/6561*step_size*y3[3]-212.0/729*step_size*y4[3]);
      double y5[4]={y[0],y[1],y[2],y[3]};

      y=f2(u1[i-1]+9017.0/3168*step_size*y1[0]-355.0/33*step_size*y2[0]+46732.0/5247*step_size*y3[0]+49.0/176*step_size*y4[0]-5103.0/18656*step_size*y5[0],u2[i-1]+9017.0/3168*step_size*y1[1]-355.0/33*step_size*y2[1]+46732.0/5247*step_size*y3[1]+49.0/176*step_size*y4[1]-5103.0/18656*step_size*y5[1],u4[i-1]+9017.0/3168*step_size*y1[2]-355.0/33*step_size*y2[2]+46732.0/5247*step_size*y3[2]+49.0/176*step_size*y4[2]-5103.0/18656*step_size*y5[2],u5[i-1]+9017.0/3168*step_size*y1[3]-355.0/33*step_size*y2[3]+46732.0/5247*step_size*y3[3]+49.0/176*step_size*y4[3]-5103.0/18656*step_size*y5[3]);
      double y6[4]={y[0],y[1],y[2],y[3]};

      y=f2(u1[i-1]+35.0/384*step_size*y1[0]+500.0/1113*step_size*y3[0]+125.0/192*step_size*y4[0]-2187.0/6784*step_size*y5[0]+11.0/84*step_size*y6[0],u2[i-1]+35.0/384*step_size*y1[1]+500.0/1113*step_size*y3[1]+125.0/192*step_size*y4[1]-2187.0/6784*step_size*y5[1]+11.0/84*step_size*y6[1],u4[i-1]+35.0/384*step_size*y1[2]+500.0/1113*step_size*y3[2]+125.0/192*step_size*y4[2]-2187.0/6784*step_size*y5[2]+11.0/84*step_size*y6[2],u5[i-1]+35.0/384*step_size*y1[3]+500.0/1113*step_size*y3[3]+125.0/192*step_size*y4[3]-2187.0/6784*step_size*y5[3]+11.0/84*step_size*y6[3]);
      double y7[4]={y[0],y[1],y[2],y[3]};

      if(order==5)
	{
	  u1[i]=u1[i-1]+35.0/383*step_size*y1[0]+500.0/1113*step_size*y3[0]+125.0/192*step_size*y4[0]-2187.0/6784*step_size*y5[0]+11.0/84*step_size*y6[0];
	  u2[i]=u2[i-1]+35.0/383*step_size*y1[1]+500.0/1113*step_size*y3[1]+125.0/192*step_size*y4[1]-2187.0/6784*step_size*y5[1]+11.0/84*step_size*y6[1];
	  u4[i]=u4[i-1]+35.0/383*step_size*y1[2]+500.0/1113*step_size*y3[2]+125.0/192*step_size*y4[2]-2187.0/6784*step_size*y5[2]+11.0/84*step_size*y6[2];
	  u5[i]=u5[i-1]+35.0/383*step_size*y1[3]+500.0/1113*step_size*y3[3]+125.0/192*step_size*y4[3]-2187.0/6784*step_size*y5[3]+11.0/84*step_size*y6[3];
	}
      else if(order==4)
	{
	  u1[i]=u1[i-1]+5179.0/57600*step_size*y1[0]+7571.0/16695*step_size*y3[0]+393.0/640*step_size*y4[0]-92097.0/339200*step_size*y5[0]+187.0/2100*step_size*y6[0]+1.0/40*step_size*y7[0];
	  u2[i]=u2[i-1]+5179.0/57600*step_size*y1[1]+7571.0/16695*step_size*y3[1]+393.0/640*step_size*y4[1]-92097.0/339200*step_size*y5[1]+187.0/2100*step_size*y6[1]+1.0/40*step_size*y7[1];
	  u4[i]=u4[i-1]+5179.0/57600*step_size*y1[2]+7571.0/16695*step_size*y3[2]+393.0/640*step_size*y4[2]-92097.0/339200*step_size*y5[2]+187.0/2100*step_size*y6[2]+1.0/40*step_size*y7[2];
	  u5[i]=u5[i-1]+5179.0/57600*step_size*y1[3]+7571.0/16695*step_size*y3[3]+393.0/640*step_size*y4[3]-92097.0/339200*step_size*y5[3]+187.0/2100*step_size*y6[3]+1.0/40*step_size*y7[3];
	}
      else
	;
    }
  /*
  if(order==4||order==5)
    {
      for(int i=0;i<=steps;i++)
	cout<<u1[i]<<' '<<u2[i]<<endl;
    }
  else
    cerr<<"Invalid input order!!"<<endl;
  */
  cout<<"the max-norm of solution error after "<<steps<<" steps is:";
  double temp1=abs(u1[steps]-u1[0]),temp2=abs(u2[steps]-u2[0]);
  if(temp1>temp2)
    cout<<temp1<<endl;
  else
    cout<<temp2<<endl;
  
}
