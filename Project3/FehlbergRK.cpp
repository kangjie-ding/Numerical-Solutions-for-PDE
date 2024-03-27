#include"generator.h"

double result1[4];

double* f1(const double x1,const double x2,const double x4,const double x5)
{
  result1[0]=x4;
  result1[1]=x5;
  result1[2]=2*x5+x1-mu*(x1+mu-1)/pow(pow(x2,2)+pow(x1+mu-1,2),1.5)-(1-mu)*(x1+mu)/pow(pow(x2,2)+pow(x1+mu,2),1.5);
  result1[3]=-2*x4+x2-mu*x2/pow(pow(x2,2)+pow(x1+mu-1,2),1.5)-(1-mu)*x2/pow(pow(x2,2)+pow(x1+mu,2),1.5);
  return result1;
}

void FehlbergRK_Method::Solver()
{
  double u1[steps+1];
  double u2[steps+1];
  double u4[steps+1];
  double u5[steps+1];
  const double step_size=period/steps;
  
  u1[0]=initial[0];u2[0]=initial[1];u4[0]=initial[2];u5[0]=initial[3];
  for(int i=1;i<=steps;i++)
    {
      double* y=f1(u1[i-1],u2[i-1],u4[i-1],u5[i-1]);
      double y1[4]={y[0],y[1],y[2],y[3]};
      
      y=f1(u1[i-1]+0.25*step_size*y1[0],u2[i-1]+0.25*step_size*y1[1],u4[i-1]+0.25*step_size*y1[2],u5[i-1]+0.25*step_size*y1[3]);
      double y2[4]={y[0],y[1],y[2],y[3]};
      
      y=f1(u1[i-1]+3.0/32*step_size*y1[0]+9.0/32*step_size*y2[0],u2[i-1]+3.0/32*step_size*y1[1]+9.0/32*step_size*y2[1],u4[i-1]+3.0/32*step_size*y1[2]+9.0/32*step_size*y2[2],u5[i-1]+3.0/32*step_size*y1[3]+9.0/32*step_size*y2[3]);
      double y3[4]={y[0],y[1],y[2],y[3]};
      
      y=f1(u1[i-1]+1932.0/2197*step_size*y1[0]-7200.0/2197*step_size*y2[0]+7296.0/2197*step_size*y3[0],u2[i-1]+1932.0/2197*step_size*y1[1]-7200.0/2197*step_size*y2[1]+7296.0/2197*step_size*y3[1],u4[i-1]+1932.0/2197*step_size*y1[2]-7200.0/2197*step_size*y2[2]+7296.0/2197*step_size*y3[2],u5[i-1]+1932.0/2197*step_size*y1[3]-7200.0/2197*step_size*y2[3]+7296.0/2197*step_size*y3[3]);
      double y4[4]={y[0],y[1],y[2],y[3]};

      y=f1(u1[i-1]+439.0/216*step_size*y1[0]-8*step_size*y2[0]+3680.0/513*step_size*y3[0]-845.0/4104*step_size*y4[0],u2[i-1]+439.0/216*step_size*y1[1]-8*step_size*y2[1]+3680.0/513*step_size*y3[1]-845.0/4104*step_size*y4[1],u4[i-1]+439.0/216*step_size*y1[2]-8*step_size*y2[2]+3680.0/513*step_size*y3[2]-845.0/4104*step_size*y4[2],u5[i-1]+439.0/216*step_size*y1[3]-8*step_size*y2[3]+3680.0/513*step_size*y3[3]-845.0/4104*step_size*y4[3]);
      double y5[4]={y[0],y[1],y[2],y[3]};

      y=f1(u1[i-1]-8.0/27*step_size*y1[0]+2*step_size*y2[0]-3544.0/2565*step_size*y3[0]+1859.0/4104*step_size*y4[0]-11.0/40*step_size*y5[0],u2[i-1]-8.0/27*step_size*y1[1]+2*step_size*y2[1]-3544.0/2565*step_size*y3[1]+1859.0/4104*step_size*y4[1]-11.0/40*step_size*y5[1],u4[i-1]-8.0/27*step_size*y1[2]+2*step_size*y2[2]-3544.0/2565*step_size*y3[2]+1859.0/4104*step_size*y4[2]-11.0/40*step_size*y5[2],u5[i-1]-8.0/27*step_size*y1[3]+2*step_size*y2[3]-3544.0/2565*step_size*y3[3]+1859.0/4104*step_size*y4[3]-11.0/40*step_size*y5[3]);
      double y6[4]={y[0],y[1],y[2],y[3]};

      if(order==4)
	{
	  u1[i]=u1[i-1]+25.0/216*step_size*y1[0]+1408.0/2565*step_size*y3[0]+2197.0/4104*step_size*y4[0]-0.2*step_size*y5[0];
	  u2[i]=u2[i-1]+25.0/216*step_size*y1[1]+1408.0/2565*step_size*y3[1]+2197.0/4104*step_size*y4[1]-0.2*step_size*y5[1];
	  u4[i]=u4[i-1]+25.0/216*step_size*y1[2]+1408.0/2565*step_size*y3[2]+2197.0/4104*step_size*y4[2]-0.2*step_size*y5[2];
	  u5[i]=u5[i-1]+25.0/216*step_size*y1[3]+1408.0/2565*step_size*y3[3]+2197.0/4104*step_size*y4[3]-0.2*step_size*y5[3];
	}
      else if(order==5)
	{
	  u1[i]=u1[i-1]+16.0/135*step_size*y1[0]+6656.0/12825*step_size*y3[0]+28561.0/56430*step_size*y4[0]-9.0/50*step_size*y5[0]+2.0/55*step_size*y6[0];
	  u2[i]=u2[i-1]+16.0/135*step_size*y1[1]+6656.0/12825*step_size*y3[1]+28561.0/56430*step_size*y4[1]-9.0/50*step_size*y5[1]+2.0/55*step_size*y6[1];
	  u4[i]=u4[i-1]+16.0/135*step_size*y1[2]+6656.0/12825*step_size*y3[2]+28561.0/56430*step_size*y4[2]-9.0/50*step_size*y5[2]+2.0/55*step_size*y6[2];
	  u5[i]=u5[i-1]+16.0/135*step_size*y1[3]+6656.0/12825*step_size*y3[3]+28561.0/56430*step_size*y4[3]-9.0/50*step_size*y5[3]+2.0/55*step_size*y6[3];
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
