#include"generator.h"

double result3[4];

double* f3(const double x1,const double x2,const double x4,const double x5)
{
  result3[0]=x4;
  result3[1]=x5;
  result3[2]=2*x5+x1-mu*(x1+mu-1)/pow(pow(x2,2)+pow(x1+mu-1,2),1.5)-(1-mu)*(x1+mu)/pow(pow(x2,2)+pow(x1+mu,2),1.5);
  result3[3]=-2*x4+x2-mu*x2/pow(pow(x2,2)+pow(x1+mu-1,2),1.5)-(1-mu)*x2/pow(pow(x2,2)+pow(x1+mu,2),1.5);
  return result3;
}

void AdMo_Method::Solver()
{
  double u1[steps+1];
  double u2[steps+1];
  double u4[steps+1];
  double u5[steps+1];
  const double step_size=period/steps;
  const double e=0.00001;
  
  u1[0]=initial[0];u2[0]=initial[1];u4[0]=initial[2];u5[0]=initial[3];
  for(int i=1;i<=steps;i++)
    {
      if(p==2)
	{
	  u1[i]=u1[i-1];u2[i]=u2[i-1];u4[i]=u4[i-1];u5[i]=u5[i-1];
	  double* y1=f3(u1[i-1],u2[i-1],u4[i-1],u5[i-1]);
	  for(int j=1;j<=1000;j++)
	    {
	      double* y2=f3(u1[i],u2[i],u4[i],u5[i]);
	      double temp[4]={u1[i-1]+0.5*step_size*y1[0]+0.5*step_size*y2[0],u2[i-1]+0.5*step_size*y1[1]+0.5*step_size*y2[1],u4[i-1]+0.5*step_size*y1[2]+0.5*step_size*y2[2],u5[i-1]+0.5*step_size*y1[3]+0.5*step_size*y2[3]};
	      
	      if(abs(temp[0]-u1[i])<e&&abs(temp[1]-u1[i])<e&&abs(temp[2]-u4[i])<e&&abs(temp[3]-u5[i])<e)
		break;
	      u1[i]=temp[0];u2[i]=temp[1];u4[i]=temp[2];u5[i]=temp[3];
	    }
	}
      else
	;
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
