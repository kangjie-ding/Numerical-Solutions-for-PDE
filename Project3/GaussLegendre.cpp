#include"generator.h"

double result5[4];

double* f5(const double x1,const double x2,const double x4,const double x5)
{
  result5[0]=x4;
  result5[1]=x5;
  result5[2]=2*x5+x1-mu*(x1+mu-1)/pow(pow(x2,2)+pow(x1+mu-1,2),1.5)-(1-mu)*(x1+mu)/pow(pow(x2,2)+pow(x1+mu,2),1.5);
  result5[3]=-2*x4+x2-mu*x2/pow(pow(x2,2)+pow(x1+mu-1,2),1.5)-(1-mu)*x2/pow(pow(x2,2)+pow(x1+mu,2),1.5);
  return result5;
}
void GLRK_Method::Solver()
{
  double u1[steps+1];
  double u2[steps+1];
  double u4[steps+1];
  double u5[steps+1];
  const double step_size=period/steps;
  const double e=0.00000001;
  u1[0]=initial[0];u2[0]=initial[1];u4[0]=initial[2];u5[0]=initial[3];
  for(int i=1;i<=steps;i++)
    {
      if(s==1)
	{
	  double y[4]={u1[i-1],u2[i-1],u4[i-1],u5[i-1]};
	  //use fixed point iteration to solve the implicit equation
	  for(int j=1;j<=1000;j++)
	    {
	      double* temp=f5(u1[i-1]+0.5*step_size*y[0],u2[i-1]+0.5*step_size*y[1],u4[i-1]+0.5*step_size*y[2],u5[i-1]+0.5*step_size*y[3]);
	      if(abs(temp[0]-y[0])<e&&abs(temp[1]-y[1])<e&&abs(temp[2]-y[2])<e&&abs(temp[3]-y[3])<e)
		break;
	      y[0]=temp[0];y[1]=temp[1];y[2]=temp[2];y[3]=temp[3];
	    }
	  u1[i]=u1[i-1]+step_size*y[0];
	  u2[i]=u2[i-1]+step_size*y[1];
	  u4[i]=u4[i-1]+step_size*y[2];
	  u5[i]=u5[i-1]+step_size*y[3];
	}
      else if(s==2)
	{
	  double y1[4]={u1[i-1],u2[i-1],u4[i-1],u5[i-1]};
	  double y2[4]={u1[i-1],u2[i-1],u4[i-1],u5[i-1]};
	  //use fixed point iteration to solve the implicit equations
	  for(int j=1;j<=1000;j++)
	    {
	      double* temp1=f5(u1[i-1]+0.25*step_size*y1[0]+(3-2*sqrt(3))/12*step_size*y2[0],u2[i-1]+0.25*step_size*y1[1]+(3-2*sqrt(3))/12*step_size*y2[1],u4[i-1]+0.25*step_size*y1[2]+(3-2*sqrt(3))/12*step_size*y2[2],u5[i-1]+0.25*step_size*y1[3]+(3-2*sqrt(3))/12*step_size*y2[3]);
	      double* temp2=f5(u1[i-1]+(3+2*sqrt(3))/12*step_size*y1[0]+0.25*step_size*y2[0],u2[i-1]+(3+2*sqrt(3))/12*step_size*y1[1]+0.25*step_size*y2[1],u4[i-1]+(3+2*sqrt(3))/12*step_size*y1[2]+0.25*step_size*y2[2],u5[i-1]+(3+2*sqrt(3))/12*step_size*y1[3]+0.25*step_size*y2[3]);
	      
	      if(abs(temp1[0]-y1[0])<e&&abs(temp1[1]-y1[1])<e&&abs(temp1[2]-y1[2])<e&&abs(temp1[3]-y1[3])<e&&abs(temp2[0]-y2[0])<e&&abs(temp2[1]-y2[1])<e&abs(temp2[2]-y2[2])<e&&abs(temp2[3]-y2[3])<e)
		break;
	      y1[0]=temp1[0];y1[1]=temp1[1];y1[2]=temp1[2];y1[3]=temp1[3];
	      y2[0]=temp2[0];y2[1]=temp2[1];y2[2]=temp2[2];y2[3]=temp2[3];
	    }
	  u1[i]=u1[i-1]+0.5*step_size*(y1[0]+y2[0]);
	  u2[i]=u2[i-1]+0.5*step_size*(y1[1]+y2[1]);
	  u4[i]=u4[i-1]+0.5*step_size*(y1[2]+y2[2]);
	  u5[i]=u5[i-1]+0.5*step_size*(y1[3]+y2[3]);
	}
      else
	{
	  double y1[4]={u1[i-1],u2[i-1],u4[i-1],u5[i-1]};
	  double y2[4]={u1[i-1],u2[i-1],u4[i-1],u5[i-1]};
	  double y3[4]={u1[i-1],u2[i-1],u4[i-1],u5[i-1]};
	  //use fixed point iteration to solve the implicit equations
	  for(int j=1;j<=1000;j++)
	    {
	      double* temp1=f5(u1[i-1]+5.0/36*step_size*y1[0]+(2.0/9-sqrt(15)/15)*step_size*y2[0]+(5.0/36-sqrt(15)/30)*step_size*y3[0],u2[i-1]+5.0/36*step_size*y1[1]+(2.0/9-sqrt(15)/15)*step_size*y2[1]+(5.0/36-sqrt(15)/30)*step_size*y3[1],u4[i-1]+5.0/36*step_size*y1[2]+(2.0/9-sqrt(15)/15)*step_size*y2[2]+(5.0/36-sqrt(15)/30)*step_size*y3[2],u5[i-1]+5.0/36*step_size*y1[3]+(2.0/9-sqrt(15)/15)*step_size*y2[3]+(5.0/36-sqrt(15)/30)*step_size*y3[3]);
	      double* temp2=f5(u1[i-1]+(5.0/36+sqrt(15)/24)*step_size*y1[0]+2.0/9*step_size*y2[0]+(5.0/36-sqrt(15)/24)*step_size*y3[0],u2[i-1]+(5.0/36+sqrt(15)/24)*step_size*y1[1]+2.0/9*step_size*y2[1]+(5.0/36-sqrt(15)/24)*step_size*y3[1],u4[i-1]+(5.0/36+sqrt(15)/24)*step_size*y1[2]+2.0/9*step_size*y2[2]+(5.0/36-sqrt(15)/24)*step_size*y3[2],u5[i-1]+(5.0/36+sqrt(15)/24)*step_size*y1[3]+2.0/9*step_size*y2[3]+(5.0/36-sqrt(15)/24)*step_size*y3[3]);
	      double* temp3=f5(u1[i-1]+(5.0/36+sqrt(15)/30)*step_size*y1[0]+(2.0/9+sqrt(15)/15)*step_size*y2[0]+5.0/36*step_size*y3[0],u2[i-1]+(5.0/36+sqrt(15)/30)*step_size*y1[1]+(2.0/9+sqrt(15)/15)*step_size*y2[1]+5.0/36*step_size*y3[1],u4[i-1]+(5.0/36+sqrt(15)/30)*step_size*y1[2]+(2.0/9+sqrt(15)/15)*step_size*y2[2]+5.0/36*step_size*y3[2],u5[i-1]+(5.0/36+sqrt(15)/30)*step_size*y1[3]+(2.0/9+sqrt(15)/15)*step_size*y2[3]+5.0/36*step_size*y3[3]);

	      if(abs(temp1[0]-y1[0])<e&&abs(temp1[1]-y1[1])<e&&abs(temp1[2]-y1[2])<e&&abs(temp1[3]-y1[3])<e&&abs(temp2[0]-y2[0])<e&&abs(temp2[1]-y2[1])<e&abs(temp2[2]-y2[2])<e&&abs(temp2[3]-y2[3])<e&&abs(temp3[0]-y3[0])<e&&abs(temp3[1]-y3[1])<e&abs(temp3[2]-y3[2])<e&&abs(temp3[3]-y3[3])<e)
		break;
	      y1[0]=temp1[0];y1[1]=temp1[1];y1[2]=temp1[2];y1[3]=temp1[3];
	      y2[0]=temp2[0];y2[1]=temp2[1];y2[2]=temp2[2];y2[3]=temp2[3];
	      y3[0]=temp3[0];y3[1]=temp3[1];y3[2]=temp3[2];y3[3]=temp3[3];
	    }
	  u1[i]=u1[i-1]+5.0/18*step_size*y1[0]+4.0/9*step_size*y2[0]+5.0/18*step_size*y3[0];
	  u2[i]=u2[i-1]+5.0/18*step_size*y1[1]+4.0/9*step_size*y2[1]+5.0/18*step_size*y3[1];
	  u4[i]=u4[i-1]+5.0/18*step_size*y1[2]+4.0/9*step_size*y2[2]+5.0/18*step_size*y3[2];
	  u5[i]=u5[i-1]+5.0/18*step_size*y1[3]+4.0/9*step_size*y2[3]+5.0/18*step_size*y3[3];
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
