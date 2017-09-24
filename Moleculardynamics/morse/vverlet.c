#include"morse.h"

void vverlet(double* lx,double* ly,double* lz,double* x,double* y,double* z,double* vx,double* vy,double* vz,double* fx,double* fy,double* fz,double mass,double dt,int ntot,double D,double ro,double alpha){
	int i;
	for(i=0;i<ntot;i++){
		x[i] = x[i] + (vx[i]*dt) + ((0.5*fx[i])/mass)*dt*dt;
		y[i] = y[i] + (vy[i]*dt) + ((0.5*fy[i])/mass)*dt*dt;
		z[i] = z[i] + (vz[i]*dt) + ((0.5*fz[i])/mass)*dt*dt;
	}
	double fxnew[ntot],fynew[ntot],fznew[ntot];
	for(i=0;i<ntot;i++){
		fxnew[i]=0.0;
		fynew[i]=0.0;
		fznew[i]=0.0;
	}
      
       morse(x,y,z,ntot,lx,ly,lz,fxnew,fynew,fznew,D,ro,alpha);	

	for(i=0;i<ntot;i++){
		vx[i]=vx[i]+dt*0.5*((fx[i]/mass)+(fxnew[i]/mass));
		vy[i]=vy[i]+dt*0.5*((fy[i]/mass)+(fynew[i]/mass));
		vz[i]=vz[i]+dt*0.5*((fz[i]/mass)+(fznew[i]/mass));
	}
}
