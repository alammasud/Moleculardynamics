#include"ljforce.h"

void vverlet(double* lx,double* ly,double* lz,double* x,double* y,double* z,double* vx,double* vy,double* vz,double* fx,double* fy,double* fz,double mass,double dt,int ntot,double E,double ro,double m,double n){
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
      
        ljforce(x,y,z,ntot,lx,ly,lz,fxnew,fynew,fznew,E,ro,m,n);

	for(i=0;i<ntot;i++){
		vx[i]=vx[i]+dt*0.5*((fx[i]/mass)+(fxnew[i]/mass));
		vy[i]=vy[i]+dt*0.5*((fy[i]/mass)+(fynew[i]/mass));
		vz[i]=vz[i]+dt*0.5*((fz[i]/mass)+(fznew[i]/mass));
	}
}
