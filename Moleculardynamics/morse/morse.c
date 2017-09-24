#include<math.h>
#include"pbc.h"


void morse(double* x,double* y, double* z,int ntot,double* lx,double* ly,double* lz,double* fx,double* fy,double* fz,double D,double ro,double alpha){
	
	int i,j,k;
	for(i=0;i<ntot;i++){
		fx[i]=0.0;
		fy[i]=0.0;
		fz[i]=0.0;
	}
	double a1,a2,a3;
	double dx,dy,dz,r,r2,sqr,sf;
	for(i=0;i<ntot;i++){
		for(j=i+1;j<ntot;j++){
			dx=x[i]-x[j];
			dy=y[i]-y[j];
			dz=z[i]-z[j];
			
			dx=dx/lx[0];
			dy=dy/ly[1];
			dz=dz/lz[2];
			
			if(fabs(dx) > 0.5) dx = pbc(dx);
			if(fabs(dy) > 0.5) dy = pbc(dy);
			if(fabs(dz) > 0.5) dz = pbc(dz);
			
			dx=dx*lx[0];
			dy=dy*ly[1];
			dz=dz*lz[2];

			r2=dx*dx+dy*dy+dz*dz;
			r=sqrt(r2);
			
			a1=exp(-alpha*(r-ro));
			a2=exp(-2*alpha*(r-ro));
			a3=(2*D*alpha)/r;
			
			fx[i]=fx[i]+a3*dx*(a2-a1);
			fy[i]=fy[i]+a3*dy*(a2-a1);
			fz[i]=fz[i]+a3*dz*(a2-a1);
			
			fx[j]=fx[j]-a3*dx*(a2-a1);
			fy[j]=fy[j]-a3*dy*(a2-a1);
			fz[j]=fz[j]-a3*dz*(a2-a1);
			
		}
	}

}



