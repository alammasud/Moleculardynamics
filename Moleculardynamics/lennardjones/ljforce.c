#include<math.h>
#include"pbc.h"


void ljforce(double* x,double* y, double* z,int ntot,double* lx,double* ly,double* lz,double* fx,double* fy,double* fz,double E,double ro,double m,double n){
	
	int i,j,k;
	double cc;
	cc=m*n/(n-m);
	for(i=0;i<ntot;i++){
		fx[i]=0.0;
		fy[i]=0.0;
		fz[i]=0.0;
	}
	
	double dx,dy,dz,r,r2,sc;

	for(i=0;i<ntot;i++){
		for(j=i+1;j<ntot;j++){
			dx=x[j]-x[i];
			dy=y[j]-y[i];
			dz=z[j]-z[i];
			
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

			sc=ro/r;

			fx[i] = fx[i]+((cc*E*sc)*(pow(sc,m-1)-pow(sc,n-1))*(dx/r2));
			fy[i] = fy[i]+((cc*E*sc)*(pow(sc,m-1)-pow(sc,n-1))*(dy/r2));
			fz[i] = fz[i]+((cc*E*sc)*(pow(sc,m-1)-pow(sc,n-1))*(dz/r2));

			fx[j] = fx[j]-((cc*E*sc)*(pow(sc,m-1)-pow(sc,n-1))*(dx/r2));
			fy[j] = fy[j]-((cc*E*sc)*(pow(sc,m-1)-pow(sc,n-1))*(dy/r2));
			fz[j] = fz[j]-((cc*E*sc)*(pow(sc,m-1)-pow(sc,n-1))*(dz/r2));
			
			
			
		}
	}

}



