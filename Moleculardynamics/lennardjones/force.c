#include<math.h>
#include"pbc.h"


void force(double latpara,double* x,double* y, double* z,int ntot,double* lx,double* ly,double* lz,double* fx,double* fy,double* fz,double E,double ro,int nx,int ny,int nz){
	
	int i,j,k;
	double m,n,c1;
	m=5.0;n=8.0;c1=m*n/(n-m);
	for(i=0;i<ntot;i++){
		fx[i]=0.0;
		fy[i]=0.0;
		fz[i]=0.0;
	}
	double a1,a2,a3;
	double dx,dy,dz,r,r2,sqr,sf;
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
			sf=ro/r;
			fx[i] = fx[i]+((c1*E*sf)*(pow(sf,m-1)-pow(sf,n-1))*(dx/r2));
			fy[i] = fy[i]+((c1*E*sf)*(pow(sf,m-1)-pow(sf,n-1))*(dy/r2));
			fz[i] = fz[i]+((c1*E*sf)*(pow(sf,m-1)-pow(sf,n-1))*(dz/r2));

			fx[j] = fx[j]-((c1*E*sf)*(pow(sf,m-1)-pow(sf,n-1))*(dx/r2));
			fy[j] = fy[j]-((c1*E*sf)*(pow(sf,m-1)-pow(sf,n-1))*(dy/r2));
			fz[j] = fz[j]-((c1*E*sf)*(pow(sf,m-1)-pow(sf,n-1))*(dz/r2));
			
			
			
		}
	}

}



