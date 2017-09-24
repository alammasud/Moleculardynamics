#include<math.h>
#include"pbc.h"


double potentialenergy(double* x,double* y, double* z,int ntot,double* lx,double* ly,double* lz,double E,double ro,double m,double n){
	
	int i,j;
	double pe=0;
	double dx,dy,dz,r,r2,sc,msc,nsc;
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
			msc=pow(sc,m);
			nsc=pow(sc,n);
			pe = pe + ((E/(n-m))*((m*nsc)-(n*msc)));
			//pe = pe + D*(1.0-exp(-sc))*(1.0-exp(-sc));
			
			
		}
	}
return pe;
}


