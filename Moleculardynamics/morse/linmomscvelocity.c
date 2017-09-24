#include<stdio.h>
#include<math.h>
#include"veltemp.h"

void linmomscvelocity(double *vx, double *vy, double *vz, int ntot, double Tdesire, double mass, double kb){
	int i;
	double sumvx, sumvy, sumvz;
	sumvx = 0;
	sumvy = 0;	
	sumvz = 0;

	for(i = 0; i<ntot; i++)	{
		sumvx = sumvx + vx[i];
		sumvy = sumvy + vy[i];
		sumvz = sumvz + vz[i];
	}

	for(i=0; i<ntot; i++){
		vx[i]= vx[i]-(sumvx/ntot);
		vy[i]= vy[i]-(sumvy/ntot);
		vz[i]= vz[i]-(sumvz/ntot);
	}

	sumvx = 0;
	sumvy = 0;
	sumvz = 0;

	for(i = 0; i<ntot; i++){
		
			sumvx = sumvx + vx[i];
			sumvy = sumvy + vy[i];
			sumvz = sumvz + vz[i];
		}

        double trand;
        trand=veltemp(vx,vy,vz,ntot, mass, kb);


for (i = 0; i <ntot; i++){
   		 vx[i] = vx[i]*sqrt(Tdesire/trand);		
    		 vy[i] = vy[i]*sqrt(Tdesire/trand);
                 vz[i] = vz[i]*sqrt(Tdesire/trand);   
  	}
	sumvx = 0;
	sumvy = 0;
	sumvz = 0;
	for(i = 0; i<ntot; i++){
		
		sumvx = sumvx + vx[i];
		sumvy = sumvy + vy[i];
		sumvz = sumvz + vz[i];
	}

}
