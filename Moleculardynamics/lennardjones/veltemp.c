#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

double veltemp(double *vx,double *vy,double *vz,int ntot,double mass,double kb){  
   int i; 	
   double Temprandom;
   double val=0.0;

   for(i=0 ;i<ntot;i++){
      val=val+vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
   }
   Temprandom=(mass*val)/(3*kb*ntot);

 return Temprandom;
}


