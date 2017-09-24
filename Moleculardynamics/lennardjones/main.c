#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include"create_fcc.h"
#include"veltemp.h"
#include"linmomscvelocity.h"
#include"ljforce.h"
#include"vverlet.h"
#include"pbc.h"
#include"potentialenergy.h"
#include"initializevel.h"
#include"kinenergy.h"

int main(){
	srand(time(NULL));
	int i,j,k,nx,ny,nz,ntot,nunit,nstep,nsave;
	//nx,ny,nz=number of duplicates along x,y,z direction,
	//ntot=total no. of atom, nunit=no. of atom in unit cell, nstep=number of steps
	double a;//lattice parameter
	double kb,mass,Tdesire,Tran,dt,ro,E,m,n,epsilon;
	//kb Boltzman constant,mass of 1 atom,desire temeprature, Tran random temeprature
	//dt time step, ro cutoff distance, D-depth of potentials, alpha -scaling paramter
	
	FILE* fp2;
	fp2=fopen("input.txt","r");
	fscanf(fp2,"%i %i %i %i %i %le %le %lf %le %le %le  %i %lf %lf",&nstep,&nx,&ny,&nz,&nunit,&a,&mass,&Tdesire,&dt,&ro,&epsilon,&nsave,&m,&n);
        
	ntot=nunit*nx*ny*nz;
	
	double *lx,*ly,*lz;//simulation cell dimension

	double *x,*y,*z; //atom position

	x=(double*) malloc(ntot*sizeof(double));
	y=(double*) malloc(ntot*sizeof(double));
	z=(double*) malloc(ntot*sizeof(double));
	
	lx=(double*) malloc(3*sizeof(double));
	ly=(double*) malloc(3*sizeof(double));
	lz=(double*) malloc(3*sizeof(double));
	
	create_fcc(nx,ny,nz,ntot,a,x,y,z,lx,ly,lz,nunit);//create fcc crystal
	
	double *vx,*vy,*vz;//velcoity of atom

	vx=(double*) malloc(ntot*sizeof(double));
	vy=(double*) malloc(ntot*sizeof(double));
	vz=(double*) malloc(ntot*sizeof(double));
	
	kb=1.38064852E-23;
	E=epsilon*kb;
	
	
	//initialize random velcotiy between -0.5 to 0.5
	for(i=0;i<ntot;i++){
		vx[i]=initializevel(-0.5, 0.5);
		vy[i]=initializevel(-0.5, 0.5);
		vz[i]=initializevel(-0.5, 0.5);
	}
	
	linmomscvelocity(vx,vy,vz,ntot,Tdesire,mass,kb);
	
	//fx-force along x
	double *fx,*fy,*fz;
	fx=(double*) malloc(ntot*sizeof(double));
	fy=(double*) malloc(ntot*sizeof(double));
	fz=(double*) malloc(ntot*sizeof(double));
        
	double time,pe,ke; //pe=potentials energy , ke=kintetic energy
	time=0;

	FILE* fp;
	fp=fopen("datates.txt","w+");
	FILE* fp3;
	fp3=fopen("energy.txt","w+");
	int cou=0;
	double tempera;
	double kin;
	
	while(time<nstep*dt){
		//call LJ force field
		 ljforce(x,y,z,ntot,lx,ly,lz,fx,fy,fz,E,ro,m,n);
		 //call velocity verlet
		 vverlet(lx,ly,lz,x,y,z,vx,vy,vz,fx,fy,fz,mass,dt,ntot,E,ro,m,n);
	       //save strcucture
		if(cou%nsave==0){
			fprintf(fp,"ITEM: TIMESTEP\n");
			fprintf(fp,"%le\n",time);
			fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
			fprintf(fp,"%i\n",ntot);
			fprintf(fp,"ITEM: BOX BOUNDS pp pp pp\n");
			fprintf(fp,"0.0 %le\n",lx[0]/a);
			fprintf(fp,"0.0 %le\n",ly[1]/a);
			fprintf(fp,"0.0 %le\n",lz[2]/a);
			fprintf(fp,"ITEM: ATOMS id type x y z\n");
			for(i=0;i<ntot;i++){
				fprintf(fp,"%i 1 %le %le %le\n",i+1,x[i]/a,y[i]/a,z[i]/a);
			}
		}
		//calculate potential energy
		pe=potentialenergy(x,y,z,ntot,lx,ly,lz,E,ro,m,n);

		kin=kineenrgy(ntot,vx,vy,vz,mass);
		if(cou%nsave==0){
			fprintf(fp3,"%le %le %le %le\n",time,pe,kin,pe+kin);
		}

		
		time=time+dt;
		cou++;
	}
	fclose(fp);
	fclose(fp2);
	fclose(fp3);
	

return 0;

}
