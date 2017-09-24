double kineenrgy(int ntot,double* vx,double* vy,double* vz,double mass){
	int i;
	double vsum=0;
	for(i=0;i<ntot;i++){
		vsum=vsum+vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
	}
	double ke=0.5*mass*vsum;
	return ke;
}

