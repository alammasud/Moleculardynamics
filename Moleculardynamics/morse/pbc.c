double pbc(double r){
	double sign;
	if(r > 0 ) sign = 1;
	if(r < 0)  sign = -1;
	r =r-sign;
	return r;
}
