#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>

//srand(time(NULL));
double initializevel(double min, double max){
//	srand(time(NULL));
	double rnd;
	double range=max-min;
	rnd = min+(double) rand()/((double) RAND_MAX/range+1);
	return rnd;
}
