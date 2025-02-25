#include "stdio.h"
#include "math.h"

class slau
{
protected:
	int N;
	int m;
	int ii[4];
	double **A;
	double *x,*b;

	int maxiter;
	double e;
public:
	slau();
	void solveSLAU();
	virtual void printx();

    double calcnormnevyaz();
	void multAx(double *temp);
	double calcnorm(double *temp);
};
