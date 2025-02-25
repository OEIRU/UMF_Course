#include "stdio.h"
#include "slau.h"
class mkr:public slau
{
private:
	int nx,ny;
	double kx,ky,hx0,hy0,x0,y0,x1,y1;
	double *f;
	double k,q;


	int typekr[4];
	double *u1,*u2,*u3,*u4;

public:
	mkr();
	int readall();
	void createSLAU();
	void printx();
};
