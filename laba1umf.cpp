#include "stdio.h"
#include "mkr.h"



int main()
{



	for(double i=0;i<=1;i+=0.2)
	{
		//printf("\n");
		for(double j=0;j<=1;j+=0.2)
			printf("\n%lf",exp(i)+exp(j));
	}
	mkr m1;
	if(!m1.readall())
		return 1;
	m1.createSLAU();
	m1.solveSLAU();
	m1.printx();

	return 0;
}
