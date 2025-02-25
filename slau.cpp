#include "slau.h"
#include "math.h"

slau::slau()
{
	N=0;
	m=0;
	ii[0]=0;
	ii[1]=0;
	ii[2]=0;
	ii[3]=0;
	A=NULL;
	x=NULL;
	b=NULL;
	maxiter=0;
	e=1000;
}

void slau::solveSLAU()
{
	int i,j,l;
	double S;
	double *X=new double[N];
	double otnnev=calcnormnevyaz()/calcnorm(b);
	for(int k=0;k<maxiter&&otnnev>e;k++)
	{
		for(i=0;i<N;i++)
		{
			S=A[i][2]*x[i];
			for(j=0;j<2;j++)
			{		
				l=ii[j];
				l+=i;
				if(l>=0&&l<N)						
					S+=A[i][j]*X[l];
			}
			for(j=2;j<4;j++)
			{
				l=ii[j];
				l+=i;
				if(l>=0&&l<N)
					S+=A[i][j+1]*x[l];
			}
			X[i]=x[i]+(b[i]-S)/A[i][2];
		}
		for(i=0;i<N;i++)
			x[i]=X[i];
		otnnev=calcnormnevyaz()/calcnorm(b);
		printf("\nnum_iteration=%d   otn nev=%lf ",k,otnnev);
	}
	printf("\n");
}


double slau::calcnormnevyaz()
{
	double *temp=new double[N];
	double *p=new double[N];

	multAx(p);

	for(int i=0;i<N;i++)
		temp[i]=b[i]-p[i];

	return calcnorm(temp);
}

double slau::calcnorm(double *temp)
{
	double S=0;
	for(int i=0;i<N;i++)
		S+=temp[i]*temp[i];
	return sqrt(S);
}

void slau::multAx(double *temp)
{
	int i,j,l;
	for(i=0;i<N;i++)
	{
		temp[i]=A[i][2]*x[i];
		for(j=0;j<2;j++)
		{
			l=ii[j];
			l+=i;
			if(l>=0&&l<N)
				temp[i]+=A[i][j]*x[l];
		}
		for(j=2;j<4;j++)
		{
			l=ii[j];
			l+=i;
			if(l>=0&&l<N)
				temp[i]+=A[i][j+1]*x[l];
		}
	}
}

void slau::printx()
{
	FILE *f=fopen("x.txt","w");
	for(int i=0;i<N;i++)
		fprintf(f,"\n%lf",x[i]);
	fclose(f);
}
