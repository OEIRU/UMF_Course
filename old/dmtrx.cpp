#include "dmtrx.h"
#include <stdio.h>
#include <math.h>
#include <conio.h>
#include <vector>
using namespace std;
vector <double> Diag, Diagn, Diagn2, Diagn1, Diagv, Diagv1, Diagv2, F;
dmtrx::~dmtrx(void)
{
}
dmtrx::dmtrx(void)
{
	L3 = new type[0x100];		U1 = new type[0x100];
	L2 = new type[0x100] ;		U2 = new type[0x100];
	L1 = new type[0x100];		U3 = new type[0x100];
	D	=	new type[0x100];	toch_resh = new type[0x100];
	f	=	new type[0x100];
}

void dmtrx::enter()
{
	int i;
	FILE *fp, *L3p, *L2p, *L1p, *Dp, *U1p, *U2p, *U3p, *ff, *tt;
	fp = fopen( "in.txt", "r");
	L3p = fopen ("inL3.txt", "r");
	L2p = fopen ("inL2.txt", "r");
	L1p = fopen ("inL1.txt", "r");
	Dp = fopen ("inD.txt", "r");
	U1p = fopen ("inU1.txt", "r");
	U2p = fopen ("inU2.txt", "r");
	U3p = fopen ("inU3.txt", "r");
	ff = fopen ("inF.txt", "r");
	tt = fopen ("toch_resh.txt", "r");
	fscanf(fp, "%d", &n);
	fscanf(fp, "%d", &m1);
	fscanf(fp, "%d", &m2);
	fscanf(fp, "%d", &maxiter);
	fscanf(fp, "%lg", &w);
	fscanf(fp, "%lg", &eps);
	n1 = n-3-m1-m2;
	n2 = n-2-m1;
	for( i = 0; i < n1; i++)
			{fscanf(L3p, "%lg", &L3[i]);Diagn2.push_back(L3[i]);}
	for( i = 0; i < n2; i++)
			{fscanf(L2p, "%lg", &L2[i]);Diagn.push_back(L2[i]);}
	for( i = 0; i < n-1; i++)
			{fscanf(L1p, "%lg", &L1[i]);Diagn1.push_back((L1[i]));}
	for( i = 0; i < n; i++)
			{fscanf(Dp, "%lg", &D[i]); Diag.push_back(D[i]);}
	for( i = 0; i < n-1; i++)
			{fscanf(U1p, "%lg", &U1[i]); Diagv1.push_back(U1[i]);}
	for( i = 0; i < n2; i++)
			{fscanf(U2p, "%lg", &U2[i]); Diagv.push_back(U2[i]);}
	for( i = 0; i < n1; i++)
			{fscanf(U3p, "%lg", &U3[i]); Diagv2.push_back(U3[i]);}
	for( i = 0; i < n; i++)
			{fscanf(ff, "%lg", &f[i]); F.push_back(f[i]);}
	for( i = 0; i < n; i++)
			fscanf(tt, "%lg", &toch_resh[i]);

	fclose( fp);
	fclose( L1p);
	fclose( L2p);
	fclose( L3p);
	fclose( Dp);
	fclose( U1p);
	fclose( U2p);
	fclose( U3p);
	fclose( ff);
	
}

void dmtrx::push( FILE *fp, int k, type pogr, type *x)
{
	fprintf ( fp, "%d\n", k);
	for ( int i = 0; i < n; i++)
	fprintf ( fp, "%.14le\n", x[i]);
	fprintf ( fp, "\n%.3le\n\n", pogr);
}

type dmtrx::otn_nevjazka(type *x)
{
	int i, k;
	type sum, norma;
	type *f_;
	f_ = new type[0x100];
	change (f_,NULL);
	sum = 0; norma = 0;
	for( i = 0; i < n; i++)
		f_[i] = D[i]*x[i];

	k = n-n1;
	for( i = 0; i < n1; i++, k++)
		{
		f_[k] += L3[i]*x[i];
		f_[i] += U3[i]*x[k];
		}

	k = n-n2;
	for( i = 0; i < n2; i++, k++)
		{
		f_[k] += L2[i]*x[i];
		f_[i] += U2[i]*x[k];
		}

	k = 1;
	for( i = 0; i < n-1; i++, k++)
		{
		f_[k] += L1[i]*x[i];
		f_[i] += U1[i]*x[k];
		}

	for( i = 0; i < n; i++)
		{	
		sum += pow(f_[i]-f[i],2);
		norma += pow(f[i],2);
		}
	norma = sqrt(sum)/sqrt(norma);
	return norma;
}

type dmtrx::pogresh(type *x)
{	
	type sum, norma;
	sum = 0; norma = 0;
	for( int i = 0; i < n; i++)
		{
		sum += pow(toch_resh[i]-x[i],2);
		norma +=pow(toch_resh[i],2);
		}
	norma = sqrt(sum)/n;///sqrt(norma);
	return norma;
}

void dmtrx::multiply( type *Xsl, type *x)
{
	int i = 0;
	type s = 0.0;
	int u1 = 1; 
	int u2 = 2+m1;
	int u3 = 3+m1+m2;
	s += D[i]*x[i]+U1[i]*x[u1]+U2[i]*x[u2]+U3[i]*x[u3];
	s = f[i]-s;
	Xsl[i] = x[i]+w*s/D[i];
	int l1 = 0;
	i = 1;
	u1 = 2;
	u2 = 3+m1;
	u3 = 4+m1+m2;
	if( n2 <= n-n2)
	{
		s = 0.0;
		for( ; i < n1; l1++, i++, u1++, u2++, u3++)
		{				 
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s += U2[i]*x[u2];
			s += U3[i]*x[u3];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		s = 0.0;
		for( ; i < n2; l1++, i++, u1++, u2++)
		{
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s += U2[i]*x[u2];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		s = 0.0;
		for( ; i < n-n2; l1++, i++, u1++)
		{
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		int l2 = 0;
		s = 0.0;
		for( ; i < n-n1; l2++, l1++, i++, u1++)
		{
			s += L2[l2]*x[l2];
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		int l3 = 0;
		s = 0.0;
		for( ; i < n-1; l3++, l2++, l1++, i++, u1++)
		{
			s += L3[l3]*x[l3];
			s += L2[l2]*x[l2];
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		s = 0.0;
		s += L3[l3]*x[l3];
		s += L2[l2]*x[l2];
		s += L1[l1]*x[l1];
		s += D[i]*x[i];
		s = f[i]-s;
		Xsl[i] = x[i]+w*s/D[i];
	}
	if( n1 <= n-n2 && n2 > n-n2)
	{
		s = 0.0;
		for( ; i < n1; l1++, i++, u1++, u2++, u3++)
		{				 
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s += U2[i]*x[u2];
			s += U3[i]*x[u3];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		s = 0.0;
		for( ; i < n-n2 ; l1++, i++, u1++, u2++)
		{
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s += U2[i]*x[u2];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		int l2 = 0;
		s = 0.0;
		for( ; i < n2; l2++, l1++, i++, u1++, u2++)
		{
			s += L2[l2]*x[l2];
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s += U2[i]*x[u2];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		s = 0.0 ;
		for( ; i < n-n1; l2++, l1++, i++, u1++)
		{
			s += L2[l2]*x[l2];
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		int l3 = 0;
		s = 0.0 ;
		for( ; i < n-1; l3++, l2++, l1++, i++, u1++)
		{
			s += L3[l3]*x[l3];
			s += L2[l2]*x[l2];
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		s = 0.0 ;
		s += L3[l3]*x[l3];
		s += L2[l2]*x[l2];
		s += L1[l1]*x[l1];
		s += D[i]*x[i];
		s = f[i]-s;
		Xsl[i] = x[i]+w*s/D[i];
	}
	if( n1 < n-n1 && n2 > n-n2)
	{
		s = 0.0;
		for( ; i < n-n2; l1++, i++, u1++, u2++, u3++)
		{				 
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s += U2[i]*x[u2];
			s += U3[i]*x[u3];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		int l2 = 0;
		s = 0.0;
		for( ; i < n1; l2++, l1++, i++, u1++, u2++, u3++)
		{
			s += L2[l2]*x[l2];
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s += U2[i]*x[u2];
			s += U3[i]*x[u3];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		s = 0.0;
		for( ; i < n-n1; l2++, l1++, i++, u1++, u2++)
		{
			s += L2[l2]*x[l2];
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s += U2[i]*x[u2];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		int l3 = 0;
		s = 0.0;
		for( ; i < n2; l3++, l2++, l1++, i++, u1++, u2++)
		{
			s += L3[l3]*x[l3];
			s += L2[l2]*x[l2];
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s += U2[i]*x[u2];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		s = 0.0;
		for( ; i < n-1; l3++, l2++, l1++, i++, u1++)
		{
			s += L3[l3]*x[l3];
			s += L2[l2]*x[l2];
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		s = 0.0;
		s += L3[l3]*x[l3];
		s += L2[l2]*x[l2];
		s += L1[l1]*x[l1];
		s += D[i]*x[i];
		s = f[i]-s;
		Xsl[i] = x[i]+w*s/D[i];
	}
	if( n1 > n-n1)
	{
		s = 0.0;
		for( ; i < n-n2; l1++, i++, u1++, u2++, u3++)
		{				 
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s += U2[i]*x[u2];
			s += U3[i]*x[u3];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		int l2 = 0;
		s = 0.0;
		for( ; i < n-n1; l2++, l1++, i++, u1++, u2++, u3++)
		{
			s += L2[l2]*x[l2];
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s += U2[i]*x[u2];
			s += U3[i]*x[u3];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		int l3 = 0;
		s = 0.0;
		for( ; i < n1; l3++, l2++, l1++, i++, u1++, u2++, u3++)
		{
			s += L3[l3]*x[l3];
			s += L2[l2]*x[l2];
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s += U2[i]*x[u2];
			s += U3[i]*x[u3];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		s = 0.0;
		for( ; i < n2; l3++, l2++, l1++, i++, u1++, u2++)
		{
			s += L3[l3]*x[l3];
			s += L2[l2]*x[l2];
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s += U2[i]*x[u2];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		s = 0.0;
		for( ; i < n-1; l3++, l2++, l1++, i++, u1++)
		{
			s += L3[l3]*x[l3];
			s += L2[l2]*x[l2];
			s += L1[l1]*x[l1];
			s += D[i]*x[i];
			s += U1[i]*x[u1];
			s = f[i]-s;
			Xsl[i] = x[i]+w*s/D[i];
			s = 0.0;
		}
		s = 0.0;
		s += L3[l3]*x[l3];
		s += L2[l2]*x[l2];
		s += L1[l1]*x[l1];
		s += D[i]*x[i];
		s = f[i]-s;
		Xsl[i] = x[i]+w*s/D[i];
	}
}

void dmtrx::change(type* a, type* b)
{
	int i;
	if (b)
		for (i=0; i<n; i++)
						a[i] = b[i];
	else
		for (i=0; i<n; i++)
						a[i] = 0.0;
}

void dmtrx::Gauss_Z()
{	   		
	FILE *fp;
	fp = fopen("out.csv","w");
	int k;
	type pogr;
	type *x, *x_;
	x_ = new type[0x100];
	x = new type[0x100];
	
	change (x,NULL);
	//change(x_,x);
	pogr = otn_nevjazka(x);
	for( k = 0; k < maxiter && pogr > eps; k++)
				{
				multiply(x,x);
				pogr = otn_nevjazka(x);
				push(fp,k,pogr,x);
				}
	fprintf( fp,"x*-x;\n");
	for( k = 0; k < n; k++)
		{
		fprintf( fp,"%5.16lg\n",fabs(toch_resh[k]-x[k]));
		}

	fprintf( fp,"\n");
	pogr = pogresh(x);
	fprintf( fp,"%5.16lg",pogr);

	fclose(fp);
}