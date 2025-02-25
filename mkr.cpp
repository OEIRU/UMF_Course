#include "mkr.h"






mkr::mkr()
{
	nx=0;
	ny=0;
	kx=1;
	ky=1;
	hx0=1;
	hy0=1;
	x0=0;
	y0=0;
	x1=0;
	y1=0;
	q=0;
	k=0;
	f=NULL;
	u1=NULL;
	u2=NULL;
	u3=NULL;
	u4=NULL;
	typekr[0]=0;
	typekr[1]=0;
	typekr[2]=0;
	typekr[3]=0;
}

int mkr::readall()
{
	FILE *fv[9];

	fv[0]=fopen("data/inputx.txt","r");
	fv[1]=fopen("data/inputy.txt","r");
	fv[2]=fopen("data/f(x).txt","r");
	fv[3]=fopen("data/k and q.txt","r");
	fv[4]=fopen("data/kr_usl1.txt","r");
	fv[5]=fopen("data/kr_usl2.txt","r");
	fv[6]=fopen("data/kr_usl3.txt","r");
	fv[7]=fopen("data/kr_usl4.txt","r");
	fv[8]=fopen("data/params.txt","r");

	for(int i=0;i<9;i++)
		if(fv[i]==NULL)
			return 0;

	if(fscanf(fv[0],"%lf %lf %lf %lf", &kx, &hx0, &x0, &x1)<4) return 0;
	if(fscanf(fv[1],"%lf %lf %lf %lf", &ky, &hy0, &y0, &y1)<4) return 0;

	if(kx!=1.0)
		nx=log(-(x1-x0)*(1.0-kx)/hx0+1.0)/log(kx)+1.0+0.1;
	else
		nx=(x1-x0)/hx0+1.0+0.1;

	if(ky!=1.0)
		ny=log(-(y1-y0)*(1.0-ky)/hy0+1.0)/log(ky)+1.0+0.1;
	else
		ny=(y1-y0)/hy0+1.0+0.1;

	N=ny*nx;

	A=new double*[N];
	for(int i=0;i<N;i++)
		A[i]=new double[5];

	b=new double[N];
	f=new double[N];

	for(int i=0;i<N;i++)
		if(fscanf(fv[2],"%lf",&f[i])==EOF) return 0;
	
	if(fscanf(fv[3],"%lf",&k)==EOF) return 0;
	if(fscanf(fv[3],"%lf",&q)==EOF) return 0;

	u1=new double[ny-2];
	if(fscanf(fv[4],"%d",&typekr[0])==EOF)return 0;
	for(int i=0;i<ny-2;i++)
		if(fscanf(fv[4],"%lf",&u1[i])==EOF) return 0;

	u2=new double[nx];
	if(fscanf(fv[5],"%d",&typekr[1])==EOF)return 0;
	for(int i=0;i<nx;i++)
		if(fscanf(fv[5],"%lf",&u2[i])==EOF) return 0;

	u3=new double[ny-2];
	if(fscanf(fv[6],"%d",&typekr[2])==EOF)return 0;
	for(int i=0;i<ny-2;i++)
		if(fscanf(fv[6],"%lf",&u3[i])==EOF) return 0;

	u4=new double[nx];
	if(fscanf(fv[7],"%d",&typekr[3])==EOF)return 0;
	for(int i=0;i<nx;i++)
		if(fscanf(fv[7],"%lf",&u4[i])==EOF) return 0;

	x=new double[N];
	for(int i=0;i<N;i++)
		if(fscanf(fv[8],"%lf",&x[i])==EOF) return 0;

	if(fscanf(fv[8],"%d %lf",&maxiter, &e)<2) return 0;

	ii[0]=-nx-1;
	ii[1]=-1;
	ii[2]=1;
	ii[3]=nx+1;

	for(int i=0;i<9;i++)
		fclose(fv[i]);

	return 1;
}




void mkr::createSLAU()
{
	double hi,him1,hj,hjm1,himax,hjmax;
	ii[0]=-nx;
	ii[1]=-1;
	ii[2]=1;
	ii[3]=nx;

	
	hjm1=hy0;
	hj=hy0*ky;

	for(int i=1;i<ny-1;i++)
	{
		him1=hx0;
		hi=hx0*kx;
		for(int l=i*nx+1;l<i*nx+nx-1;l++)
		{
			A[l][0]=-2.0*k/(hjm1*(hj+hjm1));
			A[l][1]=-2.0*k/(him1*(hi+him1));
			A[l][2]=2.0*k/(hi*him1)+2.0*k/(hj*hjm1)+q;
			A[l][3]=-2.0*k/(hi*(hi+him1));
			A[l][4]=-2.0*k/(hj*(hj+hjm1));
			b[l]=f[l];
			him1=hi;
			hi*=kx;
		}
		hjm1=hj;
		hj*=ky;
	}

	if(typekr[3]==1)
	{
		for(int i=0;i<nx;i++)
		{
			A[i][0]=0;
			A[i][1]=0;
			A[i][2]=1;
			A[i][3]=0;
			A[i][4]=0;
			b[i]=u4[i];
		}
	}
	else
	{
		for(int i=0;i<nx;i++)
		{
			A[i][0]=0;
			A[i][1]=0;
			A[i][2]=-1.0/hy0;
			A[i][3]=0;
			A[i][4]=1.0/hy0;
			b[i]=u4[i]/k;
		}
	}

	if(typekr[1]==1)
	{
		for(int i=(ny-1)*nx,j=0;i<nx*ny;i++,j++)
		{
			A[i][0]=0;
			A[i][1]=0;
			A[i][2]=1;
			A[i][3]=0;
			A[i][4]=0;
			b[i]=u2[j];
		}
	}
	else
	{
		int i=1;
		for(hjmax=hy0*kx;i<ny;i++)
		{
			hjmax*=ky;
		}
		for(int i=(ny-1)*nx,j=0;i<nx*ny;i++,j++)
		{
			A[i][0]=-1.0/hjmax;
			A[i][1]=0;
			A[i][2]=1.0/hjmax;
			A[i][3]=0;
			A[i][4]=0;
			b[i]=u2[i]/k;
		}
	}

	if(typekr[0]==1)
	{
		for(int l=1,i=0;l<ny-1;l++,i++)
		{
			A[l*nx][0]=0;
			A[l*nx][1]=0;
			A[l*nx][2]=1.0;
			A[l*nx][3]=0;
			A[l*nx][4]=0;
			b[l*nx]=u1[i];
		}
	}
	else
	{
		for(int l=1,i=0;l<ny-1;l++,i++)
		{
			A[l*nx][0]=0;
			A[l*nx][1]=0;
			A[l*nx][2]=-1.0/hx0;
			A[l*nx][3]=1.0/hx0;
			A[l*nx][4]=0;
			b[l*nx]=u1[i]/k;
		}
	}


	if(typekr[2]==1)
	{
		for(int l=2,i=0;l<ny;l++,i++)
		{
			A[l*nx-1][0]=0;
			A[l*nx-1][1]=0;
			A[l*nx-1][2]=1;
			A[l*nx-1][3]=0;
			A[l*nx-1][4]=0;
			b[l*nx-1]=u3[i];
		}
	}
	else
	{
		int i=1;
		for(himax=hx0*kx;i<nx;i++)
		{
			himax*=kx;
		}
		for(int l=2,i=0;l<ny;l++,i++)
		{
			A[l*nx-1][0]=0;
			A[l*nx-1][1]=-1.0/himax;
			A[l*nx-1][2]=1.0/himax;
			A[l*nx-1][3]=0;
			A[l*nx-1][4]=0;
			b[l*nx-1]=u3[i]/k;
		}
	}
}


void mkr::printx()
{
	FILE *f=fopen("x.txt","w");
	for(int i=0;i<ny;i=i+2)
	{
		fprintf(f,"\n");
		for(int j=i*nx;j<(i+1)*nx;j=j+2)
			fprintf(f,"%lf  ",x[j]);
	}
	fclose(f);
} 

