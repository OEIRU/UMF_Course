#pragma once
#include <stdio.h>
#include <math.h>
#include <conio.h>
typedef	double	type;

class dmtrx
{
private:
	type *L3;
		type *L2;
			type *L1;

	type *U1;
		type *U2;
			type *U3;

	type *D;
	type* f;
	type *toch_resh;

	int maxiter;	
	int n;
	int m1;
	int m2;
	int n1 ;
	int n2 ;
	type w;
	type eps;
	
public:
	void enter();
	void push(FILE *, int , type , type *);
	void  Gauss_Z();
	void change(type*, type*);
	void multiply(type*, type*);
	type otn_nevjazka(type *);
	type pogresh(type *);
	dmtrx(void);
   ~dmtrx(void);
};
