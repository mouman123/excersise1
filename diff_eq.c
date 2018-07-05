#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

int main(){
	int n=100;
	double A[101][101],B[101], h=1.0/(n), al=1.0/h;
	FILE * ptr;
	double xi=0;
	ptr=fopen("values.txt","w");
	fprintf(ptr,"#time\t #x\t #T\n");
	gsl_vector *C= gsl_vector_alloc(n-1);	
	for(int i=0;i<n-1;i++){
		
		//if(i==0)B[i]=0.0;else
		//if(i==100)B[i]=0.0;else
		B[i]=100.0;
		gsl_vector_set(C,i,B[i]);
	}
	gsl_matrix *D = gsl_matrix_alloc(n-1,n-1);
	
	for(int i=0;i<n-1;i++){
		for(int j=0;j<n-1;j++){
			A[i][j]=0.0;
			if(i==j)A[i][j]=1.0+(2.0*al);
			if(fabs(i-j)==1)A[i][j]=-al;
		
			gsl_matrix_set(D, i ,j, A[i][j]);
		}
	}
	
	
	double t=0.0;
	gsl_vector *v= gsl_vector_alloc(n-1);
	gsl_permutation * p = gsl_permutation_alloc (n-1);
	fprintf(ptr,"%.20lf\t %.20lf\t %.20lf\n ",t,0.0,0.0);
	for(int j=0;j<n-1;j++){
		xi=(j)*h;
		fprintf(ptr,"%.20lf\t %.20lf\t %.20lf\n ",t,xi,gsl_vector_get(C,j) );
	}
	fprintf(ptr,"%.20lf\t %.20lf\t %.20lf\n ",t,1.0,0.0);
	
	for(int i=0;i<n-1;i++){
		int s; t+=h;
		xi=0.0;
		//gsl_vector_set(v,0,0.0);
		//gsl_vector_set(v,100,0.0);
		gsl_linalg_LU_decomp (D, p, &s);
		gsl_linalg_LU_solve (D, p, C, v);
		fprintf(ptr,"%.20lf\t %.20lf\t %.20lf\n ",t,0.0,0.0);
		for(int j=0;j<n-1;j++){
			gsl_vector_set(C,j,gsl_vector_get(v,j));
			double Ti=gsl_vector_get(v,j);
			fprintf(ptr,"%.15lf\t %.15lf\t %.15lf\n", t, xi, Ti);	
			xi+=h;
			}
			fprintf(ptr,"%.20lf\t %.20lf\t %.20lf\n ",t,1.0,0.0);
		}
	
/*	gsl_permutation_free (p);
	gsl_vector_free (v);*/
	return(0);
}
