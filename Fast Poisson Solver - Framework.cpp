// Please check the previous code to see the use of header files.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define Gauss 1 
// #define Gauss 0 : do not Gauss elimination
double Residual(double *r, double **A, double *x, double *b, int N);
double Error(double *x, double *u, int N);
int Exact_Discretization(double **A, int N);
int Print_Matrix(double **A, int N);
int Exact_Solution(double **U, int N);
int Exact_Source(double **F, int N);
int GaussElimination(double **M, double *x, double *b, int N);
int Print_Complex_Vector(double *x, double *y, int N);
int DST2D(double **X, int N);
int iDST2D(double **X, int N);
int Transpose(double **A, int N);
int Fast_Poisson_Solver(double **F, int N);

int main()
{
	int i, j, k, N, M; 
	double **A, *x, *u, **U, *b, **F, *r;
	double y_r[16], y_i[16];
	clock_t t1, t2;
	
	// Create memory for solving Ax = b, where r = b-Ax is the residue.
	for(N=4;N<8;N*=2)      // N ¨ì 64 
	{
		// M is the total number of unknowns.
		M = (N-1)*(N-1);
	#if Gauss 
		A = (double **) malloc(M*sizeof(double*));
		A[0] = (double *) malloc(M*M*sizeof(double));
		for(i=1;i<M;++i) A[i] = A[i-1] + M;
		for(i=0;i<M*M;++i) A[0][i] = 0.0;
		x = (double *) malloc(M*sizeof(double));
		r = (double *) malloc(M*sizeof(double));
	#endif
		b = (double *) malloc(M*sizeof(double));
		// Prepare for two dimensional unknown F
		// where b is the one dimensional vector and 
		// F[i][j] = F[j+i*(N-1)];
		F = (double **) malloc((N-1)*sizeof(double*));
		F[0] = b;
		for(i=1;i<N-1;++i) F[i] = F[i-1] + (N-1);
		
		// Prepare for two dimensional unknowns U
		// where u is the one dimensional vector and
		// U[i][j] = u[j+i*(N-1)] 
		u = (double *) malloc(M*sizeof(double));
		U = (double **) malloc((N-1)*sizeof(double*));
		U[0] = u;
		for(i=1;i<N-1;++i) U[i] = U[i-1] + (N-1);
		
	#if Gauss
		Exact_Discretization(A,N);
	#endif
		Print_Matrix(A, M);
		Exact_Source(F, N);
	#if Gauss
		t1 = clock();
		GaussElimination(A,x,b,M);
		t2 = clock();
	#endif
		printf("N=%d,",N);
	#if Gauss
		printf("Gauss Elimination: %f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
		printf("%e,", Residual(r,A,x,b,M));
	#endif
		Exact_Solution(U, N);
	#if Gauss
		printf("%e\n", Error(x, u, M));
	#endif
		t1 = clock();
		Fast_Poisson_Solver(F, N-1);
		t2 = clock();
		printf("Fast Poisson Solver: %f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
		printf("%e\n", Error(b, u, M));

	#if Gauss
		free(A[0]);
		free(A);
		free(x);
		free(r);
	#endif
		free(b);
		free(u);
		free(U);
		free(F);
		//system("pause");
	}
	return 0;
}
int Exact_Discretization(double **A, int N)
{
	int i,j,k;
	for(i=0;i<N-1;++i)
	{
		for(j=0;j<N-1;++j)
		{
			k = i + j*(N-1);
			A[k][k] = -4;
			if(j>0) A[k][k-(N-1)] = 1;
			if(i>0) A[k][k-1] = 1;
			if(i<N-2) A[k][k+1] = 1;
			if(j<N-2) A[k][k+(N-1)] = 1;
		}
	}
	return 0;
}
int Print_Matrix(double **A, int N)
{
	int i, j;
	for(i=0;i<N;++i)
	{
		for(j=0;j<N;++j)
		{
			printf("%.0f ",A[i][j]);
		}
		printf("\n");
	}
}
int Exact_Solution(double **U, int N)
{
	// put the exact solution 
	int i,j,k;
	double x, y, h;
	h = 1.0/N;
	for(i=0;i<N-1;++i)
	{
		for(j=0;j<N-1;++j)
		{
			//k = j + i*(N-1);
			x = (i+1)*h;
			y = (j+1)*h;
			U[i][j] = sin(M_PI*x)*sin(2*M_PI*y);
		}
	}
	return 1;
}
int Exact_Source(double **F, int N)
{
	int i,j,k;
	double x, y, h;
	h = 1.0/N;
	for(i=0;i<N-1;++i)
	{
		for(j=0;j<N-1;++j)
		{
			//k = j + i*(N-1);
			x = (i+1)*h;
			y = (j+1)*h;
			F[i][j] = -(1.0+4.0)*h*h*M_PI*M_PI*sin(M_PI*x)*sin(2*M_PI*y);
		}
	}
	return 1;	
}
int GaussElimination(double **M, double *x, double *b, int N)
{
	// Gauss Elimination
	int i, j, k;
	double **A, v;
	
	// Copy the matrix M to A
	A = (double **) malloc(N*sizeof(double*));
	A[0] = (double *) malloc(N*N*sizeof(double));
	for(i=1;i<N;++i) A[i] = A[i-1] + N;
	for(i=0;i<N;++i)
	for(j=0;j<N;++j)
		A[i][j] = M[i][j];
	
	// put b in x 
	for(k=0;k<N;++k) x[k] = b[k];
	
	// Gauss Elimination
	for(k=0;k<N;++k)
	{
		// without pivoting
		for(i=k+1;i<N;++i)
		{
			v = A[i][k]/A[k][k];
			x[i] -= v*x[k];
			for(j=k;j<N;++j)
			{
				A[i][j] -= v*A[k][j];
			}
		}
	}
	// Backward(for upper triangular matrix) Solver
	for(i=N-1;i>=0;i--)
	{
		for(j=i+1;j<N;++j)
		{
			x[i] -= A[i][j]*x[j];
		}
		x[i] /= A[i][i];
	}
	// Free the memory of A
	free(A[0]);
	free(A);
	return 0;
}
double Residual(double *r, double **A, double *x, double *b, int N)
{
	// compute r = b - Ax and return max_i r[i]
	int i, j;
	double v = 0.0;
	for(i=0;i<N;++i)
	{
		r[i] = b[i];
		for(j=0;j<N;++j)
		{
			r[i] -= A[i][j]*x[j];
		}
		if(fabs(r[i])>v) v = fabs(r[i]);
	}
	return v;
}
double Error(double *x, double *u, int N)
{
	// return max_i |x[i] - u[i]|
	int i;
	double v = 0.0, e;
	
	for(i=0;i<N;++i)
	{
		e = fabs(x[i]-u[i]);
		if(e>v) v = e;
	}
	return v;
}
// Fast Fourier Transform in place for N = 2^p 
int FFT(double *y_r, double *y_i, int N)
{
	// Your FFT
}
int Print_Complex_Vector(double *x, double *y, int N)
{
	int n;
	for(n=0;n<N;++n)
	{
		if(y[n]>=0) printf("%d : %f +%f i\n", n, x[n], y[n]);
		else printf("%d : %f %f i\n", n, x[n], y[n]);
	}
	return 0;
}
int DST(double *x, int N)
{
	// Your DST
}
int iDST(double *x, int N)
{
	// Your iDST
}
int DST2D(double **X, int N)
{
	// Your DST for 2D
}
int iDST2D(double **X, int N)
{
	// Your iDST for 2D
}
int Transpose(double **A, int N)
{
	int i, j;
	double v;
	for(i=0;i<N;++i)
	{
		for(j=i+1;j<N;++j)
		{
			v = A[i][j];
			A[i][j] = A[j][i];
			A[j][i] = v;
		}
	}
	return 0;
}
int Fast_Poisson_Solver(double **F, int N)
{
	// Your Fast Poisson Solver
}
