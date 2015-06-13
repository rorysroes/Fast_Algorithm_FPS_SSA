#include <stdio.h>
#include <stdlib.h>
int BaseP(int *x, int N, int a, int P);
int SFT(int *X, int *x, int N, int w, int P);

int main()
{
	int a, b, i, *x, *y, *X, *Y;

	int N = 32, n = 1365, w = 10, W = 141, P = 1409; 
	printf("Input a , b : ");
	scanf("%d %d", &a, &b);
	//printf("%d x %d = %d\n", a, b, a*b);

	x = (int * ) malloc(N*sizeof(int));
	y = (int * ) malloc(N*sizeof(int));
    X = (int * ) malloc(N*sizeof(int));
    Y = (int * ) malloc(N*sizeof(int));
    
	BaseP(x, N, a, 10);
	BaseP(y, N, b, 10);
	
	
    for(i=0;i<N;++i) 
	    printf("%d ",x[i]);
    printf("\n");
    
 	SFT(X ,x ,N ,w ,P);
	for(i=0;i<N;++i) 
	    printf("%d ",X[i]);
	printf("\n");
	system("pause");
		   
	for(i=0;i<N;++i) 
	    printf("%d ",y[i]);
	printf("\n");
	  
	SFT(Y ,y ,N ,w ,P);
	for(i=0;i<N;++i) printf("%d ",Y[i]);
	printf("\n");
	system("pause");

	for(i=0;i<N;++i) 
	    X[i] = (((X[i] * Y[i]) % P)*n) % P;
	for(i=0;i<N;++i) printf("X[%d]=%d\n",i,X[i]);	    
	
	//inverse
	SFT(x ,X ,N ,W ,P );
	for(i=0;i<N;++i) printf("%d ",x[i]);  
	
	return 0;
}

int BaseP(int *x, int N, int a, int P)
{
	int i;
	for(i=0;i<N;i++)
	{
		x[i] = a % P;
		a -= x[i];
		a /= P;
	}
	return 0;
} 

int SFT(int *X, int *x, int N, int w, int P)
{
	int i,j, w0, wk;
	wk = 1;
	for(i=0;i<N;i++)
	{
		X[i] = 0; w0 = 1;
		for(j=0;j<N;++j)
		{
			//printf("%d ", w0);
			X[i] = X[i] + (w0*x[j] % P);
			w0 = (w0*wk % P);	 
		} 
		//printf("\n");
		wk = wk*w % P;
		X[i] = (X[i] % P);
	}
}
