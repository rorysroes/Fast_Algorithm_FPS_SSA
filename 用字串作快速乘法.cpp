#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define MAXCHAR 400

int ReverseOrder(char *a, int N); 
int Char2Int(char *a, int N) ;
int Multiply(char *c, char *a, int Na, char *b, int Nb);
int Find_A_Prime(int Na, int Nb, int *Nf, int *w, int *W, int *Ninv);
int Check_A_Prime(int P);
int Order(int w, int P);
int Inverse_Zp(int w, int P);
int FFT_radix_2(int *X ,int *x ,int w ,int N,int P);

int main()
{
	int i ,j ,Na ,Nb ,Nc ,P ,w ,W ,N ,Ninv ,*x ,*X ,*y ,*Y;

	char a[MAXCHAR] ,b[MAXCHAR] ,c[2*MAXCHAR] ,t;
	clock_t t1, t2 ;
	printf("Input number a and b (MAXIMUM DIGITS:%d)\n例如想計算 1234*5678 請輸入 1234 且前面補0到16位\nEX: 0000000000001234 * 0000000000005678\n",MAXCHAR);
	scanf("%s %s", a, b);
	printf("%s * %s = \n", a, b);

	Na = strlen(a); 
	Nb = strlen(b);	

	ReverseOrder(a,Na);
 	ReverseOrder(b,Nb);

	Char2Int(a,Na);
	Char2Int(b,Nb);	
	for(i=0;i<Na;++i) printf("%d",a[i]);
	printf("\n");	
	for(i=0;i<Nb;++i) printf("%d",b[i]);
	printf("\n");			

	t1 = clock();
	Nc = Multiply(c, a, Na, b, Nb);
	t2 = clock(); 
	printf("%f\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);

	ReverseOrder(c, Nc);
	for(i=0;i<Nc;++i) printf("%d",c[i]);
	printf("\n");

	P = Find_A_Prime(Na,Nb,&N,&w,&W,&Ninv);
	printf("FFT in Z[%d], with w = %d, W = %d, N = %d, Ninv = %d\n", P, w, W, N, Ninv);

	x = (int *) malloc(N*sizeof(int));
	y = (int *) malloc(N*sizeof(int));
	X = (int *) malloc(N*sizeof(int));
	Y = (int *) malloc(N*sizeof(int));
	for(i=0;i<N;++i)
	{
		x[i] = y[i] = X[i] = Y[i] = 0;
	}
	for(i=0;i<Na;++i) x[i] = a[i];
	for(i=0;i<Nb;++i) y[i] = b[i];
	t1 = clock();
	FFT_radix_2(X ,x ,w ,N ,P);	
	FFT_radix_2(Y ,y ,w ,N ,P);	
	for(i=0;i<N;++i) X[i] = ((X[i]*Y[i] % P)*Ninv) % P; 
	//for(i=0;i<N;++i) printf("%d ",X[i]);	
	FFT_radix_2(x ,X , W ,N ,P );
	//for(i=0;i<N;++i) printf("%d ",x[i]);	
	for(i=0;i<N-1;++i)
	{
		x[i+1] += x[i]/10;
		x[i] = x[i] % 10;
		//printf("%d\n",x[i]);
	}
	//for(i=0;i<N;++i) printf("%d ",x[i]);	
	for(i=0;i<Nc;++i) c[i] = x[i];
	t2 = clock();
	ReverseOrder(c, Nc);
	for(i=0;i<Nc;++i) printf("%d",c[i]);
	
	printf("\n");
	printf("%f\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
	printf("\n");		
	
	free(x);
	free(y);
	free(X);
	free(Y);
	return 0;
} 
int ReverseOrder(char *a, int N)
{
	// i <-> j
	// 0 <-> N-1
	// 1 <-> N-2
	// 2 <-> N-3
	// ...
	int i, j;
	char t;
	for(i=0;i<N/2;++i)
	{
		j = N - 1 - i;
		t = a[i];
		a[i] = a[j];
		a[j] = t;
	}
	return 0;
}
int Char2Int(char *a, int N)
{
	int i;
	for(i=0;i<N;++i)
	{
		a[i] -= 48;
	}
	return 0;
}
int Int2Char(char *a, int N)
{
	int i;
	for(i=0;i<N;++i)
	{
		a[i] += 48;
	}
	a[N] = 0;
	return 0;
}
int Multiply(char *c, char *a, int Na, char *b, int Nb)
{
	int i, j;
	for(i=0;i<Na+Nb;++i) c[i] = 0;

	for(i=0;i<Na;i++)
	{
		for(j=0;j<Nb;j++)
		{
			c[i+j] += a[i]*b[j];
			if(c[i+j] >= 10)
			{
				c[i+j+1] += (c[i+j] / 10);  // ?L?????h 
				c[i+j] = c[i+j] % 10;
			}
		}
	}
	if(c[Na+Nb-1]==0) return Na+Nb-1; // highest bit = 0 -> ????1 
	else return Na+Nb;
}
int Find_A_Prime(int Na, int Nb, int *Nf, int *w, int *W, int *Ninv)
{
	int Nm, p, n, N, find_w;
	if(Na > Nb) Nm = Na; else Nm = Nb;  // maximum of Na and Nb
	N = 1;
	while(N<2*Nm)   // Nm = D
	{
		N <<= 1;    // N = N * 2
	}
	*Nf = N;
	p = 81*(N/2)+1;  // p = 81*Nm + 1
	if (p % 2 == 0) p = p + 1;
	while(1)
	{
		if(Check_A_Prime(p)==1)
		{
			if( (p - 1) % N == 0)
			{
				find_w = 0;
				for(n=2;n<p;++n)
				{
					if(Order(n,p)==N)
					{
						find_w = 1;
						break;
					}
				}
				if(find_w == 1) break;
			}
		}
		p += 2;
	}
	*w = n;
	*W = Inverse_Zp(n,p);
	*Ninv = Inverse_Zp(N,p);
	return p;
}
int Check_A_Prime(int P)
{
	// not a prime : 0 .  A prime : 1
	if(P % 2 == 0)
	{
		return 0;
	}
	int n;
	for(n=3; n*n<=P ;n+=2)
	{
		if(P % n == 0)
		{
			return 0;
		}
	}
	return 1;
}
int Order(int w, int p)
{
	int k, w0 = w;
	for(k=1;k<p;++k)
	{
		if(w % p == 1)
			break;
		w = (w*w0) % p;
	}
	return k;
}
int Inverse_Zp(int w, int p)
{
	int k, w0 = w;
	for(k=1;k<p;++k)
	{
		if((w*w0) % p == 1)
			break;
		w = (w*w0) % p;
	}
	return w;
}
int FFT_radix_2(int *X ,int *x ,int w ,int N,int P)
{ 
	if(N==1)
	{
        X[0] = x[0]; 
	    return 0; 	
	}
	int k, n,i;
	int *even, *odd, wk;
    int *even_FT, *odd_FT;

	even = (int *) malloc(N*sizeof(int));
	odd = even + N/2; //(double *) malloc(N/2*sizeof(double));
	even_FT = (int *) malloc(N*sizeof(int));
	odd_FT = even_FT + N/2;//(double *) malloc(N/2*sizeof(double));
	 
	for(n=0;n<N/2;++n)
	{
		even[n] = x[2*n];
		odd[n] = x[2*n+1];
	} 
	
	FFT_radix_2(even_FT ,even , w * w  % P,N/2 ,P);
	FFT_radix_2(odd_FT ,odd , w * w  % P ,N/2 ,P);

	wk=1;
	
    if(N==32){
		for(k=0;k<N/2;++k)
	        {
				X[k] = (even_FT[k] + (wk * odd_FT[k]))  % P;
				//printf("wk=%d\n",wk);
				wk =( wk * (w * (w * (w * (w * (w * ( w *( w * ( w * w % P) * w % P) * w % P) * w % P)* w % P)* w % P)* w % P)* w % P) ) % P;
				//printf("wk=%d\n",wk);
				X[k+N/2] = (even_FT[k] + (wk * odd_FT[k])) % P; 
				
                if(k==14){
				    wk = (w * (w * (w * (w * (w * (w * (w * w * w % P)* w % P)* w % P) * w % P)* w % P)* w % P) * w) % P;					
				}else if(k==13){
				    wk = w * (w * (w * (w * (w * (w * w * w * w % P )* w % P) * w % P)* w % P)* w % P)* w % P;				
				}else if(k==12){
				    wk = w * (w * (w * (w * (w * (w * w * w % P) * w % P) * w % P)* w % P)* w % P)* w % P;				
				}else if(k==11){
				    wk = w * (w * (w * (w * (w * (w * w * w % P)* w % P) * w % P)* w % P)* w % P) % P;				
				}else if(k==10){
				    wk = w * (w * (w * (w * (w * w * w % P) * w % P)* w % P)* w % P)* w % P;					
				}else if(k==9){
				    wk = w * (w * (w * (w * w * w * w % P)* w % P)* w % P) * w % P;				
				}else if(k==8){
				    wk = w * (w * (w * (w * w * w % P) * w % P) * w % P) * w % P;				
				}else if(k==7){
				    wk = w * (w * (w * w * w * w % P)* w % P)* w % P;				
				}else if(k==6){
				    wk = w * (w * (w * w * w % P ) * w % P) * w % P;				
				}else if(k==5){
				    wk = w * (w * w * w * w % P)* w % P;				
				}else if(k==4){
			        wk = w * (w * w * w % P) * w % P;
				}else if(k==3){
				    wk = w * w * w * w % P;
			    }else if(k==2) {
				    wk = w * w * w % P;
				}else if(k==1) {
				    wk = w * w % P;
			    }else{
			    	wk = w % P;	
				}
			}
	}	
	else if(N==16){
		for(k=0;k<N/2;++k)
	        {
				X[k] = (even_FT[k] + (wk * odd_FT[k]))  % P;
				//printf("wk=%d\n",wk);
				wk = (wk * (w * (w * (w * w * w * w % P)* w % P)* w % P) ) % P;
				//printf("wk=%d\n",wk);
				X[k+N/2] = (even_FT[k] + (wk * odd_FT[k])) % P; 
				
				if(k==6){
				    wk = w * (w * (w * w * w % P)* w % P)* w % P;				
				}else if(k==5){
				    wk = w * (w * w * w * w % P)* w % P;				
				}else if(k==4){
			        wk = w * (w * w * w % P)* w % P;
				}else if(k==3){
				    wk = w * w * w * w % P;
			    }else if(k==2) {
				    wk = w * w * w % P;
				}else if(k==1) {
				    wk = w * w % P;
			    }else{
			    	wk = w % P;	
				}
			}
	}else if(N==8){
			for(k=0;k<N/2;++k)
			{
				X[k] = (even_FT[k] + (wk * odd_FT[k]))  % P;
				//printf("wk=%d\n",wk);
				wk = (wk * (w * w * w * w % P) )% P;
				//printf("wk=%d\n",wk);
				X[k+N/2] = (even_FT[k] + (wk * odd_FT[k])) % P; 
				
				if(k==2){
				    wk = w * w * w % P;
			    }else if(k==1) {
				    wk = w * w % P;
			    }else{
			    	wk = w % P;	
				}
		    }
	}else if(N==4){
		for(k=0;k<N/2;++k)
			{
				X[k] = (even_FT[k] + (wk * odd_FT[k]))  % P;
				//printf("wk=%d\n",wk);
				wk = wk * w * w % P;
				//printf("wk=%d\n",wk);
				X[k+N/2] = (even_FT[k] + (wk * odd_FT[k])) % P; 
				wk = w  % P; 
		    }
    }else{
		for(k=0;k<N/2;++k)
			{
				X[k] = (even_FT[k] + (wk * odd_FT[k]))  % P;
				//printf("wk=%d\n",wk);
				wk = wk * w % P;
				//printf("wk=%d\n",wk);
				X[k+N/2] = (even_FT[k] + (wk * odd_FT[k])) % P; 
		    }    	
	}
	//printf("N=%d\n",N);
	//for(i=0;i<N/2;++i) printf("X[%d]=%d\n",i,X[i]);
	//for(i=0;i<N/2;++i) printf("X[%d]=%d\n",i+N/2,X[i+N/2]);
    //printf("\n");
    
    
	free(even);
	free(even_FT);
	return 0;
}
