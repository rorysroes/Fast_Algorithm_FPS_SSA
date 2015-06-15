//引入標準輸入輸出 : printf , scanf
#include <stdio.h>
//引入標準函式(system,malloc)
#include <stdlib.h>
//引入字串處理函式
#include <string.h>
//引入時間處理函式(clock, CLOCKS_PER_SEC)
#include <time.h>
//定義在之後程式中 MAXCHAR 就是 100
#define MAXCHAR 400
int ReverseOrder(char *a, int N); 
int Char2Int(char *a, int N) ;
int  Multiply(char *c, char *a, int Na, char *b, int Nb);


int main()
{
	int i ,j ,Na ,Nb ,Nc ,P ,w ,W ,N ,Ninv ,*x ,*X ,*y ,*Y ;
	
	/*int multiplication is not working
	printf("input i , j = ");
	scanf("%d %d",&i,&j);
	printf("%d ",i*j);
	return 0;
	*/
	
	char a[MAXCHAR] ,b[MAXCHAR] ,c[2*MAXCHAR] ,t;
	/*
	printf("input a:");
	scanf("%s",a);
	printf("%s\n",a);
	//輸入 : 12345 , ascii table 電腦裡存的 49 51 52 53 
	for(i=0;i<strlen(a);i++){
		printf("%d ",a[i]);
	}
	return 0;
	*/
	/*
	a = (char *)malloc(MAXCHAR*sizeof(char));
	b = (char *)malloc(MAXCHAR*sizeof(char));
	c = (char *)malloc(2*MAXCHAR*sizeof(char));  
	*/
	clock_t t1, t2 ;
	// 用字串輸入 a,b
	printf("Input number a and b (MAXIMUM DIGITS:%d)\n",MAXCHAR);
	// 輸入字串,%s是自串的意義
	scanf("%s %s", a, b);
	// 印出字串,%s是字串的意義
	printf("%s * %s = \n", a, b);
	//for(i=0;i<MAXCHAR-1;++i) a[i] = b[i] = '9';
	//算出字串a,b的長度(strlen) 
	Na = strlen(a); 
	Nb = strlen(b);	

	//將a,b的位元順序互換. EX : a = 1234 -> 4321
	ReverseOrder(a,Na);
 	ReverseOrder(b,Nb);
	//將a,b從字元轉成0-9的數字 ,可google : ascii code, 48 = '0', ..... 
	Char2Int(a,Na);
	Char2Int(b,Nb);	
	for(i=0;i<Na;++i) printf("%d",a[i]);
	printf("\n");	
	for(i=0;i<Nb;++i) printf("%d",b[i]);
	printf("\n");			
	// 做 a , b的乘法(標準國小做法)
	t1 = clock();
	Nc = Multiply(c, a, Na, b, Nb);
	t2 = clock(); 
	printf("%f\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
	// 把c的位元順序逆序
	ReverseOrder(c, Nc);
	for(i=0;i<Nc;++i) printf("%d",c[i]);
	printf("\n");
	return 0;
	// 把c的數字變回字元
	// Int2Char(c,Nc);
	// 印出c來 
	// printf("%s\n", c); 
	//P = Find 
	
	
	
	
	 
} 

int ReverseOrder(char *a, int N)
{
	// 0 <-> N-1
	// 1 <-> N-2
	// 2 <-> N-3
	// ... 
	int i,j;
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
		a[i]-=48; 
	} 
	return 0; 
} 
int  Multiply(char *c, char *a, int Na, char *b, int Nb)
{
	int i, j;
	
	// 先把 c 歸 0 
	for(i=0;i<Na+Nb;++i) c[i] = 0;
	
	for(i=0;i<Na;i++)
	{
		for(j=0;j<Nb;j++)
		{
			// a的10^3的位置跟b的10^4放在c的10^7位置 
			c[i+j] += a[i]*b[j];
			if(c[i+j] >= 10)
			{
				c[i+j+1] += (c[i+j] / 10); // 無條件捨去 
				c[i+j] = c[i+j] % 10;
			}
		}
	}
	// c 最多有多長 
	if (c[Na+Nb-1]==0) return Na+Nb-1;  // highest bit = 0 -> 位數少 1 
	else return Na+Nb;

}







