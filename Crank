# SAUBHIK-S-CODE
*/
   Program: Write a program to solve the partial differential equation using Crank Nicolson finite difference method 
   */
   
#include<stdio.h>
#include<math.h>
#define PI 3.14159265
float x[20];
void TDMA(float a[20],float b[20],float c[20],float d[20],int m)
{
	int i;
	for(i=2;i<=m-1;i++)
	{
		b[i]=b[i]-(a[i]*c[i-1])/b[i-1];
		d[i]=d[i]-(a[i]*d[i-1])/b[i-1];
	}
	x[m-1]=d[m-1]/b[m-1];
	for(i=m-2;i>=1;i--)
	{
		x[i]=(d[i]-c[i]*x[i+1])/b[i];
	}
}
void main()
{
		int i,j,m,n,i0,j0;   
		float u[20][20],c[20],d[20],h,k,a[20],b[20],r,a1,a2,x1=0.5,t1=0.125;
		h=0.25;
		k=0.0625;
		printf("Enter the lower limit: ");
		scanf("%f",&a1);
		printf("Enter the upper limit: ");
		scanf("%f",&a2);
		n=(a2-a1)/h;
		m=(int)(a2-a1)/h;
		for(i=1;i<=m-1;i++)
		{
			u[i][0]=sin(PI*(a1+i*h));
		}
		r=k/(h*h);
		for(j=0;j<=n;j++)
		{
			u[0][j]=0;
			u[m][j]=0;
			for(i=0;i<m;i++)
	    	{
				a[i]=-r;
				b[i]=2*(1+r);
				c[i]=-r;
				d[i]=r*u[i-1][j]+2*(1-r)*u[i][j]+r*u[i+1][j];
			}
			d[1]-=a[1]*u[0][j];
			d[m-1]-=c[m-1]*u[m][j];
		    TDMA(a,b,c,d,m);
		    for(i=1;i<=m;i++)
		    {
		        u[i][j+1]=x[i];
	        }
		}
		printf("\n The Grid points are:");
		for(j=n;j>=0;j--)
		{
			for(i=0;i<=m;i++)
				printf("\n u(%d,%d)=%0.5f",i,j,u[i][j]);
			printf("\n");
		}
		i0=x1/h;
		j0=t1/k;
		printf("The value of u(%0.5f,%0.5f)=%0.5f",x1,t1,u[i0][j0]);
		printf("(correct upto 5 decimal places)");
}

/*Output
Enter the lower limit: 0 
Enter the upper limit: 1

 The Grid points are:
 u(0,4)=0.00000
 u(1,4)=0.07705
 u(2,4)=0.10704
 u(3,4)=0.07705
 u(4,4)=0.00000

 u(0,3)=0.00000
 u(1,3)=0.13703
 u(2,3)=0.20117
 u(3,3)=0.13703
 u(4,3)=0.00000

 u(0,2)=0.00000
 u(1,2)=0.26531
 u(2,2)=0.34694
 u(3,2)=0.26531
 u(4,2)=0.00000

 u(0,1)=0.00000
 u(1,1)=0.42857
 u(2,1)=0.71429
 u(3,1)=0.42857
 u(4,1)=0.00000

 u(0,0)=0.00000
 u(1,0)=0.50000
 u(2,0)=1.00000
 u(3,0)=1.50000
 u(4,0)=0.00000
The value of u(0.50000,0.12500)=0.34694(correct upto 5 decimal places) */
