#include <iostream>
#include <string>
#include <cstring>
#include <cmath>

using namespace std;

int main(void)
{
    float a[20][20],x[20],e,big,temp,relerror,sum;
    int n,i,j,maxit,itr;
    char ch;

    printf("\n\nDar la cantidad de incógnitas de la ecuación: ");
    scanf("%d",&n);

    for(i=1;i<=n;i++)
    {
        printf("\n\nDa todos los coeficientes de la ecuacion %d: \n",i);
        for(j=1;j<=n+1;j++)
          scanf("%f",&a[i][j]);
    }
    printf("\n\nDar el error relativo y el numero de iteraciones:  \n");
    scanf("%f%d",&e,&maxit);

    for(i=1;i<=n;i++)
      x[i]=0;
    for(itr=1;itr<=maxit;itr++)
    {
        big=0;
        for(i=1;i<=n;i++)
        {
            sum=0;
            for(j=1;j<=n;j++)
             {
               if(i!=j)
                sum=sum+a[i][j]*x[j];
             }
            temp=(a[i][n+1]-sum)/a[i][i];
            relerror=fabs((x[i]-temp)/temp);
            if(relerror>big)
              big=relerror;
            x[i]=temp;
        }
        if(big<=e)
        {
          printf("Converge en una solucion en %d iteraciones\n",itr);
          for(i=1;i<=n;i++)
            printf("\n%.4f\t",x[i]);
            exit(1);
        }

    }
    printf("No converge en %d iteraciones \n",maxit);

    return 0;

 }
