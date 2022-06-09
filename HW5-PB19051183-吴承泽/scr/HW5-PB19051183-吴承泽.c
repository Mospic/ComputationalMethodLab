#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi 3.1415926
#define epsilon 1e-6
#define M 50

double func(double x, int k)
{
    if(k == 0)
    {
        return log(x);
    }
    else if(k == 1)
        return sqrt(2 - pow(sin(x), 2));
    else
        exit(0);
}

double Range[2][2] = {{1, 2}, {(-1) * pi / 6, 3* pi / 4}};

void Romberg_Method(int n, int t);  //tÎªÐòºÅ

int main()
{
    Romberg_Method(1, 0);

    Romberg_Method(5, 0);

    Romberg_Method(1, 1);

    Romberg_Method(5, 1);  

    return 0;
}


void Romberg_Method(int n, int t)
{
    int i, j, k;
    double R[M][M];
    double temp = 0;
    double h = (Range[t][1] - Range[t][0]) / n;

    temp = 1.0 / 2 * (func(Range[t][0], t) + func(Range[t][1], t));
    for(i = 1;i <= n - 1; i++)
    {
        temp += func(h * i + Range[t][0], t);
    }
    R[1][1] = temp * h;
    for(k = 2; k < M; k++)
    {
        temp = 0;
        for(i = 1; i <= n * pow(2, k - 2); i++)
        {
            temp += func(Range[t][0] + (2 * i - 1) * h / pow(2, k - 1), t);
        }
        temp = temp * h / pow(2, k - 2);
        R[k][1] = (R[k - 1][1] + temp) / 2;

        for(j = 2; j <= k; j++)
        {
            R[k][j] = R[k][j - 1] + (R[k][j - 1] - R[k - 1][j - 1])/ (pow(4, j - 1) - 1);
        }
        if(fabs(R[k][k] - R[k - 1][k - 1]) < epsilon)
            break;
    }

    for(i = 1; i <= k; i++)
    {
        for(j = 1; j <= i; j++)
        {
            printf("%lf ", R[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}