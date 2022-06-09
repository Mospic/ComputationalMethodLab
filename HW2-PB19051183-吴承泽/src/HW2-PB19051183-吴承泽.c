#include <stdio.h>
#include <string.h>
#include <math.h>


//存放矩阵系数的二维数组如下所示
double A_1[5][5] = { 
                  1.0 / 9 , 1.0 / 8 , 1.0 / 7 , 1.0 / 6 , 1.0 / 5 ,
                  1.0 / 8 , 1.0 / 7 , 1.0 / 6 , 1.0 / 5 , 1.0 / 4 ,
                  1.0 / 7 , 1.0 / 6 , 1.0 / 5 , 1.0 / 4 , 1.0 / 3 ,
                  1.0 / 6 , 1.0 / 5 , 1.0 / 4 , 1.0 / 3 , 1.0 / 2 ,
                  1.0 / 5 , 1.0 / 4 , 1.0 / 3 , 1.0 / 2 , 1.0 / 1
                };
    

double b_1[5][1] = { 1 , 
                     1 ,
                     1 , 
                     1 , 
                     1 
                    };



double A_2[4][4] = {
                     7.2 , 2.3 , -4.4 , 0.5 ,
                     1.3 , 6.3 , -3.5 , 2.8 ,
                     5.6 , 0.9 , 8.1 , -1.3 ,
                     1.5 , 0.4 , 3.7 , 5.9 ,
};

double b_2[4][1] = {
                    15.1 ,
                    1.8 ,
                    16.6 ,
                    36.9
};
    /*
        用于展示系数矩阵，Display_Matrix_A1用来展示5*5的矩阵，Display_Matrix_A2用来展示4*4的矩阵，Display_Matrix用来展示n维向量
    */
void Display_Matrix_A1(double M[][5])
{
    int i, j;
    printf("\nThe matrix is :\n");
    for(i = 0; i < 5; i++)
    {
        for(j = 0 ; j < 5; j++)
        {
            printf("%lf ", M[i][j]);
        }
        printf("\n");
    }
}


void Display_Matrix_A2(double M[][4])
{
    int i, j;
    printf("\nThe matrix is :\n");
    for(i = 0; i < 4; i++)
    {
        for(j = 0 ; j < 4; j++)
        {
            printf("%lf ", M[i][j]);
        }
        printf("\n");
    }
}


void Display_Matrix(double b[][1], int n)
{
    int i;
    for(i = 0; i < n; i++)
    {
        printf(" %lf ", b[i][0]);
    }
    printf("\n");
}



    /*
      计算所得到解的误差大小,M为系数矩阵，S为计算出的解，B为常系数矩阵,返回值为该误差的2-范数
    */
double Calculate_Error_A1(double M[][5], double S[][1], double B[][1])
{
    double temp;
    double Error = 0;
    int i, j;
    for(i = 0 ; i < 5 ; i++)
    {
        temp = 0;
        for(j = 0 ; j < 5 ; j++)
        {
            temp += M[i][j] * S[j][0];
        }
        Error += pow((temp - B[i][0]),2);
    }
    return sqrt(Error);
    //展示Error数组
    
}


double Calculate_Error_A2(double M[][4], double S[][1], double B[][1])
{
    double temp;
    double Error = 0;
    int i, j;
    for(i = 0 ; i < 4 ; i++)
    {
        temp = 0;
        for(j = 0 ; j < 4 ; j++)
        {
            temp += M[i][j] * S[j][0];
        }
        Error += pow((temp - B[i][0]),2);
    }
    return sqrt(Error);
    //展示Error数组
}


    /*
      Gauss列主元消元法：输入为系数矩阵A，常数矩阵b,得到列主元消元后的上三角矩阵U，以及输出x矩阵
    */
int GaussEliminationWithPartialPivoting_A1()
{

    int i, j, k;
    double temp;
    double temp_Matrix[5][5],temp_b[5][1];
    memcpy(temp_Matrix, A_1, 5 * 5 *sizeof(double));
    memcpy(temp_b, b_1, 5 * 1 *sizeof(double));

    for(i = 0; i < 5 ; i++)
    {
        k = i;

        for(j = i + 1; j < 5; j++)
        {
            if(fabs(A_1[k][i]) < fabs(A_1[j][i]))
                k = j;
        }

        for(j = i; j < 5; j++)
        {
            temp = A_1[i][j];
            A_1[i][j] = A_1[k][j];
            A_1[k][j] = temp;
        }

        temp = b_1[i][0];
        b_1[i][0] = b_1[k][0];
        b_1[k][0] = temp;

        for(j = i + 1; j < 5 ; j++)
        {
            double lamda = A_1[j][i] / A_1[i][i];
            for(k = i; k < 5; k++)
            {
                A_1[j][k] = A_1[j][k] - lamda * A_1[i][k];
            }
            b_1[j][0] = b_1[j][0] - lamda * b_1[i][0];
        }

    }

    Display_Matrix_A1(A_1);
    //此处输出的矩阵值为上三角矩阵的值，下面计算解

    for(i = 4; i >= 0; i--)
    {
        for(j = 4; j > i ; j--)
        {
            b_1[i][0] = b_1[i][0] - A_1[i][j] * b_1[j][0];
        }
        b_1[i][0] = b_1[i][0] / A_1[i][i];
    }

    printf("\nThe Solution is :\n");
    Display_Matrix(b_1, 5);
    //此处输出为解的值


    printf("\nError is : %lf \n", Calculate_Error_A1(temp_Matrix, b_1, temp_b));
    //此处计算解的误差,用2-范数表示
}

int GaussEliminationWithPartialPivoting_A2()
{
    int i, j, k;
    double temp;
    double temp_Matrix[4][4],temp_b[4][1];
    memcpy(temp_Matrix, A_2, 4 * 4 *sizeof(double));
    memcpy(temp_b, b_2, 4 * 1 *sizeof(double));

    for(i = 0; i < 4 ; i++)
    {
        k = i;

        for(j = i + 1; j < 4; j++)
        {
            if(fabs(A_2[k][i]) < fabs(A_2[j][i]))
                k = j;
        }

        for(j = i; j < 4; j++)
        {
            temp = A_2[i][j];
            A_2[i][j] = A_2[k][j];
            A_2[k][j] = temp;
        }

        temp = b_2[i][0];
        b_2[i][0] = b_2[k][0];
        b_2[k][0] = temp;

        for(j = i + 1; j < 4 ; j++)
        {
            double lamda = A_2[j][i] / A_2[i][i];
            for(k = i; k < 4; k++)
            {
                A_2[j][k] = A_2[j][k] - lamda * A_2[i][k];
            }
            b_2[j][0] = b_2[j][0] - lamda * b_2[i][0];
        }

    }

    Display_Matrix_A2(A_2);
    //此处输出的矩阵值为上三角矩阵的值，下面计算解

    for(i = 3; i >= 0; i--)
    {
        for(j = 3; j > i ; j--)
        {
            b_2[i][0] = b_2[i][0] - A_2[i][j] * b_2[j][0];
        }
        b_2[i][0] = b_2[i][0] / A_2[i][i];
    }

    printf("\nThe Solution is :\n");
    Display_Matrix(b_2, 4);
    //此处输出为解的值

    printf("\nError is : %lf \n", Calculate_Error_A2(temp_Matrix, b_2, temp_b));
    //此处计算解的误差,用2-范数表示
}

    /*
      Doolittle分解
    */


int Doolittle_A1()
{
    int i, j, k, r;
    double temp;
    double U[5][5], L[5][5];

    for(i = 0; i < 5; i++)
    {
        for(j = 0; j < 5; j++)
        {
            L[i][j] = 0;
            U[i][j] = 0;
        }
    }
    for(i = 0; i < 5; i++)
    {
        L[i][i] = 1;
    }

    for(k = 0; k < 5; k++)
    {
        for(j = k; j < 5; j++)
        {
            temp = 0;
            for(r = 0; r < k; r++)
            {
                temp += L[k][r] * U[r][j];
            }
            U[k][j] = A_1[k][j] - temp;
        }

        for(i = k + 1; i < 5; i++)
        {
            temp = 0;
            for(r = 0; r < k; r++)
            {
                temp += L[i][r] * U[r][k];
            }
            L[i][k] = (A_1[i][k] - temp) / U[k][k];
        }
    }

    //此处计算出L、U矩阵，并打印出来
    printf("L Matrix");
    Display_Matrix_A1(L);
    printf("\n");
    printf("U Matrix");
    Display_Matrix_A1(U);
    printf("\n");

    double y[5][1], x[5][1];

    for(i = 0; i < 5; i++)
    {
        y[i][0] = 0;
        x[i][0] = 0;
    }

    for(i = 0; i < 5; i++)
    {
        temp = 0;
        for(j = 0; j < i; j++)
        {
            temp += L[i][j] * y[j][0];
        }
        y[i][0] = b_1[i][0] - temp;
    }

    for(i = 4; i >= 0; i--)
    {
        temp = 0;
        for(j = i + 1; j < 5; j++)
        {
            temp += U[i][j] * x[j][0];
        }
        x[i][0] = (y[i][0] - temp) / U[i][i];
    }

    printf("\nThe Solution is :\n");
    Display_Matrix(x, 5);
    //打印答案

    printf("\nError is : %lf \n", Calculate_Error_A1(A_1, x, b_1));
    //此处计算解的误差,用2-范数表示
}




int Doolittle_A2()
{
    int i, j, k, r;
    double temp;
    double U[4][4], L[4][4];

    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 4; j++)
        {
            L[i][j] = 0;
            U[i][j] = 0;
        }
    }
    for(i = 0; i < 4; i++)
    {
        L[i][i] = 1;
    }

    for(k = 0; k < 4; k++)
    {
        for(j = k; j < 4; j++)
        {
            temp = 0;
            for(r = 0; r < k; r++)
            {
                temp += L[k][r] * U[r][j];
            }
            U[k][j] = A_2[k][j] - temp;
        }

        for(i = k + 1; i < 4; i++)
        {
            temp = 0;
            for(r = 0; r < k; r++)
            {
                temp += L[i][r] * U[r][k];
            }
            L[i][k] = (A_2[i][k] - temp) / U[k][k];
        }
    }

    //此处计算出L、U矩阵，并打印出来
    printf("L Matrix");
    Display_Matrix_A2(L);
    printf("\n");
    printf("U Matrix");
    Display_Matrix_A2(U);
    printf("\n");

    double y[4][1], x[4][1];

    for(i = 0; i < 4; i++)
    {
        y[i][0] = 0;
        x[i][0] = 0;
    }

    for(i = 0; i < 4; i++)
    {
        temp = 0;
        for(j = 0; j < i; j++)
        {
            temp += L[i][j] * y[j][0];
        }
        y[i][0] = b_2[i][0] - temp;
    }

    for(i = 3; i >= 0; i--)
    {
        temp = 0;
        for(j = i + 1; j < 4; j++)
        {
            temp += U[i][j] * x[j][0];
        }
        x[i][0] = (y[i][0] - temp) / U[i][i];
    }

    printf("\nThe Solution is :\n");
    Display_Matrix(x, 4);
    //打印答案


    printf("\nError is : %lf \n", Calculate_Error_A2(A_2, x, b_2));
    //此处计算解的误差,用2-范数表示
}


int main()
{
    //打印Doolittle分解法产生的解
    printf("\n\nDoolittle start\n\n");
    printf("\nMatrix 1\n");
    Doolittle_A1();
    printf("\nMatrix 2\n");
    Doolittle_A2();
    printf("\n\nDoolittle over\n\n");

    //打印Gauss列主元法产生的解
    printf("\n\nGauss start\n\n");
    printf("\nMatrix 1\n");
    GaussEliminationWithPartialPivoting_A1();
    printf("\nMatrix 2\n");
    GaussEliminationWithPartialPivoting_A2();
    printf("\n\nGauss over\n\n");

    return 0;
    
}