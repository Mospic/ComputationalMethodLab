#include <stdio.h>
#include <string.h>
#include <math.h>

#define INF -65535

//存放矩阵系数的二维数组如下所示
double A_1[5][5] = { 
                  1.0 / 9 , 1.0 / 8 , 1.0 / 7 , 1.0 / 6 , 1.0 / 5 ,
                  1.0 / 8 , 1.0 / 7 , 1.0 / 6 , 1.0 / 5 , 1.0 / 4 ,
                  1.0 / 7 , 1.0 / 6 , 1.0 / 5 , 1.0 / 4 , 1.0 / 3 ,
                  1.0 / 6 , 1.0 / 5 , 1.0 / 4 , 1.0 / 3 , 1.0 / 2 ,
                  1.0 / 5 , 1.0 / 4 , 1.0 / 3 , 1.0 / 2 , 1.0 / 1
                };
    

double A_2[4][4] = {
                     4 , -1 , 1 , 3 ,
                     16 , -2 , -2 , 5 ,
                     16 , -3 , -1 , 7 ,
                     6 , -4 , 2 , 9 ,
};


void Display_Matrix(double b[][1], int n)
{
    int i;
    for(i = 0; i < n; i++)
    {
        printf(" %lf ", b[i][0]);
    }
    printf("\n");
}

//获取向量的无穷范数
double Get_Max(double X[][1], int n)
{
    int i;
    double max = INF;
    for(i = 0; i < n; i++)
    {
        if(fabs(X[i][0]) > max)
            max = fabs(X[i][0]);
    }
    return max;
}


//行列数为5的Doolittle矩阵分解函数

int Doolittle_A1(double A_1[][5], double Y[][1], double X[][1])
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
    //此处LU矩阵分解完毕
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
        y[i][0] = Y[i][0] - temp;
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

    for(i = 0; i < 5; i++)
    {
        X[i][0] = x[i][0];
    }
}



//通过反幂法迭代得到特征值与特征向量

double Get_eig_A1(double A_1[][5],double Feature_Vector[][1])
{
    int i,k = 0;
    double X[5][1] = { 1 , 
                     1 ,
                     1 , 
                     1 , 
                     1 
                    };
    double X_next[5][1] = { 1 , 
                     1 ,
                     1 , 
                     1 , 
                     1 
                    };
    double Y[5][1];
    //暂且初始化为-65535
    double eigenvalue, eigenvalue_next = INF;
    
    do
    {
        eigenvalue = eigenvalue_next;
        for(i = 0; i < 5; i++)
        {
            X[i][0] = X_next[i][0];
        }
        for(i = 0; i < 5; i++)
        {
            Y[i][0] = X[i][0] / Get_Max(X,5);
        }
        Doolittle_A1(A_1,Y,X_next);
        eigenvalue_next = Get_Max(X_next,5);

        printf("X(%d) is:\n",k);
        Display_Matrix(X,5);
        printf("Y(%d) is:\n",k);
        Display_Matrix(Y,5);
        printf("eigenvalue is : %lf\n\n", eigenvalue_next);

        k++;
    }while(fabs(eigenvalue_next - eigenvalue) >= 1e-5);

    
    for(i = 0; i < 5; i++)
    {
        Feature_Vector[i][0] = Y[i][0];
    }
    return 1 / eigenvalue;
}

//行列数为4的Doolittle矩阵分解函数

int Doolittle_A2(double A_2[][4], double Y[][1], double X[][1])
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
        y[i][0] = Y[i][0] - temp;
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
    for(i = 0; i < 4; i++)
    {
        X[i][0] = x[i][0];
    }
}

double Get_eig_A2(double A_2[][4],double Feature_Vector[][1])
{
    int i, k = 0;
    double X[4][1] = { 1 , 
                     1 ,
                     1 , 
                     1 
                    };
    double X_next[4][1] = { 1 , 
                     1 ,
                     1 , 
                     1 
                    };
    double Y[4][1];
    //暂且初始化为-65535
    double eigenvalue, eigenvalue_next = INF;
    
    do
    {
        eigenvalue = eigenvalue_next;
        for(i = 0; i < 4; i++)
        {
            X[i][0] = X_next[i][0];
        }
        for(i = 0; i < 4; i++)
        {
            Y[i][0] = X[i][0] / Get_Max(X,4);
        }
        Doolittle_A2(A_2,Y,X_next);
        eigenvalue_next = Get_Max(X_next,4);

        printf("X(%d) is:\n",k);
        Display_Matrix(X,4);
        printf("Y(%d) is:\n",k);
        Display_Matrix(Y,4);
        printf("eigenvalue is : %lf\n\n", eigenvalue_next);

        k++;
    }while(fabs(eigenvalue_next - eigenvalue) >= 1e-5);

    
    for(i = 0; i < 4; i++)
    {
        Feature_Vector[i][0] = Y[i][0];
    }
    return 1 / eigenvalue;
}

int main()
{
    int i;
    double Feature_Vector_A1[5][1];
    double min_eigenvalue_A1 = Get_eig_A1(A_1,Feature_Vector_A1);
    printf("\n\n\n The final minimum eigenvalue of A1 is :%lf\n ",min_eigenvalue_A1);
    printf("The Feature vector of A1 is:");
    Display_Matrix(Feature_Vector_A1,5);

    printf("\n\n\n");

    double Feature_Vector_A2[4][1];
    double min_eigenvalue_A2 = Get_eig_A2(A_2,Feature_Vector_A2);
    printf("\n\n\n The final minimum eigenvalue of A2 is :%lf\n ",min_eigenvalue_A2);
    printf("The Feature vector of A2 is:");
    Display_Matrix(Feature_Vector_A1,4);



    return 0;
}