#include<stdio.h>
#include<math.h>

#define ERROR -1
#define EPSILON 1E-5
#define MAXREPT 100

//函数1

double Function1(double x)
{
    return pow((x - 1), 3)  - x * x + x;
}

double Derivative_Function1(double x)
{
    return 3 * (x - 1) * (x - 1) - 2 * x + 1;
}

//函数2

double Function2(double x)
{
    return pow(sin(x), 3) + pow(cos(x), 3);
}

double Derivative_Function2(double x)
{
    return 3 * sin(x) * sin(x) * cos(x) - 3 * cos(x) * cos(x) * sin(x);
}

    /*
    二分法求解函数的根
    参数分别为函数Function调用函数指针，表征函数；a，b分别为二分法初始上下界
    若解存在，则返回值为近似解，否则返回-1
    */

double Bisection_method(double Function(double ), double a, double b)
{
    int i = 0;

    if(Function(a) * Function(b) >= 0)
        return ERROR;
    
    while(fabs(a - b) > EPSILON)
    {
        double mid = (a + b) / 2;
        i++;

        printf("iteration %d : f(%.12f) = %.12f\n", i, mid, Function(mid));
        //printf("%.12f\n",fabs(Function(mid)));

        if(fabs(Function(mid)) < EPSILON)
            return mid;

        if(Function(a) * Function(mid) < 0)
        {
            b = mid;
        }
        else if(Function(mid) * Function(b) < 0)
        {
            a = mid;
        }
    }
    return (a + b) / 2;
}

    /*
    牛顿迭代法求近似根
    参数分别为：Function调用预设的函数的计算公式；Derivative_Function调用预设的导函数；x0为初始迭代值
    若解存在，则返回值为近似解，否则返回-1
    */

double Newtons_method(double Function(double ), double Derivative_Function(double ), double x0)
{
    int k;
    for(k = 1; k < MAXREPT; k++)
    {
        double x1 = x0 - Function(x0) / Derivative_Function(x0);

        printf("iteration %d : f(%.12f) = %.12f\n", k, x1, Function(x1));
        //printf("%.12f\n",fabs(Function(x1)));
        if(fabs(x1 - x0) < EPSILON)
        {
            return x1;
        }

        x0 = x1;

    }
    return ERROR;
}


int main()
{
    printf("\nFunction1 in Bisection method in [2,3]\n");
    Bisection_method(Function1, 2, 3);
    printf("\nFunction1 in Bisection method in [2.2,3]\n");
    Bisection_method(Function1, 2.2, 3);
    printf("\nFunction1 in Bisection method in [2,2.8]\n");
    Bisection_method(Function1, 2, 2.8);

    printf("\nFunction2 in Bisection method in [2,3]\n");
    Bisection_method(Function2, 2, 3);
    printf("\nFunction2 in Bisection method in [2.2,3]\n");
    Bisection_method(Function2, 2.2, 3);
    printf("\nFunction2 in Bisection method in [2,2.8]\n");
    Bisection_method(Function2, 2, 2.8);

    printf("\nFunction1 in Newton's method with initial value 2.25\n");
    Newtons_method(Function1,Derivative_Function1,2.25);
    printf("\nFunction1 in Newton's method with initial value 2.5\n");
    Newtons_method(Function1,Derivative_Function1,2.5);
    printf("\nFunction1 in Newton's method with initial value 2.75\n");
    Newtons_method(Function1,Derivative_Function1,2.75);

    printf("\nFunction2 in Newton's method with initial value 2.25\n");
    Newtons_method(Function2,Derivative_Function2,2.25);
    printf("\nFunction2 in Newton's method with initial value 2.5\n");
    Newtons_method(Function2,Derivative_Function2,2.5);
    printf("\nFunction2 in Newton's method with initial value 2.75\n");
    Newtons_method(Function2,Derivative_Function2,2.75);


    return 0;
}
