//input an integer N and a decimal value for theta
//check that the input is a positive integer
//convert theta into radians
//calculate the trigonometric functions sin(theta) and cos(theta) using the inbuilt functions in math.h
//use functions mysine(N,theta_rad) and mycosine(theta_rad) to calculate the values of sin(theta) and cos(theta) from its Taylor series, up to N terms. e.g. sin(x) = x - x^3/6 + ....

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.14159265359 

double mysine(int N, double x) {
    double term = x, sum = 0;
    for (int k = 0; k < N; k++) {
        sum += term;
        term *= -x * x / ((2*k + 2) * (2*k + 3));
    }
    return sum;
}

double mycosine(int N, double x) {
    double term = 1, sum = 0;
    for (int k = 0; k < N; k++) {
        sum += term;
        term *= -x * x / ((2*k + 1) * (2*k + 2));
    }
    return sum;
}

double main(int argc, char **argv){
  double theta;
  theta=atof(argv[1]);
  int N;
  N=atoi(argv[2]);
  if (N<0){
    printf("enter positive value for N");
    return 0;
    }
  double theta_rad;
  theta_rad=theta*PI/180;
  double cos_taylor=cosine(theta_rad,N);
  double sin_taylor=sine(theta_rad,N);
  double cos_math=cos(theta_rad);
  double sin_math=sin(theta_rad);
  double error_sine=fabs(sin_taylor-sin_math);
  double error_cosine=fabs(cos_taylor-cos_math);
  printf(" N-value\tsine_error\tcosine_error\n");
  printf("%d\t%lf\t%lf\n",N,error_sine,error_cosine);
  return 0;
}
