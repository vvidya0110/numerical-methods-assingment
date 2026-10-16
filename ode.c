//Using different ODE solvers to find the angular velocity in the case of constant acceleration. 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef double (*ODE_RHS)(double x, double y);

/* ---------- Scalar ODE solvers: y' = f(x, y) ---------- */

double ode_euler(ODE_RHS f, double x0, double y0, double h, int num_steps) {
    double y = y0;
    for (int i = 0; i < num_steps; i++) {
        double x = x0 + i * h;
        y += h * f(x, y);
    }
    return y;
}

double ode_heun(ODE_RHS f, double x0, double y0, double h, int num_steps) {
    double y = y0;
    for (int i = 0; i < num_steps; i++) {
        double x = x0 + i * h;
        double k1 = f(x, y);
        double k2 = f(x + h, y + h * k1);
        y += 0.5 * h * (k1 + k2);
    }
    return y;
}

double ode_midpoint(ODE_RHS f, double x0, double y0, double h, int num_steps) {
    double y = y0;
    for (int i = 0; i < num_steps; i++) {
        double x = x0 + i * h;
        double k1 = f(x, y);
        double k2 = f(x + 0.5 * h, y + 0.5 * h * k1);
        y += h * k2;
    }
    return y;
}

double ode_rk4(ODE_RHS f, double x0, double y0, double h, int num_steps) {
    double y = y0;
    for (int i = 0; i < num_steps; i++) {
        double x = x0 + i * h;
        double k1 = f(x, y);
        double k2 = f(x + 0.5 * h, y + 0.5 * h * k1);
        double k3 = f(x + 0.5 * h, y + 0.5 * h * k2);
        double k4 = f(x + h, y + h * k3);
        y += (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
    }
    return y;
}


void euler_step(double *theta, double *v, double a, double R, double dt) {
    (void)R;
    *theta += *v * dt;
    *v += a * dt;
}

void heun_step(double *theta, double *v, double a, double R, double dt) {
    (void)R;
    double v_predict = *v + a * dt;
    *theta += 0.5 * (*v + v_predict) * dt;
    *v = v_predict;                     
}

void midpoint_step(double *theta, double *v, double a, double R, double dt) {
    (void)R;
    double v_mid = *v + 0.5 * a * dt;
    *theta += v_mid * dt;
    *v += a * dt;
}

void rk4_step(double *theta, double *v, double a, double R, double dt) {
    (void)R;

    double k1_theta = *v;
    double k2_theta = *v + 0.5 * dt * a;
    double k3_theta = *v + 0.5 * dt * a;
    double k4_theta = *v + dt * a;

    *theta += (dt / 6.0) * (k1_theta + 2 * k2_theta + 2 * k3_theta + k4_theta);
    *v += a * dt;                        /* RK4 on constant a reduces to this */
}



double normalize_angle(double angle) {
    double twopi = 2.0 * M_PI;
    angle = fmod(angle, twopi);
    if (angle < 0) angle += twopi;
    if (angle >= twopi) angle = 0.0;     
    return angle;
}


double trajectory_rmse(const double *t1, const double *t2, int num_points) {
    if (num_points <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < num_points; i++) {
        double d = t1[i] - t2[i];
        sum += d * d;
    }
    return sqrt(sum / num_points);
}


double trajectory_rmse_angle(const double *t1, const double *t2, int num_points) {
    if (num_points <= 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < num_points; i++) {
        double d = fmod(t1[i] - t2[i] + M_PI, 2.0 * M_PI);
        if (d < 0) d += 2.0 * M_PI;
        d -= M_PI;
        sum += d * d;
    }
    return sqrt(sum / num_points);
}



static double f_test(double x, double y) { (void)x; return y; }  

int main(void) {
    double h = 0.1;
    int n = 10;
    printf("Exact e^1   = %.8f\n", exp(1.0));
    printf("Euler       = %.8f\n", ode_euler(f_test, 0, 1, h, n));
    printf("Heun        = %.8f\n", ode_heun(f_test, 0, 1, h, n));
    printf("Midpoint    = %.8f\n", ode_midpoint(f_test, 0, 1, h, n));
    printf("RK4         = %.8f\n", ode_rk4(f_test, 0, 1, h, n));

    double theta = 0.0, v = 1.0;
    rk4_step(&theta, &v, 2.0, 0.0, 0.5);
    printf("theta=%.4f v=%.4f (exact: 0.75, 2.0)\n", theta, v);
    printf("normalize(-0.5) = %.4f\n", normalize_angle(-0.5));
    return 0;
}
