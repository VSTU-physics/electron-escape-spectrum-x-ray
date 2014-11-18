#pragma once

void spe(double *f,
         double *f_prev,
         double *z,
         int N,
         double dt,
         double pt,
         double gt,
         double a_1t,
         double b_1t,
         double c_1t,
         double a_2t,
         double b_2t,
         double c_2t);

void eval_cubic_spline(double *xs, double *ys, int M, double *x, double *y, int N);

double int_cubic_spline(double la, double lb, double *x, double *y, int N);


void test_spe();
void test_spline();
