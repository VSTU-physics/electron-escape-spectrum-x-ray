#include <stdio.h>
#include <math.h>
#include "calculations.h"

/*
	spe решает уравнение:
	\[
		\frac{\partial f}{\partial t} =
		p(z, t) \frac{\partial^2 f}{\partial z^2} + g(z, t),
	\]
	при граничных условиях:
	\[
		\a_1(t) f + \b_1(t) \frac{\partial f}{\partial z} \Bigg|_{z = z_1} = \c_1(t),
	\]
	\[
		\a_2(t) f + \b_2(t) \frac{\partial f}{\partial z} \Bigg|_{z = z_2} = \c_2(t),
	\]
	и произвольных начальных условиях неявным методом на одном шаге. Ошибка $O(h)$. A, B, C -- диагонали
	матрицы нижняя, собственно диагональ, верхняя. D -- правая часть уравнения.
*/

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
         double c_2t)
{
	double A[N], B[N], C[N], D[N], dz1, dz2, dza, dq;
	dz1 = z[1] - z[0];
	A[0] = 0;
	B[0] = a_1t - b_1t/dz1;
	C[0] = b_1t/dz1;
	D[0] = c_1t;
	for (int i = 1; i<N-1; i++){
        dz2 = dz1;
		dz1 = z[i+1] - z[i];
        dza = (dz1 + dz2)*0.5;
		A[i] = dt*pt/dz2/dza;
		B[i] = -1 - 2*pt/dz1/dz2*dt;
		C[i] = dt*pt/dz1/dza;
		D[i] = - gt*dt - f_prev[i];
	}
	A[N-1] = b_2t/(z[N-1] - z[N-2]);
	B[N-1] = a_2t - b_2t/(z[N-1] - z[N-2]);
	C[N-1] = 0;
	D[N-1] = c_2t;
	for (int i = 1; i<N; i++){
		B[i] -= A[i]/B[i-1]*C[i];
		D[i] -= A[i]/B[i-1]*D[i-1];
	}
	f[N-1] = D[N-1]/B[N-1];
	for (int i = N-2; i>=0; i--){
		f[i] = (D[i] - C[i] * f[i+1])/B[i];
	}
}

/*
	eval_cubic_spline заносит значения в точках xs[j] в массив ys[j], j = 0, 1, 2, ..., M-1 функции заданной кубиче-
	скими сплайнами для таблицы (x[i], y[i]), i = 0, 1, 2, ... , N-1. int_cubic_spline ищет интегралы, cubic_spline
	ищет коэффициенты сплайна:
	\[
		S_i(x) = a_i + b_i (x - x_i) + \frac{1}{2} c_i (x - x_i)^2 + \frac{1}{6} d_i (x - x_i)^3
		\quad x\in[x_i, x_{i+1}].
	\]
*/

void cubic_spline(double *x, double *y, int N, double *a, double *b, double *c, double *d)
{
	//double a[N], b[N], c[N], d[N];
	double A[N-2], B[N-2], C[N-2], D[N-2];
	double dx[N-1], dydx[N-1];
	c[0] = 0;
	c[N-1] = 0;
	dx[0] = x[1] - x[0];
	dydx[0] = (y[1] - y[0])/dx[0];
	for (int i = 0; i<(N-2); i++){
		dx[i+1] = x[i+2] - x[i+1];
		dydx[i+1] = (y[i+2] - y[i+1])/dx[i+1];
		A[i+1] = dx[i];
		B[i] = 2*(dx[i] + dx[i+1]);
		C[i] = dx[i];
		D[i] = 6*(dydx[i+1] - dydx[i]);
	}
	for (int i = 1; i<N-2; i++){
		B[i] -= A[i]/B[i-1]*C[i];
		D[i] -= A[i]/B[i-1]*D[i-1];
	}
	c[N-2] = D[N-3]/B[N-3];
	for (int i = N-4; i>=0; i--){
		c[i+1] = (D[i] - C[i] * c[i+2])/B[i];
	}
	for (int i = 0; i<N-1; i++){
		a[i] = y[i];
		b[i] = dydx[i] - c[i]*dx[i]/3 - c[i+1]*dx[i]/6;
		d[i] = (c[i+1] - c[i])/dx[i];
	}
	d[N-1] = 0;
	b[N-1] = 0;
	a[N-1] = y[N-1];
}

void eval_cubic_spline(double *xs, double *ys, int M, double *x, double *y, int N)
{
	double a[N], b[N], c[N], d[N];
	cubic_spline(x, y, N, a, b, c, d);
	//Пусть xs и x сортированы по возрастанию
	int j = 0;
	double h;
	for (int i = 0; i<M; i++){
		//printf("%e %e %e\n", xs[i], x[0], x[N-1]);
		if ((xs[i]>=x[0])&&(xs[i]<x[N-1])){
			while (!((xs[i]<x[j+1])&&(xs[i]>=x[j]))&&(j<N)){
				j++;
			}
			h = xs[i] - x[j];
			ys[i] = a[j] + b[j]*h + c[j]*h*h/2 + d[j]*h*h*h/6;
		}
		if (xs[i] == x[N-1]){
			ys[i] = y[N-1];
		}
	}
}

double int_cubic_spline(double la, double lb, double *x, double *y, int N)
{
	//Пусть x сортирован по возрастанию, la<lb
	if ((la>=x[0])&&(la<=x[N-1])&&(lb>=x[0])&&(lb<=x[N-1])){
		double a[N], b[N], c[N], d[N], h;
		cubic_spline(x, y, N, a, b, c, d);
		double result = 0.;
		for (int j = 1; j<N-1; j++){
			if ((la<=x[j])&&(lb>=x[j+1])) {
				h = x[j+1] - x[j];
				result+=a[j]*h + b[j]*h*h/2 + c[j]*h*h*h/6 + d[j]*h*h*h*h/24;
			}
			if ((la>x[j])&&(la<x[j+1])) {
				h = x[j+1] - la;
				result+=a[j]*h + b[j]*h*h/2 + c[j]*h*h*h/6 + d[j]*h*h*h*h/24;
			}
			if ((lb>x[j])&&(lb<x[j+1])) {
				h = lb - x[j];
				result+=a[j]*h + b[j]*h*h/2 + c[j]*h*h*h/6 + d[j]*h*h*h*h/24;
			}
		}
		return result;
	} else {
		return NAN;
	}
}


void test_spe()
{
//	int N = 1000, M = 100;
//	double t = 0., dt = 1e-2, *phi, zmax = 3.1416;
//	double dz = zmax/(N-1), z[N];
//	phi = new double[N*M];
//	for (int i = 0; i<N; i++){
//		z[i] = (i==0) ? 0 : z[i-1]+dz;
//		phi[i] = sin(z[i]);
//	}
//	for (int i = 1; i<M; i++){
//		t+=dt;
//		spe(phi+i*N, phi+(i-1)*N, z, N, t, dt);
//	}
//	FILE *file;
//	file = fopen("results.txt", "w");
//	for (int i = 0; i<N; i++){
//		fprintf(file, "%f %f %f %f\n", z[i], phi[i], *(phi+N*50+i), *(phi+N*(M-1)+i));
//	}
//	fclose(file);
}

void test_spline()
{
	int N = 1000, M = 30;
	double x[N], xx[M], y[N], yy[M], x_min = 0, x_max = 10;
	for (int i = 0; i<M; i++){
		xx[i] = (x_max - x_min)/(M - 1)*i + x_min;
		yy[i] = sin(xx[i]);
	}
	FILE *file;
	file = fopen("results.txt", "w");
	for (int i = 0; i<N; i++){
		x[i] = (x_max - x_min)/(N - 1)*i + x_min;
	}
	double pi = 3.14169;
	eval_cubic_spline(x, y, N, xx, yy, M);
	for (int i = 0; i<N; i++){
		fprintf(file, "%f %f %f\n", x[i], y[i], sin(x[i]));
	}
	fclose(file);
	printf("%f", int_cubic_spline(pi/2, 2*pi, xx, yy, M));
}


