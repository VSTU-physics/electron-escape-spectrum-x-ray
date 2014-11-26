#include <stdio.h>
#include <math.h>

#include "physics.h"
#include "calculations.h"
#include "parse.h"

double delta(double x, double dx)
{
	return (fabs(x)<dx/2) ? 1./dx : 0.;
}

int main()
{
	int Z = 32;
	//test_parse();
	double *u, *z, *up;
	double l = 1;
	int N = 100; // число точек
	subst_t s;
	auger_t a;

	u = new double[N];
	up = new double[N];
	z = new double[N];
	for (int i = 0; i<N; i++) z[i] = i*l/(N-1);

	double *E, source, *tmp, *f;
	int M = 100; // число точек в спектре
	//ltr = new double[M];
	//eps = new double[M];
	E = new double[M];
	f = new double[M];

	load_subst(Z, "data/subst.pl", s);
	load_auger(Z, "data/aug.pl", a);
    printf("Загружены данные для Z = %d:\n", Z);
    printf("Свойства материала:\n");
    printf("Зарядовое число: %d\n", s.Z);
    printf("Плотность: %f\n", s.rho);
    printf("Молярная масса: %f\n", s.M);
    printf("Макс. энергия Оже-электронов: %f\n", s.Emax);
    printf("Работа выхода: %f\n", s.U0);
    printf("\n");
    printf("Характеристики электронов:\n");
    printf("Кол-во типов: %d\n", a.N);
    printf("Энергии: [");
    for (int i = 0; i < a.N; i++) {
        printf(" %f", a.E[i]);
    }
    printf(" ]\n");
    printf("Вероятности: [");
    for (int i = 0; i < a.N; i++) {
        printf(" %f", a.P[i]);
    }
    printf(" ]\n");
	//load_ltr(ltr, E, M, "data/", s);
	//load_esharp(eps, E, M, "data/", s);

	double dE = a.E[0] / (M - 1);
	FILE *fd;
    fd = fopen("trash.out", "w");
	for (int i = 0; i<M; i++)
    {
        E[i] = i * dE;
        u[i] = a.P[0];
    }

    for (int i = M - 1; i>=0; i--)
    {
        tmp = u;
        u = up;
        up = tmp;

        source = 0;
        for (int k = 1; k<a.N; k++)
        {
            source += a.P[k]*delta(E[i] - a.E[k], dE);
        }
        spe(u, up, z, N, dE,
            - 1./(3*l_tr(s, E[i]) * eps(s, E[i])),
            source,
            I1(s, E[i])/(2 - 3 * I2(s, E[i])),
            - 1./l_tr(s, E[i])/3,
            0,
            1.,
            0.,
            0.
            );
        f[i] = l_tr(s, E[i])*3*I1(s, E[i])/(2 - 3 * I2(s, E[i]))*u[0];
    }
    fclose(fd);
    fd = fopen("data.gp", "w");
    fprintf(fd, "plot '-' with lines\n");
    for (int i = M - 1; i>=0; i--)
    {
        fprintf(fd, "%e %e\n", E[i], f[i]);
    }
	fclose(fd);
	//fd = popen("gnuplot -p data.gp", "w");
	//pclose(fd);
	return 0;
}
