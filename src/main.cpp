#include <stdio.h>
#include <math.h>
#ifdef __WIN__
    #include <windows.h>
#endif // __WIN__
#include <stdlib.h>
#include "physics.h"
#include "calculations.h"
#include "parse.h"
#include "monte-carlo.h"

void check_data(int Z, subst_t s, auger_t a, approx_t ap)
{
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
    printf("Для метода аппроксимации:\n");
    printf("Кол-во типов: %d\n", ap.N);
    printf("Энергии: [");
    for (int i = 0; i < ap.N; i++) {
        printf(" %f", ap.E[i]);
    }
    printf(" ]\n");
    printf("alpha: [");
    for (int i = 0; i < ap.N; i++) {
        printf(" %f", ap.alpha_i[i]);
    }
    printf(" ]\n");
    printf("d1: [");
    for (int i = 0; i < ap.N; i++) {
        printf(" %f", ap.d1_i[i]);
    }
    printf(" ]\n");
    printf("R0: [");
    for (int i = 0; i < ap.N; i++) {
        printf(" %f", ap.R0_i[i]);
    }
    printf(" ]\n");
}

double delta(double x, double dx)
{
    return (fabs(x) < fabs(dx/2)) ? 1./dx : 0.;
}

void solve(const char* fname, auger_t a, int N, double* z,
           int M, double* E, double* ltrs, double* epss,
           double* I1s, double* I2s, double* fs)
{
    FILE* fd = fopen(fname, "w");
    double* tmp;
    double* u = new double[N];
    double* up = new double[N];

    for (int i = 1; i < M; i++)
    {
        tmp = u;
        u = up;
        up = tmp;

        double source = 0;
        double dE = E[i] - E[i - 1];
        for (int k = 0; k < a.N; k++)
            source += a.P[k] * delta(E[i] - a.E[k], dE);

        spe(u, up, z, N, dE,
            1./3 * ltrs[i] / epss[i],
            source,
            I1s[i]/(2 - 3 * I2s[i]),
            - ltrs[i] / 3,
            0,
            1.,
            0.,
            0.
            );
        fs[i] = 3 / ltrs[i] * I1s[i]/(2 - 3 * I2s[i]) * u[0];
        if (i % 100 == 0)
        {
            for (int j = 0; j < N; j += 10)
            {
                fprintf(fd, "%e %e %e\n", E[i], z[j], u[j]);
            }
            fprintf(fd, "\n");
        }
    }
    fclose(fd);
}

void analytical(int Z, int M, int N, double l)
{
    subst_t s;
    auger_t a;
    load_subst(Z, "data/subst.pl", s);
    load_auger(Z, "data/aug.pl", a);
    double Emin = 1000, Emax = 1.1 * a.E[0];

    double* z = new double[N];
    for (int i = 0; i<N; i++)
        z[i] = i * l / (N - 1);

    double* E = new double[M];
    for (int i = 0; i < M; i++)
        E[i] = Emax + (Emin - Emax) / (M - 1) * i;

    double* I1s = new double[M];
    double* I2s = new double[M];
    double* ltrs = new double[M];
    double* epss = new double[M];
    double* rs = new double[M];
    double* bs = new double[M];
    double* fs = new double[M];

    for (int i = 0; i<M; i++)
    {
        I1s[i] = I1(s, E[i]);
        I2s[i] = I2(s, E[i]);
        ltrs[i] = l_tr(s, E[i]);
        epss[i] = eps(s, E[i]);
        bs[i] = I1s[i] / (2 - 3 * I2s[i]) / ltrs[i] * 3;
    }

    for (int i = 0; i < M-1; i++)
    {
        rs[M - 2 - i] = rs[M - 1 - i] +
            0.5 * (1 / epss[M - 1 - i] + 1 / epss[M - i - 2]) *
            (E[M - i - 2] - E[M - i - 1]);
    }

    solve("solution_a.dat", a, N, z, M, E, ltrs, epss, I1s, I2s, fs);
    FILE* fd = fopen("data_a.dat", "w");
    fprintf(fd, "# E ltr eps f r bc\n");
    for (int i = 0; i < M; i++)
    {
        fprintf(fd, "%e %e %e %e %e %e\n",
                E[i], ltrs[i], epss[i], fs[i], rs[i], bs[i]);
    }
    fclose(fd);

    delete [] I1s;
    delete [] I2s;
    delete [] ltrs;
    delete [] epss;
    delete [] rs;
    delete [] bs;
    delete [] fs;
    delete [] E;
}

void table(int Z, int M, int N, double l)
{
    subst_t s;
    auger_t a;
    load_subst(Z, "data/subst.pl", s);
    load_auger(Z, "data/aug.pl", a);
    double Emin = 1000, Emax = 1.1 * a.E[0];
    double* z = new double[N];
    for (int i = 0; i<N; i++)
        z[i] = i * l / (N - 1);

    double* E = new double[M];
    for (int i = 0; i < M; i++)
        E[i] = Emax + (Emin - Emax) / (M - 1) * i;

    double* I1s = new double[M];
    double* I2s = new double[M];
    double* ltrs = new double[M];
    double* epss = new double[M];
    double* rs = new double[M];
    double* bs = new double[M];
    double* fs = new double[M];
    load_ltr(ltrs, E, M, "data/", s);
    load_esharp(epss, E, M, "data/", s);
    for (int i = 0; i<M; i++)
    {
        I1s[i] = I1(s, E[i]);
        I2s[i] = I2(s, E[i]);
        bs[i] = I1s[i] / (2 - 3 * I2s[i]) / ltrs[i] * 3;
    }
    for (int i = 0; i < M-1; i++)
    {
        rs[M - 2 - i] = rs[M - 1 - i] +
            0.5 * (1 / epss[M - 1 - i] + 1 / epss[M - i - 2]) *
            (E[M - i - 2] - E[M - i - 1]);
    }
    solve("solution_t.dat", a, N, z, M, E, ltrs, epss, I1s, I2s, fs);
    FILE* fd = fopen("data_t.dat", "w");
    fprintf(fd, "# E ltr eps f r bc\n");
    for (int i = 0; i < M; i++)
    {
        fprintf(fd, "%e %e %e %e %e %e\n",
                E[i], ltrs[i], epss[i], fs[i], rs[i], bs[i]);
    }
    fclose(fd);
    delete [] I1s;
    delete [] I2s;
    delete [] ltrs;
    delete [] epss;
    delete [] rs;
    delete [] bs;
    delete [] fs;
    delete [] E;
}

void approximation(int Z, int M, int N, double l)
{
    subst_t s;
    auger_t a;
    approx_t ap;
    load_subst(Z, "data/subst.pl", s);
    load_auger(Z, "data/aug.pl", a);
    load_approx(Z, "data/approx.pl", ap);
    double Emin = 1000, Emax = 1.1 * a.E[0];
    double* z = new double[N];
    for (int i = 0; i<N; i++)
        z[i] = i * l / (N - 1);

    double* E = new double[M];
    for (int i = 0; i < M; i++)
        E[i] = Emax + (Emin - Emax) / (M - 1) * i;

    double* I1s = new double[M];
    double* I2s = new double[M];
    double* ltrs = new double[M];
    double* epss = new double[M];
    double* rs = new double[M];
    double* bs = new double[M];
    double* fs = new double[M];

    le_approx(ltrs, epss, E, M, ap, s);
    for (int i = 0; i<M; i++)
    {
        I1s[i] = I1(s, E[i]);
        I2s[i] = I2(s, E[i]);
        bs[i] = I1s[i] / (2 - 3 * I2s[i]) / ltrs[i] * 3;
    }

    for (int i = 0; i < M-1; i++)
    {
        rs[M - 2 - i] = rs[M - 1 - i] +
            0.5 * (1 / epss[M - 1 - i] + 1 / epss[M - i - 2]) *
            (E[M - i - 2] - E[M - i - 1]);
    }
    solve("solution_p.dat", a, N, z, M, E, ltrs, epss, I1s, I2s, fs);
    FILE* fd = fopen("data_p.dat", "w");
    fprintf(fd, "# E ltr eps f r bc\n");
    for (int i = 0; i < M; i++)
    {
        fprintf(fd, "%e %e %e %e %e %e\n",
                E[i], ltrs[i], epss[i], fs[i], rs[i], bs[i]);
    }
    fclose(fd);
    delete [] I1s;
    delete [] I2s;
    delete [] ltrs;
    delete [] epss;
    delete [] rs;
    delete [] bs;
    delete [] fs;
    delete [] E;
}
void test_all()
{
    int Z = 32;
    int N = 100;
    int M = 500;
    double l = 1;

    analytical(Z, M, N, l);
    table(Z, M, N, l);
    approximation(Z, M, N, l);
}

void test_mc()
{
    int Z = 32;
    int nparticles = 1000000;
    int ntimes = 10000;
    int N = 100;
    double Emin = 10;
    double Smax = 0.00006;
    double lmax = 0.00001;
    monte_carlo(Z, nparticles, ntimes, N, Emin, Smax, lmax);
}

int main()
{
    #ifdef __WIN__
        system("chcp 65001");
    #endif // __WIN__
    //test_all();
    test_mc();
    //test_spline();
    return 0;
}
