#include <stdio.h>
#include <math.h>
#ifdef __WIN__
    #include <windows.h>
#endif // __WIN__

#include "physics.h"
#include "calculations.h"
#include "parse.h"

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
    return (fabs(x)<dx/2) ? 1./dx : 0.;
}

void test_all()
{
    int Z = 32;
    double *u, *z, *up;
    double l = 1;
    int N = 500; // число точек
    subst_t s;
    auger_t a;
    approx_t ap;

    u = new double[N];
    up = new double[N];
    z = new double[N];
    for (int i = 0; i<N; i++) z[i] = i*l/(N-1);

    double *E, *tmp, *f, *f2, *f3, *ltrs, *epss, *I1s, *I2s, *ltrs2, *epss2, *ltrs3, *epss3,
           *Rs, *Rs2, *Rs3, *bs, *bs2, *bs3;
    int M =  5000; // число точек в спектре
    double Emax = 10000, Emin;
    ltrs = new double[M];
    epss = new double[M];
    ltrs2 = new double[M];
    epss2 = new double[M];
    ltrs3 = new double[M];
    epss3 = new double[M];
    I1s = new double[M];
    I2s = new double[M];
    E = new double[M];
    f = new double[M];
    f2 = new double[M];
    f3 = new double[M];
    Rs = new double[M];
    Rs2 = new double[M];
    Rs3 = new double[M];
    bs = new double[M];
    bs2 = new double[M];
    bs3 = new double[M];

    load_subst(Z, "data/subst.pl", s);
    load_auger(Z, "data/aug.pl", a);
    load_approx(Z, "data/approx.pl", ap);
    check_data(Z, s, a, ap);

    for (int i = 0; i<N; i++)
    {
        u[i] = 0.;
        up[i] = 0.;
    }

    Emin = 1000;
    double dE = (Emax - Emin) / (M - 1);
    for (int i = 0; i<M; i++)
    {
        E[i] = Emin + i * dE;
    }
    for (int i = 0; i<M; i++)
    {
        ltrs[i] = l_tr(s, E[i]);
        epss[i] = eps(s, E[i]);
        I1s[i] = I1(s, E[i]);
        I2s[i] = I2(s, E[i]);
    }

    FILE *fd;
    fd = fopen("solution_a.dat", "w");
    for (int i = M - 1; i>=0; i--)
    {
        tmp = u;
        u = up;
        up = tmp;

        double source = 0;
        for (int k = 0; k<a.N; k++)
        {
            source += a.P[k]*delta(E[i] - a.E[k], dE);
        }
        spe(u, up, z, N, dE,
            - 1./3 * ltrs[i] / epss[i],
            source,
            I1s[i]/(2 - 3 * I2s[i]),
            - ltrs[i] / 3,
            0,
            1.,
            0.,
            0.
            );
        f[i] = 3 / ltrs[i] * I1s[i]/(2 - 3 * I2s[i]) * u[0];
        if (i != M - 1)
            Rs[i] = Rs[i+1] + 0.5 * (1 / epss[i+1] + 1 / epss[i]) * dE;
        bs[i] = I1s[i]/(2 - 3 * I2s[i]) / ltrs[i] * 3;

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

    for (int i = 0; i<N; i++)
    {
        u[i] = 0.;
        up[i] = 0.;
    }
    //Emin = 0.;
    dE = (Emax - Emin) / (M - 1);
    for (int i = 0; i<M; i++)
    {
        E[i] = Emin + i * dE;
    }
    load_ltr(ltrs2, E, M, "data/", s );
    load_esharp(epss2, E, M, "data/", s);

    fd = fopen("solution_t.dat", "w");
    for (int i = M - 1; i>=0; i--)
    {
        tmp = u;
        u = up;
        up = tmp;

        double source = 0;
        for (int k = 0; k<a.N; k++)
        {
            source += a.P[k]*delta(E[i] - a.E[k], dE);
        }
        spe(u, up, z, N, dE,
            - 1./3 * ltrs2[i] / epss2[i],
            source,
            I1s[i]/(2 - 3 * I2s[i]),
            - ltrs2[i] / 3,
            0,
            1.,
            0.,
            0.
            );
        f2[i] = 3 / ltrs2[i] * I1s[i]/(2 - 3 * I2s[i]) * u[0];
        if (i != M - 1)
            Rs2[i] = Rs2[i+1] + 0.5 * (1 / epss2[i+1] + 1 / epss2[i]) * dE;
        bs2[i] = I1s[i]/(2 - 3 * I2s[i]) / ltrs2[i] * 3;
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
    for (int i = 0; i<N; i++)
    {
        u[i] = 0.;
        up[i] = 0.;
    }
    //Emin = 0.;
    dE = (Emax - Emin) / (M - 1);
    for (int i = 0; i<M; i++)
    {
        E[i] = Emin + i * dE;
    }

    printf("error\n");
    le_approx(ltrs3, epss3, E, M, ap, s);
    printf("error\n");
    fd = fopen("solution_p.dat", "w");
    for (int i = M - 1; i>=0; i--)
    {
        tmp = u;
        u = up;
        up = tmp;

        double source = 0;
        for (int k = 0; k<a.N; k++)
        {
            source += a.P[k]*delta(E[i] - a.E[k], dE);
        }
        spe(u, up, z, N, dE,
            - 1./3 * ltrs3[i] / epss3[i],
            source,
            I1s[i]/(2 - 3 * I2s[i]),
            - ltrs3[i] / 3,
            0,
            1.,
            0.,
            0.
            );
        f3[i] = 3 / ltrs2[i] * I1s[i]/(2 - 3 * I2s[i]) * u[0];
        if (i != M - 1)
            Rs3[i] = Rs3[i+1] + 0.5 * (1 / epss3[i+1] + 1 / epss3[i]) * dE;
        bs3[i] = I1s[i]/(2 - 3 * I2s[i]) / ltrs3[i] * 3;
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

    fd = fopen("data.dat", "w");
    for (int i = M - 1; i>=0; i--)
    {
        fprintf(fd, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
                E[i], ltrs[i], ltrs2[i], ltrs3[i],
                epss[i], epss2[i], epss3[i],
                f[i], f2[i], f3[i],
                Rs[i], Rs2[i], Rs3[i],
                bs[i], bs2[i], bs3[i]);
    }
    fclose(fd);
    fd = fopen("data.gp", "w");
    fprintf(fd, "set terminal wxt 1\n");
    fprintf(fd, "set size square\n");
    fprintf(fd, "plot 'data.dat' using 1:2 with lines title 'l_tr_A(E)',\\\n");
    fprintf(fd, "'data.dat' using 1:3 with lines title 'l_tr_T(E)',\\\n");
    fprintf(fd, "'data.dat' using 1:4 with lines title 'l_tr_P(E)'\n");
    fprintf(fd, "set terminal wxt 2\n");
    fprintf(fd, "set size square\n");
    fprintf(fd, "plot 'data.dat' using 1:5 with lines title 'eps_A(E)',\\\n");
    fprintf(fd, "'data.dat' using 1:6 with lines title 'eps_T(E)',\\\n");
    fprintf(fd, "'data.dat' using 1:7 with lines title 'eps_P(E)'\n");
    fprintf(fd, "set terminal wxt 3\n");
    fprintf(fd, "set size square\n");
    fprintf(fd, "plot 'data.dat' using 1:8 with lines title 'A',\\\n");
    fprintf(fd, "'data.dat' using 1:9 with lines title 'T',\\\n");
    fprintf(fd, "'data.dat' using 1:10 with lines title 'P'\n");
    fprintf(fd, "set terminal wxt 4\n");
    fprintf(fd, "set size square\n");
    fprintf(fd, "plot 'data.dat' using 1:11 with lines title 'R(E) A',\\\n");
    fprintf(fd, "'data.dat' using 1:12 with lines title 'R(E) T',\\\n");
    fprintf(fd, "'data.dat' using 1:13 with lines title 'R(E) P' \n");
    fprintf(fd, "set terminal wxt 5\n");
    fprintf(fd, "set size square\n");
    fprintf(fd, "plot 'data.dat' using 1:14 with lines title 'bc(E) A',\\\n");
    fprintf(fd, "'data.dat' using 1:15 with lines title 'bc(E) T',\\\n");
    fprintf(fd, "'data.dat' using 1:16 with lines title 'bc(E) P'\n");
    fclose(fd);

    fd = popen("gnuplot -p data.gp", "w");
    pclose(fd);

    delete[] E;
    delete[] f;
    delete[] ltrs;
    delete[] epss;
    delete[] f2;
    delete[] ltrs2;
    delete[] epss2;
    delete[] I1s;
    delete[] I2s;
}

void old_test()
{
        char type = 'A'; // A - аналитический расчёт коэффициентов,
                     // T - табличные значения
    int Z = 32;
    double *u, *z, *up;
    double l = 1;
    int N = 500; // число точек
    subst_t s;
    auger_t a;

    u = new double[N];
    up = new double[N];
    z = new double[N];
    for (int i = 0; i<N; i++) z[i] = i*l/(N-1);

    double *E, *tmp, *f, *ltrs, *epss, *I1s, *I2s;
    int M = 500; // число точек в спектре
    ltrs = new double[M];
    epss = new double[M];
    I1s = new double[M];
    I2s = new double[M];
    E = new double[M];
    f = new double[M];

    load_subst(Z, "data/subst.pl", s);
    load_auger(Z, "data/aug.pl", a);
    //check_data(Z, s, a);

    double dE = a.E[0] / (M - 1);
    for (int i = 0; i<M; i++)
    {
        E[i] = i * dE;
        }
    switch (type)
    {
        case 'A':
            for (int i = 0; i<M; i++)
            {
                ltrs[i] = l_tr(s, E[i]);
                epss[i] = eps(s, E[i]);
                I1s[i] = I1(s, E[i]);
                I2s[i] = I2(s, E[i]);
            }
            break;
        case 'T':
            // спецификация на эти функции: возвращают значения в массивы в
            // точках E
            //load_ltr(s, "data/", E, M, ltrs);
            //load_ebar(s, "data/", E, M, ebar);
            break;
        default:
            break;
    }

    for (int i = M - 1; i>=0; i--)
    {
        tmp = u;
        u = up;
        up = tmp;

        double source = 0;
        for (int k = 0; k<a.N; k++)
        {
            source += a.P[k]*delta(E[i] - a.E[k], dE);
        }
        spe(u, up, z, N, dE,
            - 1./3 * ltrs[i] / epss[i],
            source,
            I1s[i]/(2 - 3 * I2s[i]),
            - ltrs[i] / 3,
            0,
            1.,
            0.,
            0.
            );
        f[i] = 3 / ltrs[i] * I1s[i]/(2 - 3 * I2s[i]) * u[0];
    }
    FILE *fd;
    fd = fopen("data.gp", "w");
    fprintf(fd, "set multiplot layout 2,2\n");
    fprintf(fd, "set size square\n");
    fprintf(fd, "plot '-' with lines title 'l_tr(E)' \n");
    for (int i = M - 1; i>=0; i--)
    {
        fprintf(fd, "%e %e\n", E[i], ltrs[i]);
    }
    fprintf(fd, "end\n");
    fprintf(fd, "plot '-' with lines title 'eps(E)'\n");
    for (int i = M - 1; i>=0; i--)
    {
        fprintf(fd, "%e %e\n", E[i], epss[i]);
    }
    fprintf(fd, "end\n");
    fprintf(fd, "plot '-' with lines title 'spectrum'\n");
    for (int i = M - 1; i>=0; i--)
    {
        fprintf(fd, "%e %e\n", E[i], f[i]);
    }
    fprintf(fd, "end\n");
    //fprintf(fd, "plot '-' with lines title 'spectrum'\n");
    //for (int i = 0; i<N; i++)
    //{
        //fprintf(fd, "%e %e\n", u[i], up[i]);
    //}
    //fprintf(fd, "end\n");
    fprintf(fd, "unset multiplot\n");
    fclose(fd);
    fd = popen("gnuplot -p data.gp", "w");
    pclose(fd);
    delete[] E;
    delete[] f;
    delete[] ltrs;
    delete[] epss;
    delete[] I1s;
    delete[] I2s;
}

int main()
{
    #ifdef __WIN__
        system("chcp 65001");
    #endif // __WIN__
    test_all();
    //test_spline();
    return 0;
}
