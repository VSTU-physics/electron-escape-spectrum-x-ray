#include <stdio.h>
#include <math.h>
#ifdef __WIN__
    #include <windows.h>
#endif // __WIN32__

#include "physics.h"
#include "calculations.h"
#include "parse.h"

void check_data(int Z, subst_t s, auger_t a)
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
}

double delta(double x, double dx)
{
    return (fabs(x)<dx/2) ? 1./dx : 0.;
}

int main()
{
    #ifdef __WIN__
        system("chcp 65001");
    #endif // __WIN__
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
    check_data(Z, s, a);

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
    return 0;
}
