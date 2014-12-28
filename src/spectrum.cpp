#include <cstdio>
#include <cmath>
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

void solve(const char* fname, auger_t a, int N, double* z,
           int M, double* E, double* ltrs, double* epss,
           double* I1s, double* I2s, double* fs)
{
    FILE* fd = fopen(fname, "w");
    double* tmp;
    double* u = new double[N];
    double* up = new double[N];
    for (int i = 0; i < N; i++) u[i] = 0;

    for (int i = 1; i < M; i++)
    {
        tmp = u;
        u = up;
        up = tmp;

        double source = 0;
        double dE = E[i] - E[i - 1];

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
        for (int k = 0; k < a.N; k++)
        {
            if (fabs(E[i] - a.E[k]) < fabs(dE/2))
            {
                for (int j = 0; j < N; j++) u[j] += a.P[k];
            }
        }
        fs[i] = 3 * I1s[i] / (2 - 3 * I2s[i]) * u[0]; // *
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

void analytical(int Z, int M, int N, double l, double Emin)
{
    subst_t s;
    auger_t a;
    load_subst(Z, "data/subst.dat", s);
    load_auger(Z, "data/aug.dat", a);
    double Emax = 1.1 * a.E[0];

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

void table(int Z, int M, int N, double l, double Emin)
{
    subst_t s;
    auger_t a;
    load_subst(Z, "data/subst.dat", s);
    load_auger(Z, "data/aug.dat", a);
    double Emax = 1.1 * a.E[0];
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

void approximation(int Z, int M, int N, double l, double Emin)
{
    subst_t s;
    auger_t a;
    approx_t ap;
    load_subst(Z, "data/subst.dat", s);
    load_auger(Z, "data/aug.dat", a);
    load_approx(Z, "data/approx.dat", ap);
    double Emax = 1.1 * a.E[0];
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

void monte_carlo(int Z, int nparticles, int ntimes, int N, double Emin, double Smax, double lmax)
{
    double *E, *beta, *alpha, beta_t, alpha_t, theta_s, phi_s, theta_t, phi_t;
    E = new double[N];
    beta = new double[N];
    alpha = new double[N];

    subst_t s;
    auger_t a;
    load_subst(Z, "data/subst.dat", s);
    load_auger(Z, "data/aug.dat", a);
    double Emax = 1.1 * a.E[0];
    double dE = (Emax - Emin) / (N - 1);
    for (int i = 0; i < N; i++)
        E[i] = Emin + dE * i;

    double inv_lambda_el, *inv_lambda_in;
    inv_lambda_in = new double[N];

    load_mc_elastic(alpha, E, N, inv_lambda_el, "data/", s);
    load_mc_inelastic(beta, E, N, inv_lambda_in, "data/", s);

    double lambda;

    particle_t p;

    double x, y, z, dS;
    int k;
    dS = Smax / N;

    double *n_E, *n_S, *E_S;
    double E_exit;
    n_E = new double[N]; for (int i = 0; i<N; i++) n_E[i] = 0;
    n_S = new double[N]; for (int i = 0; i<N; i++) n_S[i] = 0;
    E_S = new double[N]; for (int i = 0; i<N; i++) E_S[i] = 0;
    int general_sum = 0;

    for (int i = 0; i<nparticles; i++)
    {
        bool stop = true;
        double S = 0;
        p.x = random(lmax/2) + lmax;
        p.y = random(lmax/2) + lmax;
        p.z = random(lmax);
        p.theta = random(M_PI);
        p.phi = random(2*M_PI);
        double rand = random(1);
        p.E = auger_source(a, rand);
        for (int j = 0; j<ntimes && stop; j++)
        {
            k = (int)((p.E - Emin) / dE);
            lambda = 1 / (inv_lambda_el + inv_lambda_in[k]);

            rand = random(1);
            x = p.x - lambda * log(rand) * sin(p.theta) * cos(p.phi);
            y = p.y - lambda * log(rand) * sin(p.theta) * sin(p.phi);
            z = p.z - lambda * log(rand) * cos(p.theta);
            S -= lambda * log(rand);
            if (z<0)
            {
                if (reflection(p.E, p.theta, s.U0))
                {
                    p.x = p.x + ( lambda * log(rand) + 2 * p.z / cos(p.theta)) * sin(p.theta) * cos(p.phi);
                    p.y = p.y + ( lambda * log(rand) + 2 * p.z / cos(p.theta)) * sin(p.theta) * sin(p.phi);
                    p.z = p.z + ( lambda * log(rand) + 2 * p.z / cos(p.theta)) * cos(p.theta);
                    p.theta = M_PI - p.theta;
                    p.phi = M_PI + p.phi;
                } else {
                    E_exit = p.E - s.U0;
                    //theta_exit = sqrt(1 - pow(sin(p.theta), 2) / (1 - s.U0 / p.E));
                    k = (int)((E_exit - Emin) / dE);
                    if (k < N) n_E[k] += 1;
                    k = (int) (S/dS);
                    //printf("%d %d %e %e\n", general_sum, k, p.E, E_S[k]);
                    if (k < N)
                    {
                        n_S[k] += 1;
                        E_S[k] += E_exit;
                    }
                    stop = 0;
                    general_sum ++;
                    //
                }
            } else {
                p.x = x;
                p.y = y;
                p.z = z;
            }

            k = (int)((p.E - Emin) / dE);
            if (random(1) < inv_lambda_in[k] * lambda)
            {
                beta_t = linterp(p.E, E[k], beta[k], E[k+1], beta[k+1]);
                p.E = p.E*(1 - 1 / beta_t * tan(M_PI * random(1) / 2));
                if (p.E < Emin)
                {
                    stop = 0;
                }
            } else {
                alpha_t = linterp(p.E, E[k], alpha[k], E[k+1], alpha[k+1]);
                theta_s = 1 / alpha_t * tan(M_PI * random(1) / 2);
                phi_s = 2*M_PI*random(1);
                phi_t = atan2(sin(p.theta)*sin(p.phi)*cos(theta_s) -
                              cos(p.phi)*sin(theta_s)*cos(phi_s) +
                              cos(p.theta)*sin(p.phi)*sin(theta_s)*sin(phi_s),
                              sin(p.theta)*cos(p.phi)*cos(theta_s) +
                              sin(p.phi)*sin(theta_s)*cos(phi_s) +
                              cos(p.theta)*cos(p.phi)*sin(theta_s)*sin(phi_s));
                theta_t = acos(cos(theta_s)*cos(p.theta) - sin(theta_s)*sin(p.theta)*sin(phi_s));
                p.theta = theta_t;
                p.phi = phi_t;
            }
        }
    }

    FILE *fd;
    fd = fopen("results_mc.dat", "w");
    for (int i = 0; i<N; i++)
    {
        fprintf(fd, "%e %e %e %e %e\n",
                E[i],
                n_E[i],
                (i+1)*dS,
                n_S[i],
                ((int) n_S[i] == 0) ? 0 : E_S[i] / n_S[i]);
    }
    fclose(fd);

    delete [] E;
    delete [] beta;
    delete [] alpha;
    delete [] n_E;
    delete [] n_S;
    delete [] E_S;
}
