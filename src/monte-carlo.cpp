#include <stdio.h>
#include <math.h>
#include <time.h>
#include "physics.h"
#include "calculations.h"
#include "parse.h"
#include <stdlib.h>

void monte_carlo(int Z, int nparticles, int ntimes, int N, double Emin, double Smax, double lmax)
{
    double *E, *beta, *alpha, beta_t, alpha_t, theta_s, phi_s, theta_t, phi_t;
    E = new double[N];
    beta = new double[N];
    alpha = new double[N];

    subst_t s;
    auger_t a;
    load_subst(Z, "data/subst.pl", s);
    load_auger(Z, "data/aug.pl", a);
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

    double rand, x, y, z, S, dS;
    int k;
    bool stop;
    dS = Smax / N;

    double *n_E, *n_S, *E_S;
    double E_exit;
    n_E = new double[N]; for (int i = 0; i<N; i++) n_E[i] = 0;
    n_S = new double[N]; for (int i = 0; i<N; i++) n_S[i] = 0;
    E_S = new double[N]; for (int i = 0; i<N; i++) E_S[i] = 0;
    int general_sum = 0;
    
    for (int i = 0; i<nparticles; i++)
    {
        stop = 1;
        S = 0;
        p.x = random(lmax/2) + lmax;
        p.y = random(lmax/2) + lmax;
        p.z = random(lmax);
        p.theta = random(M_PI);
        p.phi = random(2*M_PI);
        rand = random(1);
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

void gnuplot_mc(int &wxt)
{
    FILE *fd;
    fd = fopen("gnuplot_mc.gp", "w");
    wxt ++;
    fprintf(fd, "set terminal wxt %d\n", wxt);
    fprintf(fd, "unset key\n");
    fprintf(fd, "set title 'Спектр n(E)'\n");
    fprintf(fd, "plot 'results_mc.dat' using 1:2 lw 2 with lines\n");
    wxt ++;
    fprintf(fd, "set terminal wxt %d\n", wxt);
    fprintf(fd, "unset key\n");
    fprintf(fd, "set title 'Зависимость числа вышедших электронов от пробега n(S)'\n");
    fprintf(fd, "plot 'results_mc.dat' using 3:4 lw 2 with lines\n");
    wxt ++;
    fprintf(fd, "set terminal wxt %d\n", wxt);
    fprintf(fd, "unset key\n");
    fprintf(fd, "set title 'Зависимость средней энергии вышедших электронов от пробега E(S)'\n");
    fprintf(fd, "plot 'results_mc.dat' using 3:5 lw 2 with lines\n");
    fclose(fd);
}

