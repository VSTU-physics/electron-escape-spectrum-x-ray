#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "physics.h"
#include "calculations.h"

double J(int Z)
{
    if (Z < 10)
        return 13.6 * Z;

    return (9.76 + 58.8 * pow(Z, -1.19)) * Z;
}

double eps(subst_t s, double E)
{
    double T = (s.U0 + E) / Ee;
    return Ee * 4 * M_PI * s.rho / s.M * Na * s.Z * pow((T+1)*re, 2) /
        T / (T + 2) * log(1.166 * E / J(s.Z));
}

double eta(subst_t s, double E)
{
    double T = (s.U0 + E) / Ee;
    return .25 * pow(ALPHA * pow(s.Z, 1./3) / 0.885, 2) / T / (T + 2) *
        (1.13 + 3.76 * pow(s.Z * ALPHA * (T + 1), 2) / T / (T + 2));
}

double l_tr(subst_t s, double E)
{
    double T = (s.U0 + E) / Ee;
    double inv_l_tr = 2 * M_PI * s.rho / s.M * Na * s.Z * (s.Z + 1) *
        pow((T+1)*re/T/(T+2), 2) *
        (log1p(1 / eta(s, E)) - 1 / (eta(s, E) + 1));
    return 1 / inv_l_tr;
}

double I1(subst_t s, double E)
{
    double eps = sqrt(s.U0 / (s.U0 + E));
    double e = eps * eps;
    return (eps < 1) ? -4. / 3 / e / e * (1 - pow(1 - e, 1.5)) + 2 / e - 2 * e / 3 : 0.;
}

double I2(subst_t s, double E)
{
    double eps = sqrt(s.U0 / (s.U0 + E));
    double e = eps * eps;
    double i = 1 - e;
    return (eps < 1) ? - 8. / 7 / e / e + 8. / 5 / e - 16. * e * eps / 35 - 4 * sqrt(i) / 105 + 
                        8 * sqrt(i) / 7 / e / e  - 36 * sqrt(i) / 35 / e - 8. / 105 * e * sqrt(i)  : 0.;
}

double crsec(double theta, subst_t s, double E)
{
    double T = (s.U0 + E) / Ee;
    return pow((T + 1) * re / T / (T + 2) / (1 - cos(theta) + 2 * eta(s, E)), 2) * s.Z * (s.Z + 1);
}

double auger_source(auger_t a, double rand)
{
    double pmax = 0, p = 0;
    for (int i = 0; i<a.N; i++)
    {
        pmax += a.P[i];
    }
    for (int i = 0; i<a.N; i++)
    {
        if (p > rand)
        {
            return a.E[i-1];
        }
        p += a.P[i]/pmax;
    }
    return a.E[a.N-1];
}

bool reflection(double E, double theta, double U0)
{
    double a = sqrt(E) * cos(theta);
    double b = a * a - U0;
    if (b <= 0)
    {
        return 1;
    }
    double r = pow((a - sqrt(b)) / (a + sqrt(b)), 2);
    if (r <= random(1))
    {
        return 1;
    } else {
        return 0;
    }
};

