#include <math.h>
#include <stdio.h>
#include "physics.h"

double J(int Z)
{
    if (Z < 10)
        return 13.6 * Z;

    return (9.76 + 58.8 * pow(Z, -1.19)) * Z;
}

double eps(subst_t s, double E)
{
    double T = (s.U0 + E) / Ee;
    return Ee*4 * M_PI * s.rho / s.M * Na * s.Z * pow((T+1)*re, 2) /
        T / (T + 2) * log(1.166 * E / J(s.Z));
}

double eta(subst_t s, double E)
{
    double T = (s.U0 + E) / Ee;
    return .25 * pow(alpha * pow(s.Z, 1./3) / 0.885, 2) / T / (T + 2) *
        (1.13 + 3.76 * pow(s.Z * alpha * (T + 1), 2) / T / (T + 2));
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
    return -4. / 3 / e / e * (1 - pow(1 - e, 1.5)) + 2 / e - 2 * e / 3;
}

double I2(subst_t s, double E)
{
    double eps = sqrt(s.U0 / (s.U0 + E));
    double e = eps * eps;
    double i = 1 - e;
    return 4. / e / e * pow(i, 1.5) * (1. / 3 - i / 15 - i * i / 21) -
        8. / e / e * (1 - pow(eps, 7)) / 7. + 8. / e * (1 - pow(eps, 5)) / 5.;
}

