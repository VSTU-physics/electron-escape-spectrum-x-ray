/*!
    \file
    Файл содержит функции для всех аналитических формул в задаче.
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "physics.h"
#include "calculations.h"

/*!
    \brief Возвращает средний ионизационный потенциал

    \param[in] Z Заряд ядра атома вещества

    \f[
        J(Z) =
        \begin{cases}
            13{,}6 Z, & Z < 10; \\
            9{,}76 + 58{,}8 Z^{-1,19}Z & Z \geq 10.
        \end{cases}
    \f]
*/
double J(int Z)
{
    if (Z < 10)
        return 13.6 * Z;

    return (9.76 + 58.8 * pow(Z, -1.19)) * Z;
}

/*!
    \brief Возвращает величину \f$ \bar{\varepsilon} \f$ в аналитическом методе

    \param[in] s Параметры вещества
    \param[in] E Энергия частицы

    \f[
        \bar{\varepsilon} =
        2 \pi \frac{\rho N_AZ}{M}
        \frac{(T+1)^2 r_e^2}{T(T+2)}
        \left(
            2 \ln
            \frac{1{,}166\,T}{J/mc^2}
        \right)
    \f]

    \f$ T \f$ - кинетическая энергия электрона, отнесённая к энергии покоя.
*/
double eps(subst_t s, double E)
{
    double T = (s.U0 + E) / Ee;
    return Ee * 4 * M_PI * s.rho / s.M * Na * s.Z * pow((T+1)*re, 2) /
        T / (T + 2) * log(1.166 * E / J(s.Z));
}

/*!
    \brief Возвращает поправку \f$ \eta \f$ в аналитическом методе

    \param[in] s Параметры вещества
    \param[in] E Энергия частицы

    \f[
        \eta = \frac{1}{4}
        \Big[
            \frac{\alpha Z^{1/3}}{0{,}885}
        \Big]^2
        \frac{1}{T(T+2)} \left(1{,}13 + 3{,}76 \alpha^2 Z^2 \frac{(T+1)^2}{T(T+2)} \right);
    \f]
*/
double eta(subst_t s, double E)
{
    double T = (s.U0 + E) / Ee;
    return .25 * pow(ALPHA * pow(s.Z, 1./3) / 0.885, 2) / T / (T + 2) *
        (1.13 + 3.76 * pow(s.Z * ALPHA * (T + 1), 2) / T / (T + 2));
}

/*!
    \brief Возвращает величину транспортной длины \f$ \lambda_{tr} \f$ в аналитическом методе:

    \param[in] s Параметры вещества
    \param[in] E Энергия частицы

    \f[
        \lambda_{tr}^{-1} =
        2 \pi \frac{\rho N_AZ(Z+1)}{M}
        \frac{(T+1)^2 r_e^2}{T^2(T+2)^2}
        \left(
            \ln \left( \frac{1}{\eta} + 1 \right) - \frac{1}{\eta + 1}
        \right)
    \f]

*/
double l_tr(subst_t s, double E)
{
    double T = (s.U0 + E) / Ee;
    double inv_l_tr = 2 * M_PI * s.rho / s.M * Na * s.Z * (s.Z + 1) *
        pow((T+1)*re/T/(T+2), 2) *
        (log1p(1 / eta(s, E)) - 1 / (eta(s, E) + 1));
    return 1 / inv_l_tr;
}

/*!
    \brief Возвращает интеграл \f$ I_1 \f$ в граничном условии:

    \param[in] s Параметры вещества
    \param[in] E Энергия частицы

    \f[
        I_1 =
        -\frac{4}{3\varepsilon^4}+\frac{2}{\varepsilon^2}-\frac{2\varepsilon^2}{3}+\frac{4(1-\varepsilon^2)^\frac{3}{2}}{3\varepsilon^4},
    \f]
    где \f$ \varepsilon = \sqrt{E/U_0} \f$.
*/
double I1(subst_t s, double E)
{
    double eps = sqrt(s.U0 / (s.U0 + E));
    double e = eps * eps;
    return (eps < 1) ? -4. / 3 / e / e * (1 - pow(1 - e, 1.5)) + 2 / e - 2 * e / 3 : 0.;
}

/*!
    \brief Возвращает интеграл \f$ I_2 \f$ в граничном условии:

    \param[in] s Параметры вещества
    \param[in] E Энергия частицы

    \f[
        I_2 = \frac{4}{3}\frac{(1-\varepsilon^2)^\frac{3}{2}}{\varepsilon^4}-
        \frac{4}{21}\frac{(1-\varepsilon^2)^\frac{7}{2}}{\varepsilon^4} -
        \frac{4}{15}\frac{(1-\varepsilon^2)^\frac{5}{2}}{\varepsilon^2} -
        \frac{8(1-\varepsilon^7)}{7\varepsilon^4} +
        \frac{8(1-\varepsilon^5)}{5\varepsilon^2},
    \f]
    где \f$ \varepsilon = \sqrt{E/U_0} \f$.
*/
double I2(subst_t s, double E)
{
    double eps = sqrt(s.U0 / (s.U0 + E));
    double e = eps * eps;
    double i = 1 - e;
    return (eps < 1) ? - 8. / 7 / e / e + 8. / 5 / e - 16. * e * eps / 35 - 4 * sqrt(i) / 105 +
                        8 * sqrt(i) / 7 / e / e  - 36 * sqrt(i) / 35 / e - 8. / 105 * e * sqrt(i)  : 0.;
}

/*!
    \brief Возвращает дифференциальное сечение рассеяния в аналитическом методе:

    \param[in] theta Угол рассеяния
    \param[in] s Параметры вещества
    \param[in] E Энергия частицы

    \f[
        \frac{d\sigma_{el}}{d\Omega} =
        \frac{(T+1)^2 r_e^2}{T^2(T+2)^2} \frac{Z(Z+1)}{(1 - \cos \theta + 2 \eta)^2}
    \f]

*/
double crsec(double theta, subst_t s, double E)
{
    double T = (s.U0 + E) / Ee;
    return pow((T + 1) * re / T / (T + 2) / (1 - cos(theta) + 2 * eta(s, E)), 2) * s.Z * (s.Z + 1);
}

/*!
    \brief Возвращает энергию рождённого в Оже-процессе электрона в методе Монте-Карло:

    \param[in] s Параметры вещества
    \param[in] rand Случайное число от 0 до 1

*/
double auger_source(subst_t s, double rand)
{
    double pmax = 0, p = 0;
    for (int i = 0; i<s.N; i++)
    {
        pmax += s.P[i];
    }
    for (int i = 0; i<s.N; i++)
    {
        if (p > rand)
        {
            return s.E[i-1];
        }
        p += s.P[i]/pmax;
    }
    return s.E[s.N-1];
}

/*!
    \brief Вычисляет произошло или не произошло отражение в методе Монте-Карло:

    \param[in] E Энергия частицы
    \param[in] theta Угол между направлением движения частицы и осью \f$z\f$
    \param[in] U0 Работа выхода электрона из вещества \f$ U_0 \f$

    Вероятность отражения даёт коэффициент:
    \f[
        r = \left(\frac{\sqrt{E} \cos \theta - \sqrt{E \cos^2 \theta - U_0}}
                  {\sqrt{E} \cos \theta + \sqrt{E \cos^2 \theta - U_0}}\right)^2
    \f]

*/
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

