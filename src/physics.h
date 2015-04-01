#pragma once
#define Ee 511003.4 // эВ
#define re 2.8179403267e-13 // cм
#define Na 6.0221413e23 // моль^-1
#define ALPHA 0.00729735256

/*!
    Содержит данные о веществе.
*/
typedef struct subst
{
    int Z;          ///< Заряд ядра атома вещества.
    double rho;     ///< Плотность вещества \f$ [\text{г/см}^3] \f$.
    double M;       ///< Молярная масса вещества \f$ [\text{г/моль}] \f$.
    double U0;      ///< Работа выхода электрона из вещества \f$ [\text{эВ}] \f$.
	int N;          ///< Число Оже электронов.
	double *E;      ///< Энергии различных электронов \f$ E_i \f$ \f$ [\text{эВ}] \f$.
	double *P;      ///< Вероятности появления электронов с энергиями \f$ E_i \f$.
	char atom[30];  ///< Строка для обозначения атома.
	char shell[30]; ///< Строка для обозначения оболочек, с которых летят Оже-электроны.
} subst_t;

/*!
    Данные по сечениям упругих и неупругих столкновений в методе аппроксимации:

    \f[
        \bar{\varepsilon}(E') = \frac{E_i}{\alpha_i R_{0i}} \left( \frac{E'}{E_i} \right)^{1 - \alpha_i}
    \f]
    \f[
        \lambda_{tr}(E') = \frac{R_{0i}}{d_{1i}} \left( \frac{E'}{E_i} \right)^{\alpha_i}
    \f]
*/
typedef struct approx
{
    int N;          ///< Число промежутков, на которых проводится аппроксимация
    double *E;      ///< Энергии, означающие границы промежутков \f$ [\text{эВ}] \f$
    double *alpha_i;///< Аппроксимирующий коэффициент \f$ \alpha_i \f$
    double *d1_i;   ///< Аппроксимирующий коэффициент \f$ d_{1i} \f$
    double *R0_i;   ///< Аппроксимирующий коэффициент \f$ R_{0i} \f$
} approx_t;

/*!
    Информация о частице.
*/
typedef struct particle
{
    double x;       ///< Координата \f$ x \f$ \f$ [\text{см}] \f$
    double y;       ///< Координата \f$ y \f$ \f$ [\text{см}] \f$
    double z;       ///< Координата \f$ z \f$ \f$ [\text{см}] \f$
    double theta;   ///< Угол между направлением движения частицы и осью \f$ z \f$
    double phi;     ///< Угол между проекцией направления движения частицы на плоскость \f$ xy \f$ и  осью \f$ x \f$
    double E;       ///< Энергия частицы
} particle_t;

double l_tr(subst_t s, double E);
double eps(subst_t s, double E);
double I1(subst_t s, double E);
double I2(subst_t s, double E);
double J(int Z);
double crsec(double theta, subst_t s, double E);
double eta(subst_t s, double E);
double auger_source(subst_t s, double rand);
bool reflection(double E, double theta, double U0);
