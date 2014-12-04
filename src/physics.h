#pragma once
#define Ee 511003.4 // эВ
#define re 2.8179403267e-13 // cм
#define Na 6.0221413e23 // моль^-1
#define alpha 0.00729735256

typedef struct subst
{
    int Z;       // заряд
    double rho;  // плотность
    double M;    // молярная масса
    double Emax; // макс. энергия Оже-электронов
    double U0;   // работа выхода
} subst_t; // для описания свойств материала

typedef struct auger
{
	int N;          // число Оже электронов
	double *E;      // энергии
	double *P;      // вероятности
	char atom[30];  // строка для обозначения атома
	char shell[30]; // оболочки с которых летят Оже-электроны
} auger_t; // Оже?!

typedef struct approx
{
    int N;
    double *E;
    double *alpha_i;
    double *d1_i;
    double *R0_i;
} approx_t;

double l_tr(subst_t s, double E);
double eps(subst_t s, double E);
double I1(subst_t s, double E);
double I2(subst_t s, double E);
double J(int Z);
double crsec(double theta, subst_t s, double E);
double eta(subst_t s, double E);

