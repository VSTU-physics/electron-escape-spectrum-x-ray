#pragma once

#include <stdio.h>
#include <math.h>
#include "physics.h"
#include "calculations.h"

void test_parse();
void load_auger(int Z, const char* ch, auger_t &aug);
void load_subst(int Z, const char* ch, subst_t &subs);
void load_esharp(double *esharp, double *E, int N, const char* ch, subst_t s);
void load_ltr(double *ltr, double *E, int N, const char* ch, subst_t s);
void load_approx(int Z, const char* ch, approx_t &app);
void le_approx(double *ltr, double *eps, double *E, int N, approx_t ap, subst_t s);