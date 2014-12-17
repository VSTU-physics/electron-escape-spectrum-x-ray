#pragma once
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "physics.h"
#include "calculations.h"
#include "parse.h"
#include <stdlib.h>

void monte_carlo(int Z, int nparticles, int ntimes, int N, double Emin = 5, double Smax = 1, double lmax = 1);
void gnuplot_mc(int &wxt);

