#pragma once

#include <stdio.h>
#include <math.h>
#include "physics.h"
#include "calculations.h"

void test_parse();
void load_auger(int Z, const char* ch, auger_t &aug);
void load_subst(int Z, const char* ch, subst_t &subs);