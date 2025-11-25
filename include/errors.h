
#pragma once

#include "funzioni.h"

//blocking method
void blocking(int N_t, int num_oss);

//jackknife method
void jackknife(int N_t, int num_oss);


void err_max_blocking(int N_t, rowvec &sigma, int num_oss);
