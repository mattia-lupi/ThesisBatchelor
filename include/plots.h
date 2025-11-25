
#pragma once
#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;

//executes the plot of the lattice
void plot_reticolo();

//executes the plot of the blocking
void blocking_plot();

//executes the plot of the observables in the temperature
void plot_osservabili_T(int c_Tm);

//executes the plot of the configurations
void plot_config();
void configurations(mat &r, int iter, double L);

//configurations in temperature
void configurations(mat &r, double T, double L);
