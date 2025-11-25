
#pragma once

#include <iostream>   // Defines 'std'
#include <armadillo>  // Defines 'arma'

using namespace std;
using namespace arma;

extern int N;//number of spins
extern double magn;//value of magnation
extern int dimensione;//cluster dimension for wolff
extern int pv;//number of first nearest neighbours


//same as pow but faster
double pow1(double base, int esp);

//calculates the module of the relative distance
double mod(rowvec &r); 

//calculates the module of the relative distance
double mod(cube &r, double riga, double colonna);

//calculates the actual potential 
double potenziale_ising(mat &r, double J, int i, mat &primivicini);

//calculates the modified potential 
double modifica_potenziale(mat &r, double J, int i, mat &primivicini);

// metropolis algorithm 
void MRT2(mat &r, double *V, double T, double J, mat &primivicini, int i = -999);

//selects the spin for the cluster update
void spin_inverter(mat &r, double sp_c, mat &p_v, double prob, rowvec &spin_m);

// wolff algorithm
void Wolff(mat &r, double *V, double T, double J, mat &primivicini);
