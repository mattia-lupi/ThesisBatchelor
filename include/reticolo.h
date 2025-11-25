

#pragma once

#include "funzioni.h"
#include <fstream>
#include <cmath>

#include <iomanip>
#include <algorithm>

//prints the lattice on a file
void stampa_reticolo(mat &r);

// creates the square lattice
void crea_reticolo_quadrato(mat &r, double L);    

//needed for the creation of the exagonal lattice
double aggiunta_y(int cont_y, double H);

//needed for the creation of the exagonal lattice
double aggiunta_x(int cont_x);

// creates the exagonal lattice
void crea_reticolo_esagonale(mat &r, int L);

// creates the first neighbours getting the lattice from a file computed offline
void crea_primi_vicini(double L, string p_v_path);
