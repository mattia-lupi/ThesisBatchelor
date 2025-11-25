
#include "funzioni.h"

//same as pow but faster
double pow1(double base, int esp) {
    double ris = 1.0;
    for (int i = 0; i < esp; i++) {
        ris *= base;
    }
    return ris;
}

//calculates the module of the relative distance
double mod(rowvec &r){ 
    double mod = sqrt(pow1(r(0),2)+pow1(r(1),2));
    return mod;
}

//calculates the module of the relative distance
double mod(cube &r, double riga, double colonna){
    double mod = sqrt(pow1(r(riga,colonna,0),2)+pow1(r(riga,colonna,1),2));
    return mod;
}

//calculates the actual potential 
double potenziale_ising(mat &r, double J, int i, mat &primivicini){
    double V = 0;

    rowvec p_v = primivicini.row(i);

    for (int j = 0; j < pv; ++j){
        V += r((int)p_v(j),2) * r(i,2) * J;
    }
    return -V;
}

//calculates the modified potential 
double modifica_potenziale(mat &r, double J, int i, mat &primivicini){
    double V = 0;

    rowvec p_v = primivicini.row(i);

    for (int j = 0; j < pv; ++j){
        V += r((int)p_v(j),2);
    }
    V *= J * 2 * r(i,2);
    return V;
}

// metropolis algorithm 
void MRT2(mat &r, double *V, double T, double J, mat &primivicini, int i){

    int m;
    //find the spin for MRT2
    if (i == -999){
       // random choice for the first spin
       double s = rand() / (RAND_MAX + 1.0);
       m = (int)rint((N-1) * s);
    }else{
       //sequential metropolis algorithm
       m = i%N;
    }

    double V_mod = modifica_potenziale(r,J,m,primivicini);
    //changes the contribution to the potential for the modified spin
    
    if(V_mod>0){
        if(exp(-V_mod/T)>(rand()/((double)RAND_MAX+1.0))){
            r(m,2) *= -1;
            magn +=  2 * r(m, 2) / N;
            *V += V_mod;
        }
    }
    else{//if the change in energy is less than zero then always accept
        r(m,2) *= -1;
        magn +=  2 * r(m, 2) / N;
        *V += V_mod;
    }
    return; 
}

//selects the spin for the cluster update
void spin_inverter(mat &r, double sp_c, mat &p_v, double prob, rowvec &spin_m){
    spin_m(sp_c) = 1;

    //cycle over the first nearest neighburs
    for (int i = 0; i < pv; ++i){
        //spin to add to the cluster?
        int n = p_v(sp_c, i);
        if (r(sp_c,2) == r(n, 2) && spin_m(n) == -1){
            if(prob > (rand()/(RAND_MAX + 1.))){
                spin_inverter(r, n, p_v, prob, spin_m);
            }
        }
    }
}

// wolff algorithm
void Wolff(mat &r, double *V, double T, double J, mat &primivicini){
    //find the first spin for the cluster
    double s= rand() / (RAND_MAX + 1.0);
    int m = (int)rint((N-1) * s);
    
    rowvec spin_da_mod(N, fill::value(-1));

    double prob = 1 - exp(- 2 * J / T);
    dimensione = 0;
    
    //selects the cluster
    spin_inverter(r, m, primivicini, prob, spin_da_mod);

    for (int i = 0; i < N; ++i){
        //invert the cluster
        if(spin_da_mod(i) == 1){
            dimensione++;
            r(i, 2) *= -1;
            magn += 2 * r(i,2) / N;
        }
    }

    *V = 0;
    //recalculates the potential
    for (int i = 0; i < N; ++i){
        *V += potenziale_ising(r, J, i, primivicini);
    }
    *V/=2;
}
