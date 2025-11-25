// To execute the code you need to select the number of first nearest neighbours 
// of your lattice. This code supports 3 and 4 meaning the hexagonal
// and the square lattice.
// After selecting the number of first nearest neighbours you need to decide which 
// kind of algorithm you want the code to use (prog_t): 
// random spin choice metropolis (0), sequential spin choice metropolis (1)
// and wolff algorithm (2)
// Once the algorithm has been selected you need to choose what you are interested
// behaviour of the observables as funcion of temperature (for different L) 
// c_Tm =-1/c_Lm=-1; behaviour of the observables as funcion of T 
// for specific L: c_Tm=-1, c_Lm=0/1/2/3 or the value of the observables 
// after each iteration of the algorithm c_Tm!=-1/c_Lm!=-1.


// Check on plot_osservabili_T for the zoom on the critical temperature
// remember T_c=1.519 for pv=3 and T_c=2.269 for pv=4.
// obs does all computing for you and creates the pictures of
// the spin configurations in the folder out/config if you chose c_Tm!=-1 or 
// (c_Tm!=-1 and c_Lm!=-1)
// To see the results of your computation use plot_osservabili_T

// The following need to be used only after obs has been used
// at least once because use its results
// To see the results of the blocking and of the jackknife method use
// blocking_plot (plots only magnation as it is for example)
// To see the lattice you have been computing with, use plot_reticolo (you
// might need to zoom out if you don't see the ends of the lattice) 




#define ARMA_NO_DEBUG//no control and more speed

#include "funzioni.h"
#include "errors.h"
#include "plots.h"  
#include "reticolo.h"

/*** variabili globali ***/

int N; //number of spins
int pv = 4;//number of first nearest neighbours
double magn = 0;
int dimensione = 0;

string p_v_path = "out/primivicini.txt";
string obs_path = "out/osservabili.txt";
string dati_path = "out/dati.txt";
string risultati_path = "out/risultati.txt";

ofstream dati, osservabili, risultati,results_scaled;

void obs(rowvec L, double *T_c_vec, double J, int *data);
void obs_T(mat &r, double J, double T_c, double L, int *data);


int main() {
   
    //default seed = 1
    srand(1);

    //number of iterations of the simulation
    int N_t;

    //to be defined later in code
    double passi_eq;
    double J = 1;

    // 0 random metropolis, 1 sequential metropolis , 2 wolff
    int prog_t = 2; 

    //set -1 to have O(T), not -1 shows temperature T=T_c+2
    int c_Tm = -1;

    //set -1 to have O(T,L)
    int c_Lm = 0;

    //critical temperature value for square and exagonal lattice
    double T_c[] = {2.269, 1.519};

    //number of spins per dimension of the lattice, in total the spins are N=LxL
    //128 takes quite long, closest to infinite system
    //for exagonal use only multiples of 8
    rowvec L = {8, 16, 32, 64};

    int initial_data[] = {c_Lm, c_Tm, N_t, prog_t};

    obs(L, T_c, J, initial_data);

    //in metropolis the dimension of clusters is obviously zero
    plot_osservabili_T(c_Tm);

    // plots bloking results
    // blocking_plot();

    //execute after obs to plot the lattice
    // plot_reticolo();

    return 0;
}


void obs(rowvec L, double *T_c_vec, double J, int *data){
    int c_Lm = data[0];
    int c_Tm = data[1];
    risultati.open(risultati_path);
    results_scaled.open("out/results_scaled.txt");
    int caso_max;
    if (c_Lm == -1) {
        caso_max = L.size();
    } else {
        caso_max = c_Lm + 1;
    }
    int start;
    if(c_Lm == -1){
        start = c_Lm + 1;
    }
    else{
        start =  c_Lm;
    }
    double T_c;
    for (int caso = start; caso < caso_max; ++caso){
    //execute to have the changes in function of L
        N = pow1(L(caso), 2);
        mat r(N,3);
        if(pv==4){
            crea_reticolo_quadrato(r, L(caso));
            //creates the lattice with all spin "1"
            T_c = T_c_vec[0];
        }
        else if(pv==3){
            crea_reticolo_esagonale(r, L(caso));
            //creates the lattice with all spin "1"
            T_c = T_c_vec[1];
        }
        else{
            cout<<"Lattice not avalilable,";
            cout<<" the square lattice is going to be presented"<<endl;
            pv=4;
            crea_reticolo_quadrato(r, L(caso));
            //creates the lattice with all spin "1"
            T_c = T_c_vec[0];
        }
        crea_primi_vicini(L(caso), p_v_path);
        //creates a file with the first nearest neighbours for each spin
        if (c_Lm == -1){
            cout<<"Executing lattice of N = "<<L(caso)<<"x"<<L(caso)<<endl;
        }
        if (caso != 0 && c_Lm == -1){
            risultati << "\n\n";
            results_scaled << "\n\n";
        }
        risultati << "L=" << L(caso) << endl;
        results_scaled << "L=" << L(caso) << endl;

        // makes the simulation at prescribed temperature for the chosen lattice size
        obs_T(r, J, T_c, L(caso), data);

        if(c_Tm != -1){
        //if I want how a single observable progresses over the iterations
        // I execute it only once
            break;
        }
    }
    risultati.close();
    results_scaled.close();
}
void obs_T(mat &r, double J, double T_c, double L, int *data){
    int c_Lm = data[0];
    int c_Tm = data[1];
    int N_t = data[2];
    int prog_t = data[3];

    double T_scaled;

    mat prim_vic(N,pv);

    ifstream pr_v;
    pr_v.open(p_v_path);
    for (int i = 0; i < N; ++i){//import matrix of the first nearest neighbours
        for (int j = 0; j < pv; ++j){
            pr_v >> prim_vic(i,j);
        }
    }
    pr_v.close();

    osservabili.open(obs_path);
    dati.open(dati_path);

    for(double T = T_c-1.2; T<T_c+6; T+=0.5){
        if(T>T_c-2 && T<T_c+2){
            T-=0.25;
            if (T>T_c-0.6 && T<T_c+0.6){
                T-=0.2;
            }
        }
        if (c_Tm!=-1){
            T = T_c + 0.5;
        }
        cout<<"Computing T = "<<T<<endl;
        int passi_eq;

        if(prog_t==2){
            T_scaled = (T - T_c) / T_c;
            passi_eq = 5 * N;//sufficient to get out of transition phase
            N_t = 2e4 + passi_eq + T * 300;
            if(L==32){
                passi_eq += 9*N; 
                N_t = 5e4 + passi_eq + T * 1000;
            }
            else if(L==64){
                passi_eq += 12*N; 
                N_t = 7e4 + passi_eq + T * 2000;
            }
        }
        else{
            T_scaled = (T - T_c) / T_c;
            passi_eq = 70 * N * N;//sufficient to get out of transition phase
            if (L==64){
                N_t = 4e7 + passi_eq;
                if(T>T_c-0.6 && T<T_c+0.4){
                    N_t+=2e7;
                }
            }
            else if(L==32){
                N_t = 1e7 + passi_eq;
                if(T>T_c-0.6 && T<T_c+0.4){
                    N_t+=7e6;
                }
            }
            else if(L==16){
                N_t = 5e6 + passi_eq;
                if(T>T_c-0.6 && T<T_c+0.4){
                    N_t+=3e6;
                }
            }
            else{
                 N_t = 2e6 + passi_eq;
                 if(T>T_c-0.6 && T<T_c+0.4){
                    N_t+=2e6;
                }
            }
            N_t += T*5e4;
        }
        if (c_Tm!=-1){//to show the equilibration
            passi_eq=0;
        }

        double V = 0, V_m = 0;
        double magn_media = 0, V_quadro_m = 0;
        double M_quadro_mean = 0, M_quarta_mean = 0;
        double chi = 0, chi_ist, chi_scaled;
        double c_v = 0, c_v_ist;
        double dim_media = 0;
        double binder_ratio = 0, binder_ratio_ist;


        magn = 0;
        for (int i = 0; i < N; ++i){
        //calculation for the initial value of potential and magnation
            V += potenziale_ising(r, J, i, prim_vic);
            magn += r(i,2)/N;
        }
        V /= 2;//if not every interaction is counted twice

        ofstream dati_blocking;
        dati_blocking.open("out/dati_blocking.txt");
        
        for (int i = 0; i < N_t; ++i){//iterations
            if(i>passi_eq){//calculation of the means
                V_m = V_m * (i - passi_eq - 1.0);
                V_quadro_m = V_quadro_m * (i - passi_eq - 1.0);
                dim_media = dim_media * (i - passi_eq - 1.0);

                magn_media = magn_media * (i - passi_eq - 1.0);
                M_quadro_mean = M_quadro_mean * (i - passi_eq - 1.0);
                M_quarta_mean = M_quarta_mean * (i - passi_eq - 1.0);
            }
            
            if(prog_t==0){
                MRT2(r, &V, T, J, prim_vic);
            }
            else if(prog_t==1){
                MRT2(r, &V, T, J, prim_vic, i);
            }
            else{
                Wolff(r, &V, T, J, prim_vic);
            }
    
            if(i>passi_eq){//calculation of the means
                V_m = (V_m + V) / (i - passi_eq);
                V_quadro_m = (V_quadro_m + pow1(V,2)) / (i - passi_eq);
                dim_media = (dim_media * N + dimensione) / ((i - passi_eq) * N);

                magn_media = (magn_media + abs(magn)) / (i - passi_eq);
                M_quadro_mean = (M_quadro_mean + pow1(magn,2)) / (i - passi_eq);
                M_quarta_mean = (M_quarta_mean + pow1(magn,4)) / (i - passi_eq);

                chi_ist = N * pow1(abs(magn) - magn_media, 2) / T;
                chi = N * (M_quadro_mean - pow1(magn_media, 2)) / T;
                chi_scaled = chi / pow(L,1.75);

                c_v_ist = pow1(V - V_m, 2) / (N * pow1(T,2));
                c_v = (V_quadro_m - pow1(V_m, 2)) / (N * pow1(T,2));

                binder_ratio = 1. - M_quarta_mean / (3 * pow1(M_quadro_mean, 2));
                binder_ratio_ist = 1. - pow1(magn,4) / (3*pow1(M_quadro_mean,2));

                dati_blocking << abs(magn) << "\t" << chi_ist << "\t";
                dati_blocking << c_v_ist << "\t" << binder_ratio_ist << "\t";
                dati_blocking << dimensione << endl;
                //used to find the actual error
                
                if (c_Tm!=-1){
                    dati << i << "\t" << dim_media*N << "\t" << c_v << endl;
                    osservabili << i << "\t" << magn_media << "\t" << chi;
                    osservabili << "\t" << binder_ratio << endl;
                }            
            }
            if (pv==4 && c_Tm!=-1 && c_Lm!=-1){
                if(i==0 || i==10-1 || i==1e2-1 || i==1e3-1 || i==1e4-1){
                    configurations(r,i,L);
                }
                if(i==1e5-1 || i==1e6-1 || i==1e7-1 || i==1e8-1){
                    configurations(r,i,L);
                }
                if (i==10 || i==11 || i==12 || i==13){
                    configurations(r,i,L);
                }
            }
        }
        dati_blocking.close();
        rowvec sigma(5, fill::value(0));
        // true error for magnation, magnetic susceptibility, specific heat
        // and dimension of the clusters
        err_max_blocking(N_t-passi_eq, sigma, 5);
        results_scaled << T_scaled*L << "\t" << chi_scaled << "\t";
        results_scaled << sigma(1) / pow(L,1.75) << endl;
        risultati << T << "\t" << magn_media << "\t" << sigma(0) << "\t" << chi;
        risultati << "\t" << sigma(1) << "\t" << binder_ratio << "\t" << sigma(3);
        risultati << "\t" << L << "\t" << c_v << "\t" << sigma(2) << "\t";
        risultati << dim_media << "\t" << sigma(4)/N << endl;
        if (c_Tm!=-1){
            blocking(N_t, 5);
            jackknife(N_t,5);
            break;
        }
        if (pv==4 && c_Lm!=-1 && c_Tm==-1){
            configurations(r,T,L);
        }
    }
    dati.close();
    osservabili.close();
}
