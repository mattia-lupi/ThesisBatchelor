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



#define ARMA_NO_DEBUG //no control and more speed

#include "funzioni.h"
#include "errors.h"
#include "plots.h"
#include "reticolo.h"
#include <chrono>
#include <vector> // Added for speed over arma::mat where possible

/*** variabili globali ***/

int N; //number of spins
int pv = 4;//number of first nearest neighbours
double magn = 0;
int dimensione = 0;

string p_v_path = "out/primivicini.txt";
string obs_path = "out/osservabili.txt";
string dati_path = "out/dati.txt";
string risultati_path = "out/risultati.txt";

ofstream dati, osservabili, risultati, results_scaled;

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
    int c_Lm = -1;

    //critical temperature value for square and exagonal lattice
    double T_c[] = {2.269, 1.519};
    rowvec L = {8, 16, 32, 64};

    int initial_data[] = {c_Lm, c_Tm, N_t, prog_t};

    auto start = std::chrono::high_resolution_clock::now();
    obs(L, T_c, J, initial_data);
    auto end = std::chrono::high_resolution_clock::now();
    
    chrono::duration<double> elapsed = end - start;
    cout << "Elapsed time: " << elapsed.count() << " s\n";

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
    
    int caso_max = (c_Lm == -1) ? L.size() : c_Lm + 1;
    int start = (c_Lm == -1) ? c_Lm + 1 : c_Lm;
    
    double T_c;
    
    for (int caso = start; caso < caso_max; ++caso){
        N = (int)(L(caso) * L(caso)); 
        mat r(N,3);
        
        if(pv==4){
            crea_reticolo_quadrato(r, L(caso));
            T_c = T_c_vec[0];
        }
        else if(pv==3){
            crea_reticolo_esagonale(r, L(caso));
            T_c = T_c_vec[1];
        }
        else{
            cout << "Lattice not available, the square lattice is going to be presented" << "\n"; 
            pv=4;
            crea_reticolo_quadrato(r, L(caso));
            T_c = T_c_vec[0];
        }
        
        //creates a file with the first nearest neighbours for each spin
        crea_primi_vicini(L(caso), p_v_path);
        
        if (c_Lm == -1){
            cout << "Executing lattice of N = " << L(caso) << "x" << L(caso) << "\n";
        }
        
        if (caso != 0 && c_Lm == -1){
            risultati << "\n\n";
            results_scaled << "\n\n";
        }
        
        risultati << "L=" << L(caso) << "\n";
        results_scaled << "L=" << L(caso) << "\n";

        // makes the simulation at prescribed temperature for the chosen lattice size
        obs_T(r, J, T_c, L(caso), data);

        //if I want how a single observable progresses over the iterations
        // I execute it only once
        if(c_Tm != -1){
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

    // OPTIMIZATION: Use std::vector for faster access than arma::mat
    // If MRT2/Wolff strictly require mat, we keep mat but optimize loading.
    mat prim_vic(N, pv);
    
    ifstream pr_v;
    pr_v.open(p_v_path);
    if(pr_v.is_open()) {
        for (int i = 0; i < N; ++i){
            for (int j = 0; j < pv; ++j){
                pr_v >> prim_vic(i,j);
            }
        }
        pr_v.close();
    }

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
        
        cout << "Computing T = " << T << "\n"; // Replaced endl
        int passi_eq;

        // Optimization: Pre-calculate constants for N_t logic
        if(prog_t==2){
            T_scaled = (T - T_c) / T_c;
            passi_eq = 5 * N;
            N_t = 20000 + passi_eq + (int)(T * 300);
            
            if(L==32){
                passi_eq += 9*N; 
                N_t = 50000 + passi_eq + (int)(T * 1000);
            }
            else if(L==64){
                passi_eq += 12*N; 
                N_t = 70000 + passi_eq + (int)(T * 2000);
            }
        }
        else{
            T_scaled = (T - T_c) / T_c;
            passi_eq = 70 * N * N;
            
            if (L==64){
                N_t = 40000000 + passi_eq;
                if(T>T_c-0.6 && T<T_c+0.4) N_t += 20000000;
            }
            else if(L==32){
                N_t = 10000000 + passi_eq;
                if(T>T_c-0.6 && T<T_c+0.4) N_t += 7000000;
            }
            else if(L==16){
                N_t = 5000000 + passi_eq;
                if(T>T_c-0.6 && T<T_c+0.4) N_t += 3000000;
            }
            else{
                 N_t = 2000000 + passi_eq;
                 if(T>T_c-0.6 && T<T_c+0.4) N_t += 2000000;
            }
            N_t += (int)(T*50000);
        }
        
        //to show the equilibration
        if (c_Tm!=-1){
            passi_eq=0;
        }

        // ACCUMULATORS for speed (avoid multiply/divide every step)
        double sum_V = 0, sum_V2 = 0;
        double sum_M = 0, sum_M2 = 0, sum_M4 = 0;
        double sum_dim = 0;

        // Result holders
        double V_m = 0, V_quadro_m = 0, dim_media = 0;
        double magn_media = 0, M_quadro_mean = 0, M_quarta_mean = 0;
        double chi = 0, chi_ist, chi_scaled;
        double c_v = 0, c_v_ist;
        double binder_ratio = 0, binder_ratio_ist;
        
        // Pre-compute inverse temperature for speed
        double inv_T = 1.0 / T;
        double inv_T2 = 1.0 / (T * T);
        double inv_N = 1.0 / N;
        double L_pow_175 = pow(L, 1.75);

        double V = 0;
        magn = 0;
        
        // Initial configuration calculation
        for (int i = 0; i < N; ++i){
            V += potenziale_ising(r, J, i, prim_vic);
            magn += r(i,2); 
        }
        magn *= inv_N;
        V *= 0.5;

        ofstream dati_blocking;
        dati_blocking.open("out/dati_blocking.txt");
        // Ensure standard formatting is efficient
        dati_blocking.precision(6);
       
        for (int i = 0; i < N_t; ++i){
            
            // PHYSICS ENGINE CALLS
            if(prog_t==0){
                MRT2(r, &V, T, J, prim_vic);
            }
            else if(prog_t==1){
                MRT2(r, &V, T, J, prim_vic, i);
            }
            else{
                Wolff(r, &V, T, J, prim_vic);
            }
   
            if(i > passi_eq){
                // UPDATE ACCUMULATORS (Faster than recursive mean)
                double abs_magn = std::abs(magn);
                double magn2 = magn * magn;
                
                sum_V += V;
                sum_V2 += V * V;
                sum_dim += dimensione; // Assuming 'dimensione' is global updated by Wolff
                
                sum_M += abs_magn;
                sum_M2 += magn2;
                sum_M4 += magn2 * magn2;
                
                // CALCULATE INSTANTANEOUS MEAN
                double count = (double)(i - passi_eq);
                double inv_count = 1.0 / count;

                V_m = sum_V * inv_count;
                V_quadro_m = sum_V2 * inv_count;
                dim_media = sum_dim * inv_count * inv_N; // combined divisions

                magn_media = sum_M * inv_count;
                M_quadro_mean = sum_M2 * inv_count;
                M_quarta_mean = sum_M4 * inv_count;

                // OBSERVABLES
                double diff_magn = abs_magn - magn_media;
                chi_ist = N * (diff_magn * diff_magn) * inv_T;
                
                double diff_mag_sq = M_quadro_mean - (magn_media * magn_media);
                chi = N * diff_mag_sq * inv_T;
                
                chi_scaled = chi / L_pow_175;

                double diff_V = V - V_m;
                c_v_ist = (diff_V * diff_V) * inv_N * inv_T2;
                
                double diff_V_sq = V_quadro_m - (V_m * V_m);
                c_v = diff_V_sq * inv_N * inv_T2;

                binder_ratio = 1.0 - M_quarta_mean / (3.0 * M_quadro_mean * M_quadro_mean);
                binder_ratio_ist = 1.0 - (magn2 * magn2) / (3.0 * M_quadro_mean * M_quadro_mean);

                dati_blocking << abs_magn << "\t" << chi_ist << "\t"
                              << c_v_ist << "\t" << binder_ratio_ist << "\t"
                              << dimensione << "\n";
               
                if (c_Tm!=-1){
                    dati << i << "\t" << dim_media*N << "\t" << c_v << "\n";
                    osservabili << i << "\t" << magn_media << "\t" << chi
                                << "\t" << binder_ratio << "\n";
                }            
            }
            
            // SNAPSHOT LOGIC
            if (pv==4 && c_Tm!=-1 && c_Lm!=-1){
                // Optimized check using simple integer comparisons
                bool save = (i==0 || i==9 || i==99 || i==999 || i==9999 || 
                             i==99999 || i==999999 || i==9999999 || i==99999999 ||
                             i==10 || i==11 || i==12 || i==13);
                if(save){
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
        results_scaled << sigma(1) / L_pow_175 << "\n";
        
        risultati << T << "\t" << magn_media << "\t" << sigma(0) << "\t" << chi;
        risultati << "\t" << sigma(1) << "\t" << binder_ratio << "\t" << sigma(3);
        risultati << "\t" << L << "\t" << c_v << "\t" << sigma(2) << "\t";
        risultati << dim_media << "\t" << sigma(4)/N << "\n";

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
