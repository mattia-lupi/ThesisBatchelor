#include "errors.h"

void blocking(int N_t, int num_oss){//blocking method
    mat O(num_oss, N_t);

    ifstream dati_blocking;
    dati_blocking.open("out/dati_blocking.txt");
    for (int i = 0; i < N_t; ++i){
        for (int j = 0; j < num_oss; ++j){
            dati_blocking >> O(j, i);
        }
    }
    dati_blocking.close();

    ofstream blocking;
    blocking.open("out/blocking.txt");

    int N_B_prev=0;
    for (int B = 10; B < N_t/4; B+=10){

        int N_B = floor(N_t / B);
        
        if(N_B!=N_B_prev){
            
            mat O_mB(num_oss, N_B, fill::value(0));
            //mean observable of the block
            rowvec var_O(num_oss, fill::value(0));
            //variance of each observable
            rowvec O_mean(num_oss,fill::value(0));
            //complessive mean for each observable
            
            for (int i = 0; i < N_B; ++i){//cycle over the blocks
                for (int j = 0; j < B; ++j){
                    for (int k = 0; k < num_oss; ++k){
                        O_mB(k, i) += O(k, (i * B + j)) / B;
                    }
                }
                for (int k = 0; k < num_oss; ++k){
                    O_mean(k) += O_mB(k,i) / N_B;
                    //calculates the mean over all the blocks
                }
            }
            for (int k = 0; k < num_oss; ++k){
                for (int i = 0; i < N_B; ++i){
                    var_O(k) += pow1(O_mB(k,i) - O_mean(k),2) / N_B;
                }
            }
            blocking << B << "\t";
            for (int i = 0; i < num_oss; ++i){
                if(i==num_oss-1){
                    blocking << sqrt(var_O(i)/N_B) << endl;//error of a mean
                }
                else{
                    blocking << sqrt(var_O(i)/N_B) << "\t";
                }
            }
            N_B_prev=N_B;
        }
    }
    blocking.close();
    // blocking_plot();
}
void jackknife(int N_t, int num_oss){//jackknife method
    mat O(num_oss, N_t);

    ifstream dati_blocking;
    ofstream jackknife;
    
    jackknife.open("out/jackknife.txt");
    dati_blocking.open("out/dati_blocking.txt");
    
    for (int i = 0; i < N_t; ++i){
        for (int j = 0; j < num_oss; ++j){
            dati_blocking >> O(j, i);
        }
    }
    dati_blocking.close();

    int N_B_prev=1;//scarto il blocco singolo

    for (int B = (int)(N_t / 1e4) + 1; B < N_t/2; B+=10){
    //takes away blocks which are too small 
        int N_B = floor(N_t / B);
        
        if(N_B!=N_B_prev){
            mat O_m(num_oss, N_B, fill::value(0));
            //mean of the observables of each block
            mat O_m_jack(num_oss, N_B, fill::value(0));
            //means found by excluding a block

            rowvec var_O(num_oss, fill::value(0));
            //variance of each observable
            rowvec O_media(num_oss,fill::value(0));
            //complessive mean for each observable

            //calculate the means of each block
            for (int i = 0; i < N_B; ++i){//cycle over the blocks
                for (int j = 0; j < B; ++j){
                //cycle over the points of the block
                    for (int k = 0; k < num_oss; ++k){
                        O_m(k, i) += O(k, (i * B + j)) / B;
                    }
                }
            }

            for (int i = 0; i < N_B; ++i){
                for (int j = 0; j < N_B; ++j){
                    if(j!=i){//exclude one different block each cycle
                        for (int k = 0; k < num_oss; ++k){
                            O_m_jack(k, i) += O_m(k, j) / (N_B - 1.0);
                        }
                    }
                }
                for (int k = 0; k < num_oss; ++k){
                    O_media(k) += O_m_jack(k, i) / N_B;
                }
            }

            for (int i = 0; i < N_B; ++i){
                for (int k = 0; k < num_oss; ++k){
                    var_O(k) += pow1(O_m_jack(k,i) - O_media(k), 2);
                }
            }
            jackknife << B << "\t";
            for (int i = 0; i < num_oss; ++i){
                if (i!=num_oss-1){
                    jackknife << sqrt((N_B - 1) * var_O(i) / N_B) << "\t";
                }
                else{
                    jackknife << sqrt((N_B - 1) * var_O(i) / N_B) << endl;
                }
            }
            N_B_prev=N_B;
        }
    }

    jackknife.close();
}
void err_max_blocking(int N_t, rowvec &sigma, int num_oss){
//blocking method to find the max error to use on the errorbars for the program
    mat O(num_oss, N_t);

    ifstream dati_blocking;
    dati_blocking.open("out/dati_blocking.txt");
    for (int i = 0; i < N_t; ++i){
        for (int j = 0; j < num_oss; ++j){
            dati_blocking >> O(j, i);
        }
    }
    dati_blocking.close();

    int N_B_prev=0;
    sigma.fill(0);
    for (int B = 1; B < N_t/4; B+=10){

        int N_B = floor(N_t / B);
        
        if(N_B!=N_B_prev){
            
            mat O_mB(num_oss, N_B, fill::value(0));
            //mean observable of the block
            rowvec var_O(num_oss, fill::value(0));
            //variance of each observable
            rowvec O_mean(num_oss,fill::value(0));
            //complessive mean for each observable
            
            for (int i = 0; i < N_B; ++i){//cycle over the blocks
                for (int j = 0; j < B; ++j){
                    for (int k = 0; k < num_oss; ++k){
                        O_mB(k, i) += O(k, (i * B + j)) / B;
                    }
                }
                for (int k = 0; k < num_oss; ++k){
                    O_mean(k) += O_mB(k,i) / N_B;
                    //calculates the mean over all the blocks
                }
            }
            for (int k = 0; k < num_oss; ++k){
                for (int i = 0; i < N_B; ++i){
                    var_O(k) += pow1(O_mB(k,i) - O_mean(k),2) / N_B;
                }
            }
            for (int i = 0; i < num_oss; ++i){
                if (sigma(i)<sqrt(var_O(i)/N_B)){
                    sigma(i) = sqrt(var_O(i)/N_B);
                }
            }
            
            N_B_prev=N_B;
        }
    }
}
