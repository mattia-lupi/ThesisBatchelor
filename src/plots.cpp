#include "plots.h"

void plot_reticolo() {//executes the plot of the lattice
    string comando;
    comando = "gnuplot";
    comando += " plot_reticolo.plt";
    cout << comando << endl;
    system(comando.c_str());
    return;
}
void blocking_plot(){//executes the plot of the blocking
    //ora faccio il plot
    string comando;

    comando = "gnuplot";
    comando += " plot_blocking.plt";
    cout << comando << endl;
    system(comando.c_str());
}
void plot_osservabili_T(int c_Tm){
//executes the plot of the observables as functions of T
    if(c_Tm!=-1){
        string comando;
        comando = "gnuplot";
        comando += " plot_osservabili.plt";
        cout << comando << endl;
        system(comando.c_str());
        return;
    }
    else{
        string comando;
        comando = "gnuplot";
        comando += " plot_osservabili_T.plt";
        cout << comando << endl;
        system(comando.c_str());
        return;
    }
}

void plot_config() {//executes the plot of the configurations
    string comando;
    comando = "gnuplot";
    comando += " out/config.plt";
    cout << comando << endl;
    system(comando.c_str());
    return;
}
void configurations(mat &r, int iter, double L){
//configurations in iterations at the same temperature
    ofstream config;
    config.open("out/config.plt");
    config<<"set term png"<<endl;
    config<< "set output \"out/config/ising_2d_iter_"<<iter+1<<".png\""<<endl;
    config<<"set xrange [ 0 :   "<<L<<" ]"<<endl;
    config<<"set yrange [ 0 :   "<<L<<" ]"<<endl;
    config<<"set nokey"<<endl;
    config<<"set title \"Configuration at "<<iter+1<<" iterations for N="<<L;
    config<<"x"<<L<<"\""<<endl;
    config<<"unset tics"<<endl;
    config<<"set size ratio    1.00000"<<endl;
    for (int i = 0; i < L; ++i){
        for (int j = 0; j < L; ++j){
            if(r(i*L+j,2)==1){
                config<<"set object rectangle from "<<i<<", "<< j <<" to " <<i+1;
                config<<", "<< j+1 << "fc rgb \"white\""<<endl;
            }
            else{
                config<<"set object rectangle from "<<i<<", "<< j <<" to " <<i+1;
                config<<", "<< j+1 << "fc rgb \"black\""<<endl;
            }
        }
    }
    config<<"plot 1"<<endl;
    config<<"quit"<<endl;
    config.close();
    plot_config();
}
void configurations(mat &r, double T, double L){//configurations in temperature
    ofstream config;
    config.open("out/config.plt");
    config<<"set term png"<<endl;
    config<< "set output \"out/config/ising_2d_T_"<<T<<".png\""<<endl;
    config<<"set xrange [ 0 :   "<<L<<" ]"<<endl;
    config<<"set yrange [ 0 :   "<<L<<" ]"<<endl;
    config<<"set nokey"<<endl;
    config<<"set title \"Configuration at temperature "<<T<<" for N="<<L;
    config<<"x"<<L<<"\""<<endl;
    config<<"unset tics"<<endl;
    config<<"set size ratio    1.00000"<<endl;
    for (int i = 0; i < L; ++i){
        for (int j = 0; j < L; ++j){
            if(r(i*L+j,2)==1){
                config<<"set object rectangle from "<<i<<", "<< j <<" to " <<i+1;
                config<<", "<< j+1 << "fc rgb \"white\""<<endl;
            }
            else{
                config<<"set object rectangle from "<<i<<", "<< j <<" to " <<i+1;
                config<<", "<< j+1 << "fc rgb \"black\""<<endl;
            }
        }
    }
    config<<"plot 1"<<endl;
    config<<"quit"<<endl;
    config.close();
    plot_config();
}
