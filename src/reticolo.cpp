#include "reticolo.h"



void stampa_reticolo(mat &r){//prints the lattice on a file
    ofstream reticolo;
    reticolo.open("out/reticolo.txt");
    for (int i = 0; i < N; ++i){
        reticolo << r(i, 0) << "\t" << r(i,1) << "\t" << r(i,2) << endl;
        //0 posizione x, 1 posizione y, 2 spin
    }
    reticolo.close();
    // plot_reticolo();
    return;
}
void crea_reticolo_quadrato(mat &r, double L){// creates the square lattice    
    int n = sqrt(N);
    double L_cella = L / n;
    int cont = 0;
    for (int i = 0; i < n; ++i){
        for (int j = 0; j < n; ++j){
            r.row(cont) = {(double)i * L_cella, (double)j * L_cella, 1};
            cont++;
        }
    }
    stampa_reticolo(r);
    return;
}
double aggiunta_y(int cont_y, double H){
//needed for the creation of the exagonal lattice
    if(cont_y==0){
        return -H;
    }
    else if (cont_y==1 || cont_y==3){
        return 0;
    }
    else{
        return H;
    }
}
double aggiunta_x(int cont_x){//needed for the creation of the exagonal lattice
    if(cont_x==0){
        return 0.5;
    }
    else{
        return 1;
    }
}
void crea_reticolo_esagonale(mat &r, int L){// creates the exagonal lattice
    double H = sqrt(3.)/2.;//"height" of half a hexagon
    int cont_x = 0;
    int cont_y = 0;
    int cont_part = 1;

    double y = H, x = 0;
    do{
        r.row(cont_part-1) = {x, y, 1};
        x += aggiunta_x(cont_x);
        y += aggiunta_y(cont_y, H);
        if(cont_y==3){
            cont_y=0;
        }
        else{
            cont_y++;
        }
        if(cont_x==1){
            cont_x--;
        }
        else{ 
            cont_x++;
        }
        if(cont_part%(L)==0){//at the start of the chain
            x = 0;
            y += 2 * H;
            cont_x = 0;
            cont_y = 0;
        }
        cont_part++;
    }
    while(cont_part!=L*L+1);
    stampa_reticolo(r);
    return;
}

void crea_primi_vicini(double L, string p_v_path){
//creates a file with the first nearest neighbours for each spin
    double L_y = sqrt(3)*(L), L_x = L*3/4;
    
    mat prim_vic(N,pv), r(N,3);
    rowvec dr(2);

    ifstream reticolo;
    reticolo.open("out/reticolo.txt");

    for (int i = 0; i < N; ++i){
        reticolo >> r(i,0) >> r(i,1) >> r(i,2);
    }
    reticolo.close();

    ofstream primivicini;
    primivicini.open(p_v_path);

    for (int i = 0; i < N; ++i){
    //creates the matrix for the first nearest neighbours
        int cont = 0;
        
        for (int j = 0; j < N; ++j){
            if(i!=j){
                for (int k = 0; k < 2; ++k){
                    dr(k) = r(j,k) - r(i,k);
                    if(pv==4){
                        dr(k) -= L * rint(dr(k)/L);//moves in [-L/2,+L/2]
                    }
                }
                if(pv==3){
                    dr(0) -= L_x * rint(dr(0)/L_x);//moves in [-L/2,+L/2]
                    dr(1) -= L_y * rint(dr(1)/L_y);//moves in [-L/2,+L/2]
                }
                
                if(cont==pv){
                    break;
                }
                
                if(mod(dr) <= 1.1){
                //if they are close as much as 1 (1.1 to be sure even with
                //approximation errors) they are first nearest neighbours
                    prim_vic(i,cont) = j;
                    cont++;
                }
            }
        }
        for (int j = 0; j < pv; ++j){
            if(j!=pv-1){
                primivicini << prim_vic(i,j) << "\t";
            }
            else{
                primivicini << prim_vic(i,j) << endl;
            }
        }
    }
    primivicini.close();
}
