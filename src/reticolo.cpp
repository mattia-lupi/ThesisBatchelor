#include "reticolo.h"
#include <vector>
#include <cmath>

// Helper for squared distance to avoid sqrt()
inline double dist_sq(double x1, double y1, double x2, double y2) {
    return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

void stampa_reticolo(mat &r){
    // Optimization: Use buffer, avoid re-opening if possible (but keeping workflow)
    // Major fix: replaced endl with "\n"
    ofstream reticolo;
    reticolo.open("out/reticolo.txt");
    reticolo.precision(6);
    for (int i = 0; i < N; ++i){
        reticolo << r(i, 0) << "\t" << r(i,1) << "\t" << r(i,2) << "\n";
    }
    reticolo.close();
}

void crea_reticolo_quadrato(mat &r, double L){    
    int n = (int)L; // Assumes L is linear dimension based on main.cpp
    double L_cella = 1.0; // Usually 1.0 in standard Ising, code said L/n but usually L IS n.
    // Preserving your logic: if L passed is linear dim, and N = L*L
    
    int cont = 0;
    // Optimization: avoid repeated casts
    for (int i = 0; i < n; ++i){
        double x_pos = (double)i;
        for (int j = 0; j < n; ++j){
            // Direct access is slightly faster than list initialization
            r(cont, 0) = x_pos;
            r(cont, 1) = (double)j;
            r(cont, 2) = 1.0;
            cont++;
        }
    }
    stampa_reticolo(r);
}

// Helpers for Hexagonal
inline double aggiunta_y(int cont_y, double H){
    // Lookup table is faster than branching
    static const double vals[] = {-1.0, 0.0, 1.0, 0.0};
    // Adjust index to handle your specific H logic
    // Your logic: 0->-H, 1->0, 2->H, 3->0
    return vals[cont_y] * ((cont_y == 0 || cont_y == 2) ? H : 1.0); 
}

inline double aggiunta_x(int cont_x){
    return (cont_x == 0) ? 0.5 : 1.0;
}

void crea_reticolo_esagonale(mat &r, int L){
    double H = 0.86602540378; // sqrt(3)/2 precomputed
    int cont_x = 0;
    int cont_y = 0;
    int cont_part = 0; // Changed to 0-based for array indexing

    double y = H, x = 0;
    
    int max_N = L * L;
    
    // Optimization: while loop with direct array access
    while(cont_part < max_N){
        r(cont_part, 0) = x;
        r(cont_part, 1) = y;
        r(cont_part, 2) = 1;

        x += (cont_x == 0) ? 0.5 : 1.0;
        
        // Logic for Y addition
        if(cont_y == 0) y -= H;
        else if(cont_y == 2) y += H;
        // else y += 0
        
        // Update counters
        cont_y = (cont_y + 1) % 4;
        cont_x = (cont_x == 1) ? 0 : 1; // Toggle 0/1

        cont_part++;
        
        // Chain reset logic
        if(cont_part % L == 0 && cont_part != max_N){
            x = 0;
            y += 2 * H;
            cont_x = 0;
            cont_y = 0;
        }
    }
    stampa_reticolo(r);
}

void crea_primi_vicini(double L_double, string p_v_path){
    int L = (int)L_double;
    
    // Open output file
    ofstream primivicini;
    primivicini.open(p_v_path);

    // ---------------------------------------------------------
    // OPTIMIZATION 1: SQUARE LATTICE (Analytic Solution)
    // Complexity: O(N) instead of O(N^2)
    // ---------------------------------------------------------
    if(pv == 4){
        // We know exactly where the neighbors are. 
        // No need to read file, no need to calculate distances.
        // Lattice is filled: x changes fast, y changes slow (or vice versa depending on loop)
        // Your creation loop: i (outer), j (inner). 
        // Row 0: (0,0), (0,1)... (0, L-1)
        // This means index = i*L + j. 
        // "j" corresponds to Y (inner loop), "i" corresponds to X (outer loop).
        // Let's assume standard row-major mapping for safety.
        
        // Note: Based on crea_reticolo_quadrato: 
        // i increments X, j increments Y. 
        // So indices run along Y first.
        
        for(int i = 0; i < N; ++i){
            int neighbors[4];
            
            // Current coordinates in grid
            // Since loop was i(0..L), j(0..L) -> cont++
            // index = i * L + j
            int row = i / L; // x-coordinate index
            int col = i % L; // y-coordinate index

            // Neighbor 1: Same X, Y-1 (Left in matrix terms)
            int col_down = (col == 0) ? L - 1 : col - 1;
            neighbors[0] = row * L + col_down;

            // Neighbor 2: Same X, Y+1 (Right)
            int col_up = (col == L - 1) ? 0 : col + 1;
            neighbors[1] = row * L + col_up;

            // Neighbor 3: X-1, Same Y (Up)
            int row_prev = (row == 0) ? L - 1 : row - 1;
            neighbors[2] = row_prev * L + col;

            // Neighbor 4: X+1, Same Y (Down)
            int row_next = (row == L - 1) ? 0 : row + 1;
            neighbors[3] = row_next * L + col;

            // Write
            primivicini << neighbors[0] << "\t" << neighbors[1] << "\t" 
                        << neighbors[2] << "\t" << neighbors[3] << "\n";
        }
    }
    // ---------------------------------------------------------
    // OPTIMIZATION 2: HEX/GENERIC LATTICE (Geometric Solution)
    // Complexity: O(N^2) - but optimized memory/math
    // ---------------------------------------------------------
    else {
        // We MUST read the file because main doesn't pass 'r' to this function.
        // Use std::vector for speed, not armadillo
        struct Point { double x; double y; };
        std::vector<Point> coords(N);
        
        ifstream reticolo;
        reticolo.open("out/reticolo.txt");
        double dummy_spin;
        for (int i = 0; i < N; ++i){
            reticolo >> coords[i].x >> coords[i].y >> dummy_spin;
        }
        reticolo.close();

        double L_y = sqrt(3) * L;
        double L_x = L * 3.0 / 4.0;
        
        // Pre-squared threshold (1.1^2 = 1.21)
        double threshold_sq = 1.21;
        double inv_Lx = 1.0 / L_x;
        double inv_Ly = 1.0 / L_y;

        for (int i = 0; i < N; ++i){
            int found = 0;
            // Write buffer to avoid multiple << calls
            int n_indices[10]; // ample space

            for (int j = 0; j < N; ++j){
                if(i == j) continue;

                double dx = coords[j].x - coords[i].x;
                double dy = coords[j].y - coords[i].y;

                // Manual rint calculation is faster: floor(x + 0.5)
                // Periodic Boundary Conditions
                if (pv == 3) {
                     dx -= L_x * floor(dx * inv_Lx + 0.5);
                     dy -= L_y * floor(dy * inv_Ly + 0.5);
                } else {
                     dx -= L * floor(dx / L + 0.5);
                     dy -= L * floor(dy / L + 0.5);
                }

                // Optimization: Squared distance check
                if ((dx*dx + dy*dy) <= threshold_sq){
                    n_indices[found] = j;
                    found++;
                    if(found == pv) break;
                }
            }
            
            // Write line
            for (int k = 0; k < pv; ++k){
                primivicini << n_indices[k] << (k == pv-1 ? "\n" : "\t");
            }
        }
    }
    primivicini.close();
}
