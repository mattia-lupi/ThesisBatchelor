#ifndef FUNZIONI_H
#define FUNZIONI_H

#include <cmath>
#include <vector>
#include <cstdlib>
#include <armadillo>

using namespace std;
using namespace arma;

// Global variables (from main)
extern int N;
extern int pv;
extern double magn;
extern int dimensione;

// --- MATH HELPERS ---

// Standard x*x is significantly faster than loop-based pow or std::pow
inline double pow2(double x) {
    return x * x;
}

// Keep generic for compatibility, but don't use for squares
double pow1(double base, int esp) {
    double ris = 1.0;
    for (int i = 0; i < esp; i++) {
        ris *= base;
    }
    return ris;
}

// Optimized module calculation
double mod(rowvec &r){ 
    return sqrt(pow2(r(0)) + pow2(r(1)));
}

double mod(cube &r, double riga, double colonna){
    return sqrt(pow2(r(riga,colonna,0)) + pow2(r(riga,colonna,1)));
}

// --- PHYSICS KERNELS (Raw Pointer Optimized) ---

// Calculates potential for one spin using raw pointers for speed
// spins: pointer to the spin column (r.colptr(2))
// primivicini: raw pointer to neighbor matrix data
inline double potenziale_ising_raw(double* spins, double J, int i, double* neighbors_ptr) {
    double V = 0;
    // Armadillo stores matrices column-major, but if we treat neighbors as flat:
    // We assume primivicini is N x pv. 
    // Accessing raw arma mat: M[row + col*N_rows]
    
    // OPTIMIZATION: Manually unrolled loop for pv=4 (common case)
    if (pv == 4) {
        // Neighbors are stored in columns 0,1,2,3. 
        // Index in raw memory for (i, j) is i + j*N
        int n1 = (int)neighbors_ptr[i];         // col 0
        int n2 = (int)neighbors_ptr[i + N];     // col 1
        int n3 = (int)neighbors_ptr[i + 2*N];   // col 2
        int n4 = (int)neighbors_ptr[i + 3*N];   // col 3
        
        double s_i = spins[i];
        V += (spins[n1] + spins[n2] + spins[n3] + spins[n4]) * s_i;
    } else {
        // Fallback for hex (pv=3) or others
        double s_i = spins[i];
        double sum_neighbors = 0;
        for (int j = 0; j < pv; ++j) {
            int n_idx = (int)neighbors_ptr[i + j*N];
            sum_neighbors += spins[n_idx];
        }
        V += sum_neighbors * s_i;
    }
    
    return -V * J;
}

// Wrapper for compatibility with existing code if called elsewhere
double potenziale_ising(mat &r, double J, int i, mat &primivicini){
    return potenziale_ising_raw(r.colptr(2), J, i, primivicini.memptr());
}

// Metropolis Algorithm (Optimized)
void MRT2(mat &r, double *V, double T, double J, mat &primivicini, int i = -999){
    
    // 1. Get Raw Pointers (Avoids overhead of r(i,2) checks)
    double* S = r.colptr(2); 
    double* P = primivicini.memptr();
    
    int m;
    // Sequential vs Random selection
    if (i == -999){
        // Fast random integer generation
        m = rand() % N; 
    } else {
        m = i % N;
    }

    // 2. Calculate Energy Change Inline (Avoid function call overhead)
    // dE = 2 * J * S_m * sum(S_neighbors)
    double sum_neighbors = 0;
    if (pv == 4) {
        sum_neighbors += S[(int)P[m]] + S[(int)P[m + N]] + S[(int)P[m + 2*N]] + S[(int)P[m + 3*N]];
    } else {
        for (int j = 0; j < pv; ++j) sum_neighbors += S[(int)P[m + j*N]];
    }

    double dE = 2.0 * J * S[m] * sum_neighbors; // Cost to flip
    
    // 3. Metropolis Check
    // If dE <= 0, accept. If dE > 0, accept with prob exp(-dE/T).
    // Note: In your original code, you used V_mod. 
    // V_mod from your code calculated: V += neighbor * current. V *= J*2*current.
    // This is equivalent to dE.
    
    bool accept = false;
    if (dE <= 0) {
        accept = true;
    } else {
        // Optimization: standard rand() is faster than the floating div logic every time
        if (exp(-dE / T) > (double)rand() / (double)RAND_MAX) {
            accept = true;
        }
    }

    if (accept) {
        S[m] *= -1; // Flip spin
        magn += 2.0 * S[m] / N;
        *V += dE;
    }
}

// Wolff Algorithm (Iterative & Memory Optimized)
void Wolff(mat &r, double *V, double T, double J, mat &primivicini){
    
    // Raw pointers
    double* S = r.colptr(2);
    double* P = primivicini.memptr();

    // 1. Select Seed
    int seed = rand() % N;
    
    // 2. Setup Probabilities
    double prob = 1.0 - exp(-2.0 * J / T);
    
    // 3. Iterative Stack (Replaces Recursion and Replaces spin_da_mod array)
    // We use a static vector to avoid re-allocating memory every step.
    static std::vector<int> stack;
    stack.clear(); // Reset stack, keeps capacity (very fast)
    
    // 4. Process Seed
    // We flip the spin IMMEDIATELY. 
    // This removes the need for a "visited" array.
    // Logic: If a neighbor has a DIFFERENT spin than the current one (after flip),
    // it means they WERE parallel before check.
    double old_spin_val = S[seed];
    S[seed] *= -1; // Flip seed
    magn += 2.0 * S[seed] / N;
    dimensione = 1;
    stack.push_back(seed);

    // 5. Cluster Growth
    while (!stack.empty()) {
        int current = stack.back();
        stack.pop_back();

        // Check neighbors
        for (int j = 0; j < pv; ++j) {
            int n_idx = (int)P[current + j*N]; // Access neighbor index
            
            // If neighbor has the OLD spin value (meaning they were parallel to us before we flipped)
            if (S[n_idx] == old_spin_val) {
                if ((double)rand() / RAND_MAX < prob) {
                    S[n_idx] *= -1; // Flip neighbor
                    magn += 2.0 * S[n_idx] / N;
                    dimensione++;
                    stack.push_back(n_idx);
                }
            }
        }
    }

    // 6. Recalculate Energy
    // Since we flipped a random cluster, dE is hard to track incrementally.
    // We recalculate O(N) but using the fast raw pointer method.
    double total_V = 0;
    // Loop over all spins
    for (int i = 0; i < N; ++i) {
        // To avoid double counting, we usually sum over all pairs and divide by 2
        // Or loop over specific bonds. Here we follow your original logic:
        double local_V = 0;
        double s_i = S[i];
        
        // Only sum i < j to avoid double counting? 
        // Your original code calculated full V then divided by 2 at the end. We stick to that.
        if (pv == 4) {
             local_V += (S[(int)P[i]] + S[(int)P[i + N]] + S[(int)P[i + 2*N]] + S[(int)P[i + 3*N]]) * s_i;
        } else {
             for (int j = 0; j < pv; ++j) local_V += S[(int)P[i + j*N]] * s_i;
        }
        total_V -= local_V * J;
    }
    *V = total_V / 2.0;
}

#endif
