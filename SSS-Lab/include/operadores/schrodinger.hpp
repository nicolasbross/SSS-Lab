#include <vector>
#include <functional>
#include <Eigen/Sparse>
#include <Eigen/Dense>

// Concept 1: The Grid Class
// Represents the discrete spatial domain.
class Grid {
public:
    int N_x, N_y;       // Number of points in X and Y
    double h;           // Step size (assumed uniform for simplicity)
    double length;      // Physical length of the domain side

    // Constructor corresponds to Lua's Grid.new(256, 256)
    Grid(int width, int height, double physical_length = 1.0);

    // Helper to flatten 2D (i, j) coordinates to 1D index k
    // Necessary because matrices are 1D structures.
    int idx(int i, int j) const;
    
    // Reverse mapping (optional but useful for debugging)
    std::pair<int, int> coord(int k) const;
};

// Concept 2: The Schrödinger Operator Class
// Manages the Hamiltonian H = -Laplacian + V
class Schrodinger {
private:
    Grid grid;
    
    // Concept 3: Eigen Sparse Matrix
    // Stores the Hamiltonian. Most entries are zero (sparse), so this saves massive RAM.
    Eigen::SparseMatrix<double> hamiltonian;
    
    // Stores the Potential V(x,y)
    Eigen::VectorXd V;

    // Stores computed eigenfunctions (columns) and eigenvalues
    Eigen::MatrixXd eigenfunctions;
    Eigen::VectorXd eigenvalues;

    bool built; // Flag to check if matrix needs rebuilding

public:
    // Constructor corresponds to Schrodinger.new(domain)
    Schrodinger(Grid g);

    // Concept 4: std::function
    // Allows passing a Lua function (or C++ lambda) into the core logic.
    void set_potential(std::function<double(double, double)> func);

    // Builds the finite difference Laplacian
    void build_operator();

    // Solves H u = lambda u
    // Returns the eigenvalues for Lua to print
    std::vector<double> compute_spectrum(int num_eigenvalues);

    // Concept 5: Analytic Quantity Computation
    // Computes the doubling index: log2( int(B_2r) / int(B_r) )
    double analyze_doubling_index(double x0, double y0, double r);
};
2. The Implementation File (SchrodingerCore.cpp)
This file contains the "heavy lifting": the finite difference stencil, the integration logic for the doubling index, and the matrix construction.

C++

#include "SchrodingerCore.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

// ==========================================
// Grid Implementation
// ==========================================

Grid::Grid(int width, int height, double physical_length) 
    : N_x(width), N_y(height), length(physical_length) {
    // Determine step size h based on N and physical length
    h = length / (std::max(width, height) - 1);
}

int Grid::idx(int i, int j) const {
    // Row-major ordering: index = i * width + j
    // Boundary checks should ideally be here or in the caller.
    return i * N_y + j;
}

std::pair<int, int> Grid::coord(int k) const {
    return {k / N_y, k % N_y};
}

// ==========================================
// Schrodinger Implementation
// ==========================================

Schrodinger::Schrodinger(Grid g) : grid(g), built(false) {
    int size = grid.N_x * grid.N_y;
    V.resize(size);
    V.setZero(); // Default to free potential (V=0)
}

void Schrodinger::set_potential(std::function<double(double, double)> func) {
    // Loop over the grid and evaluate the user's function
    for (int i = 0; i < grid.N_x; ++i) {
        for (int j = 0; j < grid.N_y; ++j) {
            double x = i * grid.h - (grid.length / 2.0); // Centering coordinates
            double y = j * grid.h - (grid.length / 2.0);
            
            // Map 2D (i, j) to 1D index
            int k = grid.idx(i, j);
            V[k] = func(x, y);
        }
    }
    built = false; // Flag that matrix relies on old potential
}

void Schrodinger::build_operator() {
    if (built) return;

    int N = grid.N_x;
    int M = grid.N_y;
    int size = N * M;
    double inv_h2 = 1.0 / (grid.h * grid.h);

    // "Triplet" list is the standard way to fill Sparse Matrices in Eigen
    std::vector<Eigen::Triplet<double>> coefficients;
    coefficients.reserve(size * 5); // 5-point stencil

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            int k = grid.idx(i, j);

            // 1. Diagonal: 4/h^2 + V_{ij}
            // The Laplacian contributes 4/h^2, the Potential adds V
            double diag_val = 4.0 * inv_h2 + V[k];
            coefficients.push_back(Eigen::Triplet<double>(k, k, diag_val));

            // 2. Off-Diagonals (Neighbors): -1/h^2
            // We check bounds to apply Dirichlet Boundary Conditions (u=0 at edge)
            
            // Neighbor i+1
            if (i + 1 < N) 
                coefficients.push_back(Eigen::Triplet<double>(k, grid.idx(i + 1, j), -inv_h2));
            // Neighbor i-1
            if (i - 1 >= 0) 
                coefficients.push_back(Eigen::Triplet<double>(k, grid.idx(i - 1, j), -inv_h2));
            // Neighbor j+1
            if (j + 1 < M) 
                coefficients.push_back(Eigen::Triplet<double>(k, grid.idx(i, j + 1), -inv_h2));
            // Neighbor j-1
            if (j - 1 >= 0) 
                coefficients.push_back(Eigen::Triplet<double>(k, grid.idx(i, j - 1), -inv_h2));
        }
    }

    hamiltonian.resize(size, size);
    hamiltonian.setFromTriplets(coefficients.begin(), coefficients.end());
    built = true;
}

std::vector<double> Schrodinger::compute_spectrum(int num_eigenvalues) {
    build_operator();

    // NOTE: For a real 256x256 grid, this matrix is 65k x 65k.
    // A dense solver (SelfAdjointEigenSolver) will be slow/memory heavy.
    // In a production app, use SPECTRA (Simulated here via Eigen dense for compilation simplicity).
    
    // For this example, we proceed assuming grid is small enough or User uses Spectra.
    // Here is the "Dense" fallback for demonstration:
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(hamiltonian);
    
    eigenvalues = solver.eigenvalues().head(num_eigenvalues);
    eigenfunctions = solver.eigenvectors().leftCols(num_eigenvalues);

    // Convert to std::vector to return to Lua
    std::vector<double> evals(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());
    return evals;
}

double Schrodinger::analyze_doubling_index(double x0, double y0, double r) {
    // We analyze the ground state (first column of eigenfunctions)
    if (eigenfunctions.cols() == 0) return 0.0;
    
    Eigen::VectorXd u = eigenfunctions.col(0); // Ground state

    double integral_r = 0.0;
    double integral_2r = 0.0;
    double r_sq = r * r;
    double two_r_sq = (2 * r) * (2 * r);

    // Numerical Integration (Riemann Sum)
    for (int i = 0; i < grid.N_x; ++i) {
        for (int j = 0; j < grid.N_y; ++j) {
            double x = i * grid.h - (grid.length / 2.0);
            double y = j * grid.h - (grid.length / 2.0);
            
            double dist_sq = (x - x0)*(x - x0) + (y - y0)*(y - y0);
            
            int k = grid.idx(i, j);
            double u_val = u[k];
            double u_sq = u_val * u_val;

            if (dist_sq <= r_sq) {
                integral_r += u_sq;
            }
            if (dist_sq <= two_r_sq) {
                integral_2r += u_sq;
            }
        }
    }

    // Doubling Index formula: log_2( int(B_2r) / int(B_r) )
    if (integral_r < 1e-9) return 0.0; // Avoid division by zero
    return std::log2(integral_2r / integral_r);
}