#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm> // std::generate
#include <Eigen/Sparse>
#include <Eigen/Dense>

// Using EIGEN_PI for high precision
const double PI = EIGEN_PI;

int main() {
    // 1. Setup N values (Powers of 2: 8 to 4096)
    // std::vector helps us manage the test cases dynamically
    std::vector<int> N_vec(10); 
    int power = 3;
    
    // std::generate: Fills the vector by repeatedly calling the lambda.
    // Concise Use: Replaces manual loops for initializing algorithmic sequences.
    std::generate(N_vec.begin(), N_vec.end(), [&power]() {
        return 1 << power++; // Bitwise shift for powers of 2 (2^3, 2^4...)
    });

    std::cout << "N       Error" << std::endl;
    std::cout << "-------------" << std::endl;

    for (int N : N_vec) {
        double h = 2.0 * PI / N;

        // 2. Grid Setup
        // We replicate Python's: x = -pi + arange(1, N+1)*h
        Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N, -PI + h, PI);

        // 3. Functions u(x) and u'(x)
        // Eigen's .array() allows vectorized math (like numpy)
        Eigen::VectorXd u = x.array().sin().exp(); 
        Eigen::VectorXd u_prime = x.array().cos() * u.array();

        // 4. Construct Sparse Matrix D (4th order centered difference)
        // Stencil: 1/h * [1/12, -2/3, 0, 2/3, -1/12]
        // We use a Triplets vector to efficiently build the sparse matrix.
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(N * 4); // Reserve memory to prevent reallocations

        for (int i = 0; i < N; ++i) {
            // Periodic indices handling wraparound (modulo arithmetic)
            int i_minus_2 = (i - 2 + N) % N;
            int i_minus_1 = (i - 1 + N) % N;
            int i_plus_1  = (i + 1) % N;
            int i_plus_2  = (i + 2) % N;

            // Direct 4th-order coefficients
            // Matches (D - D.T)/h logic: 2/3 at (+1), -1/12 at (+2)
            triplets.push_back({i, i_minus_2,  1.0/12.0});
            triplets.push_back({i, i_minus_1, -2.0/3.0});
            triplets.push_back({i, i_plus_1,   2.0/3.0});
            triplets.push_back({i, i_plus_2,  -1.0/12.0});
        }

        Eigen::SparseMatrix<double> D(N, N);
        D.setFromTriplets(triplets.begin(), triplets.end());
        D *= (1.0 / h); // Scale by 1/h

        // 5. Compute Error
        Eigen::VectorXd Du = D * u;
        
        // standard L_infinity norm
        double error = (Du - u_prime).lpNorm<Eigen::Infinity>();

        std::cout << std::setw(6) << N << "  " 
                  << std::scientific << std::setprecision(4) << error << std::endl;
    }

    return 0;
}