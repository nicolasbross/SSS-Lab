#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm> // std::transform
#include <Eigen/Dense>

using namespace Eigen;
const double PI = EIGEN_PI;

// Helper: Normalized Sinc function sin(z)/z (handling z=0)
struct Sinc {
    double operator()(double z) const {
        // Use a small epsilon for stability near zero
        if (std::abs(z) < 1e-10) return 1.0;
        return std::sin(z) / z;
    }
};

int main() {
    double h = 1.0;
    double xmax = 10.0;
    
    // 1. Setup Grids
    // x: Computational grid [-10, 10] step 1
    int N_x = static_cast<int>(2 * xmax / h) + 1; // +1 to include right edge if needed
    VectorXd x = VectorXd::LinSpaced(N_x, -xmax, xmax); // -10, -9 ... 10

    // xx: Plotting grid [-10.05, 10.05] step 0.1
    // arange(-xmax-h/20, xmax+h/20, h/10)
    double start_xx = -xmax - h/20.0;
    double end_xx = xmax + h/20.0;
    double step_xx = h / 10.0;
    int N_xx = static_cast<int>((end_xx - start_xx) / step_xx); 
    // Note: LinSpaced is inclusive, arange is exclusive. We approximate count.
    VectorXd xx = VectorXd::LinSpaced(202, -10.05, 10.05); 

    // Prepare output file
    std::ofstream file("p3_data.csv");
    file << "xx,p_delta,p_square,p_hat\n";

    // Matrix to store results for 3 cases (rows=xx size, cols=3)
    MatrixXd results(xx.size(), 3);

    // 2. Process 3 cases: Delta, Square, Hat
    for (int pl = 0; pl < 3; ++pl) {
        VectorXd v(x.size());

        // Generate v based on conditions (using std::transform semantics via Eigen unaryExpr)
        if (pl == 0) { // Delta function
            v = x.unaryExpr([](double val) { return (val == 0.0) ? 1.0 : 0.0; });
        } 
        else if (pl == 1) { // Square wave
            v = x.unaryExpr([](double val) { return (std::abs(val) <= 3.0) ? 1.0 : 0.0; });
        } 
        else if (pl == 2) { // Hat function
            v = x.unaryExpr([](double val) { return std::max(0.0, 1.0 - std::abs(val) / 3.0); });
        }

        // 3. Band-limited Interpolation
        // p(xx) = sum( v[i] * sinc( pi*(xx - x[i])/h ) )
        VectorXd p = VectorXd::Zero(xx.size());

        for (int i = 0; i < x.size(); ++i) {
            double xi = x[i];
            double vi = v[i];
            
            if (vi == 0.0) continue; // Optimization: skip zero weights

            // Vectorized Sinc application
            // Argument to sinc: pi * (xx - xi) / h
            // We use unaryExpr with our Sinc functor
            VectorXd arg = (PI * (xx.array() - xi) / h);
            p += vi * arg.unaryExpr(Sinc());
        }
        
        results.col(pl) = p;
    }

    // 4. Write to CSV
    for (int i = 0; i < xx.size(); ++i) {
        file << xx[i] << "," 
             << results(i, 0) << "," 
             << results(i, 1) << "," 
             << results(i, 2) << "\n";
    }

    std::cout << "Data generated to 'p3_data.csv'.\n";
    return 0;
}