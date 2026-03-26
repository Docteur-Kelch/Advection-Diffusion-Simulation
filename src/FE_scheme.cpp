/*==========================================================
           FINITE ELEMENT SCHEME APPROXIMATION
   Project: Advection-Diffusion 2D (Stationary)
===========================================================*/

#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "mesh_file.h"
#include "global_function.h"

// Eigen shortcuts for readability
using SparseMat = Eigen::SparseMatrix<double>;
using Tripletd = Eigen::Triplet<double>;
using Vector4d = Eigen::Vector4d;

// Struct to store simulation results for one mesh level
struct SimResult {
    double h;
    double err_L2;
    double err_H1;
};

/*==========================================================
  Exact gradient of the stationary solution u = exp(eta x)
  Used for H1 error computation.
===========================================================*/
Eigen::Vector2d grad_exact_stationary(double x, double eta) {
    return Eigen::Vector2d(eta * std::exp(eta * x), 0.0);
}

/*==========================================================
  Q1 shape functions evaluation on reference cell [-1,1]^2
  Interpolates the solution u_h at local coordinates (xi, nu)
===========================================================*/
double eval_uh_local(const Vector4d& u_el, double xi, double nu) {
    // Bilinear basis functions
    double phi0 = (1.0 - xi) * (1.0 - nu);
    double phi1 = xi * (1.0 - nu);
    double phi2 = (1.0 - xi) * nu;
    double phi3 = xi * nu;

    return u_el[0] * phi0 + u_el[1] * phi1 + u_el[2] * phi2 + u_el[3] * phi3;
}

/*==========================================================
  Computes the gradient of u_h mapped to the physical cell.
  Includes the Jacobian inverse (1/h factor).
===========================================================*/
Eigen::Vector2d eval_grad_uh_local(const Vector4d& u_el, double xi, double nu, double h) {
    double dphi_dxi[4] = { -(1.0 - nu), (1.0 - nu), -nu, nu };
    double dphi_dnu[4] = { -(1.0 - xi), -xi, (1.0 - xi), xi };

    double du_dxi = 0.0;
    double du_dnu = 0.0;

    for (int i = 0; i < 4; ++i) {
        du_dxi += u_el[i] * dphi_dxi[i];
        du_dnu += u_el[i] * dphi_dnu[i];
    }
    // Transform gradients from reference to physical space
    return Eigen::Vector2d(du_dxi / h, du_dnu / h);
}

/*==========================================================
  Builds Element Matrices for Stiffness (K) and Advection (C)
  Derived using tensor products of 1D matrices.
===========================================================*/
void compute_element_matrices(double h, double eta, Eigen::Matrix4d& K_el, Eigen::Matrix4d& C_el) {
    // Diffusion Matrix (Stiffness) - Independent of h in 2D
    K_el << 4, -1, -1, -2,
        -1, 4, -2, -1,
        -1, -2, 4, -1,
        -2, -1, -1, 4;
    K_el /= 6.0;

    // Advection Matrix (Transport) - Proportional to h and eta
    C_el << -2, 2, -1, 1,
        -2, 2, -1, 1,
        -1, 1, -2, 2,
        -1, 1, -2, 2;
    C_el *= (eta * h / 12.0);
}

/*==========================================================
  Main FEM solver function.
  1. Meshes the domain.
  2. Assembles the system.
  3. Solves using SparseLU.
  4. Computes L2/H1 errors.
===========================================================*/
SimResult solve_fem_stationary(int k, double eta) {
    // 1. Mesh Generation
    Mesh mesh;
    mesh.build_Cartesian(k); // M = 2^k
    int M = mesh.M;
    double h = mesh.h;
    int N_dof = mesh.vertices.size();

    // 2. System Assembly
    Eigen::Matrix4d K_el, C_el;
    compute_element_matrices(h, eta, K_el, C_el);
    Eigen::Matrix4d A_el = K_el + C_el; // Local operator

    std::vector<Tripletd> triplets;
    Eigen::VectorXd F = Eigen::VectorXd::Zero(N_dof);
    triplets.reserve(N_dof * 9); // Reserve memory for efficiency

    // Loop over all cells to fill the triplets
    for (const auto& cell : mesh.cells) {
        int id = cell.id;
        int i = id % M;
        int j = id / M;

        // Global indices of the 4 nodes of the current cell
        int idx_BG = j * (M + 1) + i;
        int idx_BD = idx_BG + 1;
        int idx_HG = (j + 1) * (M + 1) + i;
        int idx_HD = idx_HG + 1;
        int nodes[4] = { idx_BG, idx_BD, idx_HG, idx_HD };

        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c)
                triplets.emplace_back(nodes[r], nodes[c], A_el(r, c));
    }

    SparseMat A(N_dof, N_dof);
    A.setFromTriplets(triplets.begin(), triplets.end());

    // 3. Boundary Conditions (Penalty Method)
    // Enforce Dirichlet BCs strongly by adding a large value on the diagonal
    double penalty = 1e14;
    for (int idx = 0; idx < N_dof; ++idx) {
        double x = mesh.vertices[idx].x.x;
        if (std::abs(x) < 1e-12) { // Left boundary (x=0) -> u=1
            A.coeffRef(idx, idx) += penalty;
            F(idx) = penalty * 1.0;
        }
        else if (std::abs(x - 1.0) < 1e-12) { // Right boundary (x=1) -> u=exp(eta)
            A.coeffRef(idx, idx) += penalty;
            F(idx) = penalty * exact_solution_stationary(1.0, eta);
        }
    }

    // 4. Linear Solver
    Eigen::SparseLU<SparseMat> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solver factorization failed!" << std::endl;
        return { h, -1.0, -1.0 };
    }
    Eigen::VectorXd U = solver.solve(F);

    // 5. Error Calculation using 2x2 Gauss Quadrature
    double err_L2_sq = 0.0;
    double err_H1_sq = 0.0;
    double g_pts[2] = { 0.5 - 0.5 / std::sqrt(3.0), 0.5 + 0.5 / std::sqrt(3.0) };
    double g_w[2] = { 0.5, 0.5 };

    for (const auto& cell : mesh.cells) {
        int id = cell.id;
        int i = id % M;
        int j = id / M;

        // Retrieve indices and local solution coefficients
        int idx_BG = j * (M + 1) + i;
        int idx_BD = idx_BG + 1;
        int idx_HG = (j + 1) * (M + 1) + i;
        int idx_HD = idx_HG + 1;

        Vector4d u_loc;
        u_loc << U[idx_BG], U[idx_BD], U[idx_HG], U[idx_HD];
        double x0 = i * h;

        // Loop over quadrature points
        for (int qx = 0; qx < 2; ++qx) {
            for (int qy = 0; qy < 2; ++qy) {
                double xi = g_pts[qx];
                double nu = g_pts[qy];
                double w = g_w[qx] * g_w[qy];

                double x_phy = x0 + xi * h; // Map to physical x coordinate

                // Exact values
                double u_ex = exact_solution_stationary(x_phy, eta);
                Eigen::Vector2d grad_ex = grad_exact_stationary(x_phy, eta);

                // Numerical values
                double u_h = eval_uh_local(u_loc, xi, nu);
                Eigen::Vector2d grad_h = eval_grad_uh_local(u_loc, xi, nu, h);

                // Accumulate weighted errors
                err_L2_sq += w * std::pow(u_h - u_ex, 2) * h * h;
                err_H1_sq += w * (grad_h - grad_ex).squaredNorm() * h * h;
            }
        }
    }

    return { h, std::sqrt(err_L2_sq), std::sqrt(err_H1_sq) };
}

/*==========================================================
  MAIN EXECUTION
===========================================================*/
int main() {
    // List of advection intensities to test
    std::vector<double> etas = { 1.0, 10.0, 100.0 };

    std::cout << std::fixed << std::setprecision(6);

    for (double eta : etas) {
        std::string filename = "Convergence_Eta_" + std::to_string((int)eta) + ".csv";
        std::ofstream file(filename);

        // CSV Header
        file << "k,h,Error_L2,Order_L2,Error_H1,Order_H1\n";

        std::cout << "\n========================================\n";
        std::cout << " CONVERGENCE STUDY FOR ETA = " << eta << "\n";
        std::cout << "========================================\n";
        std::cout << "k \t h \t\t Err_L2 \t Order \t Err_H1 \t Order\n";
        std::cout << "------------------------------------------------------------------\n";

        double prev_h = 0.0;
        double prev_L2 = 0.0;
        double prev_H1 = 0.0;

        // Loop over refinement levels k=2 to k=7
        for (int k = 2; k <= 7; ++k) {
            // Run simulation
            SimResult res = solve_fem_stationary(k, eta);

            // Compute experimental order of convergence (log ratio)
            double order_L2 = 0.0;
            double order_H1 = 0.0;

            if (k > 2) {
                order_L2 = std::log(res.err_L2 / prev_L2) / std::log(res.h / prev_h);
                order_H1 = std::log(res.err_H1 / prev_H1) / std::log(res.h / prev_h);
            }

            // Console output with color for orders
            std::cout << k << "\t " << res.h << "\t "
                << res.err_L2 << "\t ";
            if (k > 2) std::cout << "\033[1;32m" << order_L2 << "\033[0m";
            else std::cout << "-";

            std::cout << "\t " << res.err_H1 << "\t ";
            if (k > 2) std::cout << "\033[1;36m" << order_H1 << "\033[0m";
            else std::cout << "-";
            std::cout << "\n";

            // Save to CSV
            file << k << "," << res.h << ","
                << res.err_L2 << "," << (k > 2 ? std::to_string(order_L2) : "-") << ","
                << res.err_H1 << "," << (k > 2 ? std::to_string(order_H1) : "-") << "\n";

            // Update previous values for next iteration
            prev_h = res.h;
            prev_L2 = res.err_L2;
            prev_H1 = res.err_H1;
        }
        std::cout << "-> Results exported to : " << filename << "\n";
    }

    return 0;
}