/*==========================================================
   UNSTEADY FINITE ELEMENT SCHEME (BACKWARD EULER)
   Project: Mathematical Tools for Simulation
===========================================================*/

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "mesh_file.h"
#include "global_function.h"

// Eigen shortcuts
using SparseMat = Eigen::SparseMatrix<double>;
using Tripletd = Eigen::Triplet<double>;
using VectorXd = Eigen::VectorXd;
using Vector4d = Eigen::Vector4d;
using Matrix4d = Eigen::Matrix4d;

/*==========================================================
  HELPER FUNCTIONS
===========================================================*/

// Q1 shape functions evaluation on reference element
double eval_uh_local(const Vector4d& u_el, double xi, double nu) {
    double phi0 = (1.0 - xi) * (1.0 - nu);
    double phi1 = xi * (1.0 - nu);
    double phi2 = (1.0 - xi) * nu;
    double phi3 = xi * nu;
    return u_el[0] * phi0 + u_el[1] * phi1 + u_el[2] * phi2 + u_el[3] * phi3;
}

// Element matrices for Stiffness (K) and Advection (C)
void compute_element_matrices(double h, double eta, Eigen::Matrix4d& K_el, Eigen::Matrix4d& C_el) {
    // Diffusion 
    K_el << 4, -1, -1, -2,
        -1, 4, -2, -1,
        -1, -2, 4, -1,
        -2, -1, -1, 4;
    K_el /= 6.0;

    // Advection
    C_el << -2, 2, -1, 1,
        -2, 2, -1, 1,
        -1, 1, -2, 2,
        -1, 1, -2, 2;
    C_el *= (eta * h / 12.0);
}

// Element Mass Matrix (M) for Q1 elements
void compute_mass_matrix_elem(double h, Matrix4d& M_el) {
    M_el << 4, 2, 2, 1,
        2, 4, 1, 2,
        2, 1, 4, 2,
        1, 2, 2, 4;
    M_el *= (h * h / 36.0);
}

// Compute L2 Norm of the error between two vectors U and V
double compute_L2_error(const VectorXd& U, const VectorXd& V, const Mesh& mesh) {
    double h = mesh.h;
    int M = mesh.M;
    double err_sq = 0.0;

    // Gauss Quadrature 2x2
    double g_pts[2] = { 0.5 - 0.5 / std::sqrt(3.0), 0.5 + 0.5 / std::sqrt(3.0) };
    double g_w[2] = { 0.5, 0.5 };

    for (const auto& cell : mesh.cells) {
        int id = cell.id;
        int i = id % M;
        int j = id / M;

        int idx_BG = j * (M + 1) + i;
        int idx_BD = idx_BG + 1;
        int idx_HG = (j + 1) * (M + 1) + i;
        int idx_HD = idx_HG + 1;

        Vector4d u_loc, v_loc;
        u_loc << U[idx_BG], U[idx_BD], U[idx_HG], U[idx_HD];
        v_loc << V[idx_BG], V[idx_BD], V[idx_HG], V[idx_HD];

        for (int qx = 0; qx < 2; ++qx) {
            for (int qy = 0; qy < 2; ++qy) {
                double w = g_w[qx] * g_w[qy];
                double diff = eval_uh_local(u_loc - v_loc, g_pts[qx], g_pts[qy]);
                err_sq += w * diff * diff * h * h;
            }
        }
    }
    return std::sqrt(err_sq);
}

/*==========================================================
   MAIN SIMULATION
===========================================================*/
int main() {
    // Simulation Parameters
    double eta = 1.0;
    int k = 5; // M = 32
    double T_final = 5.0;
    double dt = 0.01;

    std::cout << "--- Unsteady Simulation (Backward Euler) ---\n";
    std::cout << "Eta = " << eta << ", k = " << k << ", dt = " << dt << "\n";

    // 1. Mesh Construction
    Mesh mesh;
    mesh.build_Cartesian(k);
    int M_grid = mesh.M;
    double h = mesh.h;
    int N_dof = mesh.vertices.size();

    std::cout << "DoFs: " << N_dof << " (h=" << h << ")\n";

    // 2. Global Assembly (Mass M, Stiffness K, Advection C)
    std::vector<Tripletd> trip_M, trip_A;
    Matrix4d M_el, K_el, C_el;

    compute_mass_matrix_elem(h, M_el);
    compute_element_matrices(h, eta, K_el, C_el);

    Matrix4d A_el = K_el + C_el; // Full spatial operator

    for (const auto& cell : mesh.cells) {
        int id = cell.id;
        int i = id % M_grid;
        int j = id / M_grid;

        int idxs[4] = {
            j * (M_grid + 1) + i,       // Bottom-Left
            j * (M_grid + 1) + i + 1,   // Bottom-Right
            (j + 1) * (M_grid + 1) + i, // Top-Left
            (j + 1) * (M_grid + 1) + i + 1 // Top-Right
        };

        for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 4; ++c) {
                trip_M.emplace_back(idxs[r], idxs[c], M_el(r, c));
                trip_A.emplace_back(idxs[r], idxs[c], A_el(r, c));
            }
        }
    }

    SparseMat MatM(N_dof, N_dof), MatA(N_dof, N_dof);
    MatM.setFromTriplets(trip_M.begin(), trip_M.end());
    MatA.setFromTriplets(trip_A.begin(), trip_A.end());

    // 3. Compute Discrete Steady State (Reference for E_stab)
    // Solve A * U = F with penalty for Dirichlet BCs
    SparseMat SysStat = MatA;
    VectorXd F_stat = VectorXd::Zero(N_dof);
    double penalty = 1e14;

    for (int idx = 0; idx < N_dof; ++idx) {
        double x = mesh.vertices[idx].x.x;
        if (std::abs(x) < 1e-12) {
            SysStat.coeffRef(idx, idx) += penalty;
            F_stat(idx) = penalty * 1.0;
        }
        else if (std::abs(x - 1.0) < 1e-12) {
            SysStat.coeffRef(idx, idx) += penalty;
            F_stat(idx) = penalty * std::exp(eta);
        }
    }

    Eigen::SparseLU<SparseMat> solverStat;
    solverStat.compute(SysStat);
    VectorXd U_steady = solverStat.solve(F_stat);
    std::cout << "Steady state computed.\n";

    // 4. Initialization for Unsteady Problem
    // Interpolated exact solution (Reference for E_tot)
    VectorXd U_exact_interp(N_dof);
    for (int i = 0; i < N_dof; ++i) {
        U_exact_interp[i] = exact_solution_stationary(mesh.vertices[i].x.x, eta);
    }

    // Initial Condition u0
    VectorXd U_n(N_dof);
    for (int i = 0; i < N_dof; ++i) {
        double x = mesh.vertices[i].x.x;
        // u0 from Eq (2) at t=0
        U_n[i] = std::exp(eta * x) + std::exp(0.5 * eta * x) * std::sin(M_PI * x);
    }

    // Setup Linear System for Backward Euler
    // (1/dt * M + A) U^{n} = RHS
    SparseMat LHS = MatM * (1.0 / dt) + MatA;

    // Apply penalty to LHS (Constant matrix)
    for (int idx = 0; idx < N_dof; ++idx) {
        double x = mesh.vertices[idx].x.x;
        if (std::abs(x) < 1e-12 || std::abs(x - 1.0) < 1e-12) {
            LHS.coeffRef(idx, idx) += penalty;
        }
    }

    Eigen::SparseLU<SparseMat> solverUnsteady;
    solverUnsteady.compute(LHS);

    // Output file
    std::ofstream file("Resultats_Q26.csv");
    file << "n,time,Err_vs_Steady,Err_vs_Exact\n";

    double t = 0.0;
    int n_step = 0;
    VectorXd U_prev = U_n;

    std::cout << "Starting time loop...\n";

    // 5. Time Loop
    while (t < T_final) {
        n_step++;
        t += dt;

        // Construct RHS: (1/dt) * M * U^{n-1}
        VectorXd rhs = (MatM * U_prev) * (1.0 / dt);

        // Apply Dirichlet BCs to RHS
        for (int idx = 0; idx < N_dof; ++idx) {
            double x = mesh.vertices[idx].x.x;
            if (std::abs(x) < 1e-12) {
                rhs(idx) = penalty * 1.0; // u(0,y) = 1
            }
            else if (std::abs(x - 1.0) < 1e-12) {
                rhs(idx) = penalty * std::exp(eta); // u(1,y) = exp(eta)
            }
        }

        // Solve
        U_n = solverUnsteady.solve(rhs);

        // Compute Errors
        double err_steady = compute_L2_error(U_n, U_steady, mesh);
        double err_exact = compute_L2_error(U_n, U_exact_interp, mesh);

        // Save
        file << n_step << "," << t << "," << err_steady << "," << err_exact << "\n";

        if (n_step % 10 == 0) {
            std::cout << "Iter " << n_step << " t=" << t
                << " | E_Stab=" << err_steady
                << " | E_Tot=" << err_exact << "\n";
        }
        U_prev = U_n;
    }

    std::cout << "Done. Results saved to Resultats_Q26.csv\n";
    return 0;
}