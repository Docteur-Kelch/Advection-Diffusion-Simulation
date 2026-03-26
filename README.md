# 2D Advection-Diffusion Simulation Project

## Overview
This project provides a comprehensive C++ implementation for solving the **2D Advection-Diffusion equation** on a unit square domain $[0,1] \times [0,1]$.

The project explores two major numerical methods:
1.  **Finite Volume Method (FVM):** Focusing on conservation laws, implementing both Classical (Central/Upwind) and Scharfetter-Gummel schemes.
2.  **Finite Element Method (FEM):** Using bilinear ($Q_1$) elements for both stationary and unsteady problems, leveraging the **Eigen3** library for linear algebra.

This codebase was developed for educational purposes to analyze convergence rates ($L^2$ and $H^1$ errors) and numerical stability in diffusion-dominated vs. advection-dominated regimes.

---

## Prerequisites

Before compiling, ensure you have the following installed:
* **C++ Compiler:** GCC (`g++`) supporting C++11 or higher.
* **Eigen3 Library:** A header-only C++ library for linear algebra (required for the FEM part).
    * *Note:* The instructions below assume Eigen is located at `/mingw64/include/eigen3`. Please adjust the include path (`-I`) to match your local installation.

---

## 1. Finite Volume Part (VF)

The Finite Volume implementation is contained in a single source file (`VF_scheme.cpp`). Different exercises and solvers are activated using **preprocessor flags** (`-D`) during compilation.

**Source File:** `VF_scheme.cpp`

### Available Simulations

#### A. Classical Scheme Convergence (Q7)
Solves the transient problem using a standard approach:
* **Diffusion:** Central Finite Difference approximation.
* **Advection:** Upwind scheme (stable for transport).
* **Goal:** Validates the order of convergence in the $L^2$ norm.

* **Compilation:**
    ```bash
    g++ -o exec_q7 VF_scheme.cpp -DRUN_Q7
    ```
* **Execution:**
    ```bash
    ./exec_q7
    ```
* **Output:** `results_q7.csv`

#### B. Scharfetter-Gummel Scheme (Q12)
Solves the problem using the Scharfetter-Gummel flux approximation. This scheme is particularly robust for high Péclet numbers (advection-dominated flows) as it respects the exponential nature of the exact solution.

* **Compilation:**
    ```bash
    g++ -o exec_q12 VF_scheme.cpp -DRUN_Q12
    ```
* **Execution:**
    ```bash
    ./exec_q12
    ```
* **Output:** `results_SG.csv`

#### C. Steady-State Convergence (Q13)
Simulates the transient evolution over a long time period ($T=5.0$) to verify convergence towards the theoretical stationary solution. It compares the numerical steady state against the exact stationary solution.

* **Compilation:**
    ```bash
    g++ -o exec_q13 VF_scheme.cpp -DRUN_Q13
    ```
* **Execution:**
    ```bash
    ./exec_q13
    ```
* **Output:**
    * `results_Q13_summary.csv` (Convergence summary)
    * `evolution_k4.csv` (Time evolution of error for mesh $k=4$)

---

## 2. Finite Element Part (FEM)

The Finite Element implementation uses the **Eigen3** library to assemble and solve sparse linear systems. The code is split into stationary and unsteady solvers.

### A. Stationary Simulation (Q25)
Solves the steady-state equation:
$$ - \Delta u + \mathbf{V} \cdot \nabla u = 0 $$
It computes both $L^2$ and $H^1$ errors for different mesh refinement levels ($k=2$ to $7$) and different advection intensities ($\eta = 1, 10, 100$).

**Source File:** `FE_scheme.cpp`

* **Compilation:**
    ```bash
    g++ -I /mingw64/include/eigen3 FE_scheme.cpp -o FE_scheme
    ```
* **Execution:**
    ```bash
    ./FE_scheme
    ```
* **Output:** Generates `Convergence_Eta_1.csv`, `Convergence_Eta_10.csv`, etc.

### B. Unsteady Simulation (Q26 / Bonus)
Solves the time-dependent problem using an implicit **Backward Euler** time-stepping scheme.
$$ \frac{\partial u}{\partial t} - \Delta u + \mathbf{V} \cdot \nabla u = 0 $$
This simulation monitors the error against the steady-state solution and the exact transient solution over time.

**Source File:** `FE_scheme_Bonus.cpp`

* **Compilation:**
    ```bash
    g++ -I /mingw64/include/eigen3 FE_scheme_Bonus.cpp -o FE_scheme_B
    ```
* **Execution:**
    ```bash
    ./FE_scheme_B
    ```
* **Output:** `Resultats_Q26.csv`

---

## Project Structure

| File | Description |
| :--- | :--- |
| `VF_scheme.cpp` | Main entry point for all Finite Volume exercises (Q7, Q12, Q13). |
| `FE_scheme.cpp` | Main entry point for FEM Stationary solver (Q25). |
| `FE_scheme_Bonus.cpp` | Main entry point for FEM Unsteady solver (Q26). |
| `mesh_file.h` | Header-only library defining mesh structures (`Cell`, `Edge`, `Vertex`) and Cartesian mesh generation. |
| `global_function.h` | Helper functions for exact solutions and common math utilities. |
| `*.csv` | Generated results files containing error norms and convergence orders. |

---

## Cleaning Up

To remove all generated executables and clean the directory, run:

**On Linux / Mac / Git Bash:**
```bash
rm -f exec_q7 exec_q12 exec_q13 FE_scheme FE_scheme_B