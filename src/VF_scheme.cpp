#include <fstream>
#include <chrono>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "global_function.h"
#include "mesh_file.h"

// ============================================================================
//   SÉLECTEUR D'EXERCICE (Préprocesseur)
//   Permet de compiler un seul exercice à la fois sans modifier le code.
//   Utilisation via terminal : g++ ... -DRUN_Q7
// ============================================================================

// Décommenter une ligne ci-dessous si vous ne compilez pas via le terminal :
//#define RUN_Q7
//#define RUN_Q12
//#define RUN_Q13

// ============================================================================

/*-----------------------------------------------------------------------------
 * PARTIE 1 : SOLUTIONS EXACTES ET FONCTIONS UTILITAIRES
 *-----------------------------------------------------------------------------*/


// --- Fonction de Bernoulli (Flux Scharfetter-Gummel) ---
// B(z) = z / (e^z - 1). Gère la singularité en z=0 pour éviter la division par zéro.
double B_Bernoulli(double z) {
    if (std::abs(z) < 1e-12) {
        return 1.0; // Limite en 0
    }
    return z / (std::exp(z) - 1.0);
}

/*-----------------------------------------------------------------------------
 * PARTIE 2 : SOLVEURS VOLUMES FINIS
 *-----------------------------------------------------------------------------*/

 // ------ SOLVEUR CLASSIQUE (Q5 & Q7) ------
 /* Implémente un schéma centré pour la diffusion et Upwind pour l'advection.
  * Ce schéma est conditionnellement stable (CFL diffusive et convective). */
void update_solution(Mesh& mesh, const std::vector<double>& u_n,
    std::vector<double>& u_next, double dt, double eta) {

    // Boucle sur toutes les cellules du maillage
    for (size_t i = 0; i < mesh.cells.size(); ++i) {
        Cell& K = mesh.cells[i];
        double flux_total = 0.0;
        double u_center = u_n[i];

        // Somme des flux sur les 4 faces de la cellule
        for (int k = 0; k < 4; ++k) {
            Edge& s = *(K.E[k]);
            double orient = K.w[k]; // +1 sortant, -1 entrant

            // Calcul de la vitesse normale V.n
            // V = (eta, 0), donc V.n dépend de nx uniquement
            double nx = s.n_sigma.x * orient;
            double Vn = eta * nx;

            // --- Gestion des Conditions aux Limites (BC) ---

            // 1. Neumann Homogène (bnd == 2) -> Flux nul (ex: parois haut/bas)
            if (s.bnd == 2) continue;

            double u_other;      // Valeur "voisine" (u_j ou bord)
            double d_effective;  // Distance pour le gradient (h ou h/2)

            if (s.bnd == 0) { // Face Intérieure
                Cell* neigh = (s.T[0]->id == K.id) ? s.T[1] : s.T[0];
                u_other = u_n[neigh->id];
                d_effective = s.d_sigma; // = h
            }
            else { // Face Dirichlet (bnd == 1)
                // Conditions imposées sur les bords verticaux
                if (s.x_sigma.x < 1e-10) {
                    u_other = 1.0;          // x=0 -> u=1
                }
                else {
                    u_other = std::exp(eta); // x=1 -> u=exp(eta)
                }
                d_effective = s.d_sigma; // = h/2
            }

            // --- Calcul des Flux Numériques ---

            // A. Flux Diffusif (Centré) : -|sigma| * (u_voisin - u_centre) / d
            double flux_diff = -s.length * (u_other - u_center) / d_effective;

            // B. Flux Advectif (Upwind) : Décentrage amont selon le signe de Vn
            // Si Vn > 0 (flux sortant), on prend u_center. Sinon u_other.
            double flux_adv;
            if (Vn >= 0.0) {
                flux_adv = s.length * Vn * u_center;
            }
            else {
                flux_adv = s.length * Vn * u_other;
            }

            flux_total += flux_diff + flux_adv;
        }

        // Mise à jour explicite : u^{n+1} = u^n - (dt / Aire) * Somme(Flux)
        u_next[i] = u_center - (dt / K.area) * flux_total;
    }
}

// ------ SOLVEUR SCHARFETTER-GUMMEL (Q9-12) ------
/* Implémente le flux SG qui est exact pour les solutions stationnaires 1D.
 * Plus robuste pour les régimes dominés par l'advection (Péclet élevé). */
void update_solution_SG(Mesh& mesh, const std::vector<double>& u_n,
    std::vector<double>& u_next, double dt, double eta) {

    for (size_t i = 0; i < mesh.cells.size(); ++i) {
        Cell& K = mesh.cells[i];
        double flux_total = 0.0;

        for (int k = 0; k < 4; ++k) {
            Edge& s = *(K.E[k]);
            double orient = K.w[k];
            double Vn = eta * (s.n_sigma.x * orient);

            if (s.bnd == 2) continue; // Neumann

            double u_other;
            if (s.bnd == 0) {
                Cell* neigh = (s.T[0]->id == K.id) ? s.T[1] : s.T[0];
                u_other = u_n[neigh->id];
            }
            else { // Dirichlet
                u_other = (s.x_sigma.x < 0.5) ? 1.0 : std::exp(eta);
            }

            // --- Calcul du Flux SG ---
            // Formule : J = (1/d) * [ B(-z)u_i - B(z)u_j ]
            // où z est le nombre de Péclet local de l'arête
            double z = Vn * s.d_sigma;

            // Note: B_Bernoulli(-z) correspond au coefficient diffusif modifié + upwind implicite
            double flux = (s.length / s.d_sigma) * (B_Bernoulli(-z) * u_n[i] - B_Bernoulli(z) * u_other);

            flux_total += flux;
        }
        u_next[i] = u_n[i] - (dt / K.area) * flux_total;
    }
}

/*-----------------------------------------------------------------------------
 * PARTIE 3 : PROGRAMMES PRINCIPAUX (MAIN)
 *-----------------------------------------------------------------------------*/

#ifdef RUN_Q7
 // ======================= MAIN (Q7) - Étude de Convergence (Schéma Classique) =======================
int main() {
    std::cout << ">>> Lancement de l'exercice Q7 (Convergence Classique) <<<" << std::endl;
    double T_final = 0.5;
    double eta = 1.0;
    // Liste des raffinements : k=1 (grossier) -> k=6 (fin)
    std::vector<int> k_values = { 1, 2, 3, 4, 5, 6 };
    double prev_err = 0.0, prev_h = 0.0;

    std::ofstream csv("results_q7.csv");
    csv << "k,M,h,Error_L2,Order,Time(s)\n";

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "k\tM\th\t\tError L2\t\tOrder\t\tTime(s)\n";
    std::cout << "--------------------------------------------------------\n";

    for (int k : k_values) {
        Mesh mesh;
        mesh.build_Cartesian(k); // Construction du maillage 2^k * 2^k
        int M = mesh.M;
        double h = mesh.h;

        std::vector<double> u_n(mesh.cells.size()), u_next(mesh.cells.size());

        // Initialisation (t=0)
        for (size_t i = 0; i < u_n.size(); ++i)
            u_n[i] = exact_solution(mesh.cells[i].x_ij.x, 0.0, eta);

        // Condition CFL (Q6) : Pas de temps adaptatif
        double dt_diff = h * h / 4.0;      // Stabilité diffusive
        double dt_adv = h / std::abs(eta); // Stabilité convective
        double dt = 0.4 * std::min(dt_diff, dt_adv); // Coeff sécurité 0.4

        double t = 0.0;

        // --- Boucle Temporelle ---
        auto start_time = std::chrono::high_resolution_clock::now();
        while (t < T_final) {
            if (t + dt > T_final) dt = T_final - t; // Ajustement dernier pas
            update_solution(mesh, u_n, u_next, dt, eta);
            std::swap(u_n, u_next); // Échange rapide des vecteurs
            t += dt;
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        double duration_sec = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / 1000.0;

        // Calcul de l'erreur L2
        double err = 0.0;
        for (size_t i = 0; i < mesh.cells.size(); ++i) {
            double diff = u_n[i] - exact_solution(mesh.cells[i].x_ij.x, T_final, eta);
            err += mesh.cells[i].area * diff * diff;
        }
        err = std::sqrt(err);

        // Calcul de l'ordre de convergence expérimental
        double order = (prev_err > 0.0) ? std::log(err / prev_err) / std::log(h / prev_h) : 0.0;

        // Affichage console
        std::cout << k << "\t" << M << "\t" << h << "\t" << err << "\t"
            << std::fixed << std::setprecision(2) << order
            << std::scientific << "\t" << duration_sec << std::endl;

        // Sauvegarde CSV
        csv << k << "," << M << "," << h << "," << err << "," << order << "," << duration_sec << "\n";

        prev_err = err; prev_h = h;
    }
    csv.close();
    std::cout << "\nRésultats sauvegardés dans 'results_q7.csv'\n";
    return 0;
}
#endif

#ifdef RUN_Q12
// ======================= MAIN (Q12) - Étude Scharfetter-Gummel (SG) =======================
int main() {
    std::cout << ">>> Lancement de l'exercice Q12 (Convergence SG) <<<" << std::endl;
    double T_final = 0.5;
    double eta = 1.0;
    std::vector<int> k_values = { 1, 2, 3, 4, 5, 6 };
    double prev_err = 0.0, prev_h = 0.0;

    std::ofstream csv("results_SG.csv");
    csv << "k,M,h,Error_L2,Order,Time(s)\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "k\tM\th\t\tError L2\t\tOrder\tTime(s)\n";

    for (int k : k_values) {
        Mesh mesh;
        mesh.build_Cartesian(k);
        int M = mesh.M;
        double h = mesh.h;

        std::vector<double> u_n(mesh.cells.size()), u_next(mesh.cells.size());
        for (size_t i = 0; i < u_n.size(); ++i)
            u_n[i] = exact_solution(mesh.cells[i].x_ij.x, 0.0, eta);

        // Condition CFL améliorée pour SG (Q11)
        double dt = 0.4 * (h * h) / (4.0 + 2.0 * h * std::abs(eta));
        double t = 0.0;

        auto start_time = std::chrono::high_resolution_clock::now();
        while (t < T_final) {
            if (t + dt > T_final) dt = T_final - t;
            update_solution_SG(mesh, u_n, u_next, dt, eta); // Appel du solveur SG
            u_n = u_next; // Note: u_n = u_next est plus lent que swap, mais correct
            t += dt;
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        double duration_sec = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / 1000.0;

        double err = 0.0;
        for (size_t i = 0; i < mesh.cells.size(); ++i) {
            double diff = u_n[i] - exact_solution(mesh.cells[i].x_ij.x, T_final, eta);
            err += mesh.cells[i].area * diff * diff;
        }
        err = std::sqrt(err);
        double order = (prev_err > 0.0) ? std::log(err / prev_err) / std::log(h / prev_h) : 0.0;

        std::cout << k << "\t" << M << "\t" << h << "\t" << err << "\t"
            << std::fixed << std::setprecision(2) << order
            << std::scientific << "\t" << duration_sec << std::endl;

        csv << k << "," << M << "," << h << "," << err << "," << order << "," << duration_sec << "\n";
        prev_err = err; prev_h = h;
    }
    csv.close();
    std::cout << "\nRésultats sauvegardés dans 'results_SG.csv'\n";
    return 0;
}
#endif

#ifdef RUN_Q13
// ======================= MAIN (Q13) - Convergence vers l'État Stationnaire =======================
int main() {
    std::cout << ">>> Lancement de l'exercice Q13 (Stationnaire) <<<" << std::endl;
    double T_final = 5.0; // Temps long
    double eta = 1.0;
    std::vector<int> k_values = { 1, 2, 3, 4, 5, 6 };
    double prev_err = 0.0, prev_h = 0.0;

    std::ofstream csv_summary("results_Q13_summary.csv");
    csv_summary << "k,M,h,Error_L2,Order,n_inf,Time_s\n";

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "k\tM\th\t\tError_L2\tOrder\tn_inf\tTime(s)" << std::endl;

    for (int k : k_values) {
        Mesh mesh;
        mesh.build_Cartesian(k);
        int M = mesh.M;
        double h = mesh.h;

        std::vector<double> u_n(mesh.cells.size()), u_next(mesh.cells.size());

        // Initialisation avec la solution exacte instationnaire à t=0
        for (size_t i = 0; i < u_n.size(); ++i)
            u_n[i] = exact_solution(mesh.cells[i].x_ij.x, 0.0, eta);

        double dt = 0.4 * (h * h) / (4.0 + 2.0 * h * std::abs(eta));
        double t = 0.0;
        int n_inf = 0; // Compteur d'itérations pour atteindre le régime permanent

        std::ofstream evolution;
        if (k == 4) evolution.open("evolution_k4.csv");

        auto start_time = std::chrono::high_resolution_clock::now();

        // --- Boucle Temporelle avec Critère d'Arrêt ---
        while (t < T_final) {
            if (t + dt > T_final) dt = T_final - t;

            // CHOIX DU SOLVEUR (Décommenter celui désiré)
            update_solution(mesh, u_n, u_next, dt, eta);       // Schéma Classique
            // update_solution_SG(mesh, u_n, u_next, dt, eta); // Schéma SG (Optionnel pour Q13)

            // Calcul de la variation L2 entre deux pas de temps (Critère de Cauchy)
            double diff_sum = 0.0;
            for (size_t i = 0; i < u_n.size(); ++i) {
                double delta = u_next[i] - u_n[i];
                diff_sum += mesh.cells[i].area * (delta * delta);
            }
            double diff_norm_L2 = std::sqrt(diff_sum);

            // Sauvegarde de l'historique de convergence pour k=4
            if (k == 4 && n_inf % 10 == 0) {
                double err_to_inf = 0.0;
                for (size_t i = 0; i < u_n.size(); ++i) {
                    double d = u_next[i] - exact_solution_stationary(mesh.cells[i].x_ij.x, eta);
                    err_to_inf += mesh.cells[i].area * d * d;
                }
                evolution << n_inf << "," << std::sqrt(err_to_inf) << "," << diff_norm_L2 << "\n";
            }

            std::swap(u_n, u_next);
            t += dt;
            n_inf++;

            // Arrêt si la solution ne bouge plus (tolérance 1e-13)
            if (diff_norm_L2 < 1e-13) break;
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        double duration_sec = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / 1000.0;

        // --- Comparaison avec la Solution Stationnaire Théorique ---
        double err = 0.0;
        for (size_t i = 0; i < mesh.cells.size(); ++i) {
            double diff = u_n[i] - exact_solution_stationary(mesh.cells[i].x_ij.x, eta);
            err += mesh.cells[i].area * diff * diff;
        }
        err = std::sqrt(err);

        double order = (prev_err > 0.0) ? std::log(err / prev_err) / std::log(h / prev_h) : 0.0;

        std::cout << std::fixed << std::setprecision(0) << k << "\t" << M << "\t"
            << std::scientific << std::setprecision(2) << h << "\t"
            << err << "\t" << std::fixed << std::setprecision(2) << order << "\t"
            << n_inf << "\t" << duration_sec << std::endl;

        csv_summary << k << "," << M << "," << h << "," << err << "," << order << "," << n_inf << "," << duration_sec << "\n";

        prev_err = err; prev_h = h;
        if (evolution.is_open()) evolution.close();
    }
    csv_summary.close();
    std::cout << "\nRésultats sauvegardés dans 'results_Q13_summary.csv'\n";
    std::cout << "Évolution temporelle pour k=4 sauvegardée dans 'evolution_k4.csv'\n";
    return 0;
}
#endif