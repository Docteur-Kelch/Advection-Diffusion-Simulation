// ------------------ STRUCTURES DE DONNÉES (Q2 - Q7) ------------------

#ifndef MAILLAGE
#define MAILLAGE
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
#include <algorithm>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Vector2d {
    double x, y;

      Vector2d(double x_ = 0.0, double y_ = 0.0) : x(x_), y(y_) {}
};

// Classe Vertex : représente les nœuds du maillage
class Vertex {
public:
    Vector2d x; // Coordonnées géométriques du sommet
    Vertex(double x_ = 0.0, double y_ = 0.0) : x(x_, y_) {}
};

class Cell; // Forward declaration

// Classe Edge : interface entre deux cellules
class Edge {
public:
    Vector2d x_sigma;        // Centre de l'arête
    Vector2d n_sigma;        // Normale de référence
    std::array<Vertex*, 2> V; // Les deux sommets extrémités
    std::array<Cell*, 2> T;  // Cellules adjacentes (T[0] gauche/bas, T[1] droite/haut)
    int bnd;                 // 0: intérieur, 1: Dirichlet, 2: Neumann
    double length;           // Longueur |sigma|
    double d_sigma;          // Distance entre centres (ou centre-bord)

    Edge() : bnd(0), length(0.0), d_sigma(0.0) {
        V[0] = V[1] = nullptr;
        T[0] = T[1] = nullptr;
    }
};

// Classe Cell : volume de contrôle
class Cell {
public:
    int id;
    Vector2d x_ij;           // Barycentre
    std::array<Edge*, 4> E;  // Arêtes [gauche, droite, bas, haut]
    std::array<double, 4> w; // Orientation du flux
    double area;             // Surface |K|

    Cell() : id(-1), area(0.0) {
        for (int i = 0; i < 4; ++i) {
            E[i] = nullptr;
            w[i] = 0.0;
        }
    }
};

// ------------------ MAILLAGE (Q2) ------------------

class Mesh {
public:
    std::vector<Vertex> vertices; // Liste des sommets
    std::vector<Edge> edges;
    std::vector<Cell> cells;
    double h;
    int M;

    void build_Cartesian(int k) {
        M = std::pow(2, k); // M = 2^k subdivisions
        h = 1.0 / M;
        cells.resize(M * M);

        // 1. Création des sommets (Grille de (M+1)x(M+1))
        vertices.resize((M + 1) * (M + 1));
        for (int j = 0; j <= M; ++j) {
            for (int i = 0; i <= M; ++i) {
                vertices[j * (M + 1) + i].x = Vector2d(i * h, j * h);
            }
        }

        // 2. Création des cellules
        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < M; ++i) {
                int id = i + j * M;
                cells[id].id = id;
                cells[id].x_ij = Vector2d((i + 0.5) * h, (j + 0.5) * h);
                cells[id].area = h * h;
            }
        }

        // 3. Création des arêtes (Verticales puis Horizontales)
        int N_vert = (M + 1) * M;
        int N_horiz = M * (M + 1);
        edges.resize(N_vert + N_horiz);

        // Arêtes Verticales
        for (int j = 0; j < M; ++j) {
            for (int i = 0; i <= M; ++i) {
                int eid = j * (M + 1) + i;
                Edge& e = edges[eid];
                e.length = h;
                e.x_sigma = Vector2d(i * h, (j + 0.5) * h);
                e.n_sigma = Vector2d(1.0, 0.0);
                e.V[0] = &vertices[j * (M + 1) + i];       // Sommet bas
                e.V[1] = &vertices[(j + 1) * (M + 1) + i]; // Sommet haut

                if (i == 0) { // Bord gauche Dirichlet
                    e.bnd = 1;
                    e.d_sigma = h / 2.0;
                    e.T[0] = nullptr;
                    e.T[1] = &cells[j * M];
                }
                else if (i == M) { // Bord droit Dirichlet
                    e.bnd = 1;
                    e.d_sigma = h / 2.0;
                    e.T[0] = &cells[j * M + i - 1];
                    e.T[1] = nullptr;
                }
                else { // Intérieur
                    e.bnd = 0;
                    e.d_sigma = h;
                    e.T[0] = &cells[j * M + i - 1];
                    e.T[1] = &cells[j * M + i];
                }
            }
        }

        // Arêtes Horizontales
        int offset = N_vert;
        for (int j = 0; j <= M; ++j) {
            for (int i = 0; i < M; ++i) {
                int eid = offset + j * M + i;
                Edge& e = edges[eid];
                e.length = h;
                e.x_sigma = Vector2d((i + 0.5) * h, j * h);
                e.n_sigma = Vector2d(0.0, 1.0);
                e.V[0] = &vertices[j * (M + 1) + i];     // Sommet gauche
                e.V[1] = &vertices[j * (M + 1) + i + 1]; // Sommet droit

                if (j == 0) { // Bord bas Neumann
                    e.bnd = 2;
                    e.d_sigma = h / 2.0;
                    e.T[0] = nullptr;
                    e.T[1] = &cells[i];
                }
                else if (j == M) { // Bord haut Neumann
                    e.bnd = 2;
                    e.d_sigma = h / 2.0;
                    e.T[0] = &cells[(j - 1) * M + i];
                    e.T[1] = nullptr;
                }
                else { // Intérieur
                    e.bnd = 0;
                    e.d_sigma = h;
                    e.T[0] = &cells[(j - 1) * M + i];
                    e.T[1] = &cells[j * M + i];
                }
            }
        }

        // 4. Connectivité Cellules -> Arêtes
        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < M; ++i) {
                Cell& C = cells[i + j * M];
                C.E[0] = &edges[j * (M + 1) + i];           C.w[0] = -1.0; // Gauche
                C.E[1] = &edges[j * (M + 1) + i + 1];       C.w[1] = 1.0;  // Droite
                C.E[2] = &edges[offset + j * M + i];        C.w[2] = -1.0; // Bas
                C.E[3] = &edges[offset + (j + 1) * M + i];  C.w[3] = 1.0;  // Haut
            }
        }
    }
};

#endif // !MAILLAGE
