#pragma once
#include "Grid.h"
#include "LU_solve.h"
#include <iostream>
class FixedPointIteration
{
public:
    FixedPointIteration(Grid& grid) {
        this->n = grid.getN();
        q.resize(n);
        q0.resize(n);
        product.resize(n);
        l.resize(n - 1);
        d.resize(n);
        r.resize(n - 1);
        b.resize(n);
        L.resize(n - 1);
        D.resize(n);
        U.resize(n - 1);
        read_init_vector();
        make_q0(grid);
    }

    // w - параметр релаксации
    void fixed_point_iteration(Grid& grid, double w, double t1, double t2) {
        M_matrix_assembly(grid, t1, t2);
        b_vector_assembly(grid, t1, t2);
        G_matrix_assembly(grid, t2);
        boundary_condition(grid, t2);
        static double b_norm = vector_norm(b);
        double prod_norm;
        int iter = 0;
        do {
            for (int i = 0; i < n; i++)
                product[i] = q[i];
            calculation_lu(L, D, U, l, d, r, n);
            calculation_y(L, b, q, n);
            calucaltion_x(D, U, q, n);
            update_matrix(grid, t1, t2);
            for (int i = 0; i < n; i++)
                product[i] *= (1 - w);
            for (int i = 0; i < n; i++)
                q[i] *= w;
            for (int i = 0; i < n; i++)
                q[i] += product[i];
            multMatrVect();
            for (int i = 0; i < n; i++)
                product[i] -= b[i];
            prod_norm = vector_norm(product);
            iter++;
            //cout << iter << endl;
        } while (prod_norm / b_norm > eps && iter <= maxiter);
        for (int i = 0; i < n - 1; i++) {
            d[i] = 0;
            l[i] = 0;
            r[i] = 0;
            b[i] = 0;
        }
        d[n - 1] = 0;
        b[n - 1] = 0;
        q0 = q;
    }

    void print_vector(ofstream& fout, double t2) {
        fout << "t = " << t2 << endl;
        for (int i = 0; i < n; i++)
            fout << q[i] << endl;
        fout << "-----------------------" << endl;
    }

private:
    vector <vector <double>> G = { {0, 0},
                                   {0, 0} };
    vector <vector<double>> G1 = { {1.0 / 2, -1.0 / 2},
                                  {-1.0 / 2, 1.0 / 2} };

    vector <vector <double>> M = { {0, 0},
                                   {0, 0} };
    vector <vector<double>> C = { {1.0 / 3, 1.0 / 6},
                                  {1.0 / 6, 1.0 / 3} };

    int n;  // Размерность вектора d
    double eps = 1E-13;
    int maxiter = 100;
    vector <double> product;
    vector <double> q;
    vector <double> b;
    // Глобальная матрица
    vector <double> l;
    vector <double> d;
    vector <double> r;
    static const int m = 1;  // Полуширина ленты
    // Матрица LU разложения
    vector <double> L;
    vector <double> U;
    vector <double> D;
    // Вектор разложения функции начального условия u0
    vector <double> q0;

    void read_init_vector() {
        ifstream init_vect;
        init_vect.open("init_vect.txt");
        for (int i = 0; i < n; i++)
            init_vect >> q[i];
        init_vect.close();
    }

    void make_q0(Grid& grid) {
        for (int i = 0; i < n; i++)
            q0[i] = grid.u0(grid.getX(i));
    }

    // Занесение локальной матрицы a в глобальную матрицу
    // k1, k2, k3 - глобальные номера базисных функций
    void add_local_matrix(vector<vector <double>>& a, int k) {
        d[k] += a[0][0];
        d[k + 1] += a[1][1];
        l[k] += a[1][0];
        r[k] += a[0][1];
    }


    // Умножение глобальной матрицы{l, d, r} на вектор q, результат в product
    void multMatrVect() {
        if (n >= 2)
            product[0] = d[0] * q[0] + r[0] * q[1];
        for (int i = 1; i < n - 1; i++) {
            product[i] = l[i - 1] * q[i - 1] +
                d[i] * q[i] + r[i] * q[i + 1];
        }
        product[n - 1] = l[n - 2] * q[n - 2] + d[n - 1] * q[n - 1];
    }


    void G_matrix_assembly(Grid& grid , double t2) {
        double lambda1, lambda2;
        double h;
        double x1, x2;
        double g;
        int wi;  // Номер подобласти
        for (int i = 0; i < n - 1; i++) {
            x2 = grid.getX(i + 1);
            x1 = grid.getX(i);
            h = x2 - x1;
            wi = grid.inSubArea(i);
            lambda1 = grid.lambda(wi, q[i], t2);
            lambda2 = grid.lambda(wi, q[i + 1], t2);
            g = (lambda1 + lambda2) / h;
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    G[i][j] = g * G1[i][j];
            add_local_matrix(G, i);
        }
    }


    void M_matrix_assembly(Grid& grid, double t1, double t2) {
        double delta_t = t2 - t1;
        double h;
        double x1, x2;
        double g;
        int wi;  // Номер подобласти
        for (int i = 0; i < n - 1; i++) {
            x2 = grid.getX(i + 1);
            x1 = grid.getX(i);
            h = x2 - x1;
            wi = grid.inSubArea(i);
            g = h * grid.sigma(wi) / delta_t;
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    M[i][j] = g * C[i][j];
            add_local_matrix(M, i);
        }
    }


    void b_vector_assembly(Grid& grid, double t1, double t2) {
        double delta_t = t2 - t1;
        double h;
        double x1, x2;
        double g;
        int wi;  // Номер подобласти
        for (int i = 0; i < n - 1; i++) {
            x2 = grid.getX(i + 1);
            x1 = grid.getX(i);
            h = x2 - x1;
            wi = grid.inSubArea(i);
            g = grid.sigma(wi) / delta_t * h / 6;
            // ТУТ ДОЛЖНО БЫТЬ ВРЕМЯ(t)
            b[i] += h / 6 * (2 * grid.f(wi, x1, t2) + grid.f(wi, x2, t2));
            b[i + 1] += h / 6 * (grid.f(wi, x1, t2) + 2 * grid.f(wi, x2, t2));

            b[i] +=  g * (2 * q0[i] + q0[i + 1]);
            b[i + 1] += g * (q0[i] + 2 * q0[i + 1]);
        }
    }

    double u1(double t) {
        return 55 + 2 * t;
    }

    double u2(double t) {
        return 125 + 2 * t;
    }

    void boundary_condition(Grid& grid, double t2) {
        double u1_ = u1(t2);
        double u2_ = u2(t2);
        r[0] = 0;
        d[0] = 1;
        b[0] = u1_;

        l[n - 2] = 0;
        d[n - 1] = 1;
        b[n - 1] = u2_;
    }


    double vector_norm(vector <double>& v) {
        double norm = 0;
        for (int i = 0; i < n; i++)
            norm += v[i] * v[i];
        return sqrt(norm);
    }

    // ЖУТКО НЕЭФФЕКТИВНО
    void update_matrix(Grid& grid ,double t1, double t2) {
        for (int i = 0; i < n - 1; i++) { // костыль
            d[i] = 0;
            l[i] = 0;
            r[i] = 0;
        }
        d[n - 1] = 0;
        M_matrix_assembly(grid, t1, t2); // костыль
        G_matrix_assembly(grid, t2);
        boundary_condition(grid, t2);  // костыль
    }
};