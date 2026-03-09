#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <cfloat>
#include <iomanip>
#include <algorithm>
#include <stdexcept>


/*
  Решение тестовой задачи (вариант 1) методом баланса + аналитическое решение.
  Вариант 1:
    xi = 0.4, mu1 = 0, mu2 = 1
    k1(x) = x+1  => k1* = xi + 1
    k2(x) = x    => k2* = xi
    q1(x) = x    => q1* = xi
    q2(x) = x^2  => q2* = xi^2
    f1(x) = x    => f1* = xi
    f2(x) = e^{-x} => f2* = e^{-xi}
*/

// Коэффициенты общего решения
// u(x) = A e^{λx} + B e^{-λx} + f/q
// (см. решение уравнения -k u'' + q u = f)
struct Coeffs {
    double A1, B1, A2, B2;
};

// Решаем коэффициенты A1,B1,A2,B2 (для случая q>0 на обоих интерв.)
Coeffs solve_constants_analytic(double xi,
    double k1, double k2,
    double q1, double q2,
    double f1, double f2,
    double mu1, double mu2)
{
    // λ = sqrt(q/k) — корни характеристического уравнения
    double lambda1 = sqrt(q1 / k1);
    double lambda2 = sqrt(q2 / k2);

    // Частные решения u_p = f/q
    double up1 = (q1 != 0.0) ? f1 / q1 : 0.0;
    double up2 = (q2 != 0.0) ? f2 / q2 : 0.0;

    // Экспоненты для компактности записи
    double e1p = exp(lambda1 * xi);
    double e1m = exp(-lambda1 * xi);
    double e2p = exp(lambda2 * xi);
    double e2m = exp(-lambda2 * xi);

    double e2p1 = exp(lambda2 * 1.0);
    double e2m1 = exp(-lambda2 * 1.0);

    // Матрица системы (4 уравнения из условий выше)
    std::vector<std::vector<double>> M(4, std::vector<double>(4, 0.0));
    std::vector<double> rhs(4, 0.0);

    // (1) u1(0)=mu1 → A1 + B1 = mu1 - up1
    M[0][0] = 1.0; M[0][1] = 1.0; M[0][2] = 0.0; M[0][3] = 0.0;
    rhs[0] = mu1 - up1;

    // (2) u2(1)=mu2 → A2 e^{λ2} + B2 e^{-λ2} = mu2 - up2
    M[1][0] = 0.0; M[1][1] = 0.0; M[1][2] = e2p1; M[1][3] = e2m1;
    rhs[1] = mu2 - up2;

    // u1(ksi)=u2(ksi) — непрерывность температуры
    M[2][0] = e1p; M[2][1] = e1m; M[2][2] = -e2p; M[2][3] = -e2m;
    rhs[2] = up2 - up1;

    // (4) k1 u1'(ksi) = k2 u2'(ksi) — непрерывность теплового потока
    // w = -k u'
    M[3][0] = k1 * lambda1 * e1p;
    M[3][1] = -k1 * lambda1 * e1m;
    M[3][2] = -k2 * lambda2 * e2p;
    M[3][3] = k2 * lambda2 * e2m;
    rhs[3] = 0.0;

    // Решение СЛАУ методом Гаусса (размер 4x4)
    Coeffs C;
    std::vector<std::vector<double>> Mc = M;
    std::vector<double> rc = rhs;
    int n = 4;
    for (int i = 0; i < n; ++i) {
        int piv = i;
        for (int j = i + 1; j < n; ++j)
            if (fabs(Mc[j][i]) > fabs(Mc[piv][i])) piv = j;
        if (piv != i) {
            std::swap(Mc[i], Mc[piv]);
            std::swap(rc[i], rc[piv]);
        }
        double diag = Mc[i][i];
        if (fabs(diag) < 1e-18) { std::cerr << "Singular system while solving constants\n"; break; }
        for (int j = i; j < n; ++j) Mc[i][j] /= diag;
        rc[i] /= diag;
        for (int r = 0; r < n; ++r) if (r != i) {
            double m = Mc[r][i];
            for (int c = i; c < n; ++c) Mc[r][c] -= m * Mc[i][c];
            rc[r] -= m * rc[i];
        }
    }
    C.A1 = rc[0]; C.B1 = rc[1]; C.A2 = rc[2]; C.B2 = rc[3];
    return C;
}

// Аналитическое u(x) (использует заранее найденные коэффициенты)
struct AnalyticSolver {
    double xi;
    double k1, k2, q1, q2, f1, f2;
    double mu1, mu2;
    double lambda1, lambda2;
    double up1, up2;
    Coeffs C;

    AnalyticSolver(double xi_,
        double k1_, double k2_,
        double q1_, double q2_,
        double f1_, double f2_,
        double mu1_, double mu2_)
        : xi(xi_), k1(k1_), k2(k2_), q1(q1_), q2(q2_), f1(f1_), f2(f2_), mu1(mu1_), mu2(mu2_)
    {
        lambda1 = sqrt(q1 / k1);
        lambda2 = sqrt(q2 / k2);
        up1 = (q1 != 0.0) ? f1 / q1 : 0.0;
        up2 = (q2 != 0.0) ? f2 / q2 : 0.0;
        C = solve_constants_analytic(xi, k1, k2, q1, q2, f1, f2, mu1, mu2);
    }

    double u(double x) const {
        if (x <= xi + 1e-14) {
            double e1 = exp(lambda1 * x);
            double e2 = exp(-lambda1 * x);
            return C.A1 * e1 + C.B1 * e2 + up1;
        }
        else {
            double e1 = exp(lambda2 * x);
            double e2 = exp(-lambda2 * x);
            return C.A2 * e1 + C.B2 * e2 + up2;
        }
    }
};

// ----------------------------- метод баланса -----------------------------

std::vector<double> balance_method(int n, int type,
    double xi,
    double k1s, double k2s,
    double q1s, double q2s,
    double f1s, double f2s,
    double mu1, double mu2)

{

    double h = 1.0 / n; // шаг сетки h = (b-a)/n
    int N = n + 1; // число узлов

    // Коэффициенты трехдиагональной матрицы
    std::vector<double> a(N, 0.0), b(N, 0.0), c(N, 0.0), d(N, 0.0), v(N, 0.0);

    // Граничные условия u0 = μ1, uN = μ2
    b[0] = 1.0; c[0] = 0.0; d[0] = mu1;
    a[n] = 0.0; b[n] = 1.0; d[n] = mu2;

    for (int i = 1; i < N - 1; ++i) {
        // Границы контрольного объема [x_{i-1/2}, x_{i+1/2}]
        double xL = i * h - h / 2.0;
        double xR = i * h + h / 2.0;

        double kL, kR, qL, qR, fL, fR;

        if (type == 0) {
            // Кусочно-заданные коэффициенты (формулы из условия задачи)
            kL = (xL < xi ? k1s : k2s);
            kR = (xR < xi ? k1s : k2s);
            qL = (xL < xi ? q1s : q2s);
            qR = (xR < xi ? q1s : q2s);
            fL = (xL < xi ? f1s : f2s);
            fR = (xR < xi ? f1s : f2s);
        }
        else if (type == 1) {
            // Кусочно-заданные функции
            kL = (xL < xi ? xL + 1 : xL);
            kR = (xR < xi ? xR + 1 : xR);

            qL = (xL < xi ? xL : xL * xL);
            qR = (xR < xi ? xR : xR * xR);

            fL = (xL < xi ? xL : exp(-xL));
            fR = (xR < xi ? xR : exp(-xR));
        }
        else {
            throw std::runtime_error("Unknown type in balance_method!");
        }

        // Ai = - a_i / h 
        // Bi = - a_{i+1} / h
        // Ci = Ai + Bi + d_i * h
        // фиi = (f_L + f_R)/2 * h 
        // Формируем коэффициенты для трёхдиагональной системы
        a[i] = -kL / h; //Ai 
        c[i] = -kR / h; //Bi 
        b[i] = kL / h + kR / h + (qL + qR) / 2.0 * h; //Ci 
        d[i] = (fL + fR) / 2.0 * h; // -фиi 
    }


    // Прогонка (трёхдиаг)

    int nEq = N;
    std::vector<double> P(nEq, 0.0), Q(nEq, 0.0);
    // Прямой ход
    P[0] = -c[0] / b[0]; // alpha_i
    Q[0] = d[0] / b[0]; // betta_i
    for (int i = 1; i < nEq; ++i) {
        double denom = b[i] + a[i] * P[i - 1];
        P[i] = (i == nEq - 1 ? 0.0 : -c[i] / denom); // подсчет alpha_i
        Q[i] = (d[i] - a[i] * Q[i - 1]) / denom; // Подсчет betta_i
    }
    v[nEq - 1] = Q[nEq - 1]; // Обратный ход: v_i = α_i v_{i+1} + β_i
    for (int i = nEq - 2; i >= 0; --i)
        v[i] = P[i] * v[i + 1] + Q[i];

    return v;
}

