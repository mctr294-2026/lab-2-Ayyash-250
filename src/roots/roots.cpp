#include "roots.hpp"
#include <cmath>
#include <iostream>

const double TOLERANCE = 1e-6;
const int MAX_ITERATIONS = 1e6;

bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root) {
    if (f(a) * f(b) >= 0) 
        return false;
    for (int i = 0; i < MAX_ITERATIONS; ++i) {
        double c = (a + b) / 2.0;
            if (std::abs(f(c)) < TOLERANCE || std::abs(b-a) < TOLERANCE) { // Checks if we close enough to the root
                *root = c; // Writes the answer to the pointer
                return true;
            }
            if (f(c) * f(a) > 0){
                a = c; //  Move the left wall (a) to the middle (c)
            }
            else {
                b = c; // Move the right wall (b) to the middle (c)
            }
    }
    return false;
    }


bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root) {
    if (f(a) * f(b) >= 0)
        return false;
    for (int i = 0; i < MAX_ITERATIONS; ++i) {
        double c = a - (f(a) * (b - a)) / (f(b) - f(a));
        if (std::abs(f(c)) < TOLERANCE) { // Checks if we close enough to the root
            *root = c; // Writes the answer to the pointer
            return true;
        }
        if (f(a) * f(c) > 0) {
            a = c; //  Move the left wall (a) to the middle (c)
        }
        else {
            b = c; // Move the right wall (b) to the middle (c)
            }
    }
    return false;
}

bool newton_raphson(std::function<double(double)> f,
    std::function<double(double)> g,
    double a, double b, double c,
    double *root) {
    double x_n = c;
    double x_n1 = 0;
    for (int i = 0; i < MAX_ITERATIONS; ++i) {
        if (std::abs(g(x_n)) < 1e-12) { // Checks if the formula's denominator is not zero
            return false;
        }
            x_n1 = x_n - (f(x_n) / g(x_n));
            if (std::abs(x_n1 - x_n) < TOLERANCE) { // Checks if we close enough to the root
                *root = x_n1; // Writes the answer to the pointer
                return true;
            }
            else {
                x_n = x_n1;
            }
            if (x_n1 < a || x_n > b) { // Checks if we go beyond the boundaries
                return false;
            }
    }
    return false;
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root) {
    double x_0 = c;
    double x_1 = c + 1e-6;
    for (int i = 0; i < MAX_ITERATIONS; ++i) {
        if (std::abs(f(x_1) - f(x_0)) < 1e-12) { // Checks if the formula's denominator is not zero
            return false;
        }
        double x_2 = x_1 - (f(x_1) * ((x_1 - x_0) / (f(x_1) - f(x_0)))); 
            if (std::abs(f(x_2) - f(x_1)) < TOLERANCE) { // Checks if we close enough to the root
                *root = x_2; // Writes the answer to the pointer
                return true;
        }
            else {
                x_0 = x_1;
                x_1 = x_2;
            }
            if (x_0 < a || x_1 > b) { // Checks if we go beyond the boundaries
                return false;
            }
    }

    return false;
}

