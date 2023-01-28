#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <math.h>
#include <random>
#include <string>
#include "constants.h"

bool iszero(double x);
bool iszero(double x, double err);
double rd(double x, double y, double z);
bool file_exists(std::string fileName);
std::string time_to_string();
std::string time_to_filename();

template <typename T>
int gaussian_random(int n, T& random_num, double sigma=1, double avg=0){
    std::default_random_engine generator;
    generator.seed(rand());
    std::normal_distribution<double> distribution(avg,sigma);
    for(int i=0; i<n; ++i) random_num[i] = distribution(generator);
    return 0;
}

template <typename T>
int gaussian_random_adjust(int n, T& random_num, double sigma, double avg=0) {
    double mean = 0;
    for(int i=0; i<n; ++i) mean += random_num[i];
    mean /= n;

    double sigma_calc = 0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+ : sigma_calc)
    #endif // _OPENMP
    for(int i=0; i<n; ++i) {
        random_num[i] -= mean;
        sigma_calc += random_num[i]*random_num[i];
    }
    sigma_calc = sqrt(sigma_calc/n);

    double adjust_width = sigma/sigma_calc;
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for(int i=0; i<n; ++i) {
        random_num[i] *= adjust_width;
        random_num[i] += avg;
    }
    return 0;
}

template <typename T>
int uniform_random(int n, T& random_num, double r_min, double r_max){
    std::default_random_engine generator;
    generator.seed(rand());
    std::uniform_real_distribution<double> uniform_dis(r_min,r_max);
    for(int i=0; i<n; ++i) random_num[i] = uniform_dis(generator);
    return 0;
}

template <typename T>
int uniform_random_adjust(int n, T& random_num, double avg) {
    double mean = 0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+ : mean)
    #endif // _OPENMP
    for(unsigned int i=0; i<n; ++i) mean += random_num[i];
    mean /= n;
    double adjust = avg - mean;
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for(unsigned int i=0; i<n; ++i) random_num[i] += adjust;
    return 0;
}

template <typename T>
int uniform_random_in_circle(int n, double radius, T& x, T& y){
    uniform_random(n, x, 0, 1);
    uniform_random(n, y, 0, 2*k_pi);

    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for(int i=0; i<n; ++i) {
        double r = radius*sqrt(x.at(i));
        double u = r * cos(y.at(i));
        x.at(i) = u;
        y.at(i) = sqrt(r*r - u*u);
    }
    return 0;
}

template <typename T>
int uniform_random_in_hollow_circle(int n, double radius_out, double radius_in, T& x, T& y){
    uniform_random(n, x, 0, 1);
    uniform_random(n, y, 0, 2*k_pi);

    double out_r2 = radius_out * radius_out;
    double in_r2 = radius_in * radius_in;
    out_r2 -= in_r2;

    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for(int i=0; i<n; ++i) {
        double r = sqrt(x.at(i)*out_r2+in_r2);
        double u = r * cos(y.at(i));
        x.at(i) = u;
        y.at(i) = sqrt(r*r - u*u);
    }
    return 0;
}

template <typename T>
int uniform_random_in_ellipse(int n, double a, double b, T& x, T& y){
    uniform_random(n, x, 0, 1);
    std::vector<double> u(n);
    uniform_random(n, u, 0, 1);
    uniform_random(n, y, 0, 1);

    double c = b/a;
    double ab = a*b;
    double bb = b*b;
    double aa = a*a;
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for(int i=0; i<n; ++i) {
        double theta = atan(c*tan(k_pi*x.at(i)/2));
        if (u.at(i)<0.25) x.at(i) = theta;
        else if (u.at(i)<0.5) x.at(i) = k_pi - theta;
        else if (u.at(i)<0.75) x.at(i) = k_pi + theta;
        else x.at(i) = -theta;
        double cs = cos(x.at(i));
        double cc = cs*cs;
        double ss = 1-cc;
        double r = ab*sqrt(y.at(i)/(bb*cc+aa*ss));
        cs *= r;
        x.at(i) = cs;
        y.at(i) = sqrt(r*r-cs*cs);
    }
    return 0;
}

template <typename T>
double rms(int n, T& v) {
    double sum = 0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+ : sum)
    #endif // _OPENMP
    for(int i=0; i<n; ++i) {
        sum += v[i];
    }
    double avg = sum/n;
    sum = 0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+ : sum)
    #endif // _OPENMP
    for(int i=0; i<n; ++i) {
       double adj = v[i] - avg;
       sum += adj*adj;
    }
    return sqrt(sum/n);
}

#endif // FUNCTIONS_H
