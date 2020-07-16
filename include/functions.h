#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <random>
#include <string>

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
    for(int i=0; i<n; ++i) {
        random_num[i] -= mean;
        sigma_calc += random_num[i]*random_num[i];
    }
    sigma_calc = sqrt(sigma_calc/n);

    double adjust_width = sigma/sigma_calc;
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
    for(unsigned int i=0; i<n; ++i) mean += random_num[i];
    mean /= n;
    double adjust = avg - mean;
    for(unsigned int i=0; i<n; ++i) random_num[i] += adjust;
    return 0;
}

template <typename T>
double rms(int n, T& v) {
    double sum = 0;
    for(int i=0; i<n; ++i) {
        sum += v[i];
    }
    double avg = sum/n;
    sum = 0;
    for(int i=0; i<n; ++i) {
       double adj = v[i] - avg;
       sum += adj*adj;
    }
    return sqrt(sum/n);
}

#endif // FUNCTIONS_H
