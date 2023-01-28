
#include "beam.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include "arbitrary_electron_beam.h"
#include "functions.h"
#include "ions.h"

Beam::Beam(int charge_number, double mass_number, double kinetic_energy, double emit_nx, double emit_ny, double dp_p,
           double sigma_s, double n_particle): charge_number_(charge_number), mass_number_(mass_number),
           kinetic_energy_(kinetic_energy), emit_nx_(emit_nx), emit_ny_(emit_ny), dp_p_(dp_p), sigma_s_(sigma_s),
           particle_number_(n_particle) {
    mass_ = mass_number*k_u;
    gamma_ = 1+kinetic_energy_/mass_;
    beta_ = sqrt(gamma_*gamma_-1)/gamma_;
    r_ = k_ke*charge_number_*charge_number_*k_e*1e-6/mass_;
    bunched_ = (sigma_s_>0)?true:false;
    emit_x_ = emit_nx_/(beta_*gamma_);
    emit_y_ = emit_ny_/(beta_*gamma_);
    energy_spread_ = beta_*beta_*dp_p_;
    dv_v_ = dp_p_/(gamma_*gamma_);
    p0_ = gamma_*mass_*1e6*k_e*beta_/k_c;
}

Beam::Beam(int charge_number, double mass_number, double kinetic_energy, double emit_nx, double emit_ny, double dp_p,
           double n_particle): charge_number_(charge_number), mass_number_(mass_number), kinetic_energy_(kinetic_energy),
           emit_nx_(emit_nx), emit_ny_(emit_ny), dp_p_(dp_p), particle_number_(n_particle) {
    mass_ = mass_number_*k_u;
    gamma_ = 1+kinetic_energy_/mass_;
    beta_ = sqrt(gamma_*gamma_-1)/gamma_;
    r_ = k_ke*charge_number_*charge_number_*k_e*1e-6/mass_;
    bunched_ = false;
    sigma_s_ = -1;
    emit_x_ = emit_nx_/(beta_*gamma_);
    emit_y_ = emit_ny_/(beta_*gamma_);
    energy_spread_ = beta_*beta_*dp_p_;
    dv_v_ = dp_p_/(gamma_*gamma_);
    p0_ = gamma_*mass_*1e6*k_e*beta_/k_c;
}

int Beam::set_center(int i, double x) {
    if(i<3) {
        center_[i] = x;
        return 0;
    }
    else {
        perror("Error index for electron beam center!");
        return 1;
    }
}

vector<double>& EBeam::get_v(EBeamV v) {
    switch(v) {
    case EBeamV::TPR_TR: return tpr_t;
    case EBeamV::TPR_L: return tpr_l;
    case EBeamV::V_RMS_TR: return v_rms_t;
    case EBeamV::V_RMS_L: return v_rms_l;
    case EBeamV::V_AVG_X: return v_avg_x;
    case EBeamV::V_AVG_Y: return v_avg_y;
    case EBeamV::V_AVG_L: return v_avg_l;
    default:perror("Wrong phase coordinate selected!");return tpr_t;
    }
}

void EBeam::set_tpr(double tpr_tr, double tpr_long) {
    tpr_t.resize(1);
    tpr_l.resize(1);
    tpr_t.at(0) =  tpr_tr;
    tpr_l.at(0) = tpr_long;

    v_rms_t.resize(1);
    v_rms_l.resize(1);
    v_rms_t.at(0) = sqrt(tpr_tr/k_me)*0.001*k_c;
    v_rms_l.at(0) = sqrt(tpr_long/k_me)*0.001*k_c;
}

void EBeam::set_v_rms(double v_rms_tr, double v_rms_long) {
    v_rms_t.resize(1);
    v_rms_l.resize(1);
    v_rms_t.at(0) = v_rms_tr;
    v_rms_l.at(0) = v_rms_long;
    tpr_t.resize(1);
    tpr_l.resize(1);
    tpr_t.at(0) = v_rms_tr*v_rms_tr*k_me*1e6/(k_c*k_c);
    tpr_l.at(0) = v_rms_long*v_rms_long*k_me*1e6/(k_c*k_c);
}

void EBeam::set_v_avg(double v_avg_tx, double v_avg_ty, double v_avg_long) {
    v_avg_x.resize(1);
    v_avg_y.resize(1);
    v_avg_l.resize(1);
    v_avg_x.at(0) = v_avg_tx;
    v_avg_y.at(0) = v_avg_ty;
    v_avg_l.at(0) = v_avg_long;
}

void GaussianBunch::set_angles(double sigma_xp, double sigma_yp, double sigma_dpp) {
    sigma_xp_ = sigma_xp;
    sigma_yp_ = sigma_yp;
    sigma_dpp_ = sigma_dpp;
    double tpr_t = (sigma_xp*sigma_xp + sigma_yp*sigma_yp)*beta_*beta_*gamma_*gamma_*k_me*1e6/2;
    double tpr_l = sigma_dpp*sigma_dpp*beta_*beta_*k_me*1e6;
    set_tpr(tpr_t, tpr_l);
}

void UniformCylinder::density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n_particle) {
    int nq = this->charge_number();
    if (nq<0) nq *= -1;
    double r2 = radius_*radius_;
    double density = current_/(k_pi*r2*nq*k_e*this->beta()*k_c);
    std::fill(ne.begin(), ne.end(), 0);
    for(int i=0; i<n_particle; ++i){
        if(x.at(i)*x.at(i)+y.at(i)*y.at(i)<=r2) ne.at(i) = density;
    }
}

void UniformCylinder::density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n_particle,
                double cx, double cy, double cz) {
    int nq = this->charge_number();
    if (nq<0) nq *= -1;
    double r2 = radius_*radius_;
    double density = current_/(k_pi*r2*nq*k_e*beta_*k_c);
    std::fill(ne.begin(), ne.end(), 0);
    //ion_center - electron_center
    cx -= center_[0];
    cy -= center_[1];
    cz -= center_[2];
    for(int i=0; i<n_particle; ++i){
        if((x.at(i)+cx)*(x.at(i)+cx)+(y.at(i)+cy)*(y.at(i)+cy)<=r2) ne.at(i) = density;
    }
}

void UniformHollow::density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n_particle) {
    int nq = this->charge_number();
    double out_r2 = out_radius_*out_radius_;
    double in_r2 = in_radius_*in_radius_;
    double area = k_pi*(out_r2-in_r2);
    double density;
    if (area!=0)
        density = current_/(area*nq*k_e*this->beta()*k_c);
    else
        density = 0;
    if (density<0) density *= -1;
    std::fill(ne.begin(),ne.end(),0);

    for(int i=0; i<n_particle; ++i){
        double r2 = x[i]*x[i]+y[i]*y[i];
        if(r2<=out_r2 && r2>=in_r2) ne[i] = density;
    }
}

void UniformHollow::density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n_particle,
                double cx, double cy, double cz) {
    int nq = this->charge_number();
    double out_r2 = out_radius_*out_radius_;
    double in_r2 = in_radius_*in_radius_;
    double area = k_pi*(out_r2-in_r2);
    double density;
    if (area!=0)
        density = current_/(area*nq*k_e*this->beta()*k_c);
    else
        density = 0;
    if (density<0) density *= -1;
    std::fill(ne.begin(),ne.end(),0);
     //ion_center - electron_center
    cx -= this->center(0);
    cy -= this->center(1);
    cz -= this->center(2);
    for(int i=0; i<n_particle; ++i){
        double r2 = (x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy);
        if(r2<=out_r2 && r2>=in_r2) ne[i] = density;
    }
}


void UniformHollowBunch::density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n_particle) {

    int nq = this->charge_number();
    double out_r2 = out_radius_*out_radius_;
    double in_r2 = in_radius_*in_radius_;
    double area = k_pi*(out_r2-in_r2);
    double density;
    if (area!=0)
        density = current_/(area*nq*k_e*this->beta()*k_c);
    else
        density = 0;
    if (density<0) density *= -1;
    std::fill(ne.begin(),ne.end(),0);

    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    for(int i=0; i<n_particle; ++i){
        double r2 = x[i]*x[i]+y[i]*y[i];
        if(z[i]<=right_end && z[i]>=left_end &&r2<=out_r2 && r2>=in_r2) ne[i] = density;
    }
}

void UniformHollowBunch::density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n_particle,
                double cx, double cy, double cz) {

    int nq = this->charge_number();
    double out_r2 = out_radius_*out_radius_;
    double in_r2 = in_radius_*in_radius_;
    double area = k_pi*(out_r2-in_r2);
    double density;
    if (area!=0)
        density = current_/(area*nq*k_e*this->beta()*k_c);
    else
        density = 0;
    if (density<0) density *= -1;
    //ion_center - electron_center
    cx -= this->center(0);
    cy -= this->center(1);
    cz -= this->center(2);
    std::fill(ne.begin(),ne.end(),0);

    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    for(int i=0; i<n_particle; ++i){
        double r2 = (x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy);
        double z_shifted = z[i]+cz;
        if(z_shifted<=right_end && z_shifted>=left_end &&r2<=out_r2 && r2>=in_r2) ne[i] = density;
    }
}

double UniformHollowBunch::n_electron(){
    int nq = this->charge_number();
    double out_r2 = out_radius_*out_radius_;
    double in_r2 = in_radius_*in_radius_;
    double area = k_pi*(out_r2-in_r2);
    double density;
    if (area!=0)
        density = current_/(area*nq*k_e*this->beta()*k_c);
    else
        density = 0;
    if (density<0) density *= -1;

    return density*area;
}
//
//void UniformHollowBunch::create_particle_location(){
//    int n_sample = samples->n_sample();
//
//    uniform_random(n_sample, samples->z, -length_/2, length_/2);
//    uniform_random(n_sample, samples->x, 0, 1);
//    uniform_random(n_sample, samples->y, 0, 2*k_pi);
//    double out_r2 = out_radius_*out_radius_;
//    double in_r2 = in_radius_*in_radius_;
//    out_r2 -= in_r2;
//    for(int i=0; i<n_sample; ++i) {
//        double r = sqrt(samples->x.at(i)*out_r2+in_r2);
//        double x = r * cos(samples->y.at(i));
//        samples->x.at(i) = x;
//        samples->y.at(i) = sqrt(r*r - x*x);
//    }
//}


void UniformHollowBunch::create_particle_location(){
    int n_sample = samples->n_sample();
    uniform_random(n_sample, samples->z, -length_/2, length_/2);
    uniform_random_in_hollow_circle(n_sample, out_radius_, in_radius_, samples->x, samples->y);
}


double UniformBunch::n_electron() {
    int nq = this->charge_number();
    if (nq<0) nq *= -1;
    double r2 = radius_*radius_;
    double density = current_/(k_pi*r2*nq*k_e*this->beta()*k_c);
    double volume = k_pi*r2*length_;
    return density*volume;
}

//void UniformBunch::create_particle_location() {
//    int n_sample = samples->n_sample();
//
//    uniform_random(n_sample, samples->z, -length_/2, length_/2);
//    uniform_random(n_sample, samples->x, 0, 1);
//    uniform_random(n_sample, samples->y, 0, 2*k_pi);
//    for(int i=0; i<n_sample; ++i) {
//        double r = radius_*sqrt(samples->x.at(i));
//        double x = r * cos(samples->y.at(i));
//        samples->x.at(i) = x;
//        samples->y.at(i) = sqrt(r*r - x*x);
//    }
//}

void UniformBunch::create_particle_location() {
    int n_sample = samples->n_sample();
    uniform_random(n_sample, samples->z, -length_/2, length_/2);
    uniform_random_in_circle(n_sample, radius_, samples->x, samples->y);
}


void UniformBunch::density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n_particle) {

    int nq = this->charge_number();
    if (nq<0) nq *= -1;
    double r2 = radius_*radius_;
    double density = current_/(k_pi*r2*nq*k_e*this->beta()*k_c);
    std::fill(ne.begin(),ne.end(),0);
    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    for(int i=0; i<n_particle; ++i){
        if(z[i]<=right_end && z[i]>= left_end && x[i]*x[i]+y[i]*y[i]<=r2)
            ne[i] = density;
    }
}

void UniformBunch::density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n_particle,
                double cx, double cy, double cz) {
    int nq = this->charge_number();
    if (nq<0) nq *= -1;
    double r2 = radius_*radius_;
    double density = current_/(k_pi*r2*nq*k_e*this->beta()*k_c);

    //ion_center - electron_center
    cx -= this->center(0);
    cy -= this->center(1);
    cz -= this->center(2);
    std::fill(ne.begin(),ne.end(),0);
    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    for(int i=0; i<n_particle; ++i){
        if((z[i]+cz)<=right_end && (z[i]+cz)>=left_end && (x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy)<=r2)
            ne[i] = density;
    }
}

double EllipticUniformBunch::n_electron() {
    int nq = this->charge_number();
    if (nq<0) nq *= -1;
    double area = k_pi*rh_*rv_;
    double density = current_/(area*nq*k_e*this->beta()*k_c);
    double volume = area*length_;
    return density*volume;
}

void EllipticUniformBunch::create_particle_location() {
    int n_sample = samples->n_sample();
    uniform_random(n_sample, samples->z, -length_/2, length_/2);
    uniform_random_in_ellipse(n_sample, rh_, rv_, samples->x, samples->y);
}


void EllipticUniformBunch::density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n_particle) {
    int nq = this->charge_number();
    if (nq<0) nq *= -1;
    double density = current_/(k_pi*rh_*rv_*nq*k_e*this->beta()*k_c);
    std::fill(ne.begin(),ne.end(),0);
    double inv_rh2 = 1.0/(rh_*rh_);
    double inv_rv2 = 1.0/(rv_*rv_);
    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    for(int i=0; i<n_particle; ++i){
        if(z[i]<=right_end && z[i]>=left_end && inv_rh2*x[i]*x[i]+inv_rv2*y[i]*y[i]<=1)
            ne[i] = density;
    }
}

void EllipticUniformBunch::density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n_particle,
                double cx, double cy, double cz) {
    int nq = this->charge_number();
    if (nq<0) nq *= -1;
    double density = current_/(k_pi*rh_*rv_*nq*k_e*this->beta()*k_c);

    //ion_center - electron_center
    cx -= this->center(0);
    cy -= this->center(1);
    cz -= this->center(2);
    std::fill(ne.begin(),ne.end(),0);
    double inv_rh2 = 1.0/(rh_*rh_);
    double inv_rv2 = 1.0/(rv_*rv_);
    double left_end = -0.5*length_;
    double right_end = 0.5*length_;
    for(int i=0; i<n_particle; ++i){
        if((z[i]+cz)<=right_end && (z[i]+cz)>=left_end &&
           inv_rh2*(x[i]+cx)*(x[i]+cx)+inv_rv2*(y[i]+cy)*(y[i]+cy)<=1)
            ne[i] = density;
    }
}

void GaussianBunch::create_particle_location() {
    int n_sample = samples->n_sample();
    gaussian_random(n_sample, samples->x, sigma_x_);
    gaussian_random(n_sample, samples->y, sigma_y_);
    gaussian_random(n_sample, samples->z, sigma_s_);
    gaussian_random_adjust(n_sample, samples->x, sigma_x_);
    gaussian_random_adjust(n_sample, samples->y, sigma_y_);
    gaussian_random_adjust(n_sample, samples->z, sigma_s_);
}

void GaussianBunch::density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n_particle) {
    double amp = n_electron_/(sqrt(8*k_pi*k_pi*k_pi)*sigma_x_*sigma_y_*sigma_s_);
    double sigma_x2 = -1/(2*sigma_x_*sigma_x_);
    double sigma_y2 = -1/(2*sigma_y_*sigma_y_);
    double sigma_s2 = -1/(2*sigma_s_*sigma_s_);
    for(int i=0; i<n_particle; ++i){
        ne[i] = amp*exp(x[i]*x[i]*sigma_x2+y[i]*y[i]*sigma_y2+z[i]*z[i]*sigma_s2);
    }
}

void GaussianBunch::density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n_particle,
                double cx, double cy, double cz) {
    double amp = n_electron_/(sqrt(8*k_pi*k_pi*k_pi)*sigma_x_*sigma_y_*sigma_s_);
    double sigma_x2 = -1/(2*sigma_x_*sigma_x_);
    double sigma_y2 = -1/(2*sigma_y_*sigma_y_);
    double sigma_s2 = -1/(2*sigma_s_*sigma_s_);
    //ion_center - electron_center
    cx -= this->center(0);
    cy -= this->center(1);
    cz -= this->center(2);
    for(int i=0; i<n_particle; ++i){
        ne[i] = amp*exp((x[i]+cx)*(x[i]+cx)*sigma_x2+(y[i]+cy)*(y[i]+cy)*sigma_y2+(z[i]+cz)*(z[i]+cz)*sigma_s2);
    }
}

void ParticleBunch::load_particle(long int n) {
    if(n>=0) n_ = load_electrons(x, y, z, vx, vy, vz, filename_, n, line_skip_, binary_, buffer_);
    create_e_tree(x, y, z, n_, s_, tree_, list_e_);
    if(length_==0) {
        auto itr = z.begin();
        double z_max = *itr;
        double z_min = *itr;
        ++itr;
        for(; itr!=z.end(); ++itr) {
            if(*itr>z_max) z_max = *itr;
            if(*itr<z_min) z_min = *itr;
        }
        length_ = z_max - z_min;
    }
}

void ParticleBunch::load_particle() {
    load_particle(0);
}


void ParticleBunch::density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n) {
    double rate = n_electron_/n_;
    std::vector<int> list_i;
    int idx_out;
    create_ion_tree(x, y, z, n, tree_, list_i, idx_out);
    if (v_x_corr_)
        {::density(tree_, list_e_, vx, vy, vz, n_, list_i, idx_out, n, ne, v_avg_l, v_rms_t, v_rms_l);}
    else
        {::density(tree_, list_e_, vx, vy, vz, n_, list_i, idx_out, n, ne, v_rms_t, v_rms_l);}
    for(int i=0; i<n; ++i) ne[i] *= rate;
    return ;
}

void ParticleBunch::density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n, double cx, double cy,
                       double cz) {
    double rate = n_electron_/n_;
    cx -= this->center(0);
    cy -= this->center(1);
    cz -= this->center(2);
    for(int i=0; i<n; ++i) {
        x[i] += cx;
        y[i] += cy;
        z[i] += cz;
    }
    std::vector<int> list_i;
    int idx_out;
    create_ion_tree(x, y, z, n, tree_, list_i, idx_out);
    if (v_x_corr_){
        velocity_ = Velocity::VARY_Z;
        ::density(tree_, list_e_, vx, vy, vz, n_, list_i, idx_out, n, ne, v_avg_l, v_rms_t, v_rms_l);
    }
    else{
        velocity_ = Velocity::CONST;
        ::density(tree_, list_e_, vx, vy, vz, n_, list_i, idx_out, n, ne, v_rms_t, v_rms_l);
    }
    for(int i=0; i<n; ++i) {
        x[i] -= cx;
        y[i] -= cy;
        z[i] -= cz;
    }
    for(int i=0; i<n; ++i) ne[i] *= rate;
}

void ParticleBunch::create_particle_location() {
    std::copy(samples->x.begin(), samples->x.end(), x.begin());
    std::copy(samples->y.begin(), samples->y.end(), y.begin());
    std::copy(samples->z.begin(), samples->z.end(), z.begin());
}

void EBeam::multi_density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n) {
    vector<double> d(n);
    for(int i=0; i<n_; ++i) {
        density(x, y, z, d, n, -cx_.at(i), -cy_.at(i), -cz_.at(i));
        for(int j=0; j<n; ++j) ne.at(j) += d.at(j);
    }
}

void EBeam::multi_density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n, double cx, double cy,
                       double cz) {
    vector<double> d(n);
    for(int i=0; i<n_; ++i) {
        density(x, y, z, d, n, cx-cx_.at(i), cy-cy_.at(i), cz-cz_.at(i));
        for(int j=0; j<n; ++j) ne.at(j) += d.at(j);
    }
}

void UniformBunch::edge_field(Cooler& cooler, vector<double>&x, vector<double>& y, vector<double>&z,
                             vector<double>& field, int n, double cx, double cy, double cz){
    double g = 1 + 2*log(cooler.pipe_radius()/radius_);
    double v = beta_*k_c;
    double dz_rising = t_rising_*v;
    double dz_falling = -t_falling_*v;
    double fld_rising = 0;
    double fld_falling = 0;
    double coef = -g*k_ke*current_/(v*gamma_*gamma_);
    if(t_rising_>0) fld_rising =coef/dz_rising;
    if(t_falling_>0) fld_falling = coef/dz_falling;
    double left_end = -0.5*length_;
    double falling_end = left_end + dz_falling;
    double right_end = 0.5*length_;
    double rising_end = right_end + dz_rising;
    std::fill(field.begin(), field.end(), 0);
    double r2 = radius_*radius_;
    //ion_center - electron_center
    cx -= this->center(0);
    cy -= this->center(1);
    cz -= this->center(2);
    for(int i=0; i<n; ++i) {
        if((x[i]+cx)*(x[i]+cx)+(y[i]+cy)*(y[i]+cy)<=r2) {
            if((z[i]+cz)>falling_end && (z[i]+cz)<left_end) field.at(i) = fld_falling;
            else if((z[i]+cz)>right_end && (z[i]+cz)<rising_end) field.at(i) = fld_rising;
        }
    }
}

void UniformBunch::edge_field(Cooler& cooler, vector<double>&x, vector<double>& y, vector<double>&z,
                             vector<double>& field, int n){
    double g = 1 + 2*log(cooler.pipe_radius()/radius_);
    double v = beta_*k_c;
    double dz_rising = t_rising_*v;
    double dz_falling = -t_falling_*v;
    double fld_rising = 0;
    double fld_falling = 0;
    double coef = -g*k_ke*current_/(v*gamma_*gamma_);
    if(t_rising_>0) fld_rising =coef/dz_rising;
    if(t_falling_>0) fld_falling = coef/dz_falling;
    double left_end = -0.5*length_;
    double falling_end = left_end + dz_falling;
    double right_end = 0.5*length_;
    double rising_end = right_end + dz_rising;
    std::fill(field.begin(), field.end(), 0);
    double r2 = radius_*radius_;
    for(int i=0; i<n; ++i) {
        if(x[i]*x[i]+y[i]*y[i]<=r2) {
            if(z[i]>falling_end && z[i]<left_end) field.at(i) = fld_falling;
            else if(z[i]>right_end && z[i]<rising_end) field.at(i) = fld_rising;
        }
    }
}

void EBeam::multi_edge_field(Cooler& cooler, vector<double>&x, vector<double>& y, vector<double>&z,
                             vector<double>& field, int n) {
    vector<double> d(n);
    for(int i=0; i<n_; ++i) {
        edge_field(cooler, x, y, z, d, n, -cx_.at(i), -cy_.at(i), -cz_.at(i));
        for(int j=0; j<n; ++j) field.at(j) += d.at(j);
    }
}
void EBeam::multi_edge_field(Cooler& cooler, vector<double>&x, vector<double>& y, vector<double>&z,
                             vector<double>& field, int n, double cx, double cy, double cz) {
    vector<double> d(n);
    for(int i=0; i<n_; ++i) {
        edge_field(cooler, x, y, z, d, n, cx-cx_.at(i), cy-cy_.at(i), cz-cz_.at(i));
        for(int j=0; j<n; ++j) field.at(j) += d.at(j);
    }
}

void EBeam::create_samples(int n_sample) {
    std::string filename = "no_such_file";
    double length = this->length();
    samples.reset(new ParticleBunch(n_electron(), filename, length));
    samples->set_gamma(gamma_);
    samples->set_n_sample(n_sample);
    samples->set_corr();
    create_particle_location();
    if(shape()==Shape::PARTICLE_BUNCH) {
        ParticleBunch* this_ptr = dynamic_cast<ParticleBunch*>(this);
        std::copy(samples->vx.begin(), samples->vx.end(), this_ptr->vx.begin());
        std::copy(samples->vy.begin(), samples->vy.end(), this_ptr->vy.begin());
        std::copy(samples->vz.begin(), samples->vz.end(), this_ptr->vz.begin());
    }
    else {
        create_particle_velocity();
    }
    adjust_particle_location();
    samples->load_particle(-1);
}

void EBeam::create_particle_velocity(){
    int n_sample = samples->n_sample();
    gaussian_random(n_sample, samples->vx, v_rms_t.at(0));
    gaussian_random(n_sample, samples->vz, v_rms_l.at(0));
    gaussian_random_adjust(n_sample, samples->vx, v_rms_t.at(0));
    gaussian_random_adjust(n_sample, samples->vz, v_rms_l.at(0));
    uniform_random(n_sample, samples->vy, 0, 2*k_pi);
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for(int i=0; i<n_sample; ++i) {
        double vy = samples->vx.at(i)*sin(samples->vy.at(i));
        samples->vy.at(i) = vy;
        samples->vx.at(i) = sqrt(samples->vx.at(i)*samples->vx.at(i) - vy*vy);
    }
}

void EBeam::adjust_particle_location() {
    int n_sample = samples->n_sample();
    if (!iszero(center_[0])) {
        for(auto& x: samples->x) x += center_[0];
    }
    if (!iszero(center_[1])) {
        for(auto& y: samples->y) y += center_[1];
    }
    if (!iszero(center_[2])) {
        for(auto& z: samples->z) z += center_[2];
    }
    if (disp_) {
        std::vector<double> dp_p(n_sample);
        double k = 1/(beta_*k_c);
        #ifdef _OPENMP
            #pragma omp parallel for
        #endif // _OPENMP
        for(int i=0; i<n_sample; ++i) {
            dp_p.at(i) = k*samples->vz.at(i);
        }
        if (!iszero(dx_)) adjust_disp(dx_, samples->x, dp_p, samples->x, n_sample);
        if (!iszero(dy_)) adjust_disp(dy_, samples->y, dp_p, samples->y, n_sample);
    }
}


