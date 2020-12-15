#include "ions.h"
#include <assert.h>
#include <cmath>
#include <fstream>

#include "arbitrary_electron_beam.h"
#include "functions.h"
#include <chrono>



void Ions::set_twiss(Twiss& t) {
    twiss.bet_x = t.bet_x;
    twiss.bet_y = t.bet_y;
    twiss.alf_x = t.alf_x;
    twiss.alf_y = t.alf_y;
    twiss.disp_x = t.disp_x;
    twiss.disp_y = t.disp_y;
    twiss.disp_dx = t.disp_dx;
    twiss.disp_dy = t.disp_dy;
}

void Ions::set_twiss(Cooler& cooler) {
    twiss.bet_x = cooler.beta_h();
    twiss.bet_y = cooler.beta_v();
    twiss.alf_x = cooler.alpha_h();
    twiss.alf_y = cooler.alpha_v();
    twiss.disp_x = cooler.disp_h();;
    twiss.disp_y = cooler.disp_v();;
    twiss.disp_dx = cooler.der_disp_h();
    twiss.disp_dy = cooler.der_disp_v();
}

const vector<double>& Ions::get_cdnt(Phase p) const {
    switch(p) {
        case Phase::X_BET: {
            return x_bet;
        }
        case Phase::XP_BET: {
            return xp_bet;
        }
        case Phase::Y_BET: {
            return y_bet;
        }
        case Phase::YP_BET: {
            return yp_bet;
        }
        case Phase::X: {
            return x;
        }
        case Phase::Y: {
            return y;
        }
        case Phase::XP: {
            return xp;
        }
        case Phase::DS: {
            return ds;
        }
        case Phase::DP_P: {
            return dp_p;
        }
        default: {
            assert(false&&"Wrong phase coordinate selected!");
            return x;
        }
    }
}

vector<double>& Ions::cdnt(Phase p) {
    switch(p) {
        case Phase::X_BET: {
            return x_bet;
        }
        case Phase::XP_BET: {
            return xp_bet;
        }
        case Phase::Y_BET: {
            return y_bet;
        }
        case Phase::YP_BET: {
            return yp_bet;
        }
        case Phase::X: {
            return x;
        }
        case Phase::Y: {
            return y;
        }
        case Phase::XP: {
            return xp;
        }
        case Phase::YP: {
            return yp;
        }
        case Phase::DS: {
            return ds;
        }
        case Phase::DP_P: {
            return dp_p;
        }
        default: {
            assert(false&&"Wrong phase coordinate selected!");
            return x;
        }
    }
}

void Ions::save_ions_sdds(string filename) const {
    using std::endl;
    std::ofstream output_particles;
    output_particles.open(filename);
    output_particles<<"SDDS1"<<endl;
    output_particles<<"! Define colums:"<<endl
        <<"&column name=x, type=double, units=m, description=NULL, &end"<<endl
        <<"&column name=xp, type=double, units=NULL, description=NULL, &end"<<endl
        <<"&column name=y, type=double, units=m, description=NULL, &end"<<endl
        <<"&column name=yp, type=double, units=NULL, description=NULL, &end"<<endl
        <<"&column name=ds, type=double, units=m, description=NULL, &end"<<endl
        <<"&column name=dp/p, type=double, units=NULL, description=NULL, &end"<<endl
        <<"!Declare ASCII data and end the header"<<endl
        <<"&data mode=ascii, &end"<<endl
        <<n_<<endl;
    output_particles.precision(10);
    output_particles<<std::showpos;
    output_particles<<std::scientific;
    for(int i=0; i<n_; ++i) {
        output_particles<<x[i]<<' '<<xp[i]<<' '<<y[i]<<' '<<yp[i]<<' '<<ds[i]<<' '<<dp_p[i]<<std::endl;
    }
    output_particles.close();
}

//Calculate the transverse emittance statistically
double emit(vector<double>& x, vector<double>&xp, int n) {
    double emit, x_mean, xp_mean, dlt2_x, dlt2_xp, dlt_xxp;
    x_mean = 0;
    xp_mean = 0;

    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:x_mean,xp_mean)
    #endif // _OPENMP
    for(int i=0; i<n; ++i){
        x_mean += x[i];
       xp_mean += xp[i];
    }
    x_mean /= n;
    xp_mean /= n;

    dlt2_x = 0;
    dlt2_xp = 0;
    dlt_xxp = 0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:dlt2_x,dlt2_xp,dlt_xxp)
    #endif // _OPENMP
    for(int i=0; i<n; ++i){
        double x_adj = x[i]-x_mean;
        double xp_adj = xp[i]-xp_mean;
        dlt2_x += x_adj*x_adj;
        dlt2_xp += xp_adj*xp_adj;
        dlt_xxp += x_adj*xp_adj;
    }
    emit = sqrt(dlt2_x*dlt2_xp-dlt_xxp*dlt_xxp)/n;
    return emit;
}

//Calculate the longitudinal emittance as (dp/p)^2/n
double emit_p(vector<double>& dp_p, int n){
    double emit_p = 0;
    double dp_p_mean = 0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:dp_p_mean)
    #endif // _OPENMP
    for(int i=0; i<n; ++i){
        dp_p_mean += dp_p[i];
    }
    dp_p_mean /= n;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:emit_p)
    #endif // _OPENMP
    for(int i=0; i<n; ++i){
        double dp_p_adj = dp_p[i] - dp_p_mean;
        emit_p += dp_p_adj*dp_p_adj;
    }
    emit_p /= n;
    return emit_p;
}

void adjust_disp_inv(double dx, vector<double>& x_bet, vector<double>& dp_p, vector<double>& x, int n) {
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for(int i=0; i<n; ++i) x_bet[i] = x[i]-dx*dp_p[i];
}

void adjust_disp(double dx, vector<double>& x_bet, vector<double>& dp_p, vector<double>& x, int n){
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for(int i=0; i<n; ++i) x[i] = x_bet[i]+dx*dp_p[i];
}

void Ions_MonteCarlo::emit(vector<double>& x_bet, vector<double>& xp_bet, vector<double>& y_bet, vector<double>& yp_bet,
                      vector<double>& dp_p, vector<double>& ds, double& emit_x, double& emit_y, double& emit_s) {
    int n_sample = n_;
    emit_x = ::emit(x_bet, xp_bet, n_sample);
    emit_y = ::emit(y_bet, yp_bet, n_sample);
    emit_s = ::emit_p(dp_p, n_sample);
    if(bunched_)
        emit_s += emit_p(ds, n_sample)/(beta_s_*beta_s_);
}

void Ions_MonteCarlo::emit(double& emit_x, double& emit_y, double& emit_s){
    emit(x_bet, xp_bet, y_bet, yp_bet, dp_p, ds, emit_x, emit_y, emit_s);
}

void Ions_SingleParticle::emit(vector<double>& x_bet, vector<double>& xp_bet, vector<double>& y_bet, vector<double>& yp_bet,
                      vector<double>& dp_p, vector<double>& ds, double& emit_x, double& emit_y, double& emit_s) {
    double alf_x = twiss.alf_x;
    double alf_y = twiss.alf_y;
    double beta_x = twiss.bet_x;
    double beta_y = twiss.bet_y;
    double gamma_x = (1+alf_x*alf_x)/beta_x;
    double gamma_y = (1+alf_y*alf_y)/beta_y;

    emit_x = 0;
    emit_y = 0;
    emit_s = 0;
    int n_sample = n_;
    double inv_beta_s2 = 0;
    if(bunched_) inv_beta_s2 = 1/(beta_s_*beta_s_);
    for(int i=0; i<n_sample; ++i) {
        emit_x += beta_x*xp_bet[i]*xp_bet[i]+2*alf_x*x_bet[i]*xp_bet[i]+gamma_x*x_bet[i]*x_bet[i];
        emit_y += beta_y*yp_bet[i]*yp_bet[i]+2*alf_y*y_bet[i]*yp_bet[i]+gamma_y*y_bet[i]*y_bet[i];
        emit_s += dp_p[i]*dp_p[i];
        if(bunched_) emit_s += ds[i]*ds[i]*inv_beta_s2;
    }
    emit_x /= 2*n_sample;
    emit_y /= 2*n_sample;
    emit_s /= n_sample;
}

void Ions_SingleParticle::emit(double& emit_x, double& emit_y, double& emit_s) {
    emit(x_bet, xp_bet, y_bet, yp_bet, dp_p, ds, emit_x, emit_y, emit_s);
}

//Generate Gaussian random number in S frame with given Twiss parameters
//First, rotate to O frame where alf = 0;
//Second, Generate x and xp with Gaussian random number in O frame
//Third, rotate back to S frame
int gaussian_bet_cod(double beta_xs, double alf_xs, double emit_x, vector<double>& x_bet, vector<double>& xp_bet, int n){

    double gamma_xs = (1+alf_xs*alf_xs)/beta_xs;
    double theta = atan(2*alf_xs/(gamma_xs-beta_xs))/2;     //rotation angle between O frame and S frame

    //Transfer matrix between O and S frames
    double matrix_os[2][2], matrix_so[2][2];
    matrix_os[0][0] = cos(theta);
    matrix_os[0][1] = -sin(theta);
    matrix_os[1][0] = sin(theta);
    matrix_os[1][1] = cos(theta);
    matrix_so[0][0] = matrix_os[0][0];
    matrix_so[0][1] = -matrix_os[0][1] ;
    matrix_so[1][0] = -matrix_os[1][0];
    matrix_so[1][1] = matrix_os[1][1];

    //Calculate beta and sigma in O frame
    double beta_xo = matrix_so[0][0]*matrix_so[0][0] * beta_xs-2*matrix_so[0][0]*matrix_so[0][1]*alf_xs+
                     matrix_so[0][1]*matrix_so[0][1]*gamma_xs;
    double sigma_xo = sqrt(emit_x*beta_xo);
    double sigma_xpo = sqrt(emit_x/beta_xo);
    //Generate x and xp in O frame
    gaussian_random(n, x_bet, sigma_xo);
    gaussian_random(n, xp_bet, sigma_xpo);
    gaussian_random_adjust(n, x_bet, sigma_xo);
    gaussian_random_adjust(n, xp_bet, sigma_xpo);

    //Rotate back to S frame
    for(int i=0; i<n;++i){
        double x = matrix_os[0][0]*x_bet[i]+matrix_os[0][1]*xp_bet[i];
        double xp = matrix_os[1][0]*x_bet[i]+matrix_os[1][1]*xp_bet[i];
        x_bet[i] = x;
        xp_bet[i] = xp;
    }
    return 0;
}

void Ions::adjust_disp(){
    ::adjust_disp(twiss.disp_x, x_bet, dp_p, x, n_);
    ::adjust_disp(twiss.disp_y, y_bet, dp_p, y, n_);
    ::adjust_disp(twiss.disp_dx, xp_bet, dp_p, xp, n_);
    ::adjust_disp(twiss.disp_dy, yp_bet, dp_p, yp, n_);
}

void Ions::adjust_disp_inv(){
    ::adjust_disp_inv(twiss.disp_x, x_bet, dp_p, x, n_);
    ::adjust_disp_inv(twiss.disp_y, y_bet, dp_p, y, n_);
    ::adjust_disp_inv(twiss.disp_dx, xp_bet, dp_p, xp, n_);
    ::adjust_disp_inv(twiss.disp_dy, yp_bet, dp_p, yp, n_);
}

Ions_MonteCarlo::Ions_MonteCarlo(int n_sample){
    n_=n_sample;
    x_bet.resize(n_sample,0);
    y_bet.resize(n_sample,0);
    xp_bet.resize(n_sample,0);
    yp_bet.resize(n_sample,0);
    ds.resize(n_sample,0);
    dp_p.resize(n_sample,0);
    x.resize(n_sample,0);
    y.resize(n_sample,0);
    xp.resize(n_sample,0);
    yp.resize(n_sample,0);
}

Ions_MonteCarlo::Ions_MonteCarlo(std::string filename, int n, int skip, bool binary, int n_buffer) {
    auto n_loaded = load_electrons(x, xp, y, yp, ds, dp_p, filename, n, skip, binary,n_buffer);
    if (n_loaded!=n_) n_ = n_loaded;
}

void Ions_MonteCarlo::create_samples(Beam& ion) {
    double beta_xs = twiss.bet_x;
    double beta_ys = twiss.bet_y;
    double alf_xs = twiss.alf_x;
    double alf_ys = twiss.alf_y;
    double emit_x = ion.emit_x();
    double emit_y = ion.emit_y();

    int n_sample = n_;

    gaussian_bet_cod(beta_xs, alf_xs, emit_x, x_bet, xp_bet, n_sample);
    gaussian_bet_cod(beta_ys, alf_ys, emit_y, y_bet, yp_bet, n_sample);

    double sigma_p = ion.dp_p();
    gaussian_random(n_sample, dp_p, sigma_p);
    gaussian_random_adjust(n_sample, dp_p, sigma_p);

    //longitudinal sampling
    if(ion.bunched()) {
        double sigma_s = ion.sigma_s();
        gaussian_random(n_sample, ds, sigma_s);
        gaussian_random_adjust(n_sample, ds, sigma_s);
    }

    double dx = twiss.disp_x;
    double dy = twiss.disp_y;
    double dpx = twiss.disp_dx;
    double dpy = twiss.disp_dy;
    ::adjust_disp(dx, x_bet, dp_p, x, n_sample);
    ::adjust_disp(dy, y_bet, dp_p, y, n_sample);
    ::adjust_disp(dpx, xp_bet, dp_p, xp, n_sample);
    ::adjust_disp(dpy, yp_bet, dp_p, yp, n_sample);

    bunched_ = ion.bunched();
    if(bunched_) beta_s_ = ion.sigma_s()/ion.dp_p();
}

Ions_SingleParticle::Ions_SingleParticle(int n_tr, int n_l):n_tr_(n_tr),n_l_(n_l) {
    n_=n_tr*n_tr*n_l;
    x_spl.resize(n_tr_);
    y_spl.resize(n_tr_);
    xp_spl.resize(n_tr_);
    yp_spl.resize(n_tr_);
    dp_p_spl.resize(n_l_);
    ds_spl.resize(n_l_);

    x_bet.resize(n_);
    y_bet.resize(n_);
    xp_bet.resize(n_);
    yp_bet.resize(n_);
    ds.resize(n_);
    dp_p.resize(n_);
    x.resize(n_);
    y.resize(n_);
    xp.resize(n_);
    yp.resize(n_);

};

void Ions_SingleParticle::single_particle_grid(Beam &ion){


    double alf_x = twiss.alf_x;
    double alf_y = twiss.alf_y;
    double dphi = 2.0*k_pi/n_tr_;
    double phi = 0;
    for(int i=0; i<n_tr_; ++i){
        x_spl[i] = sin(phi);
        y_spl[i] = sin(phi);
        xp_spl[i] = cos(phi)-alf_x*sin(phi);
        yp_spl[i] = cos(phi)-alf_y*sin(phi);
        phi += dphi;
    }

    if(ion.bunched()){
        phi = 0;
        dphi = 2.0*k_pi/n_l_;
        for(int i=0; i<n_l_; ++i){
            ds_spl[i] = sin(phi);
            dp_p_spl[i] = cos(phi);
            phi += dphi;
        }
    }
}

void Ions_SingleParticle::create_samples(Beam& ion) {

    double emit_x = ion.emit_x();
    double emit_y = ion.emit_y();
    double sigma_p = ion.dp_p();
    double beta_x = twiss.bet_x;
    double beta_y = twiss.bet_y;
    double dx = twiss.disp_x;
    double dy = twiss.disp_y;
    double dpx = twiss.disp_dx;
    double dpy = twiss.disp_dy;

    double y_amp = sqrt(2.0*emit_y*beta_y);
    double yp_amp = sqrt(2.0*emit_y/beta_y);
    double x_amp = sqrt(2.0*emit_x*beta_x);
    double xp_amp = sqrt(2.0*emit_x/beta_x);

    double ds_amp, dp_amp;
    if(ion.bunched()){  //bunched beam
        double sigma_s = ion.sigma_s();
        ds_amp = sqrt(2.0)*sigma_s;
        dp_amp = sqrt(2.0)*sigma_p;
    }

    int cnt = 0;
    for(int i=0; i<n_tr_; ++i){
        double y_spl_tmp = y_amp*y_spl[i];
        double yp_spl_tmp = yp_amp*yp_spl[i];
        for(int j=0; j<n_tr_; ++j){
            double x_spl_tmp = x_amp*x_spl[j];
            double xp_spl_tmp = xp_amp*xp_spl[j];
            if(ion.bunched()){  //bunched beam
                for(int k=0; k<n_l_; ++k){
                    double ds_spl_tmp = ds_amp*ds_spl[k];
                    double dp_spl_tmp = dp_amp*dp_p_spl[k];
                    x_bet[cnt] = x_spl_tmp;
                    xp_bet[cnt] = xp_spl_tmp;
                    y_bet[cnt] = y_spl_tmp;
                    yp_bet[cnt] = yp_spl_tmp;
                    ds[cnt] = ds_spl_tmp;
                    dp_p[cnt] = dp_spl_tmp;
                    ++cnt;
                }
            }
            else{   //coasting beam, ds=s-s0 is set to be zero!
                for(int k=-1; k<3; k +=2){
                    x_bet[cnt] = x_spl_tmp;
                    xp_bet[cnt] = xp_spl_tmp;
                    y_bet[cnt] = y_spl_tmp;
                    yp_bet[cnt] = yp_spl_tmp;
                    dp_p[cnt] = k*sigma_p;
                    ++cnt;
                }
            }
        }
    }

    ::adjust_disp(dx, x_bet, dp_p, x, cnt);
    ::adjust_disp(dy, y_bet, dp_p, y, cnt);
    ::adjust_disp(dpx, xp_bet, dp_p, xp, cnt);
    ::adjust_disp(dpy, yp_bet, dp_p, yp, cnt);

    bunched_ = ion.bunched();
    if(bunched_) beta_s_ = ion.sigma_s()/ion.dp_p();
}
