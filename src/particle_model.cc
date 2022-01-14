#include "particle_model.h"
#include <chrono>
#include <cmath>
#include "constants.h"
#include "ecooling.h"
#include "functions.h"

void ParticleModel::update_ibeam(Beam& ion, Ions& ion_sample, Ring& ring, EBeam& ebeam, Cooler& cooler, ECoolRate* ecool_solver) {
    vector<double>& dp_p = ion_sample.cdnt(Phase::DP_P);
    if(ecool) {
        double freq = k_c*ion.beta()/ring.circ()*cooler.section_number();
        ecool_solver->adjust_rate(ion, ebeam, {&freq});
        apply_cooling_kick(freq, ion, ion_sample, ecool_solver);
    }
    if(ibs) {
        apply_ibs_kick(ion, ion_sample);
    }

//    if(edge_effect) {
//        apply_edge_kick(cooler, ebeam, ion, ion_sample, ecool_solver);
//    }

    move_particles(ion, ion_sample, ring);
    update_beam_parameters(ion, ion_sample);

    if(fixed_bunch_length && ion.bunched()) {
        ring.update_bet_s();
        ring.rf.v = ring.calc_rf_voltage();
    }
}

void ParticleModel::apply_cooling_kick(double freq, Beam& ion, Ions& ion_sample, ECoolRate* ecool_solver) {
    vector<double>& xp = ion_sample.cdnt(Phase::XP);
    vector<double>& yp = ion_sample.cdnt(Phase::YP);
    vector<double>& dp_p = ion_sample.cdnt(Phase::DP_P);
    vector<double>& force_x = ecool_solver->scratch(ECoolRateScratch::FORCE_X);
    vector<double>& force_y = ecool_solver->scratch(ECoolRateScratch::FORCE_Y);
    vector<double>& force_z = ecool_solver->scratch(ECoolRateScratch::FORCE_Z);
    double p0 = ion.p0_SI();
    double t_cooler = ecool_solver->t_cooler();
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for(int i=0; i<ion_sample.n_sample(); ++i) {
        xp[i] = !iszero(xp[i])?xp[i]*exp(force_x[i]*t_cooler*dt*freq/(xp[i]*p0)):xp[i];
        yp[i] = !iszero(yp[i])?yp[i]*exp(force_y[i]*t_cooler*dt*freq/(yp[i]*p0)):yp[i];
        double dp = force_z[i]*t_cooler*dt*freq/(dp_p[i]*p0);
        if(!iszero(dp_p[i],1e-7)) {
            dp_p[i] = dp>0.15?dp_p[i]*(1+dp):dp_p[i]*exp(dp);
        }
//        dp_p[i] = !iszero(dp_p[i],1e-7)?dp_p[i]*exp(force_z[i]*t_cooler*dt*freq/(dp_p[i]*p0)):dp_p[i];
    }
}

void ParticleModel::apply_ibs_kick(Beam& ion, Ions& ion_sample) {
    auto twiss = ion_sample.get_twiss();
    assert(twiss.bet_x>0&& twiss.bet_y>0
           &&"TWISS parameters for the reference point not defined! Define twiss_ref.");
    ibs_kick(ion_sample.n_sample(), r_ibs.at(0), twiss.bet_x, ion.emit_x(), ion_sample.cdnt(Phase::XP));
    ibs_kick(ion_sample.n_sample(), r_ibs.at(1), twiss.bet_y, ion.emit_y(), ion_sample.cdnt(Phase::YP));
    if (ion.bunched()) ibs_kick(ion_sample.n_sample(), r_ibs.at(2), 1, ion.dp_p()*ion.dp_p(), ion_sample.cdnt(Phase::DP_P));
    else ibs_kick(ion_sample.n_sample(), r_ibs.at(2), 2, ion.dp_p()*ion.dp_p(), ion_sample.cdnt(Phase::DP_P));
}

void ParticleModel::ibs_kick(int n_sample, double rate, double twiss, double emit, vector<double>& v) {
    if(rate>0) {
        double theta = sqrt(2*rate*dt*emit/twiss);
        gaussian_random(n_sample, rdn, 1, 0);
        #ifdef _OPENMP
            #pragma omp parallel for
        #endif // _OPENMP
        for(int i=0; i<n_sample; ++i) v[i] += theta*rdn[i];
    }
    else {
        double k = exp(rate*dt);
        #ifdef _OPENMP
            #pragma omp parallel for
        #endif // _OPENMP
        for(int i=0; i<n_sample; ++i) v[i] *= k;
    }
}

void ParticleModel::move_particles(Beam& ion, Ions& ion_sample, Ring& ring) {
    //New betatron oscillation coordinates
    auto twiss = ion_sample.get_twiss();
    vector<double>& x_bet = ion_sample.cdnt(Phase::X_BET);
    vector<double>& xp_bet = ion_sample.cdnt(Phase::XP_BET);
    vector<double>& y_bet = ion_sample.cdnt(Phase::Y_BET);
    vector<double>& yp_bet = ion_sample.cdnt(Phase::YP_BET);
    vector<double>& dp_p = ion_sample.cdnt(Phase::DP_P);

    int n_sample = ion_sample.n_sample();
    ion_sample.adjust_disp_inv();

     //random phase advance
    double alf_x = twiss.alf_x;
    double alf_y = twiss.alf_y;
    double beta_x = twiss.bet_x;
    double beta_y = twiss.bet_y;

    double gamma_x = (1+alf_x*alf_x)/beta_x;
    double gamma_y = (1+alf_y*alf_y)/beta_y;

    uniform_random(n_sample, rdn, -1, 1);
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for(int i=0; i<n_sample; ++i){
        double I = beta_x*xp_bet[i]*xp_bet[i]+2*alf_x*x_bet[i]*xp_bet[i]+gamma_x*x_bet[i]*x_bet[i];
        double phi = k_pi*rdn[i];
        x_bet[i] = sqrt(I*beta_x)*sin(phi);
        xp_bet[i] = sqrt(I/beta_x)*(cos(phi)-alf_x*sin(phi));
    }
    uniform_random(n_sample, rdn, -1, 1);
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for(int i=0; i<n_sample; ++i){
        double I = beta_y*yp_bet[i]*yp_bet[i]+2*alf_y*y_bet[i]*yp_bet[i]+gamma_y*y_bet[i]*y_bet[i];
        double phi = k_pi*rdn[i];
        y_bet[i] = sqrt(I*beta_y)*sin(phi);
        yp_bet[i] = sqrt(I/beta_y)*(cos(phi)-alf_y*sin(phi));
    }

    if(ion.bunched()){
        uniform_random(n_sample, rdn, -1, 1);
        double beta_s = ring.beta_s();
        if(fixed_bunch_length) beta_s =  ion.sigma_s()/rms(n_sample, dp_p);
        double beta_s2_inv = 1/(beta_s*beta_s);
        vector<double>& ds = ion_sample.cdnt(Phase::DS);
        #ifdef _OPENMP
            #pragma omp parallel for
        #endif // _OPENMP
        for(int i=0; i<n_sample; ++i){
            double I = ds[i]*ds[i]*beta_s2_inv+dp_p[i]*dp_p[i];
            I = sqrt(I);
            double phi = k_pi*rdn[i];
            dp_p[i] = I*sin(phi);
            ds[i] = I*beta_s*cos(phi);
        }
        if(fixed_bunch_length) {
            gaussian_random_adjust(n_sample, ds, ion.sigma_s());
        }
    }
    ion_sample.adjust_disp();

}

void ParticleModel::update_beam_parameters(Beam &ion, Ions& ion_sample) {
    double emit_x, emit_y, emit_z;
    ion_sample.emit(emit_x, emit_y, emit_z);
    ion.set_emit_x(emit_x);
    ion.set_emit_y(emit_y);

    if(ion.bunched()) {
        if(fixed_bunch_length) {
            ion.set_dp_p(rms(ion_sample.n_sample(), ion_sample.cdnt(Phase::DP_P)));
            ion_sample.update_bet_s(ion);
        }
        else {
            ion.set_sigma_s(rms(ion_sample.n_sample(), ion_sample.cdnt(Phase::DS)));
            ion.set_dp_p(rms(ion_sample.n_sample(), ion_sample.cdnt(Phase::DP_P)));
        }
    }
    else {
        ion.set_dp_p(sqrt(emit_z));
    }
}

void ParticleModel::save_ions(int i, Ions& ion_sample) {
    if (ion_save_itvl>0 && i%ion_save_itvl==0) {
            std::size_t found = outfilename.find_last_of(".");
            if (found == string::npos)
                ion_sample.save_ions_sdds(outfilename+"_ions"+std::to_string(i)+".txt");
            else
                ion_sample.save_ions_sdds(outfilename.substr(0,found)+"_ions_"+std::to_string(i)+".txt");
    }
}

