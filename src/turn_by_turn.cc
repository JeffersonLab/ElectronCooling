#include "turn_by_turn.h"
#include <chrono>
#include <cmath>
#include "beam.h"
#include "constants.h"
#include "ecooling.h"
#include "functions.h"
#include "particle_model.h"
#include "ring.h"

void TurnByTurnModel::move_particles(Beam& ion, Ions& ion_sample, Ring& ring) {
    //Transverse
    //New betatron oscillation coordinates
    auto twiss = ion_sample.get_twiss();
//    double dx = twiss.disp_x;
//    double dpx = twiss.disp_dx;
//    double dy = twiss.disp_y;
//    double dpy = twiss.disp_dy;

    int n_sample = ion_sample.n_sample();
    ion_sample.adjust_disp_inv();

    //Transverse motion by tunes
    assert(ring.tunes.qx>0&&ring.tunes.qy>0&&"Transverse tunes are needed for Turn_by_turn model");
    double Qx = ring.tunes.qx;
    double Qy = ring.tunes.qy;
    auto x_bet = ion_sample.cdnt(Phase::X_BET);
    auto xp_bet = ion_sample.cdnt(Phase::XP_BET);
    auto y_bet = ion_sample.cdnt(Phase::Y_BET);
    auto yp_bet = ion_sample.cdnt(Phase::YP_BET);
    double bet_x = twiss.bet_x;
    double bet_y = twiss.bet_y;
    double alf_x = twiss.alf_x;
    double alf_y = twiss.alf_y;
    double gamma_x = (1+alf_x*alf_x)/bet_x;
    double gamma_y = (1+alf_y*alf_y)/bet_y;
    for (int i=0; i<n_sample; ++i) {
        double phi = 2*k_pi*Qx;
        double x_bet_0 = x_bet[i];
        double xp_bet_0 = xp_bet[i];
        x_bet[i] = (cos(phi)+alf_x*sin(phi))*x_bet_0 + bet_x*sin(phi)*xp_bet_0;
        xp_bet[i] = -gamma_x*sin(phi)*x_bet_0 + (cos(phi)-alf_x*sin(phi))*xp_bet_0;
        phi = 2*k_pi*Qy;
        double y_bet_0 = y_bet[i];
        double yp_bet_0 = yp_bet[i];
        y_bet[i] = (cos(phi)+alf_y*sin(phi))*y_bet_0 + bet_y*sin(phi)*yp_bet_0;
        yp_bet[i] = -gamma_y*sin(phi)*y_bet_0 + (cos(phi)-alf_y*sin(phi))*yp_bet_0;
    }

    //Longitudinal motion.
    auto dp_p = ion_sample.cdnt(Phase::DP_P);
    auto ds = ion_sample.cdnt(Phase::DS);
    if (ring.tunes.qs>0||ring.rf.v>0) {    //RF, synchrotron oscillation.
//        assert(ring.tunes->qs>0||ring.rf->v>0&&"Longitudinal tune or RF cavity needed for Turn_by_turn model");

        if(ring.rf.v>0) { //Longitudinal motion by RF.
            double circ = ring.circ();
            double beta2 = ion.beta()*ion.beta();
            double beta2_inv = 1/beta2;
            double total_energy = ion.gamma()*ion.mass(); //ion total energy [MeV/c^2]
            double total_energy_inv = 1/total_energy;
            double adj_dp2dE = beta2*total_energy;

            double volt = ring.rf.v;
            double phi_s = ring.rf.phi;
//            double phi_0 = ring.rf->phi_0();
            double h = ring.rf.h;
//            double s_s = phi_s*circ/(h*2*k_pi);
            double half_phase = h*k_pi;
            double total_phase = h*2*k_pi;
            double adj_s2phi = total_phase/circ;
            double adj_phi2s = 1/adj_s2phi;
            double adj_dE = ion.charge_number()*volt*1e-6; // [MeV/c^2]
    //                double adj_dE = ion.charge_number()*ring.rf_->volt()*1e-6; // [MeV/c^2]
            double sin_phi_s = sin(phi_s);
//            double sin_phi_s = sin(phi_s+phi_0);
            double eta = 1/(ring.rf.gamma_tr*ring.rf.gamma_tr) - 1/(ion.gamma()*ion.gamma()); //phase slip factor
            double adj_dE2dphi = total_phase*eta*beta2_inv*total_energy_inv;
            for(int i = 0; i < n_sample; ++i) {
                dp_p[i] *= adj_dp2dE; //dp/p -> dE/E -> dE in [MeV/c^2]
//                ds[i] += s_s;  //s = ds + s_s: adjust ds to be measured from the start of the ring
                ds[i] *= adj_s2phi;  //phi = s*h*2*pi/circ: s -> phi
//                dp_p[i] += adj_dE*(sin(ds[i]+phi_0)-sin_phi_s); //dE_n+1 = dE_n + q*v*[sin(phi)-sin(phi_s)]
                dp_p[i] += adj_dE*(sin(ds[i])-sin_phi_s); //dE_n+1 = dE_n + q*v*[sin(phi)-sin(phi_s)]
                ds[i] += adj_dE2dphi*dp_p[i]; //phi_n+1 = phi_n + 2*pi*h*eta*dE_n+1 / (beta^2 * E)
                ds[i] += half_phase;    //phi in [0, total phase]
                ds[i] = fmod(ds[i], total_phase);  //phi_n+1 module h*2*pi
                if (ds[i]<0) ds[i] += total_phase; //adjust phi is phi is less than zero
                ds[i] -= half_phase;    //phi in [-half_phase, half_phase]
                ds[i] *= adj_phi2s;  //phi -> s
//                ds[i] -= s_s;        // ds = s - s_s: adjust ds back to be centered about s_s
                dp_p[i]  = dp_p[i]*total_energy_inv*beta2_inv; //dE -> dE/E -> dp/p = beta*beta*dE/E;
            }
        }
        else if(ring.tunes.qs>0) {//Longitudinal motion by tune
            double phi = 2*k_pi*ring.tunes.qs;
            double inv_beta_s = 1/ring.beta_s();
            double beta_s = ring.beta_s();
            for (int i=0; i<n_sample; ++i) {
                double dp_p_0 = dp_p[i];
                double ds_0 = ds[i];
                dp_p[i] = cos(phi)*dp_p_0 - sin(phi)*ds_0*inv_beta_s;
                ds[i] = sin(phi)*dp_p_0*beta_s + cos(phi)*ds_0;
            }
        }
    }
    else {  //No RF.
        double gamma_0 = ion.gamma();
        double beta_0 = ion.beta();
        double half_length = 0.5*ring.circ();
        for(int i=0; i<n_sample; ++i) {
            double gamma2 = 1+(1+dp_p[i])*(1+dp_p[i])*(gamma_0*gamma_0-1);
            double beta = sqrt(1-1/gamma2);
            double s = (beta/beta_0-1)*2*half_length;
            ds[i] += half_length;    //s in [0, 2*half_length]
            ds[i] += s;
            ds[i] = fmod(ds[i], 2*half_length);
            if(ds[i]<0) ds[i] += 2*half_length;
            ds[i] -= half_length;   //s in [-half_length, half_length]
        }
    }

    ion_sample.adjust_disp();
}
