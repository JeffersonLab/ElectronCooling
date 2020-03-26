#include "force.h"
#include <cmath>
#include <cstdio>

#include <iostream>
#include <fstream>

void rho_lamor_dlt2_eff_e(double v2_eff_e, double mag_field, vector<double>& v_rms_l, vector<double>& v_rms_t, Temperature tpr,
                          int ion_number, vector<double>& dlt2_eff_e, vector<double>& rho_lamor) {
    switch(tpr) {
        case Temperature::CONST: {
            dlt2_eff_e.push_back(v_rms_l[0]*v_rms_l[0]+v2_eff_e);
            rho_lamor.push_back(k_me*1e6*v_rms_t[0]/(mag_field*k_c*k_c));
            break;
        }
        case Temperature::VARY: {
            dlt2_eff_e.resize(ion_number);
            rho_lamor.resize(ion_number);
            for(int i=0; i<ion_number; ++i) {
                dlt2_eff_e.at(i) = v_rms_l[i]*v_rms_l[i]+v2_eff_e;
                rho_lamor.at(i) = k_me*1e6*v_rms_t[i]/(mag_field*k_c*k_c);
            }
            break;
        }
        case Temperature::USER_DEFINE: {
            break;
        }
        case Temperature::SPACE_CHARGE: {
            break;
        }
        default: {
            break;
        }
    }

}

double dlt(Temperature tpr,  double v2, vector<double>& dlt2_eff_e, int i) {
    double dlt = 0;
    switch(tpr) {
        case Temperature::CONST: {
            dlt = v2+dlt2_eff_e[0];
            break;
        }
        case Temperature::VARY: {
            dlt = v2+dlt2_eff_e[i];
            break;
        }
        case Temperature::USER_DEFINE: {
            break;
        }
        case Temperature::SPACE_CHARGE: {
            break;
        }
        default: {
            break;
        }
    }
    return dlt;
}

double lc(Temperature tpr, double rho_max, double rho_min, vector<double>& rho_lamor, int i) {
    double lc = 0;
    switch(tpr) {
        case Temperature::CONST: {
            lc = log((rho_max+rho_min+rho_lamor[0])/(rho_min+rho_lamor[0]));
            break;
        }
        case Temperature::VARY: {
            lc = log((rho_max+rho_min+rho_lamor[i])/(rho_min+rho_lamor[i]));
            break;
        }
        case Temperature::USER_DEFINE: {
            break;
        }
        case Temperature::SPACE_CHARGE: {
            break;
        }
        default: {
            break;
        }
    }
    return lc;
}

void Force_Park::friction_force(int charge_number, int ion_number, vector<double>& v_tr, vector<double>& v_long,vector<double>& density_e,
                EBeam& ebeam, vector<double>& force_tr, vector<double>& force_long) {
    double f_const = -4*charge_number*charge_number*k_c*k_c*k_ke*k_ke*k_e*k_e*k_e/(k_me*1e6);
    double v2_eff_e = t_eff*k_c*k_c/(k_me*1e6);

    double wp_const = 4*k_pi*k_c*k_c*k_e*k_ke/(k_me*1e6);
    double rho_min_const = charge_number*k_e*k_ke*k_c*k_c/(k_me*1e6);

    vector<double> dlt2_eff_e;
    vector<double> rho_lamor;
    auto tpr = ebeam.temperature();
    rho_lamor_dlt2_eff_e(v2_eff_e, mag_field, ebeam.get_v(EBeamV::V_RMS_L), ebeam.get_v(EBeamV::V_RMS_TR), tpr, ion_number,
        dlt2_eff_e, rho_lamor);

    force_tr.resize(ion_number);
    force_long.resize(ion_number);

    for(int i=0; i<ion_number; ++i){
        double v2 = v_tr[i]*v_tr[i]+v_long[i]*v_long[i];
        if(v2>0){
            double dlt = ::dlt(tpr, v2, dlt2_eff_e, i);

            //Calculate rho_min
            double rho_min = rho_min_const/dlt;
            dlt = sqrt(dlt);
            double wp = sqrt(wp_const*density_e[i]);

            //Calculate rho_max
            double rho_max = dlt/wp;
            double rho_max_2 = pow(3*charge_number/density_e[i], 1.0/3);
            if(rho_max<rho_max_2) rho_max = rho_max_2;
            double rho_max_3 = dlt*time_cooler;
            if(rho_max>rho_max_3) rho_max = rho_max_3;

            double lc = ::lc(tpr, rho_max, rho_min, rho_lamor, i);   //Coulomb Logarithm
            //Calculate friction force
            double f = f_const*density_e[i]*lc/(dlt*dlt*dlt);
            force_tr[i] = f*v_tr[i];
            force_long[i] = f*v_long[i];
        }
        else{
            force_tr[i] = 0;
            force_long[i] = 0;
        }
    }
}
