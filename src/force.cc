#include "force.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>



void ForcePark::rho_lamor_dlt2_eff_e(double v2_eff_e, double mag_field, vector<double>& v_rms_l, vector<double>& v_rms_t, Temperature tpr,
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

double ForcePark::dlt(Temperature tpr,  double v2, vector<double>& dlt2_eff_e, int i) {
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

double ForcePark::lc(Temperature tpr, double rho_max, double rho_min, vector<double>& rho_lamor, int i) {
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

void ForcePark::friction_force(int charge_number, int ion_number, vector<double>& v_tr, vector<double>& v_long,vector<double>& density_e,
                EBeam& ebeam, vector<double>& force_tr, vector<double>& force_long) {
    double f_const = charge_number*charge_number*k_f;
    double v2_eff_e = v_eff*v_eff;

    double wp_const = k_wp;
    double rho_min_const = charge_number*k_rho_min;

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
            double dlt = this->dlt(tpr, v2, dlt2_eff_e, i);

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

            double lc = this->lc(tpr, rho_max, rho_min, rho_lamor, i);   //Coulomb Logarithm
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

double ForceNonMag::rho_max(int charge_number, double v2, double ve2, double ne) {
    double rho_max = 0;
    if(smooth_rho_max) {
        rho_max = sqrt((v2+ve2)/(k_wp*ne));
    }
    else {
        if(v2>ve2) rho_max = sqrt(v2/(k_wp*ne));
        else rho_max = sqrt(ve2/(k_wp*ne));
    }
    double rho_max_1 = this->rho_max_1(charge_number, ne);
    double rho_max_2 = this->rho_max_2(sqrt(v2));
    if(rho_max < rho_max_1) rho_max = rho_max_1;
    if(rho_max > rho_max_2) rho_max = rho_max_2;
    return rho_max;
}

void ForceNonMag::friction_force(int charge_number, int ion_number, vector<double>& v_tr, vector<double>& v_long,vector<double>& density_e,
                EBeam& ebeam, vector<double>& force_tr, vector<double>& force_long) {
    init();
    double f_const = this->f_const(charge_number);
    double rho_min_const = this->rho_min_const(charge_number);
    auto tpr = ebeam.temperature();

    vector<double>& ve_rms_l = ebeam.get_v(EBeamV::V_RMS_L);
    vector<double>& ve_rms_tr = ebeam.get_v(EBeamV::V_RMS_TR);

    force_tr.resize(ion_number);
    force_long.resize(ion_number);

    switch(tpr) {
    case Temperature::CONST: {
        double ve_l = ve_rms_l.at(0);
        double ve_tr = ve_rms_l.at(0);
        double ve2 = ve_l*ve_l + ve_tr*ve_tr;

        for(int i=0; i<ion_number; ++i) {
            double v2 = v_tr[i]*v_tr[i]+v_long[i]*v_long[i];
            double v = sqrt(v2);
            if(v2>0) {
                 force(v, v_tr[i], v_long[i], v2, ve_tr, ve_l, ve2, f_const, rho_min_const,  charge_number,
                       density_e.at(i), force_tr[i], force_long[i]);
            }
            else {
                force_tr[i] = 0;
                force_long[i] = 0;
            }
        }
        break;
    }
    case Temperature::VARY: {
        for(int i=0; i<ion_number; ++i) {
            double v2 = v_tr[i]*v_tr[i]+v_long[i]*v_long[i];
            double v = sqrt(v2);
            if(v2>0) {
                double ve_l = ve_rms_l.at(i);
                double ve_tr = ve_rms_l.at(i);
                double ve2 = ve_l*ve_l + ve_tr*ve_tr;
                force(v, v_tr[i], v_long[i], v2, ve_tr, ve_l, ve2, f_const, rho_min_const, charge_number,
                      density_e.at(i), force_tr[i], force_long[i]);
            }
            else {
                force_tr[i] = 0;
                force_long[i] = 0;
            }
        }
        break;
    }
    default: {
        break;
    }
    }
    fin();
}


void ForceNonMagMeshkov::force(double v, double v_tr, double v_l, double v2, double ve_tr, double ve_l, double ve2,
                               double f_const, double rho_min_const, int charge_number, double ne,
                               double& force_tr, double& force_l) {
    double rho_min = rho_min_const/(v2+ve2);
    double rho_max = this->rho_max(charge_number, v2, ve2, ne);
    double lc = rho_max>rho_min?log(rho_max/rho_min):0;
    if(v >= ve_tr) {
        double coef = f_const*ne*lc/(v*v*v);
        force_tr = coef*v_tr;
        force_l = coef*v_l;
    }
    else if(v>=ve_l) {
        double coef = f_const*ne*lc/(ve_tr*ve_tr);
        force_tr = coef*v_tr/ve_tr;
        force_l = v_l>0?coef:-coef;
    }
    else {
        force_tr = 0;
        force_l = f_const*ne*lc*v_l/(ve_l*ve_tr*ve_tr);
    }
}

double ForceNonMagDerbenev::rho_max_ve_tr(double ve_tr, double ne) {
    return sqrt(ve_tr*ve_tr/(k_wp*ne));
}

void ForceNonMagDerbenev::force(double v, double v_tr, double v_l, double v2, double ve_tr, double ve_l, double ve2,
                               double f_const, double rho_min_const, int charge_number, double ne,
                               double& force_tr, double& force_l) {
    double rho_min = rho_min_const/(v2+ve2);
    double rho_max = this->rho_max(charge_number, v2, ve2, ne);
    double rho_max_ve_tr = this->rho_max_ve_tr(ve_tr, ne);
    double lc = rho_max>rho_min?log(rho_max/rho_min):0;
    double lc_ve_tr = rho_max_ve_tr>rho_min?log(rho_max_ve_tr/rho_min):0;
    if(v >= ve_tr) {
        double coef = f_const*ne*lc/(v*v*v);
        force_tr = coef*v_tr;
        force_l = coef*v_l - f_const*sqrt(2/k_pi)*lc_ve_tr/(v_l*v_l);
    }
    else if(v>=ve_l) {
        double coef = f_const*ne/(ve_tr*ve_tr);
        force_tr = coef*lc*v_tr/ve_tr;
        force_l = coef*v_l*(lc/sqrt(v_l*v_l+ve_l*ve_l) - sqrt(2/k_pi)*lc_ve_tr/ve_tr);
    }
    else {
        force_tr = 0;
        force_l = f_const*ne*lc*v_l/(sqrt(v_l*v_l+ve_l*ve_l)*ve_tr*ve_tr);
    }
}

double ForceNonMagNumeric1D::b(double q, void* params) {
    double v_tr = ((P*)params)->v_tr;
    double v_l = ((P*)params)->v_l;
    double ve_tr = ((P*)params)->ve_tr;
    double ve_l = ((P*)params)->ve_l;
    int flag =  ((P*)params)->flag;

    double ve2_tr = ve_tr*ve_tr;
    double bot = (ve_l*ve_l)/ve2_tr+q;
    double integrand = exp(-v_tr*v_tr/(2*ve2_tr*(1+q))-v_l*v_l/(2*ve2_tr*bot))/((1+q)*sqrt(bot));
    if(flag==0) return integrand/(1+q); //B_tr.
    else return integrand/bot;  //B_l;
}

void ForceNonMagNumeric1D::force(double v, double v_tr, double v_l, double v2, double ve_tr, double ve_l, double ve2,
                               double f_const, double rho_min_const, int charge_number, double ne,
                               double& force_tr, double& force_l) {
    double rho_min = rho_min_const/(v2+ve2);
    double rho_max = this->rho_max(charge_number, v2, ve2, ne);
    double lc = rho_max>rho_min?log(rho_max/rho_min):0;

    p.ve_tr = ve_tr;
    p.ve_l = ve_l;
    p.v_tr = v_tr;
    p.v_l = v_l;
    p.flag = 0;

    ForceNonMagNumeric1D* ptr2 = this;
    auto ptr = [=](double x)->double{return ptr2->b(x,&p);};
    gsl_function_pp<decltype(ptr)> fp(ptr);
    gsl_function *f = static_cast<gsl_function*>(&fp);

    double b_tr;
    double error;
    gsl_integration_qagiu(f, 0, espabs, esprel, limit, gw, &b_tr, &error);

    double b_l;
    p.flag = 1;
    gsl_integration_qagiu(f, 0, espabs, esprel, limit, gw, &b_l, &error);

    double ve3_tr = ve_tr*ve_tr*ve_tr;
    force_tr = f_const*ne*lc*v_tr*b_tr/ve3_tr;
    force_l = f_const*ne*lc*v_l*b_l/ve3_tr;
}

double ForceNonMagNumeric3D::inner_norm_integrand(double vl, void*params) {
    double vtr = ((P*)params)->vtr;
    double ve_tr = ((P*)params)->ve_tr;
    double ve_l = ((P*)params)->ve_l;
    return exp(-vtr*vtr/(2*ve_tr*ve_tr)-vl*vl/(2*ve_l*ve_l))*vtr;
}

double ForceNonMagNumeric3D::outter_norm_integrand(double vtr, void*params) {
    double result;
    double error;

    P* p =(P*)params;
    p->vtr = vtr;

    ForceNonMagNumeric3D* ptr2 = this;
    auto ptr = [=](double x)->double{return ptr2->inner_norm_integrand(x,params);};
    gsl_function_pp<decltype(ptr)> fp(ptr);
    gsl_function *f = static_cast<gsl_function*>(&fp);

    gsl_integration_qagi(f, espabs, esprel, limit, giw, &result, &error);
    return result;
}

double ForceNonMagNumeric3D::inner_integrand(double phi, void* params) {
    double vtr = ((P*)params)->vtr;
    double vl = ((P*)params)->vl;
    double ve_tr = ((P*)params)->ve_tr;
    double ve_l = ((P*)params)->ve_l;
    double v_tr = ((P*)params)->v_tr;
    double v_l = ((P*)params)->v_l;
    double rho_max = ((P*)params)->rho_max;
    int flag = ((P*)params)->flag;
    int charge_number = ((P*)params)->charge_number;

    double sub_vl = v_l - vl;
    double sub_vtr = v_tr - vtr*cos(phi);
    double f_inv_bot = 1/(sub_vl*sub_vl + sub_vtr*sub_vtr + vtr*vtr*sin(phi)*sin(phi));
    double rho_min = rho_min_const(charge_number)*f_inv_bot;
    double lc = rho_max>rho_min?log(rho_max/rho_min):0;
    double f = exp(-vtr*vtr/(2*ve_tr*ve_tr)-vl*vl/(2*ve_l*ve_l))*vtr*f_inv_bot;
    if(flag==0) return f*sub_vtr;
    else return f*sub_vl;
}

double ForceNonMagNumeric3D::middle_integrand(double vl, void* params) {

    P* p = (P*)params;
    p->vl = vl;

    ForceNonMagNumeric3D* class_ptr = this;
    auto ptr = [=](double x)->double{return class_ptr->inner_integrand(x,params);};
    gsl_function_pp<decltype(ptr)> fp(ptr);
    gsl_function *f = static_cast<gsl_function*>(&fp);

    int key = 1;
    double result;
    double error;
    gsl_integration_qag(f, 0, k_pi, espabs, esprel, limit, key, giw, &result, &error);
    return result;
}

double ForceNonMagNumeric3D::outter_integrand(double vtr, void* params) {

    P* p = (P*)params;
    p->vtr = vtr;

    ForceNonMagNumeric3D* class_ptr = this;
    auto ptr = [=](double x)->double{return class_ptr->middle_integrand(x,params);};
    gsl_function_pp<decltype(ptr)> fp(ptr);
    gsl_function *f = static_cast<gsl_function*>(&fp);

    double result;
    double error;
    gsl_integration_qagi(f, espabs, esprel, limit, gmw, &result, &error);
    return result;
}

void ForceNonMagNumeric3D::force(double v, double v_tr, double v_l, double v2, double ve_tr, double ve_l, double ve2,
                               double f_const, double rho_min_const, int charge_number, double ne,
                               double& force_tr, double& force_l) {
    double rho_max = this->rho_max(charge_number, v2, ve2, ne);


    p.ve_tr = ve_tr;
    p.ve_l = ve_l;

    ForceNonMagNumeric3D* class_ptr = this;
    auto norm_ptr = [=](double x)->double{return class_ptr->outter_norm_integrand(x,&p);};
    gsl_function_pp<decltype(norm_ptr)> norm_fp(norm_ptr);
    gsl_function *f = static_cast<gsl_function*>(&norm_fp);

    double norm;
    double error;
    gsl_integration_qagiu(f, 0, espabs, esprel, limit, gow, &norm, &error);
    norm *= k_pi;
    double inv_norm = 1/norm;

    p.v_tr = v_tr;
    p.v_l = v_l;
    p.charge_number = charge_number;
    p.rho_max = rho_max;

    auto ptr = [=](double x)->double{return class_ptr->outter_integrand(x,&p);};
    gsl_function_pp<decltype(ptr)> fp(ptr);
    f = static_cast<gsl_function*>(&fp);

    double f_tr;
    p.flag = 0;
    gsl_integration_qagiu(f, 0, espabs, esprel, limit, gow, &f_tr, &error);
    double f_l;
    p.flag = 1;
    gsl_integration_qagiu(f, 0, espabs, esprel, limit, gow, &f_l, &error);

    double ve3_tr = ve_tr*ve_tr*ve_tr;
    force_tr = f_const*ne*f_tr*inv_norm;
    force_l = f_const*ne*f_l*inv_norm;
}

