#include "force.h"
#include <cmath>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <iostream>
#ifdef _OPENMP
    #include <omp.h>
#endif // _OPENMP

#include "functions.h"

double FrictionForceSolver::max_impact_factor(double v_dlt, int charge_number,double density_e){

    //double wp_const = 4 * k_pi * k_c*k_c * k_e * k_ke / (k_me*1e6);
//    double wp = sqrt(4 * k_pi * k_c*k_c * k_re * density_e);
    double wp = sqrt(k_wp*density_e);
    double rho_max = v_dlt / wp; //The shelding radius, rho_sh
    double rho_max_2 = pow(3 * charge_number / density_e, 1.0/3);
    if(rho_max<rho_max_2) rho_max = rho_max_2;
    double rho_max_3 = v_dlt * time_cooler;
    if(rho_max>rho_max_3) rho_max = rho_max_3;

    return rho_max;
}

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
    #pragma omp parallel for
    for(int i=0; i<ion_number; ++i){
        if(iszero(density_e.at(i))) {
            force_tr[i] = 0;
            force_long[i] = 0;
            continue;
        }
        double v2 = v_tr[i]*v_tr[i]+v_long[i]*v_long[i];
        if(v2>0){
            double dlt = this->dlt(tpr, v2, dlt2_eff_e, i);

            //Calculate rho_min
            double rho_min = rho_min_const/dlt;
            dlt = sqrt(dlt);

            double rho_max = max_impact_factor(dlt, charge_number, density_e[i]);
//            double wp = sqrt(wp_const*density_e[i]);
//
//            //Calculate rho_max
//            double rho_max = dlt/wp;
//            double rho_max_2 = pow(3*charge_number/density_e[i], 1.0/3);
//            if(rho_max<rho_max_2) rho_max = rho_max_2;
//            double rho_max_3 = dlt*time_cooler;
//            if(rho_max>rho_max_3) rho_max = rho_max_3;

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
    init(ebeam);
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
        double ve_tr = ve_rms_tr.at(0);
        double ve2 = ve_l*ve_l + ve_tr*ve_tr;

        #ifdef _OPENMP
            #pragma omp parallel for
        #endif // _OPENMP
        for(int i=0; i<ion_number; ++i) {
            if(iszero(density_e.at(i))) {
                force_tr[i] = 0;
                force_long[i] = 0;
                continue;
            }
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
        #ifdef _OPENMP
            #pragma omp parallel for
        #endif // _OPENMP
        for(int i=0; i<ion_number; ++i) {
            if(iszero(density_e.at(i))) {
                force_tr[i] = 0;
                force_long[i] = 0;
                continue;
            }
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

    #ifdef _OPENMP
    int i = omp_get_thread_num();
    P* p = &(this->p.at(i));
    gsl_integration_workspace* gw = this->gw.at(i);
    #else
    P* p = &(this->p);
    gsl_integration_workspace* gw = this->gw;
    #endif // _OPENMP

    p->ve_tr = ve_tr;
    p->ve_l = ve_l;
    p->v_tr = v_tr;
    p->v_l = v_l;
    p->flag = 0;

    ForceNonMagNumeric1D* ptr2 = this;
    auto ptr = [=](double x)->double{return ptr2->b(x,p);};
    gsl_function_pp<decltype(ptr)> fp(ptr);
    gsl_function *f = static_cast<gsl_function*>(&fp);

    double b_tr = 0;
    double error = 0;
    int status =  gsl_integration_qagiu(f, 0, espabs, esprel, limit, gw, &b_tr, &error);
    if(status == GSL_EDIVERGE){
        status = gsl_integration_qagiu(f, 0, 1e-10, 1e-10, 10*limit, gw, &b_tr, &error);
        if(status == GSL_EDIVERGE) std::cout<<"GSL integration qagui error for transverse friction force!"<<std::endl;
    }

    double b_l = 0;
    p->flag = 1;
    status = gsl_integration_qagiu(f, 0, espabs, esprel, limit, gw, &b_l, &error);
    if(status == GSL_EDIVERGE){
        status = gsl_integration_qagiu(f, 0, 1e-10, 1e-10, 10*limit, gw, &b_l, &error);
        if(status == GSL_EDIVERGE) std::cout<<"GSL integration qagui error for longitudinal friction force!"<<std::endl;
    }

    double ve3_tr = ve_tr*ve_tr*ve_tr;
    double ff = -f_const*ne*lc/ve3_tr;
    force_tr = ff*v_tr*b_tr;
    force_l = ff*v_l*b_l;
}

ForceNonMagNumeric1D::ForceNonMagNumeric1D(int n):limit(n){
    gsl_set_error_handler_off();
    #ifdef _OPENMP
    #pragma omp parallel
    {
    int n = omp_get_num_threads();
    }
    gw.resize(n);
    p.resize(n);
    for(int i=0; i<n; ++i) {
        gw.at(i) = gsl_integration_workspace_alloc(limit);
    }
    #else
    gw = gsl_integration_workspace_alloc(limit);
    #endif // _OPENMP
}
ForceNonMagNumeric1D::~ForceNonMagNumeric1D(){
    #ifdef _OPENMP
    #pragma omp parallel for
    for(int i=0; i<omp_get_num_threads(); ++i) {
        gsl_integration_workspace_free(gw.at(i));
    }
    #else
    gsl_integration_workspace_free(gw);
    #endif // _OPENMP
}

ForceNonMagNumeric3D::ForceNonMagNumeric3D(int n):limit(n){
    #ifdef _OPENMP
    #pragma omp parallel
    {
    int n = omp_get_num_threads();
    }
    giw.resize(n);
    gmw.resize(n);
    gow.resize(n);
    p.resize(n);

    for(int i=0; i<n; ++i) {
        giw.at(i) = gsl_integration_workspace_alloc(limit);
        gmw.at(i) = gsl_integration_workspace_alloc(limit);
        gow.at(i) = gsl_integration_workspace_alloc(limit);
    }

    #else
    giw = gsl_integration_workspace_alloc(limit);
    gmw = gsl_integration_workspace_alloc(limit);
    gow = gsl_integration_workspace_alloc(limit);
    #endif // _OPENMP
}
ForceNonMagNumeric3D::~ForceNonMagNumeric3D(){
    #ifdef _OPENMP
    if(use_gsl) {
        #pragma omp parallel for
        for(int i=0; i<omp_get_num_threads(); ++i) {
            gsl_integration_workspace_free(giw.at(i));
            gsl_integration_workspace_free(gmw.at(i));
            gsl_integration_workspace_free(gow.at(i));
        }
    }
    #else
    if(use_gsl) {
        gsl_integration_workspace_free(giw);
        gsl_integration_workspace_free(gmw);
        gsl_integration_workspace_free(gow);
    }
    #endif // _OPENMP
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
    #ifdef _OPENMP
    int i = omp_get_thread_num();
    gsl_integration_qagi(f, espabs, esprel, limit, giw.at(i), &result, &error);
    #else
    gsl_integration_qagi(f, espabs, esprel, limit, giw, &result, &error);
    #endif // _OPENMP

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
    double lc = 0;

    double k = 1.41421356237;  //sqrt(2)
    double sub_vl = v_l - k*ve_l*vl;
    double sub_vtr = v_tr - k*ve_tr*vtr*cos(phi);
    double f_bot = sub_vl*sub_vl + sub_vtr*sub_vtr + 2*ve_tr*ve_tr*vtr*vtr*sin(phi)*sin(phi);
    double f_inv_bot = 1/f_bot;
    double f = exp(-vtr*vtr-vl*vl)*vtr*f_inv_bot;
    if(use_mean_rho_min) lc = mean_lc;
    else {
        double rho_min = rho_min_const(charge_number)*f_inv_bot;
        lc = rho_max>rho_min?log(rho_max/rho_min):0;
    }

    f *= lc*sqrt(f_inv_bot);
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
    #ifdef _OPENMP
    int i = omp_get_thread_num();
    gsl_integration_qag(f, 0, k_pi, espabs, esprel, limit, key, giw.at(i), &result, &error);
    #else
    gsl_integration_qag(f, 0, k_pi, espabs, esprel, limit, key, giw, &result, &error);
    #endif // _OPENMP
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
    #ifdef _OPENMP
    int i = omp_get_thread_num();
    gsl_integration_qagi(f, espabs, esprel, limit, gmw.at(i), &result, &error);
    #else
    gsl_integration_qagi(f, espabs, esprel, limit, gmw, &result, &error);
    #endif // _OPENMP
    return result;
}

void ForceNonMagNumeric3D::init(EBeam& ebeam) {
    auto tpr = ebeam.temperature();
    if(tpr==Temperature::CONST) {
        const_tmpr = true;
    }
    else {
        const_tmpr = false;
    }
}

void ForceNonMagNumeric3D::pre_int(double sgm_vtr, double sgm_vl) {
    hlf_v2tr.resize(n_tr);
    hlf_v2l.resize(n_l);
    vtr_cos.resize(n_phi,vector<double>(n_tr));
    vl.resize(n_l);
    vtr.resize(n_tr);
    v2tr_sin2.resize(n_phi,vector<double>(n_tr));

    vector<double> phi(n_phi);

    double d_phi = k_pi/n_phi;
    phi.at(0) = d_phi/2;
    for(int i=1; i<n_phi; ++i) phi.at(i) = phi.at(i-1) + d_phi;

    double d_vtr = 3*sgm_vtr/n_tr;
    vtr.at(0) = d_vtr/2;
    for(int i=1; i<n_tr; ++i) vtr.at(i) = vtr.at(i-1) + d_vtr;

    double d_vl = 6*sgm_vl/n_l;
    vl.at(0) = -3*sgm_vl + d_vl/2;
    for(int i=1; i<n_l; ++i) vl.at(i) = vl.at(i-1) + d_vl;

    d = d_phi*d_vtr*d_vl;

    for(int i=0; i<n_tr; ++i) hlf_v2tr.at(i) = vtr.at(i)*vtr.at(i);
    for(int i=0; i<n_l; ++i) hlf_v2l.at(i) = -vl.at(i)*vl.at(i)/2;
    for(auto& e: vl) e*= -1;
    for(auto& e: phi) e = cos(e);   //cos(phi)
    for(int i=0; i<n_phi; ++i) {
        for(int j=0; j<n_tr; ++j) {
            vtr_cos.at(i).at(j) = -phi.at(i)*vtr.at(j);
        }
    }
    for(auto& e: phi) e = 1 - e*e;   //sin(phi)*sin(phi)
    for(int i=0; i<n_phi; ++i) {
        for(int j=0; j<n_tr; ++j) {
            v2tr_sin2.at(i).at(j) = hlf_v2tr.at(j)*phi.at(i);
        }
    }

    for(auto& e: hlf_v2tr) e /= -2;
}

void ForceNonMagNumeric3D::calc_exp_vtr(double sgm_vtr, double sgm_vl) {
    exp_vtr.resize(n_l, vector<double>(n_tr));
    double inv_ve2_tr = 1/(sgm_vtr*sgm_vtr);
    double inv_ve2_l = 1/(sgm_vl*sgm_vl);
    for(int i=0; i<n_l; ++i) {
        for(int j=0; j<n_tr; ++j) {
            exp_vtr.at(i).at(j) = exp(hlf_v2tr.at(j)*inv_ve2_tr + hlf_v2l.at(i)*inv_ve2_l)*vtr.at(j);
        }
    }
    f_inv_norm = 0;
    for(auto&v:exp_vtr)
        for(auto&e:v)
            f_inv_norm+=e;
    f_inv_norm = 1/(n_phi*f_inv_norm);
}

#ifdef _OPENMP
bool ForceNonMagNumeric3D::first_run = true;
vector<vector<double>> ForceNonMagNumeric3D::exp_vtr;
vector<double> ForceNonMagNumeric3D::hlf_v2tr;
vector<double> ForceNonMagNumeric3D::hlf_v2l;
vector<vector<double>> ForceNonMagNumeric3D::vtr_cos;
vector<double> ForceNonMagNumeric3D::vl;
vector<double> ForceNonMagNumeric3D::vtr;
vector<vector<double>> ForceNonMagNumeric3D::v2tr_sin2;
#endif // _OPENMP
void ForceNonMagNumeric3D::force_grid(double v, double v_tr, double v_l, double v2, double ve_tr, double ve_l, double ve2,
                               double f_const, double rho_min_const, int charge_number, double ne,
                               double& force_tr, double& force_l) {
    if(first_run) {
        pre_int(ve_tr, ve_l);
        calc_exp_vtr(ve_tr, ve_l);
        first_run = false;
    }
    else if(!const_tmpr) {
        calc_exp_vtr(ve_tr, ve_l);
    }
    double f_tr = 0;
    double f_l = 0;

    double rho_max = this->rho_max(charge_number, v2, ve2, ne);
    if(use_mean_rho_min) {
        mean_rho_min = this->rho_min_const(charge_number)/(v2+ve2);
        mean_lc = rho_max>mean_rho_min?log(rho_max/mean_rho_min):0;
    }
    for(int i=0; i<n_tr; ++i) {
        for(int j=0; j<n_l; ++j) {
            for(int k=0; k<n_phi; ++k) {
                double sub_vl = v_l + vl.at(j);
                double sub_vtr = v_tr + vtr_cos.at(k).at(i);
                double f_bot = sub_vl*sub_vl + sub_vtr*sub_vtr + v2tr_sin2.at(k).at(i);
                double f_inv_bot = 1/f_bot;
                double rho_min = this->rho_min_const(charge_number)*f_inv_bot;
                double lc = rho_max>rho_min?log(rho_max/rho_min):0;
                double f = exp_vtr.at(j).at(i)*lc*(f_inv_bot*sqrt(f_inv_bot));
                f_tr += sub_vtr*f;
                f_l += sub_vl*f;
            }
        }
    }

    double ff = f_const*ne*f_inv_norm;
    force_tr = ff*f_tr;
    force_l = ff*f_l;
}


#ifdef _OPENMP
double ForceNonMagNumeric3D::mean_rho_min = 0;
double ForceNonMagNumeric3D::mean_lc = 0;
#endif // _OPENMP

void ForceNonMagNumeric3D::force_gsl(double v, double v_tr, double v_l, double v2, double ve_tr, double ve_l, double ve2,
                               double f_const, double rho_min_const, int charge_number, double ne,
                               double& force_tr, double& force_l) {
    double rho_max = this->rho_max(charge_number, v2, ve2, ne);
    if(use_mean_rho_min) {
        mean_rho_min = this->rho_min_const(charge_number)/(v2+ve2);
        mean_lc = rho_max>mean_rho_min?log(rho_max/mean_rho_min):0;
    }

    #ifdef _OPENMP
    int i = omp_get_thread_num();
    P* p = &(this->p.at(i));
    gsl_integration_workspace* gow = this->gow.at(i);
    #else
    P* p = &(this->p);
    #endif // _OPENMP

    p->ve_tr = ve_tr;
    p->ve_l = ve_l;

    ForceNonMagNumeric3D* class_ptr = this;

    double error;
    double inv_norm = 0.35917424443382906;  //1/norm with norm = 0.886226925*k_pi

    p->v_tr = v_tr;
    p->v_l = v_l;
    p->charge_number = charge_number;
    p->rho_max = rho_max;

    auto ptr = [=](double x)->double{return class_ptr->outter_integrand(x,const_cast<ForceNonMagNumeric3D::P*>(p));};
    gsl_function_pp<decltype(ptr)> fp(ptr);
    gsl_function *f = static_cast<gsl_function*>(&fp);

    double f_tr = 0;
    p->flag = 0;
    gsl_integration_qagiu(f, 0, espabs, esprel, limit, gow, &f_tr, &error);
    double f_l = 0;
    p->flag = 1;
    gsl_integration_qagiu(f, 0, espabs, esprel, limit, gow, &f_l, &error);

    double ff = f_const*ne*inv_norm;
    force_tr = ff*f_tr;
    force_l = ff*f_l;
}

void ForceNonMagNumeric3D::force(double v, double v_tr, double v_l, double v2, double ve_tr, double ve_l, double ve2,
                               double f_const, double rho_min_const, int charge_number, double ne,
                               double& force_tr, double& force_l) {
    if(use_gsl) {
        force_gsl(v, v_tr, v_l, v2, ve_tr, ve_l, ve2, f_const, rho_min_const, charge_number, ne, force_tr, force_l);
    }
    else {
        force_grid(v, v_tr, v_l, v2, ve_tr, ve_l, ve2, f_const, rho_min_const, charge_number, ne, force_tr, force_l);
    }

}


void ForceMeshkov::force(double ve_tr, double ve_l, double ve2_tr, double ve2_l, double v_tr, double v_l,
                         double v2, double rho_min_const, int charge_number,  double density, double f_const,
                        double& force_tr,double& force_l) {
        double rho_L    = k_me_kg * ve_tr / ( mag_field * k_e ); //SI units, in m
        double wp = sqrt(k_wp*density);

        //The Number of multiple adiabatic collisions with a single electron
        double v_sh = sqrt(v2 + ve2_l);
        double N_col = 1 + ve_tr/(k_pi*v_sh);
//        //dynamic shielding radius
//        double rho_sh = sqrt(v2 + ve2_l) / wp;
        //minimum impact parameter
        double rho_min = rho_min_const / (v2 + ve2_l);  //in units of m
        //intermediate impact parameter
        double rho_F = rho_L * v_sh / ve_tr;
        double rho_max = max_impact_factor(v_sh,charge_number,density);
        double rho_min_mag = k*rho_L;
        if(rho_min_mag<rho_min) rho_min_mag = rho_min;
        //Coulomb Logarithms
        double L_M = rho_max>rho_min_mag? log(rho_max/rho_min_mag): 0;
        double L_A = k*rho_L>rho_F? log(k*rho_L/rho_F): 0;
        double L_F = rho_F>rho_min? log(rho_F/rho_min): 0;

        double ellipse = (v_tr*v_tr)/ve2_tr + (v_l*v_l)/ve2_l;
        double v3 = v2*sqrt(v2);

        double f_tr = 0;
        double f_l = 0;

        if(v2 > ve2_tr) { //Region I
           double k1 = 1-3*v_l*v_l/v2;  //(v2_tr-2*v2_l)/v2
           double k2 = 2 + k1;         //3*v2_tr/v2

           f_tr = 2*L_F + L_M*k1;
           f_tr /= v3;

           f_l = 2 + 2*L_F + L_M*k2;
           f_l /= v3;
        }
        else if(v2 < ve2_l) { //Region III
            double ve3_tr = ve2_tr*ve_tr;
            double ve3_l = ve2_l*ve_l;
            f_tr = 2*(L_F + N_col*L_A)/ve3_tr + L_M/ve3_l;
            f_l = 2*(L_F + N_col*L_A)/(ve2_tr*ve_l) + L_M/ve3_l;
          }

        //Region II (a or b)
        else{ //This constrains to a donut shape (d_paral_e < sqrt(v2) < d_perp_e).
            f_tr = (1-3*v_l*v_l/v2) * (L_M/v3);
            f_tr += (2/(ve2_tr*ve_tr)) * (L_F + N_col*L_A);

            if( ellipse <= 1.0 ){ // Region IIb
                  //This is the same as result_long in Region 3
                  f_l = 2*(L_F + N_col*L_A)/(ve2_tr*ve_l) + L_M/(ve2_l*ve_l);
            }
            else{ //It must be Region IIa
                  f_l = (((3*v_tr*v_tr*L_M)/v2) + 2) / v3;
                  f_l += 2 * (L_F + N_col*L_A)/(ve2_tr * v_l);
            }
        }

        force_tr = f_const * density * v_tr * f_tr;
        force_l = f_const * density * v_l * f_l;
}

void ForceMeshkov::friction_force(int charge_number, int ion_number,
            vector<double>& v_tr, vector<double>& v_l, vector<double>& density,
            EBeam& ebeam, vector<double>& force_tr, vector<double>& force_long) {

    double rho_min_const = charge_number * k_rho_min;
    double f_const = k_f*charge_number*charge_number/2;
    force_tr.resize(ion_number);
    force_long.resize(ion_number);
    auto tpr = ebeam.temperature();
    vector<double>& ve_rms_l = ebeam.get_v(EBeamV::V_RMS_L);
    vector<double>& ve_rms_tr =  ebeam.get_v(EBeamV::V_RMS_TR);

    switch(tpr) {
    case Temperature::CONST: {
        double ve_l = ve_rms_l.at(0);
        double ve_tr = ve_rms_tr.at(0);
        double ve2_tr = ve_tr*ve_tr;
        double ve2_l = ve_l*ve_l;
        double ve2 = ve2_tr + ve2_l;
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif // _OPENMP
        for(int i=0; i<ion_number; ++i) {
            if(iszero(density.at(i))) {
                force_tr[i] = 0;
                force_long[i] = 0;
                continue;
            }
            double v2 = v_tr[i]*v_tr[i]+v_l[i]*v_l[i];
            double v = sqrt(v2);
            if(v2>0) {
                force(ve_tr, ve_l, ve2_tr, ve2_l, v_tr[i], v_l[i],v2, rho_min_const, charge_number, density[i], f_const, force_tr[i], force_long[i]);
            }
        }
        break;
    }
    case Temperature::VARY: {
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif // _OPENMP
        for(int i=0; i<ion_number; ++i) {
            if(iszero(density.at(i))) {
                force_tr[i] = 0;
                force_long[i] = 0;
                continue;
            }
            double v2 = v_tr[i]*v_tr[i]+v_l[i]*v_l[i];
            double v = sqrt(v2);
            if(v2>0) {
                double ve_l = ve_rms_l.at(i);
                double ve_tr = ve_rms_tr.at(i);
                double ve2_tr = ve_tr*ve_tr;
                double ve2_l = ve_l*ve_l;
                double ve2 = ve2_tr + ve2_l;
                force(ve_tr, ve_l, ve2_tr, ve2_l, v_tr[i], v_l[i],v2, density[i], rho_min_const, charge_number, f_const, force_tr[i], force_long[i]);
            }
        }
        break;
    }
    }
}

#ifdef _OPENMP
bool ForceDSM::first_run = true;
vector<double> ForceDSM::a;
vector<double> ForceDSM::cos_a;
vector<double> ForceDSM::tan_a;
vector<double> ForceDSM::t2;
vector<double> ForceDSM::ve;
vector<double> ForceDSM::exp_ve2;
bool ForceDSM::first_run_fa = true;
double ForceDSM::f_inv_norm;
vector<vector<double>> ForceDSM::exp_vtr;
vector<double> ForceDSM::hlf_v2tr;
vector<double> ForceDSM::hlf_v2l;
vector<vector<double>> ForceDSM::vtr_cos;
vector<double> ForceDSM::vl;
vector<double> ForceDSM::vtr;
vector<vector<double>> ForceDSM::v2tr_sin2;
#endif // _OPENMP

void ForceDSM::calc_alpha() {
//    sin_a.resize(n_a);
    cos_a.resize(n_a);
    tan_a.resize(n_a);
    double da = k_pi/n_a;
    double a = (-k_pi+da)/2;
    for(int i=0; i<n_a; ++i) {
//        sin_a.at(i) = sin(a);
        cos_a.at(i) = cos(a);
        tan_a.at(i) = tan(a);
        a += da;
    }
}

void ForceDSM::calc_ve() {
    t2.resize(n_ve);
    ve.resize(n_ve);
    double dt = 1.0/n_ve;
    double t = dt/2;
    for(int i=0; i<n_ve; ++i) {
        ve.at(i) = (1-t)/t;
        t2.at(i) = t*t;
        t += dt;
    }
}

void ForceDSM::calc_exp_ve2(double ve2_l, vector<double>& t2) {
    exp_ve2.resize(n_ve);
    double coef = -1/(2*ve2_l);
    for(int i=0; i<n_ve; ++i) {
        exp_ve2.at(i) = exp(coef*ve.at(i)*ve.at(i))/t2.at(i);
    }
}

void ForceDSM::init(EBeam& ebeam) {
    auto tpr = ebeam.temperature();
    if(tpr==Temperature::CONST) {
        const_tpr = true;
    }
    else {
        const_tpr = false;
    }
}


void ForceDSM::pre_int(double sgm_vtr, double sgm_vl) {
    hlf_v2tr.resize(n_tr);
    hlf_v2l.resize(n_l);
    vtr_cos.resize(n_phi,vector<double>(n_tr));
    vl.resize(n_l);
    vtr.resize(n_tr);
    v2tr_sin2.resize(n_phi,vector<double>(n_tr));

    vector<double> phi(n_phi);

    double d_phi = k_pi/n_phi;
    phi.at(0) = d_phi/2;
    for(int i=1; i<n_phi; ++i) phi.at(i) = phi.at(i-1) + d_phi;

    double d_vtr = 3*sgm_vtr/n_tr;
    vtr.at(0) = d_vtr/2;
    for(int i=1; i<n_tr; ++i) vtr.at(i) = vtr.at(i-1) + d_vtr;

    double d_vl = 6*sgm_vl/n_l;
    vl.at(0) = -3*sgm_vl + d_vl/2;
    for(int i=1; i<n_l; ++i) vl.at(i) = vl.at(i-1) + d_vl;

//    d = d_phi*d_vtr*d_vl;

    for(int i=0; i<n_tr; ++i) hlf_v2tr.at(i) = vtr.at(i)*vtr.at(i);
    for(int i=0; i<n_l; ++i) hlf_v2l.at(i) = -vl.at(i)*vl.at(i)/2;
    for(auto& e: vl) e*= -1;
    for(auto& e: phi) e = cos(e);   //cos(phi)
    for(int i=0; i<n_phi; ++i) {
        for(int j=0; j<n_tr; ++j) {
            vtr_cos.at(i).at(j) = -phi.at(i)*vtr.at(j);
        }
    }
    for(auto& e: phi) e = 1 - e*e;   //sin(phi)*sin(phi)
    for(int i=0; i<n_phi; ++i) {
        for(int j=0; j<n_tr; ++j) {
            v2tr_sin2.at(i).at(j) = hlf_v2tr.at(j)*phi.at(i);
        }
    }

    for(auto& e: hlf_v2tr) e /= -2;
}

void ForceDSM::calc_exp_vtr(double sgm_vtr, double sgm_vl) {
    exp_vtr.resize(n_l, vector<double>(n_tr));
    double inv_ve2_tr = 1/(sgm_vtr*sgm_vtr);
    double inv_ve2_l = 1/(sgm_vl*sgm_vl);
    for(int i=0; i<n_l; ++i) {
        for(int j=0; j<n_tr; ++j) {
            exp_vtr.at(i).at(j) = exp(hlf_v2tr.at(j)*inv_ve2_tr + hlf_v2l.at(i)*inv_ve2_l)*vtr.at(j);
        }
    }
    f_inv_norm = 0;
    for(auto&v:exp_vtr)
        for(auto&e:v)
            f_inv_norm+=e;
    f_inv_norm = 1/(n_phi*f_inv_norm);
}



void ForceDSM::force(double ve_tr, double ve_l, double ve2_tr, double ve2_l, double v_tr, double v_l,
                         double v2, double rho_min_const, int charge_number,  double density, double f_const,
                        double& force_tr,double& force_l) {
        double rho_L    = k_me_kg * ve_tr / ( mag_field * k_e ); //SI units, in m
        double wp = sqrt(k_wp*density);

        //The Number of multiple adiabatic collisions with a single electron
        double v_sh = sqrt(v2 + ve2_l);
        double N_col = 1 + ve_tr/(k_pi*v_sh);
//        //dynamic shielding radius
//        double rho_sh = sqrt(v2 + ve2_l) / wp;
        //minimum impact parameter
        double rho_min = rho_min_const / (v2 + ve2_l);  //in units of m
        //intermediate impact parameter
        double rho_F = rho_L * v_sh / ve_tr;
        double rho_max = max_impact_factor(v_sh,charge_number,density);
        double rho_min_mag = k*rho_L;
        if(rho_min_mag<rho_min) rho_min_mag = rho_min;
        //Coulomb Logarithms
        double L_M = rho_max>rho_min_mag? log(rho_max/rho_min_mag): 0;
        double L_A = k*rho_L>rho_F? log(k*rho_L/rho_F): 0;
        double L_F = rho_F>rho_min? log(rho_F/rho_min): 0;
        if(rho_F>rho_max) L_F = rho_max>rho_min? log(rho_max/rho_min): 0;
        double f_tr = 0;
        double f_l = 0;
        double f_mag_tr = 0;
        double f_mag_l = 0;
        double f_fa_tr = 0;
        double f_fa_l = 0;

        force_l = 0;
        force_tr = 0;

        if(L_M>0) {
            if(first_run) {
                calc_alpha();
                calc_ve();
                calc_exp_ve2(ve2_l, t2);
                first_run = false;
            }

            double y = v_l/ve_l;
            double z = v_tr/ve_l;
            int sgn = 1;
            if(z<0) {
                sgn = -1;
                z *= -1;
            }

            for(int i=0; i<n_a; ++i) {
                double yztan = y + z*tan_a.at(i);
                double exp_yztan2 = exp(-yztan*yztan/2);
                double ycos_zsin = yztan*cos_a.at(i);
                double fi_l = ycos_zsin*exp_yztan2;
                f_tr += tan_a.at(i)*fi_l;
                f_l += fi_l;
            }
            double da = k_pi/n_a;
            f_tr *= da;
            f_l *= da;
            if(iszero(z,1e-8)) {
                f_l = 2*v_l/ve_l*exp(-v_l*v_l/(2*ve2_l));
            }
            double coef = f_const*density*L_M/(sqrt(2*k_pi)*ve2_l);
            f_mag_tr += sgn*coef*f_tr;
            f_mag_l += coef*f_l;

            f_l = 0;
            if(v2>ve2_l) {  //Non-logarithm contribution
                if(!const_tpr) calc_exp_ve2(ve2_l,t2);
                for(int i=0; i<n_ve; ++i) {
                    double vm = v_l - ve.at(i);
                    double vp = v_l + ve.at(i);
                    double bm = v_tr*v_tr + vm*vm;
                    bm = sqrt(bm*bm*bm);
                    double bp = v_tr*v_tr + vp*vp;
                    bp = sqrt(bp*bp*bp);
                    f_l += (vm/bm+vp/bp)*exp_ve2.at(i);
                }
                double dt = 1.0/n_ve;
                f_l *= dt;
                coef = 2*coef/L_M;
                f_mag_l += coef*f_l;
            }

            force_tr += f_mag_tr;
            force_l += f_mag_l;
        }

        if(mag_only) return;
        if(L_A>0 || L_F>0) {    //Fast collision & adiabatic collision
            if(first_run_fa) {
                pre_int(ve_tr, ve_l);
                calc_exp_vtr(ve_tr, ve_l);
                first_run_fa = false;
            }
            else if(!const_tpr) {
                calc_exp_vtr(ve_tr, ve_l);
            }
            f_tr = 0;
            f_l = 0;

            for(int i=0; i<n_tr; ++i) {
                for(int j=0; j<n_l; ++j) {
                    for(int k=0; k<n_phi; ++k) {
                        double sub_vl = v_l + vl.at(j);
                        double sub_vtr = v_tr + vtr_cos.at(k).at(i);
                        double f_bot = sub_vl*sub_vl + sub_vtr*sub_vtr + v2tr_sin2.at(k).at(i);
                        double f_inv_bot = 1/f_bot;
                        double rho_min = rho_min_const*f_inv_bot;
                        L_F = rho_F>rho_min? log(rho_F/rho_min): 0;
                        if(rho_F>rho_max) L_F = rho_max>rho_min? log(rho_max/rho_min): 0;
                        double f = exp_vtr.at(j).at(i)*(L_F+N_col*L_A)*(f_inv_bot*sqrt(f_inv_bot));
                        f_tr += sub_vtr*f;
                        f_l += sub_vl*f;
                    }
                }
            }

            double ff = 2*f_const*density*f_inv_norm;
            f_fa_tr = ff*f_tr;
            f_fa_l = ff*f_l;
        }
        force_tr += f_fa_tr;
        force_l += f_fa_l;

}

void ForceDSM::friction_force(int charge_number, int ion_number,
            vector<double>& v_tr, vector<double>& v_l, vector<double>& density,
            EBeam& ebeam, vector<double>& force_tr, vector<double>& force_long) {

    init(ebeam);

    double rho_min_const = charge_number * k_rho_min;
    double f_const_non_mag = k_f*charge_number*charge_number;
    double f_const = f_const_non_mag/2;
    force_tr.resize(ion_number);
    force_long.resize(ion_number);
    auto tpr = ebeam.temperature();
    vector<double>& ve_rms_l = ebeam.get_v(EBeamV::V_RMS_L);
    vector<double>& ve_rms_tr =  ebeam.get_v(EBeamV::V_RMS_TR);

    switch(tpr) {
    case Temperature::CONST: {
        double ve_l = ve_rms_l.at(0);
        double ve_tr = ve_rms_tr.at(0);
        double ve2_tr = ve_tr*ve_tr;
        double ve2_l = ve_l*ve_l;
        double ve2 = ve2_tr + ve2_l;
        auto start_time = std::chrono::high_resolution_clock::now();
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif // _OPENMP
        for(int i=0; i<ion_number; ++i) {
            if(iszero(density.at(i))) {
                force_tr[i] = 0;
                force_long[i] = 0;
                continue;
            }
            double v2 = v_tr[i]*v_tr[i]+v_l[i]*v_l[i];
            double v = sqrt(v2);
            if(v2>0) {
                force(ve_tr, ve_l, ve2_tr, ve2_l, v_tr[i], v_l[i],v2, rho_min_const, charge_number, density[i], f_const, force_tr[i], force_long[i]);
            }
        }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto time = end_time - start_time;
    std::cout << "Force mag DSM took " << time/std::chrono::microseconds(1) << " us to run.\n";
        break;
    }
    case Temperature::VARY: {

        for(int i=0; i<ion_number; ++i) {
            if(iszero(density.at(i))) {
                force_tr[i] = 0;
                force_long[i] = 0;
                continue;
            }
            double v2 = v_tr[i]*v_tr[i]+v_l[i]*v_l[i];
            double v = sqrt(v2);
            if(v2>0) {
                double ve_l = ve_rms_l.at(i);
                double ve_tr = ve_rms_tr.at(i);
                double ve2_tr = ve_tr*ve_tr;
                double ve2_l = ve_l*ve_l;
                double ve2 = ve2_tr + ve2_l;
                force(ve_tr, ve_l, ve2_tr, ve2_l, v_tr[i], v_l[i],v2, density[i], rho_min_const, charge_number, f_const, force_tr[i], force_long[i]);
            }
        }
        break;
    }
    }

}
