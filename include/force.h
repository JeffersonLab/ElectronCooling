#ifndef FORCE_H
#define FORCE_H

#include <functional>
#include <gsl/gsl_integration.h>
#include <vector>
#include "beam.h"
#include "constants.h"

using std::vector;
enum class ForceFormula {PARKHOMCHUK, NONMAG_DERBENEV, NONMAG_MESHKOV, NONMAG_NUM1D, NONMAG_NUM3D, MESHKOV, DSM};

class FrictionForceSolver{
protected:
    double time_cooler;
    double mag_field = 0;
    const double k_f = -4*k_pi*k_c*k_c*k_ke*k_ke*k_e*k_e*k_e/(k_me*1e6);
    const double k_wp = 4*k_pi*k_c*k_c*k_e*k_ke/(k_me*1e6);
    const double k_rho_min = k_e*k_ke*k_c*k_c/(k_me*1e6);
    virtual void init(EBeam& ebeam){};
    virtual void fin(){};
    double max_impact_factor(double v_dlt, int charge_number,double density_e);
public:
    void set_time_cooler(double t){time_cooler = t;}
    void set_mag_field(double x){mag_field = x;}
    double t_cooler(){return time_cooler;}
    virtual void friction_force(int charge_number, int ion_number,
            vector<double>& v_tr, vector<double>& v_l, vector<double>& density,
            EBeam& ebeam, vector<double>& force_tr, vector<double>& force_long) = 0;
};

class ForcePark: public FrictionForceSolver {
private:
    const double k_f = -4*k_c*k_c*k_ke*k_ke*k_e*k_e*k_e/(k_me*1e6);
    double t_eff = 0; //Effective temperature.
    double v_eff = 0; //Effective velocity.
    void rho_lamor_dlt2_eff_e(double v2_eff_e, double mag_field, vector<double>& v_rms_l, vector<double>& v_rms_t, Temperature tpr,
                          int ion_number, vector<double>& dlt2_eff_e, vector<double>& rho_lamor);
    double dlt(Temperature tpr,  double v2, vector<double>& dlt2_eff_e, int i);
    double lc(Temperature tpr, double rho_max, double rho_min, vector<double>& rho_lamor, int i);
public:
    void set_t_eff(double x){t_eff = x; v_eff = sqrt(t_eff*k_c*k_c/(k_me*1e6));}
    void set_v_eff(double v){v_eff = v; t_eff = v_eff*v_eff*k_me*1e6/(k_c*k_c);}
    virtual void friction_force(int charge_number, int ion_number,
            vector<double>& v_tr, vector<double>& v_l, vector<double>& density,
            EBeam& ebeam, vector<double>& force_tr, vector<double>& force_long);
};

class ForceNonMag: public FrictionForceSolver {
protected:
    bool smooth_rho_max = false;
    virtual double f_const(int charge_number){return charge_number*charge_number*k_f;}
    double rho_min_const(int charge_number) {return charge_number*k_rho_min;}
    double rho_max_1(int charge_number, double density_e){return pow(3*charge_number/density_e, 1.0/3);}
    double rho_max_2(double v){return v*time_cooler;};
    double rho_max(int charge_number, double v2, double ve2, double ne);
    virtual void force(double v, double v_tr, double v_l, double v2, double ve_tr, double ve_l, double ve2,
                       double f_const,double rho_min_const, int charge_number, double ne,
                       double& force_tr, double& force_l) = 0;

public:
    void set_smooth_rho_max(bool b){smooth_rho_max = b;}
    virtual void friction_force(int charge_number, int ion_number,
            vector<double>& v_tr, vector<double>& v_l, vector<double>& density,
            EBeam& ebeam, vector<double>& force_tr, vector<double>& force_long);
};

class ForceNonMagDerbenev: public ForceNonMag {
private:
    double rho_max_ve_tr(double ve_tr, double ne);
    void force(double v, double v_tr, double v_l, double v2, double ve_tr, double ve_l, double ve2,
                               double f_const, double rho_min_const, int charge_number, double ne,
                               double& force_tr, double& force_l);
public:

};

class ForceNonMagMeshkov: public ForceNonMag {
private:
    void force(double v, double v_tr, double v_l, double v2, double ve_tr, double ve_l, double ve2,
                               double f_const, double rho_min_const, int charge_number, double ne,
                               double& force_tr, double& force_l);
public:

};

class ForceNonMagNumeric1D: public ForceNonMag {
private:
    double f_const(int charge_number){return 0.28209479177387814*charge_number*charge_number*k_f;} //Coef - 1/sqrt(4*pi)
    const double k_f = 2*sqrt(2*k_pi)*k_pi*k_c*k_c*k_ke*k_ke*k_e*k_e*k_e/(k_me*1e6);
    gsl_integration_workspace *gw = nullptr;
    size_t limit = 100;
    double espabs = 1e-6;
    double esprel = 1e-6;
    struct P{
        double v_tr;
        double v_l;
        double ve_tr;
        double ve_l;
        int flag;   //0: calculate B_tr; else: calculate B_l;
    }p;
    double b(double q, void* params);
    void force(double v, double v_tr, double v_l, double v2, double ve_tr, double ve_l, double ve2,
                               double f_const, double rho_min_const, int charge_number, double ne,
                               double& force_tr, double& force_l);
public:
    void set_espabs(double x){espabs = x;}
    void set_esprel(double x){esprel = x;}
    ForceNonMagNumeric1D(int n=100):limit(n){gw = gsl_integration_workspace_alloc(limit);}
    ~ForceNonMagNumeric1D(){gsl_integration_workspace_free(gw);}
};

//
class ForceNonMagNumeric3D: public ForceNonMag {
private:
    gsl_integration_workspace *giw;
    gsl_integration_workspace *gmw;
    gsl_integration_workspace *gow;

    size_t limit = 100;
    double espabs = 1e-5;
    double esprel = 1e-3;
    struct P{
        double v_tr;
        double v_l;
        double ve_tr;
        double ve_l;
        double vtr;
        double vl;
        double rho_max;
        int charge_number;
        int flag;   //0: calculate B_tr; else: calculate B_l;
    }p;

    bool use_gsl = false;
    bool use_mean_rho_min = false;
    double mean_rho_min = 0;
    double mean_lc = 0;

    bool first_run = true;
    bool const_tmpr = true;
    int n_tr = 20;
    int n_l = 10;
    int n_phi = 10;
    double d;
    double f_inv_norm;
    vector<vector<double>> exp_vtr;
    vector<double> hlf_v2tr;
    vector<double> hlf_v2l;
    vector<vector<double>> vtr_cos;
    vector<double> vl;
    vector<double> vtr;
    vector<vector<double>> v2tr_sin2;

    void pre_int(double sgm_vtr, double sgm_vl);
    void calc_exp_vtr(double sgm_vtr, double sgm_vl);



    void init(EBeam& ebeam);
    double inner_integrand(double phi, void* params);
    double middle_integrand(double vl, void* params);
    double outter_integrand(double vtr, void* params);
    double inner_norm_integrand(double vl, void* params);
    double outter_norm_integrand(double vtr, void* params);
    void force_grid(double v, double v_tr, double v_l, double v2, double ve_tr, double ve_l, double ve2,
                               double f_const, double rho_min_const, int charge_number, double ne,
                               double& force_tr, double& force_l);
    void force_gsl(double v, double v_tr, double v_l, double v2, double ve_tr, double ve_l, double ve2,
                               double f_const, double rho_min_const, int charge_number, double ne,
                               double& force_tr, double& force_l);
    void force(double v, double v_tr, double v_l, double v2, double ve_tr, double ve_l, double ve2,
                               double f_const, double rho_min_const, int charge_number, double ne,
                               double& force_tr, double& force_l);
public:
    void set_espabs(double x){espabs = x;}
    void set_esprel(double x){esprel = x;}
    void set_gsl(bool b) {use_gsl = b;}
    void set_mean_rho_min(bool b) {use_mean_rho_min = b;}
    void set_grid(int ntr, int nl, int nphi){n_tr = ntr; n_l = nl; n_phi = nphi; first_run = true;}
    ForceNonMagNumeric3D(int n=100):limit(n){giw = gsl_integration_workspace_alloc(limit); gmw = gsl_integration_workspace_alloc(limit);
        gow = gsl_integration_workspace_alloc(limit);}
    ~ForceNonMagNumeric3D(){
        if(use_gsl) {
            gsl_integration_workspace_free(giw);
            gsl_integration_workspace_free(gmw);
            gsl_integration_workspace_free(gow);
        }
    }
};

//gsl function wrapper for member functions in class.
template< typename F >
    class gsl_function_pp : public gsl_function {
    public:
    gsl_function_pp(const F& func) : _func(func) {
    function = &gsl_function_pp::invoke;
    params=this;
    }
    private:
    const F& _func;
    static double invoke(double x, void *params) {
    return static_cast<gsl_function_pp*>(params)->_func(x);
    }
};


class ForceMeshkov: public FrictionForceSolver {
protected:
    //A fudge factor to smooth the friction force shape. Using the
        // "Classical" default definition here from the BETACOOL documentation
    double k = 2;
    void force(double ve_tr, double ve_l, double ve2_tr, double ve2_l, double v_tr, double v_l, double v2,
               double rho_min_const, int charge_number,  double density, double f_const, double& force_tr,double& force_l);
public:
    void set_smooth_factor(double x){k = x;}
    virtual void friction_force(int charge_number, int ion_number,
            vector<double>& v_tr, vector<double>& v_l, vector<double>& density,
            EBeam& ebeam, vector<double>& force_tr, vector<double>& force_long);
};

class ForceDSM: public FrictionForceSolver {    //Derbenev-Skrinsky-Meshkov formula for magnetized cooling.
    bool mag_only = false;

     //A fudge factor to smooth the friction force shape. Using the
    // "Classical" default definition here from the BETACOOL documentation
    double k = 2;
    int n_a = 100;
    int n_ve = 100;
    bool first_run = true;
    bool const_tpr = true;
    vector<double> a;
//    vector<double> sin_a;
    vector<double> cos_a;
    vector<double> tan_a;
    vector<double> t2;
    vector<double> ve;
    vector<double> exp_ve2;

    bool first_run_fa = true;
    int n_tr = 20;
    int n_l = 10;
    int n_phi = 10;
    double d;
    double f_inv_norm;
    vector<vector<double>> exp_vtr;
    vector<double> hlf_v2tr;
    vector<double> hlf_v2l;
    vector<vector<double>> vtr_cos;
    vector<double> vl;
    vector<double> vtr;
    vector<vector<double>> v2tr_sin2;

protected:
    void init(EBeam& ebeam);
    void pre_int(double sgm_vtr, double sgm_vl);
    void calc_exp_vtr(double sgm_vtr, double sgm_vl);
    void calc_alpha();
    void calc_ve();
    void calc_exp_ve2(double ve2_l, vector<double>& t2);
    void force(double ve_tr, double ve_l, double ve2_tr, double ve2_l, double v_tr, double v_l, double v2,
               double rho_min_const, int charge_number,  double density, double f_const, double& force_tr,double& force_l);
public:
    void set_smooth_factor(double x){k = x;}
    void set_steps(int i) {n_a = i; n_ve = i; first_run = true;}
    void set_grid(int ntr, int nl, int nphi){n_tr = ntr; n_l = nl; n_phi = nphi; first_run_fa = true;}
    void set_mag_only(bool b){mag_only = b;}
    virtual void friction_force(int charge_number, int ion_number,
            vector<double>& v_tr, vector<double>& v_l, vector<double>& density,
            EBeam& ebeam, vector<double>& force_tr, vector<double>& force_long);
};
#endif // FORCE_H
