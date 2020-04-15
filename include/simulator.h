#ifndef SIMULATOR_H_INCLUDED
#define SIMULATOR_H_INCLUDED

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "cooler.h"
#include "beam.h"
#include "ecooling.h"
#include "ions.h"
#include "ibs.h"
#include "luminosity.h"
#include "ring.h"

using std::string;
using std::ofstream;

enum class DynamicModel {RMS, PARTICLE, MODEL_BEAM = PARTICLE, TURN_BY_TURN};

//extern double vl_emit_nx, vl_emit_ny, vl_dp_p, vl_sigma_s, vl_rx_ibs, vl_ry_ibs, vl_rs_ibs,
//    vl_rx_ecool, vl_ry_ecool, vl_rs_ecool, vl_rx_total, vl_ry_total, vl_rs_total, vl_t;

class Simulator{
 protected:
    double t = 0;
    int n_step;
    double dt;
    double t0 = 0;
    bool ibs = true;
    bool ecool = true;
    bool fixed_bunch_length = false;
    bool reset_time = true;
    bool overwrite = false;
    bool calc_luminosity = false;
    int output_itvl = 1;
    int ion_save_itvl = -1;
    ofstream outfile;
    string outfilename = "output_dynamic.txt";
    vector<double> r_ibs = {0,0,0};
    vector<double> r_ecool = {0,0,0};
    vector<double> r = {0,0,0};
    vector<double> emit = {0,0,0,0};
    void output_sddshead();
    void output_to_file();
    void output(bool bunched=true, double v_rf=0, double lum=0);
    virtual void update_ibeam(Beam& ion, Ions& ion_sample, Ring& ring, EBeam& ebeam, Cooler& cooler, ECoolRate* ecool_solver)=0;
    virtual void adjust_rf_voltage(Ring& ring) = 0;
    virtual void save_ions(int i, Ions& ion_sample) = 0;
    virtual void precondition(Ions& ion_sample){};

 public:
    Simulator(double time, int n):t(time),n_step(n){dt = time/n;}
    void set_ibs(bool b){ibs=b;}
    void set_ecool(bool b){ecool=b;}
    void set_ion_save(int x){ion_save_itvl = x; }
    void set_output_file(string filename){outfilename = filename; }
    void set_output_intvl(int x){output_itvl = x; }
    void set_fixed_bunch_length(bool b){fixed_bunch_length = b; }
    void set_ini_time(double t){t0 = t;}
    void set_reset_time(bool b){reset_time = b;}
    void set_overwrite(bool b) {overwrite = b; }
    void set_calc_lum(bool b) {calc_luminosity = b; }

    virtual void run(Beam& ion, Ions& ion_sample, Cooler& cooler, EBeam& ebeam,
                     Ring& ring, IBSSolver* ibs_solver, ECoolRate* ecool_solver,
                     FrictionForceSolver* force_solver, LuminositySolver* lum_solver);
};

#endif // SIMULATOR_H_INCLUDED
