#ifndef ECOOLING_H
#define ECOOLING_H

#include <initializer_list>
#include <vector>

#include "beam.h"
#include "cooler.h"
#include "force.h"
#include "ions.h"
#include "ring.h"

using std::vector;
using std::initializer_list;

enum class IonSample {SINGLE_PARTICLE, MONTE_CARLO};
enum class ECoolRateScratch {XP_BET, YP_BET, XP, YP, DP_P, V_TR, V_LONG, FORCE_X, FORCE_Y, FORCE_Z};

class ECoolRate{
    bool p_shift_ = false;             //Position shift. false: ion center and e- center overlap, true: there's a shift between the beam
    bool v_shift_ = false;             //Vecocity shift.
    double bunch_separate_ = 0;
    double t_cooler_ = 0;
    int n_long_sample_ = 50;
    int scratch_size = 0;
    vector<double> ne;
    vector<double> xp_bet, yp_bet, xp, yp, dp_p, x, y, x_bet, y_bet;
    vector<double> v_tr, v_long;
    vector<double> force_x, force_y, force_z;
    void electron_density(Ions& ion_sample, EBeam &ebeam);
    void init_scratch(int n_sample);
    void space_to_dynamic(int n_sample, Beam &ion, Ions &ion_sample);
    void beam_frame(int n_sample, double gamma_e);
    void force(int n_sample, Beam &ion, EBeam &ebeam, Cooler &cooler, FrictionForceSolver &force_solver);
    void restore_velocity(int n_sample, EBeam &ebeam);
    void bunched_to_coasting(Beam &ion, Ions& ion_sample, EBeam &ebeam, Cooler &cooler, FrictionForceSolver &force_solver);
    void lab_frame(int n_sample, double gamma_e);
    void force_distribute(int n_sample, Beam &ion, Ions &ion_sample);
    void apply_kick(int n_sample, Beam &ion, Ions& ion_sample);
public:
    void adjust_rate(Beam &ion, EBeam &ebeam, initializer_list<double*> func);
    vector<double>& scratch(ECoolRateScratch s);
    double t_cooler(){return t_cooler_;}
    void set_n_long_sample(int n){n_long_sample_ = n;}
    void set_p_shift(bool b){p_shift_ = b;}
    void set_v_shift(bool b){v_shift_ = b;}
    void ecool_rate(FrictionForceSolver &force, Beam &ion, Ions &ptcl, Cooler &cooler, EBeam &ebeam,
                  Ring &ring, double &rate_x, double &rate_y, double &rate_s);
};

#endif // ECOOLING_H
