#ifndef FORCE_H
#define FORCE_H

#include <vector>
#include "constants.h"
#include "beam.h"

using std::vector;
enum class ForceFormula {PARKHOMCHUK};

class FrictionForceSolver{
protected:
    double time_cooler;
    double mag_field = 0;
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
    double t_eff = 0; //Effective temperature.
    double v_eff = 0; //Effective velocity.
public:
    void set_t_eff(double x){t_eff = x; v_eff = sqrt(t_eff*k_c*k_c/(k_me*1e6));}
    void set_v_eff(double v){v_eff = v; t_eff = v_eff*v_eff*k_me*1e6/(k_c*k_c);}
    virtual void friction_force(int charge_number, int ion_number,
            vector<double>& v_tr, vector<double>& v_l, vector<double>& density,
            EBeam& ebeam, vector<double>& force_tr, vector<double>& force_long);
};

#endif // FORCE_H
