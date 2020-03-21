#ifndef BEAM_H
#define BEAM_H

#include "constants.h"
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>
#include "arbitrary_electron_beam.h"

class Beam{
    int charge_number_;   //Number of charges
    double mass_number_;       //mass = A*u [MeV/c^2]
    double mass_;    //unit in MeV/c^2
    double r_;       //classical radius, in m
    double kinetic_energy_;      //kinetic energy, in MeV
    double beta_;    //Lorentz factors
    double gamma_;   //Lorentz factors
    double emit_nx_; //normalized horizontal emittance, in m
    double emit_ny_; //normalized vertical emittance, in m
    double emit_x_;  //geometrical horizontal emittance, in m
    double emit_y_;  //geometrical vertical emittance, in m
    double dp_p_;     //momentum spread dp/p
    double energy_spread_;       // dE/E
    double sigma_s_; //RMS bunch length. set it to -1 for coasting beam, in m
    double particle_number_; //number of particles
    double p0_; //momentum in kg*m/s
    bool bunched_;   //Return true if beam is bunched.
    double center_[3] = {0,0,0};

public:
    int set_emit_nx(double x){emit_nx_ = x; emit_x_ = emit_nx_/(beta_*gamma_); return 0;}
    int set_emit_ny(double x){emit_ny_ = x; emit_y_ = emit_ny_/(beta_*gamma_); return 0;}
    int set_emit_x(double x){emit_x_ = x; emit_nx_ = beta_*gamma_*emit_x_; return 0;}
    int set_emit_y(double x){emit_y_ = x; emit_ny_ = beta_*gamma_*emit_y_; return 0;}
    int set_dp_p(double x){dp_p_ = x; energy_spread_ = beta_*beta_*dp_p_; return 0;}
    int set_sigma_s(double x){sigma_s_ = x; return 0;}
    int set_center(double cx, double cy, double cz){center_[0] = cx; center_[1] = cy; center_[2] = cz; return 0;}
    int set_center(int i, double x);
    int charge_number() const {return charge_number_;}
    double mass() const {return mass_;}
    double kinetic_energy() const {return kinetic_energy_;}
    double beta() const {return beta_;}
    double gamma() const {return gamma_;}
    double emit_nx() const {return emit_nx_;}
    double emit_ny() const {return emit_ny_;}
    double emit_x() const {return emit_x_;}
    double emit_y() const {return emit_y_;}
    double dp_p() const {return dp_p_;}
    double energy_spread() const {return energy_spread_;}
    double sigma_s() const {return sigma_s_;}
    double r() const {return r_;}
    double particle_number() const {return particle_number_;}
    double mass_number() const {return mass_number_;}
    double mass_SI() const {return mass_*1e6*k_e;}
    double p0_SI() const{return p0_;}
    double p0() const{return beta_*gamma_*mass_;}  //Momentum in [MeV/c]
    bool bunched()const {return bunched_;}
    int center(double &cx, double &cy, double &cz){cx = center_[0]; cy = center_[1]; cz = center_[2]; return 0;}
    double center(int i){ if (i<3) return center_[i]; else perror("Error index for electron beam center!"); return 1.0;}
    Beam(int charge_number, double mass_number, double kinetic_energy, double emit_nx, double emit_ny, double dp_p,
        double sigma_s, double n_particle);
    Beam(int charge_number, double mass_number, double kinetic_energy, double emit_nx, double emit_ny, double dp_p,
        double n_particle);
};

enum class Shape {UNIFORM_CYLINDER, GAUSSIAN_BUNCH, UNIFORM_BUNCH, GAUSSIAN_CYLINDER, ELLIPTIC_UNIFORM_BUNCH,
    UNIFORM_HOLLOW, UNIFORM_HOLLOW_BUNCH, PARTICLE_BUNCH};

enum class Velocity {CONST, USER_DEFINE, SPACE_CHARGE, VARY, VARY_X, VARY_Y, VARY_Z}  ;
enum class Temperature {CONST, USER_DEFINE, SPACE_CHARGE, VARY, VARY_X, VARY_Y, VARY_Z}  ;
enum class EBeamV {TPR_TR, TPR_L, V_RMS_TR, V_RMS_L, V_AVG_X, V_AVG_Y, V_AVG_L};

class EBeam {
protected:
    double kinetic_energy_ = 0;
    double gamma_ = 1;
    double beta_ = 0;
    bool bunched_ = true;
    double center_[3] = {0,0,0};
    double neutralisation_ = 2;
    Velocity velocity_ = Velocity::CONST;
    Temperature temperature_ = Temperature::CONST;
    vector<double> tpr_t;
    vector<double> tpr_l;
    vector<double> v_rms_t;
    vector<double> v_rms_l;
    vector<double> v_avg_x;
    vector<double> v_avg_y;
    vector<double> v_avg_l;
public:
    virtual ~EBeam(){};
    Velocity velocity() const {return velocity_;}
    Temperature temperature() const {return temperature_;}
    int charge_number() const {return -1;}
    double mass() const {return k_me;}
    double mass_number() const {return k_me/k_u;}
    double mass_SI() const {return k_me*1e6*k_e;}
    double kinetic_energy() const {return kinetic_energy_;}
    double gamma() const {return gamma_;}
    double beta() const {return beta_;}
    bool bunched()const {return bunched_;}
    virtual Shape shape() const = 0;
    virtual double length() const = 0;
    double neutral() const {return neutralisation_;}
    void set_kinetic_energy(double ke){kinetic_energy_ = ke; gamma_ = 1 + ke/k_me;
            beta_ = sqrt(gamma_*gamma_-1)/gamma_;}
    void set_gamma(double g){gamma_ = g; beta_ = sqrt(g*g-1)/g; kinetic_energy_ = (g-1)*k_me;}
    void set_beta(double b){beta_ = b; gamma_ = 1/sqrt(1-b*b); kinetic_energy_ = (gamma_-1)*k_me;}
    int set_center(double cx, double cy, double cz){center_[0] = cx; center_[1] = cy; center_[2] = cz; return 0;}
    void center(double &cx, double &cy, double &cz) const {cx = center_[0]; cy = center_[1]; cz = center_[2];}
    double center(int i) const { if (i<3&&i>-1) return center_[i]; else perror("Error index for electron beam center!"); return 1.0;}
    void set_center(int i, double x){if(i<3&&i>-1) center_[i]=x; else perror("Error index for electron beam center!");}
    void set_tpr(double tpr_tr, double trp_long);
    void set_v_rms(double v_rms_tr, double v_rms_long);
    void set_v_avg(double v_avg_tx, double v_avg_ty, double v_avg_long);
    void set_neutral(double x){neutralisation_ = x;}
    vector<double>& get_v(EBeamV v);
    virtual void density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n) = 0;
    virtual void density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n,
                         double cx, double cy, double cz) = 0;
};

class UniformCylinder: public EBeam{
    double current_;                   //Current of the beam in A
    double radius_;              //Radius of the beam in meter
 public:
    void density (vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n);
    void density (vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n,
                double cx, double cy, double cz);
    double current() const {return current_;}
    double radius() const {return radius_;}
    Shape shape() const {return Shape::UNIFORM_CYLINDER;}
    double length() const {perror("length() not defined for UniformCylinder, which is coasting"); return 0;}
    UniformCylinder(double current, double radius, double neutralisation=2):current_(current),radius_(radius)
                    {bunched_ = false;};
};

class UniformHollow: public EBeam {
    double current_;    //Peak current, the current as if the beam is coasting.
    double in_radius_;
    double out_radius_;
 public:
    void density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n);
    void density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n,
                double cx, double cy, double cz);
    double current() const {return current_;}
    double out_radius() const {return out_radius_;}
    double in_radius() const {return in_radius_;}
    Shape shape() const {return Shape::UNIFORM_HOLLOW;}
    double length() const {perror("length() not defined for UniformHollow, which is coasting"); return 0;}
    UniformHollow(double current, double in_radius, double out_radius):current_(current),
        in_radius_(in_radius), out_radius_(out_radius){bunched_ = false;};
};


class UniformHollowBunch: public EBeam {
    double current_;
    double in_radius_;
    double out_radius_;
    double length_;
 public:
    void density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n);
    void density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n,
                double cx, double cy, double cz);
    double current() const {return current_;}
    double out_radius() const {return out_radius_;}
    double in_radius() const {return in_radius_;}
    Shape shape() const {return Shape::UNIFORM_HOLLOW_BUNCH;}
    double length() const {return length_;}
    UniformHollowBunch(double current, double in_radius, double out_radius, double length):current_(current),
        in_radius_(in_radius), out_radius_(out_radius), length_(length) {}
};


class UniformBunch: public EBeam{
    double current_;                   //Current of the beam in A, assuming the beam is DC.
    double radius_;              //Radius of the beam in meter
    double length_;
public:
    //Calculate the charge density for a given position (x,y,z) in Lab frame.
    void density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n);
    void density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n,
                double cx, double cy, double cz);
    Shape shape() const {return Shape::UNIFORM_BUNCH;}
    double length() const {return length_;}
    double current() const {return current_;}
    double radius() const {return radius_;}
    UniformBunch(double current, double radius, double length):current_(current),radius_(radius),
            length_(length){};

};


class EllipticUniformBunch: public EBeam{
    double current_;
    double rh_;         //half horizontal axis
    double rv_;         //half vertical axis
    double length_;     //bunch length
public:
    //Calculate the charge density for a given position (x,y,z) in Lab frame.
    void density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n);
    void density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n,
                double cx, double cy, double cz);
    Shape shape() const {return Shape::ELLIPTIC_UNIFORM_BUNCH;}
    double length() const {return length_;}
    EllipticUniformBunch(double current, double rh, double rv, double length):current_(current),
            rh_(rh),rv_(rv),length_(length){};
};


class GaussianBunch: public EBeam{
    double n_electron_;
    double sigma_x_;
    double sigma_y_;
    double sigma_s_;
 public:
    void density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n);
    void density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n,
                double cx, double cy, double cz);
    Shape shape() const {return Shape::GAUSSIAN_BUNCH;}
    double length() const {return 6*sigma_s_;}
    GaussianBunch(double n_electron, double sigma_x, double sigma_y, double sigma_s):n_electron_(n_electron),
                sigma_x_(sigma_x),sigma_y_(sigma_y),sigma_s_(sigma_s){};

};

class ParticleBunch: public EBeam {
    double n_electron_;
    std::string filename_;
    long int n_ = 0;
    double length_ = 0;
    bool v_x_corr_ = false;    //Velocity position correlation
    int line_skip_ = 0;
    vector<Box> tree_;
    vector<long int> list_e_;
    int s_ = 100;
    bool binary_ = false;
    int buffer_ = 1000;
public:
    std::vector<double> x, y, z, vx, vy, vz;  //Electron phase space coordinates
    //Calculate the charge density for a given position (x,y,z) in Lab frame.
    void density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n);
    void density(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n,
                double cx, double cy, double cz);
    Shape shape() const {return Shape::PARTICLE_BUNCH;}
    double length() const {return length_;}
    bool bunched(){return true;}
    bool corr() const {return v_x_corr_;}
    void set_corr(bool corr = true){v_x_corr_ = corr;}
    void set_buffer(int n) {buffer_ = n;}
    void set_s(int s) {s_ = s;}
    void set_binary(bool b) {binary_ = b;}
    void set_skip(int n) {line_skip_ = n;}

    ParticleBunch(double n_electron, std::string filename, double length):n_electron_(n_electron),
        filename_(filename),length_(length){temperature_ = Temperature::VARY;};
    ParticleBunch(double n_electron, std::string filename):n_electron_(n_electron),
        filename_(filename){temperature_ = Temperature::VARY;};
    void load_particle(long int n);
    void load_particle();

};

class MultiBunches: public EBeam{
    int n_; //Number of bunches
    vector<double> cx_;     //List of cxs.
    vector<double> cy_;     //List of cys.
    vector<double> cz_;     //List of czs.
    EBeam* bunches_;
 public:
    EBeam* bunch(){return bunches_;}
    vector<double>& cx(){return cx_;}
    vector<double>& cy(){return cy_;}
    vector<double>& cz(){return cz_;}
    MultiBunches(int n):n_(n){cx_.resize(n); cy_.resize(n); cz_.resize(n);}
    ~MultiBunches(){delete bunches_;}
    MultiBunches(const MultiBunches& obj) = delete;
    MultiBunches& operator=(const MultiBunches& obj) = delete;
    Shape shape() const {return bunches_->shape();}
    double length() const {return bunches_->length();}
    void density (vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n);
    void density (vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n,
                double cx, double cy, double cz);
};

#endif // BEAM_H
