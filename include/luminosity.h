#ifndef LUMINOSITY_H
#define LUMINOSITY_H

#include <vector>

struct CollidingBeam {
    double np = 0;
    double sigma_x = 0;
    double sigma_y = 0;
    double geo_emit_x = 0;
    double geo_emit_y = 0;
    double bet_x_star = 0;
    double bet_y_star = 0;
    double bet_x_max = 0;
    double bet_y_max = 0;
    double focus_length_x = 0;
    double focus_length_y = 0;
    double aper_x = 0;
    double aper_y = 0;
    bool adjust_bet = false;
};

class LuminositySolver {
    std::vector<CollidingBeam> beam = std::vector<CollidingBeam>(2);
    double dx_ = 0;
    double dy_ = 0;
    double freq_ = 1;
    double aper_ratio_ = 10;
    bool match_= false;
    bool use_ion_emit_ = false;
public:
    void set_distance(double dx, double dy){dx_=dx; dy_=dy;}
    void set_freq(double f){freq_=f;}
    void set_match(bool b){match_ = b;}
    void set_geo_emit(double emit_x, double emit_y, int i);
    void set_beam_size(double sigma_x, double sigma_y, int i);
    void set_particle_number(double n, int i){beam.at(i).np=n;}
    void set_bet_star(double bet_x, double bet_y, int i){beam.at(i).bet_x_star=bet_x; beam.at(i).bet_y_star=bet_y;}
    void set_bet_max(double bet_x, double bet_y, int i){beam.at(i).bet_x_max=bet_x; beam.at(i).bet_y_max=bet_y;}
    bool match(){return match_;}
    void set_focus_length(double x, double y, int i){beam.at(i).focus_length_x=x; beam.at(i).focus_length_y=y;}
    void set_aperture(double x, double y, int i){beam.at(i).aper_x=x; beam.at(i).aper_y=y;};
    void calc_bet_max(int i);
    void adjust_bet(int i);
    void set_adjust_bet(bool b, int i){beam.at(i).adjust_bet=b;}
    void set_aper_ratio(double i){aper_ratio_ = i;}
    void set_use_ion_emit(bool b){use_ion_emit_ = b;}
    bool use_ion_emit(){return use_ion_emit_;}
    void match(int i); //match the size of beam i to the other beam.
    double bet_x_star(int i){return beam.at(i).bet_x_star;}
    double bet_y_star(int i){return beam.at(i).bet_y_star;}
    double bet_x_max(int i){return beam.at(i).bet_x_max;}
    double bet_y_max(int i){return beam.at(i).bet_y_max;}
    double luminosity();
};

#endif // LUMINOSITY_H
