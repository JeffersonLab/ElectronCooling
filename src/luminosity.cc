#include "luminosity.h"
#include "math.h"
#include "constants.h"

double LuminositySolver::luminosity() {
    double sig_x2 = beam.at(0).sigma_x*beam.at(0).sigma_x + beam.at(1).sigma_x*beam.at(1).sigma_x;
    double sig_y2 = beam.at(0).sigma_y*beam.at(0).sigma_y + beam.at(1).sigma_y*beam.at(1).sigma_y;
    double lum = beam.at(0).np*beam.at(1).np*freq_/(2*k_pi*sqrt(sig_x2*sig_y2))*exp(-dx_*dx_/(2*sig_x2)-dy_*dy_/(2*sig_y2));
    return lum;
}

void LuminositySolver::set_geo_emit(double emit_x, double emit_y, int i) {
    beam.at(i).geo_emit_x=emit_x;
    beam.at(i).geo_emit_y=emit_y;
    beam.at(i).sigma_x = sqrt(emit_x*beam.at(i).bet_x_star);
    beam.at(i).sigma_y = sqrt(emit_y*beam.at(i).bet_y_star);
}

void LuminositySolver::set_beam_size(double sigma_x, double sigma_y, int i) {
    beam.at(i).sigma_x=sigma_x;
    beam.at(i).sigma_y=sigma_y;
    beam.at(i).geo_emit_x = sigma_x*sigma_x/beam.at(i).bet_x_star;
    beam.at(i).geo_emit_y = sigma_y*sigma_y/beam.at(i).bet_y_star;
};

