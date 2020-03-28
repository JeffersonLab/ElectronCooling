#include "other_effects.h"
#include <vector>

using std::vector;

void edge_effect(EBeam& ebeam, Beam& ion, Ions& ion_sample, Cooler& cooler, double dt) {
    vector<double> ez;
    int n_sample = ion_sample.n_sample();
    ez.resize(n_sample);
    ebeam.edge_field(cooler, ion_sample.cdnt(Phase::X), ion_sample.cdnt(Phase::Y), ion_sample.cdnt(Phase::DS), ez, n_sample);
    vector<double>& dp_p = ion_sample.cdnt(Phase::DP_P);
    double p0 = ion.p0_SI();
    double inv_p0 = 1/p0;
    for(int i=0; i<n_sample; ++i) {
        dp_p.at(i) = (dp_p.at(i)*p0 + ez.at(i)*dt)*inv_p0;
    }
}
