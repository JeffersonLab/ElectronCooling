#include "../include/rms_dynamic.h"
#include <chrono>
#include <cmath>
#include "../include/constants.h"
#include "../include/functions.h"

void RMSModel::update_ibeam(Beam& ion, Ions& ion_sample, Ring& ring, EBeam& ebeam, Cooler& cooler, ECoolRate* ecool_solver) {
    double emit_nx = ion.emit_nx();
    double emit_ny = ion.emit_ny();
    double dp = ion.dp_p();
    emit_nx *= exp(r.at(0)*dt);
    emit_ny *= exp(r.at(1)*dt);
    dp *= dp*exp(r.at(2)*dt);
    dp = sqrt(dp);

    ion.set_emit_nx(emit_nx);
    ion.set_emit_ny(emit_ny);
    ion.set_dp_p(dp);

    if(ion.bunched()) {
        if(fixed_bunch_length) {
            ring.rf.v = ring.calc_rf_voltage();
        }
        else {
            ion.set_sigma_s(ring.beta_s()*dp);
        }
    }

    if(ecool) ion_sample.create_samples(ion);
}

