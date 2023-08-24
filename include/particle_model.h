#ifndef PARTICLE_MODEL_HPP
#define PARTICLE_MODEL_HPP

#include "simulator.h"
#include <vector>
#include "beam.h"
#include "cooler.h"
#include "ring.h"

class ParticleModel: public Simulator {
 protected:
    void update_ibeam(Beam& ion, Ions& ion_sample, Ring& ring, EBeam& ebeam, Cooler& cooler, ECoolRate* ecool_solver);
    void apply_cooling_kick(double freq, Beam& ion, Ions& ion_model, ECoolRate* ecool_solver);
    void apply_ibs_kick(Beam& ion, Ions& ion_sample);
    void ibs_kick(int n_sample, double rate, double twiss, double emit, vector<double>& v);
    virtual void move_particles(Beam& ion, Ions& ion_sample, Ring& ring);
    virtual void apply_edge_kick(Cooler& cooler, EBeam& ebeam, Beam& ion, Ions& ion_sample, ECoolRate* ecool_solver){};
    void update_beam_parameters(Beam &ion, Ions& ion_sample);
    void adjust_rf_voltage(Ring& ring){};
    void save_ions(int i, Ions& ion_sample);
    vector<double> rdn;
    void precondition(Ions& ion_sample){resize_rdn(ion_sample.n_sample());};
 public:
//    using Simulator::Simulator;
    void resize_rdn(int n_sample){rdn.resize(n_sample);}
    ParticleModel(double time, int n):Simulator(time, n){}
};

#endif // IBS_HPP
