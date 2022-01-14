#ifndef TURN_BY_TURN_H_INCLUDED
#define TURN_BY_TURN_H_INCLUDED

#include <fstream>
#include <string>

#include <vector>
#include "beam.h"
#include "cooler.h"
#include "ring.h"
#include "particle_model.h"

void initialize_turn_by_turn_model(Beam &ion, Ring &ring);
void turn_by_turn_move_particles(Beam &ion, Ring &ring, Cooler &cooler);

class TurnByTurnModel: public ParticleModel {
protected:
    void move_particles(Beam& ion, Ions& ion_sample, Ring& ring);
    void apply_edge_kick(Cooler& cooler, EBeam& ebeam, Beam& ion, Ions& ion_sample, ECoolRate* ecool_solver);
    ofstream out_single_particle;
    string filename_single_particle = "single_particle_";
    int idx = -1;

public:
    using ParticleModel::ParticleModel;
};

#endif // TURN_BY_TURN_H_INCLUDED
