#ifndef OTHER_EFFECTS_H
#define OTHER_EFFECTS_H

#include "beam.h"
#include "cooler.h"
#include "ions.h"

void edge_effect(EBeam& ebeam,Beam& ion, Ions& ion_sample, Cooler& cooler, double dt);
#endif // OTHER_EFFECTS_H
