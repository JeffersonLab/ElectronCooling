#ifndef RMS_DYNAMIC_H_INCLUDED
#define RMS_DYNAMIC_H_INCLUDED

#include "simulator.h"
#include <vector>



using std::vector;

class RMSModel:public Simulator {
 private:
    void update_ibeam(Beam& ion, Ions& ion_sample, Ring& ring, EBeam& ebeam, Cooler& cooler, ECoolRate* ecool_solver);
    void adjust_rf_voltage(Ring& ring){if(ring.rf.gamma_tr>0) ring.rf.v = ring.calc_rf_voltage();};
    void save_ions(int i, Ions& ion_sample){};
 public:
    using Simulator::Simulator;
};



#endif // RMS_DYNAMIC_H_INCLUDED
