#include "../include/simulator.h"
#include "../include/functions.h"
#include "../include/luminosity.h"
#include "../include/ui.h"

extern Record uircd;

//ECoolRate* ecool_solver = nullptr;
//FrictionForceSolver* force_solver = nullptr;
//LuminositySolver* lum_solver = nullptr;
//Simulator* simulator = nullptr;
//extern IBSSolver* ibs_solver;



void Simulator::output_sddshead() {
    using std::endl;
    outfile<<"SDDS1"<<endl;
    outfile<<"! Define colums:"<<endl
        <<"&column name=t, type=double, units=s, description=time, &end"<<endl
        <<"&column name=emit_x, type=double, units=m*rad, description=\"normalized horizontal emittance\", &end"<<endl
        <<"&column name=emit_y, type=double, units=m*rad, description=\"normalized vertical emittance\", &end"<<endl
        <<"&column name=dp/p, type=double, units=NULL, description=\"momentum spread\", &end"<<endl
        <<"&column name=sigma_s, type=double, units=m, description=\"RMS bunch length\", &end"<<endl
        <<"&column name=rx, type=double, units=1/s, description=\"horizontal expansion rate\", &end"<<endl
        <<"&column name=ry, type=double, units=1/s, description=\"vertical expansion rate\", &end"<<endl
        <<"&column name=rs, type=double, units=1/s, description=\"longitudinal expansion rate\", &end"<<endl
        <<"&column name=rx_ibs, type=double, units=1/s, description=\"horizontal IBS expansion rate\", &end"<<endl
        <<"&column name=ry_ibs, type=double, units=1/s, description=\"vertical IBS expansion rate\", &end"<<endl
        <<"&column name=rs_ibs, type=double, units=1/s, description=\"longitudinal IBS expansion rate\", &end"<<endl
        <<"&column name=rx_ecool, type=double, units=1/s, description=\"horizontal electron cooling rate\", &end"<<endl
        <<"&column name=ry_ecool, type=double, units=1/s, description=\"vertical electron cooling rate\", &end"<<endl
        <<"&column name=rs_ecool, type=double, units=1/s, description=\"longitudinal electron cooling rate\", &end"<<endl
        <<"&column name=rf_voltage, type=double, units=V, description=\"Voltage of the RF cavity\", &end"<<endl
        <<"&column name=luminosity, type=double, units=1/s*1/cm^2, description=\"Instant luminosity\", &end"<<endl
        <<"!Declare ASCII data and end the header"<<endl
        <<"&data mode=ascii, &end"<<endl
        <<n_step+1<<endl;
}

void Simulator::output_to_file() {
    if(!overwrite && file_exists(outfilename)) {
        string filename = outfilename;
        int i = 1;
        do {
            filename = std::to_string(i)+'_'+filename;
            ++i;
        } while(file_exists(filename));
        outfile.open(filename);
    }
    else {
        outfile.open(outfilename);
    }
    output_sddshead();
    outfile.precision(10);
    outfile<<std::showpos;
    outfile<<std::scientific;
}

void Simulator::output(bool bunched, double v_rf, double lum) {
    outfile<<t<<' '<<emit.at(0)<<' '<<emit.at(1)<<' '<<emit.at(2)<<' ';
    if(bunched) outfile<<emit.at(3)<<' ';
    else outfile<<0<<' ';
    outfile<<r.at(0)<<' '<<r.at(1)<<' '<<r.at(2)<<' ';
    outfile<<r_ibs.at(0)<<' '<<r_ibs.at(1)<<' '<<r_ibs.at(2)<<' ';
    outfile<<r_ecool.at(0)<<' '<<r_ecool.at(1)<<' '<<r_ecool.at(2)<<' ';
    outfile<<v_rf<<' '<<lum<<' ';
    outfile<<std::endl;
}

void Simulator::run(Beam& ion, Ions& ion_sample, Cooler& cooler, EBeam& ebeam,
                    Ring& ring, IBSSolver* ibs_solver, ECoolRate* ecool_solver,
                    FrictionForceSolver* force_solver, LuminositySolver* lum_solver) {
    precondition(ion_sample);
    //Set time for new simulation or continued simulation.
    if(reset_time) t = t0;
    else t = uircd.t;
    output_to_file();

    adjust_rf_voltage(ring);

for(int i=0; i<n_step+1; ++i) {
        save_ions(i, ion_sample);
        //record
        emit.at(0) = ion.emit_nx();
        emit.at(1) = ion.emit_ny();
        emit.at(2) = ion.dp_p();
        if (ion.bunched()) emit.at(3) = ion.sigma_s();

        if(ibs) {
            ibs_solver->rate(ring.lattice(), ion, r_ibs.at(0), r_ibs.at(1), r_ibs.at(2));
        }

        if(ecool) {
            ecool_solver->ecool_rate(*force_solver, ion, ion_sample, cooler, ebeam, ring, r_ecool.at(0), r_ecool.at(1),
                                     r_ecool.at(2));
        }

        for(int i=0; i<3; ++i) r.at(i) = r_ibs.at(i) + r_ecool.at(i);

        if (output_itvl==1 || i%output_itvl==0) {
            double lum = 0;
            if(calc_luminosity) {
                if(lum_solver->use_ion_emit()) lum_solver->set_geo_emit(ion.emit_x(), ion.emit_y(), 0);
                lum = lum_solver->luminosity();
            }
            output(ion.bunched(), ring.rf.v, lum*10000);
        }

        update_ibeam(ion, ion_sample, ring, ebeam, cooler, ecool_solver);

        t += dt;
        std::cout<<i<<std::endl;
    }
    save_ions(n_step, ion_sample);

    uircd.emit_nx = emit.at(0);
    uircd.emit_ny = emit.at(1);
    uircd.dp_p = emit.at(2);
    if(ion.bunched()) uircd.sigma_s = ion.sigma_s();
    else uircd.sigma_s = 0;
    uircd.rx_ecool = r_ecool.at(0);
    uircd.ry_ecool = r_ecool.at(1);
    uircd.rs_ecool = r_ecool.at(2);
    uircd.rx_ibs = r_ibs.at(0);
    uircd.ry_ibs = r_ibs.at(1);
    uircd.rs_ibs = r_ibs.at(2);
    uircd.rx_total = r.at(0);
    uircd.ry_total = r.at(1);
    uircd.rs_total = r.at(2);
    uircd.t = t;

    //Need to add: save the results before quiting the simulation.
    outfile.close();
    std::cout<<"Finished dynamic simulation."<<std::endl;
}
