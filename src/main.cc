#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <map>


#include "../include/beam.h"
#include "../include/constants.h"
#include "../include/cooler.h"
#include "../include/ecooling.h"
#include "../include/force.h"
#include "../include/ibs.h"
#include "../include/ions.h"
#include "../include/luminosity.h"
#include "../include/particle_model.h"
#include "../include/ring.h"
#include "../include/rms_dynamic.h"
#include "../include/simulator.h"
#include "../include/turn_by_turn.h"
#include "../include/ui.h"

using std::string;

extern std::map<std::string, Section> sections;
extern muParserHandle_t math_parser;

enum class Test {IBS, ECOOL, BOTH,DYNAMICBOTHBUNCHED, MATH_PARSER};

int main(int argc, char** argv) {
    srand(time(NULL));

    if(argc>1) {
        std::ifstream input_file(argv[1]);
        string line;
        Section sec_flag = Section::NONE;
        Set_ptrs ptrs;
        if (math_parser == NULL) {
            initialize_parser(math_parser);
        }
        while (std::getline(input_file,line)) {
            if(!line.empty() && line[line.size()-1] == '\r') line.erase(line.size()-1);
            if (!line.empty()) {
                line = remove_comments(line);
                line = trim_blank(line);
                line = trim_tab(line);
                string line_orgn = line;
                str_toupper(line);
                if (!line.empty()) {
                    if (sections.find(line)!=sections.end()) {
                        sec_flag = sections[line];
                        switch (sec_flag) {
                            case Section::SECTION_ION: {
                                if (ptrs.ion_ptr.get() == nullptr) ptrs.ion_ptr.reset(new Set_ion());
                                break;
                            }
                            case Section::SECTION_RING: {
                                if (ptrs.ring_ptr.get() == nullptr) ptrs.ring_ptr.reset(new Set_ring());
                                break;
                            }
                            case Section::SECTION_IBS: {
                                if (ptrs.ibs_ptr.get() == nullptr) ptrs.ibs_ptr.reset(new Set_ibs());
                                break;
                            }
                            case Section::SECTION_COOLER: {
                                if (ptrs.cooler_ptr.get() == nullptr) ptrs.cooler_ptr.reset(new Set_cooler());
                                break;
                            }
                            case Section::SECTION_E_BEAM: {
                                if (ptrs.e_beam_ptr.get() == nullptr) ptrs.e_beam_ptr.reset(new Set_e_beam());
                                break;
                            }
                            case Section::SECTION_ECOOL: {
                                if (ptrs.ecool_ptr.get() == nullptr) ptrs.ecool_ptr.reset(new Set_ecool());
                                break;
                            }
                            case Section::SECTION_LUMINOSITY: {
                                if (ptrs.luminosity_ptr.get() == nullptr) ptrs.luminosity_ptr.reset(new Set_luminosity());
                                break;
                            }
                            case Section::SECTION_SIMULATION: {
                                if (ptrs.dynamic_ptr.get() == nullptr) ptrs.dynamic_ptr.reset(new Set_dynamic());
                                break;
                            }
                            case Section::SECTION_RUN: {
                                set_section_run(ptrs);
                                break;
                            }
                            case Section::SECTION_SCRATCH: {
                                break;
                            }
                            default : {
                                assert(false && "WRONG SECTION NAME!");
                                break;
                            }
                        }
                    }
                    else {
                        assert(sec_flag>Section::NONE && "Script file must start with a section!");
                        switch (sec_flag) {
                            case Section::SECTION_ION: {
                                define_ion_beam(line, ptrs.ion_ptr.get());
                                break;
                            }
                            case Section::SECTION_RING: {
                                define_ring(line_orgn, ptrs.ring_ptr.get());
                                break;
                            }
                            case Section::SECTION_IBS: {
                                set_ibs(line, ptrs.ibs_ptr.get());
                                break;
                            }
                            case Section::SECTION_COOLER: {
                                define_cooler(line, ptrs.cooler_ptr.get());
                                break;
                            }
                            case Section::SECTION_E_BEAM: {
                                define_e_beam(line_orgn, ptrs.e_beam_ptr.get());
                                break;
                            }
                            case Section::SECTION_ECOOL: {
                                set_ecool(line, ptrs.ecool_ptr.get());
                                break;
                            }
                            case Section::SECTION_LUMINOSITY: {
                                set_luminosity(line, ptrs.luminosity_ptr.get());
                                break;
                            }
                            case Section::SECTION_SIMULATION: {
                                set_simulation(line, ptrs.dynamic_ptr.get());
                                break;
                            }
                            case Section::SECTION_RUN: {
                                run(line, ptrs);
                                break;
                            }
                            case Section::SECTION_SCRATCH: {
                                parse(line, math_parser);
                                break;
                            }
                        }

                    }
                }
            }
        }


    }
    else {
        Test test = Test::MATH_PARSER;
        switch (test) {
            case Test::MATH_PARSER: {
                break;
            }
            case Test::IBS: {
                //********************************
                // Test IBS rate
                //********************************
                // define proton beam;
                double m0, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl;
                int Z;
                Z = 1;
                m0 = 938.272;
                KE = 8000;
                emit_nx0 = 2.2e-6;
                emit_ny0 = 2.2e-6;
                dp_p0 = 0.001*0.6;
                N_ptcl = 6.58E11;
                sigma_s0 = 7;
                Beam p_beam(Z,m0/k_u, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl);

                 // define the lattice of the proton ring
                std::string filename = "test.tfs";
                Lattice lattice(filename);

                //Define the ring
                Ring ring(lattice, p_beam);

                //Set IBS parameters.
                int nu = 100;
                int nv = 100;
                int nz = 40;
                double log_c = 39.2/2;
                IBSSolver* ibs_solver;
                ibs_solver = new IBSSolver_Martini(nu, nv, nz, log_c, 0);
        //
                double rx_ibs, ry_ibs, rz_ibs;

                ibs_solver->rate(lattice, p_beam, rx_ibs, ry_ibs, rz_ibs);
                std::cout<<rx_ibs<<' '<<ry_ibs<<' '<<rz_ibs<<std::endl;

                ibs_solver->set_k(0.2);
                ibs_solver->rate(lattice, p_beam, rx_ibs, ry_ibs, rz_ibs);
                std::cout<<rx_ibs<<' '<<ry_ibs<<' '<<rz_ibs<<std::endl;
                break;
            }
            case Test::ECOOL: {
                //********************************
                // Test Electron Cooling Rate
                //********************************

                //define coasting 12C6+ beam
                int n_charge = 6;
                double n_mass = 12;
                double kinetic_energy = 30*n_mass;
                double gamma = 1+kinetic_energy/(n_mass*k_u);
                double beta = sqrt(1-1/(gamma*gamma));
                double emit_nx0 = beta*gamma*5e-6;
                double emit_ny0 = beta*gamma*5e-6;
                double dp_p0 = 0.0004;
                double n_ptcl = 5E8;
                Beam c_beam(n_charge, n_mass, kinetic_energy, emit_nx0, emit_ny0, dp_p0, n_ptcl);

                //define the ring
                std::string filename = "csrm.tfs";
                Lattice lattice(filename);
                Ring ring(lattice, c_beam);

                //define the cooler
                double cooler_length = 3.4;
                double n_section = 1;
                double magnetic_field = 0.039;
                double beta_h = 10;
                double beta_v = 17;
                double dis_h = 0;
                double dis_v = 0;
                Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);

                //define electron beam
                double current = 0.03;
                double radius = 0.025;
                double radius_in = 0.005;
    //                UniformCylinder uniform_cylinder(current, radius);
                UniformHollow uni_hl(current, radius_in, radius);
                double gamma_e = c_beam.gamma();
                double tmp_tr = 0.05;
                double tmp_long = 0.1;

                uni_hl.set_gamma(gamma_e);
                uni_hl.set_tpr(tmp_tr, tmp_long);

                //define cooling model
                int n_sample = 10000;
                double rate_x, rate_y, rate_s;
                Ions_MonteCarlo c_sample(n_sample);
                c_sample.set_twiss(cooler);
                c_sample.create_samples(c_beam);
                ECoolRate ecool;
                Force_Park force_park;

                ecool.ecool_rate(force_park, c_beam, c_sample, cooler, uni_hl, ring, rate_x, rate_y, rate_s);
                std::cout<<"rate_x = "<<rate_x<<" rate_y = "<<rate_y<<" rate_s = "<<rate_s<<std::endl;

                break;
            }
            case Test::DYNAMICBOTHBUNCHED: {
                    // define proton beam;
                    double m0, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl;
                    int Z;
                    Z = 1;
                    m0 = 938.272;
                    KE = 3e4;
                    emit_nx0 = 0.4962695094e-6;
                    emit_ny0 = 0.4962695094e-6;
                    dp_p0 = 4e-4;
                    sigma_s0 = 1.994525702e-2;
                    N_ptcl = 6.56E9;
                    Beam p_beam(Z,m0/k_u, KE, emit_nx0, emit_ny0, dp_p0, sigma_s0, N_ptcl);

                     // define the lattice of the proton ring
                    std::string filename = "test.tfs";
                    Lattice lattice(filename);

                    //Define the ring
                    Ring ring(lattice, p_beam);

                    //define the cooler
                    double cooler_length = 60;
                    double n_section = 1;
                    double magnetic_field = 1;
                    double beta_h = 100;
                    double beta_v = 100;
                    double dis_h = 0.5;
                    double dis_v = 0;
                    Cooler cooler(cooler_length,n_section,magnetic_field,beta_h,beta_v,dis_h, dis_v);

                    //define electron beam
                    double sigma_x = 0.0012;
                    double sigma_y = 0.0012;
                    double sigma_s = 0.025;
                    double ne = 1.25E10;
                    GaussianBunch e_beam(ne, sigma_x, sigma_y, sigma_s);
                    double gamma_e = p_beam.gamma();
                    double tmp_tr = 0.5;
                    double tmp_long = 0.1;
                    e_beam.set_gamma(gamma_e);
                    e_beam.set_tpr(tmp_tr, tmp_long);
                    FrictionForceSolver* force_solver = nullptr;
                    ECoolRate* ecool_solver = nullptr;
                    IBSSolver* ibs_solver = nullptr;
                    LuminositySolver* lum_solver = nullptr;
                    force_solver = new Force_Park();
                    ecool_solver = new ECoolRate();

                    //Set IBS parameters.
                    int nu = 200;
                    int nv = 200;
                    int nz = 40;
                    double log_c = 39.9/2;
                    ibs_solver = new IBSSolver_Martini(nu, nv, nz, log_c, 0);

                    double time = 60;
                    int n_step = 600;
                    ParticleModel simulator(time, n_step);
                    simulator.set_ion_save(100);
    //                RMSModel simulator(time, n_step);
                    simulator.set_output_file("test_SIMULATION_particle_BUNCHED_both.txt");
                    simulator.set_ibs(true);
                    simulator.set_ecool(true);

                    Ions_MonteCarlo ion_sample(10000);
                    ion_sample.set_twiss(cooler);
                    ion_sample.create_samples(p_beam);

                    simulator.resize_rdn(ion_sample.n_sample());
                    simulator.run(p_beam, ion_sample, cooler, e_beam, ring, ibs_solver, ecool_solver, force_solver, lum_solver);


                    break;

                }
            default: {
                assert(false&&"Wrong test type selected!");
            }
        }
    }

    return 0;
}
