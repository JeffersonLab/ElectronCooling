#ifndef ARBITRARY_ELECTRON_BEAM_H
#define ARBITRARY_ELECTRON_BEAM_H

#ifndef BOX_H
#define BOX_H

#include <iostream>
#include <vector>

using std::vector;
using std::cout;
using std::endl;

typedef struct Box{
	double center[3] = {0,0,0};
	long int parent = 0;
	long int child[8] = {0};
	long int first_ptcl = 0;
	int n_ptcl = 0;
	int n_child = 0;
	double box_size = 0;
	long int first_ion = 0;
	int n_ion = 0;
} Box;

typedef struct Colleague{
	long int clg[28] = {0}; //At most 27 colleagues, the last 0 means the end.
} Colleague;

long int load_electrons(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector<double>& vx,
                                 std::vector<double>& vy, std::vector<double>& vz, std::string filename, long int n,
                                int skip = 0, bool binary = false, int n_buffer = 1000);

int create_e_tree(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, const unsigned long int n,
                  const unsigned int s, vector<Box> &tree, std::vector<unsigned long int>& list);

int create_ion_tree(double * x, double * y, double * z, const unsigned int n, vector<Box> &tree,
                std::vector<unsigned int>& list, unsigned int &idx_out);

std::ostream& operator<<(std::ostream& os, Box& box);
std::ostream& operator<<(std::ostream& os, Colleague& clg);

void create_e_tree(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, const long int n,
                  const int s, vector<Box> &tree, std::vector<long int>& list);

void create_ion_tree(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, const int n, std::vector<Box> &tree,
                std::vector<int>& list, int &idx_out);

void density(vector<Box> &tree, std::vector<long int>& list_e, std::vector<double>& vx, std::vector<double>& vy,
             std::vector<double>& vz, const long int ne,  std::vector<int>& list_i,
             int idx_out, const int ni, std::vector<double>& density_e,
             std::vector<double>& v_rms_t, std::vector<double>& v_rms_l) ;

void density(vector<Box> &tree, std::vector<long int>& list_e, std::vector<double>& vx, std::vector<double>& vy,
             std::vector<double>& vz, const long int ne,  std::vector<int>& list_i, int idx_out,
             const int ni, std::vector<double>& density_e, std::vector<double>& v_avg_z, std::vector<double>& v_rms_t,
             std::vector<double>& v_rms_l);


#endif

#endif // ARBITRARY_ELECTRON_BEAM_H
