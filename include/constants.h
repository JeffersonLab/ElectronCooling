#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>

const double k_ke = 8.9875517873681764E9;   //Coulomb's constant, ke= 1/(4*pi*epsilon_0), in N*m^2/C^2
const double k_c = 299792458.0;             //speed of light, in m/s
const double k_e = 1.602176565E-19;         //Elementary charge, in C
const double k_pi = 3.1415926535897932384626;
const double k_u = 931.49406121;            //Atomic mass unit, in MeV/c^2
const double k_me = 0.510998928;            //electron mass, in MeV/c^2
const std::string k_whitespace = " \n\r\t\v\f"; //\n newline, \r carriage return, \t horizontal tab, \v vertical tab, \f form feed.
const double k_re = 2.8179403227E-15;		//Classical electron radius, k_re = k_ke * (k_e^2)/(k_me_kg * k_c * k_c)
const double k_me_kg = 9.1938356e-31;        //electron mass in kg
const double k_N_eVm = 6.242e18;            //Convert Newtons to eV/m


#endif  //CONSTANTS_H
