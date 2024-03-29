/**
 * @file isb.h
 * @brief Calculates IBS expansion rate.
 * @details Defines the IBSSolver class to handle IBS rate calculation.
 * @author He Zhang
 * @email hezhang@jlab.org
 */

#ifndef IBS_HPP
#define IBS_HPP

#include <array>
#include <assert.h>
#include <memory>
#include <tuple>
#include <vector>
#include "beam.h"
#include "ring.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <force.h>

class Lattice;
class Beam;

/**
 * @enum IBSModel
 * @brief Choose different formulas for IBS rate calculation.
 * @see IBSSolver_Martini, IBSSolver_BM, IBSSolver_BMZ, and IBSSolver_BM_Complete.
 */
enum class IBSModel {MARTINI, BM, BMC, BMZ};

/**
 * @brief Base class for IBS rate calculators.
 */
class IBSSolver {
protected:
    double log_c_ = 0.0;     //Coulomb logarithm.
    double k_ = 0.0;          //Coupling rate in transverse directions.
    bool cache_invalid = true;
    bool ibs_by_element = false; //Calculate and output the ibs rate contribution element by element.

    void ibs_coupling(double &rx, double &ry, double k, double emit_x, double emit_y);
    void ibs_by_element_sddshead(std::ofstream& outfile, int n_element);
public:
    double log_c() const { return log_c_; } ///< Get the value of the Coulomb logarithm.
    double k() const { return k_; } ///< Get the value of the transverse coupling rate.
    void set_k(double x) { k_ = x; } ///< Set the value of the transverse coupling rate.
    void set_log_c(double x) { log_c_ = x; } ///< Set the value of the Coulomb logarithm.
    /**
     * @brief Choose whether to calculate the contribution to IBS from each individual element.
     */
    void set_ibs_by_element(bool b) {ibs_by_element = b;}
    /**
     * @brief Call it if the lattice or the beam energy changed. Cached intermediate results will be calculated again.
     */
    void invalidate_cache() { cache_invalid = true; }

    /**
     * @brief Constructor of IBSSolver.
     * \param[in] log_c Coulomb logarithm
     * \param[in] k The transverse coupling rate, from 0 to 1.
     */
    IBSSolver(double log_c, double k);

    /**
     * @brief Calculate the IBS expansion rate.
     * This is the virtual function in the base class to calculate the IBS expansion rate of the ion beam.
     * Should be overloaded in the derived classes.
     * \param[in] lattice The lattice of the ion ring.
     * \param[in] beam The ion beam.
     * \param[out] rx The horizontal IBS rate.
     * \param[out] ry The vertical IBS rate.
     * \param[out] rs The longitudinal IBS rate.
     */
    virtual void rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs) = 0;

    /**
     * @brief Calculate the IBS expansion rate.
     * This is the virtual function in the base class to calculate the IBS expansion rate of the ion beam.
     * Should be overloaded in the derived classes.
     * \param[in] lattice The lattice of the ion ring.
     * \param[in] beam The ion beam.
     * \return The horizontal, vertical, and longitudinal IBS rates.
     */
    virtual std::tuple<double, double, double> rate(const Lattice &lattice, const Beam &beam) = 0;
};

/**
 * @brief IBS rate calculation using the Martini model.
 * Ref: M. Martini, Intrabeam scattering in the ACOL-AA machines, CERN-PS-84-9 (AA), Geneva, Switzerland, May 1984.
 */
class IBSSolver_Martini : public IBSSolver {
private:
    struct TrigonometryStorageUV {
        double sin_u2_cos_v2;
        double g1;
        double g2_1;
        double g2_2;
    };
    struct TrigonometryStorageV {
        double sin_v;
        double cos_v;
    };
    struct TrigonometryStorageU {
        double sin_u;
        double sin_u2;
        double cos_u2;
        double g3;
        std::vector<TrigonometryStorageUV> uv;
    };
    struct OpticalStorage {
        double a;
        double b2;
        double c2;
        double d2;
        double dtld;
        double k1;
        double k2;
        double k3;
    };

    int nu_ = 0;                //Grid number in u direction.
    int nv_ = 0;                //Grid number in v direction.
    int nz_ = 0;                //Grid number in z direction.

    // Scratch variables for IBS calculation (Martini model)
    std::vector<double> sigma_xbet, sigma_xbetp, sigma_y, sigma_yp;
    std::vector<TrigonometryStorageU> storage_u;
    std::vector<TrigonometryStorageV> storage_v;
    std::vector<OpticalStorage> storage_opt;
    std::vector<double> f1, f2, f3;

    void bunch_size(const Lattice &lattice, const Beam &beam);
    void abcdk(const Lattice &lattice, const Beam &beam);
    void coef_f();
    void f();
    double coef_a(const Lattice &lattice, const Beam &beam) const;
public:
    int nu() const { return nu_; } ///< Grid number of 3D numerical integration.
    int nv() const { return nv_; } ///< Grid number of 3D numerical integration.
    int nz() const { return nz_; } ///< Grid number of 3D numerical integration.

    /**
     * @brief Set grid number of 3D numerical integration.
     * \param[in] nu Grid number of 3D numerical integration.
     */
    void set_nu(int nu) { assert(nu>0&&"Wrong value of nu in IBS parameters!"); nu_ = nu; invalidate_cache(); }
    /**
     * @brief Set grid number of 3D numerical integration.
     * \param[in] nv Grid number of 3D numerical integration.
     */
    void set_nv(int nv) { assert(nv>0&&"Wrong value of nv in IBS parameters!"); nv_ = nv; invalidate_cache(); }
    /**
     * @brief Set grid number of 3D numerical integration.
     * \param[in] nz Grid number of 3D numerical integration.
     */
    void set_nz(int nz) { assert(nz>0&&"Wrong value of nz in IBS parameters!"); nz_ = nz; invalidate_cache(); }
    /**
     * @brief Constructor of IBS Solver using Martini model.
     * \param[in] nu Grid number of 3D numerical integration.
     * \param[in] nv Grid number of 3D numerical integration.
     * \param[in] nz Grid number of 3D numerical integration.
     * \param[in] log_c Coulomb logarithm
     * \param[in] k The transverse coupling rate, from 0 to 1.
     */
    IBSSolver_Martini(int nu, int nv, int nz, double log_c, double k);
    /**
     * @brief Calculate the IBS expansion rate  using Martini model.
     * This is the virtual function in the base class to calculate the IBS expansion rate of the ion beam.
     * Should be overloaded in the derived classes.
     * \param[in] lattice The lattice of the ion ring.
     * \param[in] beam The ion beam.
     * \param[out] rx The horizontal IBS rate.
     * \param[out] ry The vertical IBS rate.
     * \param[out] rs The longitudinal IBS rate.
     */
    virtual void rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs);
    /**
     * @brief Calculate the IBS expansion rate using Martini model.
     * \param[in] lattice The lattice of the ion ring.
     * \param[in] beam The ion beam.
     * \return The horizontal, vertical, and longitudinal IBS rates.
     */
    virtual std::tuple<double, double, double> rate(const Lattice &lattice, const Beam &beam) {
        double rx, ry, rs;
        rate(lattice, beam, rx, ry, rs);
        return std::make_tuple(rx, ry, rs);
    }
};

/**
 * @brief IBS rate calculation using the Bjorken-Mtingwa model and Nagaitsev's fast numerical evaluation.
 * Ref: J.D.Bjorken,S.K.Mtingwa,��IntrabeamScattering,��Part.Acc.Vol.13,pp.115�C143(1983).
 * Ref: S. Nagaitsev, Intrabeam scattering formulas for fast numerical evaluation, Phys. Rev. ST Accel. Beams 8, 064403,
 * 30 June 2005.
 */
class IBSSolver_BM : public IBSSolver {
 private:
     struct OpticalStorage { //variables only depends on the TWISS parameters and the energy.
         double phi;
         double dx2; //D_x * D_x
         double dx_betax_phi_2; // D_x * D_x / (beta_x * beta_x) + phi * phi
         double sqrt_betay; // sqrt(beta_y)
         double gamma_phi_2; // gamma * gamma * phi * phi
     };
     struct Kernels {
         double  psi;
         double sx;
         double sp;
         double sxp;
         double inv_sigma;
     };

     // Scratch variables for IBS calculation (Bjorken-Mtingwa model using Sergei Nagitsev's formula)
     std::vector<OpticalStorage> optical_strage;
     std::vector<Kernels> kernels;
     void init_fixed_var(const Lattice &lattice, const Beam &beam);
     void calc_kernels(const Lattice &lattice, const Beam &beam);
     double coef_bm(const Lattice &lattice, const Beam &beam) const;
 public:
     /**
     * @brief Constructor of IBS Solver using Bjorken-Mtingwa model and Nagaitsev's fast numerical evaluation.
     * \param[in] log_c Coulomb logarithm
     * \param[in] k The transverse coupling rate, from 0 to 1.
     */
     IBSSolver_BM(double log_c, double k);
     /**
     * @brief Calculate the IBS expansion rate  using Bjorken-Mtingwa model and Nagaitsev's fast numerical evaluation.
     * \param[in] lattice The lattice of the ion ring.
     * \param[in] beam The ion beam.
     * \param[out] rx The horizontal IBS rate.
     * \param[out] ry The vertical IBS rate.
     * \param[out] rs The longitudinal IBS rate.
     */
     virtual void rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs);
     /**
     * @brief Calculate the IBS expansion rate using Bjorken-Mtingwa model and Nagaitsev's fast numerical evaluation.
     * \param[in] lattice The lattice of the ion ring.
     * \param[in] beam The ion beam.
     * \return The horizontal, vertical, and longitudinal IBS rates.
     */
     virtual std::tuple<double, double, double> rate(const Lattice &lattice, const Beam &beam) {
        double rx, ry, rs;
        rate(lattice, beam, rx, ry, rs);
        return std::make_tuple(rx, ry, rs);
     }
};

/**
 * @brief IBS rate calculation using the Bjorken-Mtingwa model with vertical dispersion and non-ultrarelativistic terms.
 * Ref: F. Antoniou & F. Zimmermann, Revision of intrabeam scattering with non-ultrarelativistic corrections and
 * vertical dispersion for MAD-X, CERN-ATS-2012-066, Geneva, Switzerland, May 8, 2012.
 */
class IBSSolver_BMZ : public IBSSolver {
private:
    int nt_;     //Number of steps for integration.
    struct optcl{
        double phi_x2;
        double phi_y2;
        double hx;
        double hy;
        double dx_2_over_beta_x;
        double dy_2_over_beta_y;
        double beta_x_over_hx;
        double hy_beta_x_over_hx;
        double beta_phi_x2;
        double beta_phi_y2;
        double hy_over_beta_y;
    };
    std::vector<optcl> optical;
    double factor = 3;
    void init_optical(const Lattice &lattice);
    double calc_abc(const Lattice &lattice, const Beam& beam, int i, double& a, double& b, double& c,
                               double& ax, double& bx, double& ay, double& by,double& al, double& bl);
    double coef(const Lattice &lattice, const Beam &beam) const;
    void calc_integral(double a, double b, double c, double ax, double bx, double ay, double by, double al,
                                  double bl, double& ix, double& iy, double& is, int nt, double u);
public:
    /**
     * @brief Constructor of IBS Solver using Bjorken-Mtingwa model with vertical dispersion and non-ultrarelativistic terms.
     * \param[in] nt Number of steps to carry out the numerical integration.
     * \param[in] log_c Coulomb logarithm
     * \param[in] k The transverse coupling rate, from 0 to 1.
     */
    IBSSolver_BMZ(int nt, double log_c, double k);

    /**
    * @brief Set the number of  steps to carry out the numerical integration.
    * \param[in] n The number of steps
    */
    void set_nt(int n){assert(n>0&&"Wrong value of nt in IBS parameters!"); nt_ = n; invalidate_cache();}
    /**
     * @brief Calculate the IBS expansion rate  using Bjorken-Mtingwa model with vertical dispersion and non-ultrarelativistic terms.
     * \param[in] lattice The lattice of the ion ring.
     * \param[in] beam The ion beam.
     * \param[out] rx The horizontal IBS rate.
     * \param[out] ry The vertical IBS rate.
     * \param[out] rs The longitudinal IBS rate.
     */
    virtual void rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs);
    /**
    * @brief Set the upper bound to carry out the numerical integration.
    * \param[in] x This number is in proportion to the upper bound of the numerical integration.
    */
    void set_factor(double x){factor = x;}
    /**
     * @brief Calculate the IBS expansion rate using Bjorken-Mtingwa model with vertical dispersion and non-ultrarelativistic terms.
     * \param[in] lattice The lattice of the ion ring.
     * \param[in] beam The ion beam.
     * \return The horizontal, vertical, and longitudinal IBS rates.
     */
    virtual std::tuple<double, double, double> rate(const Lattice &lattice, const Beam &beam) {
        double rx, ry, rs;
        rate(lattice, beam, rx, ry, rs);
        return std::make_tuple(rx, ry, rs);
    }
};

/**
 * @brief IBS rate calculation using the Bjorken-Mtingwa model with vertical dispersion as in BETACOOL.
 * Ref: I. Meshkov, A. Sidorin, A. Smirnov, G. Trubnikov, P. Pivin, BETACOOL Physics Guide for simulation of long term
 * beam dynamics in ion storage rings (since 1995), Joint Institute for Nuclear Research, Joliot Curie, 6, Dubna,
 * 141980 Russian Federation, 2007.
 */
class IBSSolver_BM_Complete : public IBSSolver {
private:
    int nt_;
    double inv_ex;
    double inv_ey;
    double inv_dp2;
    double gamma;
    double gamma2;
    double factor = 3;
    struct Itgrl{
        double lambda;
        double lambda_sqrt;
        double ct;   // ct = 1/(1-t)^2
    };
    struct Optc {
        double phix;
        double phiy;
        double hx;
        double hy;
    };

    std::vector<Optc> optc;
    std::vector<Itgrl> itgrl;
    void init_optc(const Lattice &lattice);
    double det(std::array<std::array<double, 3>,3>& l);
    double inv(std::array<std::array<double, 3>,3>& l, std::array<std::array<double, 3>,3>& v);
    double trace(std::array<std::array<double, 3>,3>& l){return l[0][0]+l[1][1]+l[2][2];}
    void calc_l(const Lattice& lattice, int i, std::array<std::array<double, 3>,3>& lh,
                std::array<std::array<double, 3>,3>& lv, std::array<std::array<double, 3>,3>& ls);
    void calc_itgl(int i, std::array<std::array<double, 3>,3>& ii, std::array<std::array<double, 3>,3>& l,
                  std::array<std::array<double, 3>,3>& ll, std::array<std::array<double, 3>,3>& lh,
                  std::array<std::array<double, 3>,3>& lv, std::array<std::array<double, 3>,3>& ls);
    void calc_beam_const(const Beam& beam);
    double coef(const Lattice &lattice, const Beam &beam) const;
public:
    /**
     * @brief Constructor of IBS Solver using Bjorken-Mtingwa model with vertical dispersion vertical dispersion as in BETACOOL.
     * \param[in] nt Number of steps to carry out the numerical integration.
     * \param[in] log_c Coulomb logarithm
     * \param[in] k The transverse coupling rate, from 0 to 1.
     */
     IBSSolver_BM_Complete(int nt, double log_c, double k);
     /**
    * @brief Set the upper bound to carry out the numerical integration.
    * \param[in] x This number is in proportion to the upper bound of the numerical integration.
    */
     void set_factor(double x){factor = x;}
     /**
     * @brief Calculate the IBS expansion rate  using Bjorken-Mtingwa model with vertical dispersion as in BETACOOL.
     * \param[in] lattice The lattice of the ion ring.
     * \param[in] beam The ion beam.
     * \param[out] rx The horizontal IBS rate.
     * \param[out] ry The vertical IBS rate.
     * \param[out] rs The longitudinal IBS rate.
     */
     virtual void rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs);
      /**
     * @brief Calculate the IBS expansion rate using Bjorken-Mtingwa model with vertical dispersion as in BETACOOL.
     * \param[in] lattice The lattice of the ion ring.
     * \param[in] beam The ion beam.
     * \return The horizontal, vertical, and longitudinal IBS rates.
     */
     virtual std::tuple<double, double, double> rate(const Lattice &lattice, const Beam &beam) {
        double rx, ry, rs;
        rate(lattice, beam, rx, ry, rs);
        return std::make_tuple(rx, ry, rs);
     }

};
#endif
