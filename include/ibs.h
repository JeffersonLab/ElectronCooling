#ifndef IBS_HPP
#define IBS_HPP

#include <assert.h>
#include <memory>
#include <vector>

class Lattice;
class Beam;

enum class IBSModel {MARTINI, BM, BMZ};

class IBSSolver {
protected:
    double log_c_ = 0.0;     //Coulomb logarithm.
    double k_ = 0.0;          //Coupling rate in transverse directions.
    bool cache_invalid = true;
    bool ibs_by_element = false; //Calculate and output the ibs rate contribution element by element.

    void ibs_coupling(double &rx, double &ry, double k, double emit_x, double emit_y);
public:
    double log_c() const { return log_c_; }
    double k() const { return k_; }
    void set_k(double x) { k_ = x; }
    void set_log_c(double x) { log_c_ = x; }
    void set_ibs_by_element(bool b) {ibs_by_element = b;}
    void invalidate_cache() { cache_invalid = true; }

    IBSSolver(double log_c, double k);

    virtual void rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs) = 0;
};

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
    int nu() const { return nu_; }
    int nv() const { return nv_; }
    int nz() const { return nz_; }
    void set_nu(int nu) { assert(nu>0&&"Wrong value of nu in IBS parameters!"); nu_ = nu; invalidate_cache(); }
    void set_nv(int nv) { assert(nv>0&&"Wrong value of nv in IBS parameters!"); nv_ = nv; invalidate_cache(); }
    void set_nz(int nz) { assert(nz>0&&"Wrong value of nz in IBS parameters!"); nz_ = nz; invalidate_cache(); }
    IBSSolver_Martini(int nu, int nv, int nz, double log_c, double k);
    virtual void rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs);
};

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
     IBSSolver_BM(double log_c, double k);
     virtual void rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs);

};

class IBSSolver_BMZ : public IBSSolver {
private:
    int nt_;     //Number of steps for integration.
    struct itgrl{
        double lambda;
//        double lambda_2;
//        double lambda_3;
        double lambda_sqrt;
        double ct;   // ct = 1/(1-t)^2
    };
    double gamma_2;
    double gamma_2_inv;
    double gamma_4;
    double emit_x_inv;
    double emit_y_inv;
    double emit_x2_inv;
    double emit_y2_inv;

    struct optcl{
        double phi_x2;
        double phi_y2;
        double hx;
        double hy;
        double dx_2_over_beta_x;
        double dy_2_over_beta_y;
        double beta_xy;
        double beta_x_over_hx;
        double hy_beta_x_over_hx;
        double beta_phi_x2;
        double beta_phi_y2;
        double hy_over_beta_y;
        double hx_hy_over_beta_y;
        double beta_x_hy_over_beta_y;
    };
    std::vector<itgrl> integral;
    std::vector<optcl> optical;
    void init_integral(int n);
    void init_optical(const Lattice &lattice);
    double calc_abc(const Lattice &lattice, const Beam& beam, int i, double& a, double& b, double& c,
                               double& ax, double& bx, double& ay, double& by,double& al, double& bl);
    double coef(const Lattice &lattice, const Beam &beam) const;
    void calc_integral(double a, double b, double c, double ax, double bx, double ay, double by, double al,
                                  double bl, double ix, double iy, double is, int nt, std::vector<itgrl>& g);
public:
    IBSSolver_BMZ(int nt, double log_c, double k);
    set_nt(int n){assert(n>0&&"Wrong value of nt in IBS parameters!"); nt_ = n; invalidate_cache();}
    virtual void rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs);
};
#endif
