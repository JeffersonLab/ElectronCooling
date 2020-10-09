#include <fstream>
#include <memory>
#include <cmath>
#include <cstring>
#include "functions.h"
#include "ibs.h"
#include "ring.h"
#include "beam.h"

IBSSolver::IBSSolver(double log_c, double k)
    : log_c_(log_c), k_(k)
{
}

void IBSSolver::ibs_coupling(double &rx, double &ry, double k, double emit_x, double emit_y)
{
    double rxc = 0.5*(rx*(2-k)+ry*k*emit_y/emit_x);
    double ryc = 0.5*(ry*(2-k)+rx*k*emit_x/emit_y);
    rx = rxc;
    ry = ryc;
}

void IBSSolver::ibs_by_element_sddshead(std::ofstream& outfile, int n_element) {
    using std::endl;
    outfile<<"SDDS1"<<endl;
    outfile<<"! Define colums:"<<endl
        <<"&column name=s, type=double, units=m, description=\"element position\", &end"<<endl
        <<"&column name=beta_x, type=double, units=m, description=\"TWISS parameter beta x\", &end"<<endl
        <<"&column name=beta_y, type=double, units=m, description=\"TWISS parameter beta y\", &end"<<endl
        <<"&column name=alfa_x, type=double, units=NULL, description=\"TWISS parameter alfa x\", &end"<<endl
        <<"&column name=alfa_y, type=double, units=NULL, description=\"TWISS parameter alfa y\", &end"<<endl
        <<"&column name=dx, type=double, units=m, description=\"horizontal dispersion function\", &end"<<endl
        <<"&column name=dy, type=double, units=m, description=\"vertical dispersion function\", &end"<<endl
        <<"&column name=rx_i, type=double, units=1/s, description=\"horizontal IBS expansion rate contribution by the element\", &end"<<endl
        <<"&column name=ry_i, type=double, units=1/s, description=\"vertical IBS expansion rate contribution by the element\", &end"<<endl
        <<"&column name=rs_i, type=double, units=1/s, description=\"longitudinal IBS expansion rate contribution by the element\", &end"<<endl
        <<"&column name=rx, type=double, units=1/s, description=\"accumulated horizontal IBS expansion rate\", &end"<<endl
        <<"&column name=ry, type=double, units=1/s, description=\"accumulated vertical IBS expansion rate\", &end"<<endl
        <<"&column name=rs, type=double, units=1/s, description=\"accumulated longitudinal IBS expansion rate\", &end"<<endl
        <<"!Declare ASCII data and end the header"<<endl
        <<"&data mode=ascii, &end"<<endl
        <<n_element-1<<endl;
}


IBSSolver_Martini::IBSSolver_Martini(int nu, int nv, int nz, double log_c, double k)
    : IBSSolver(log_c, k), nu_(nu), nv_(nv), nz_(nz)
{
}

//Calculate sigma_xbet, sigma_xbetp, sigma_y, sigma_yp
void IBSSolver_Martini::bunch_size(const Lattice &lattice, const Beam &beam)
{
    int n = lattice.n_element();

    sigma_xbet.resize(n);
    sigma_xbetp.resize(n);
    sigma_y.resize(n);
    sigma_yp.resize(n);

    double emit_x = beam.emit_x();
    double emit_y = beam.emit_y();
    for(int i=0; i<n; ++i) {
        sigma_xbet[i] = sqrt(lattice.betx(i)*emit_x);
        sigma_y[i] = sqrt(lattice.bety(i)*emit_y);
        double alf2 = lattice.alfx(i);
        alf2 *= alf2;
        sigma_xbetp[i] = sqrt((1+alf2)*emit_x/lattice.betx(i));
        alf2 = lattice.alfy(i);
        alf2 *= alf2;
        sigma_yp[i] = sqrt((1+alf2)*emit_y/lattice.bety(i));
    }
}

//Calculate a, b, c, d and dtld
//Call bunch_size() before calling this one
void IBSSolver_Martini::abcdk(const Lattice &lattice, const Beam &beam)
{
    double d_tld, q, sigma_x, sigma_tmp;
    const int n = lattice.n_element();
    storage_opt.resize(n);

    const double dp_p = beam.dp_p();
    const double beta = beam.beta();
    const double gamma = beam.gamma();
    const double r = beam.r();
    for(int i=0; i<n; ++i){
        const double betx = lattice.betx(i);
        const double alfx = lattice.alfx(i);
        const double dx = lattice.dx(i);
        const double dpx = lattice.dpx(i);
        const double alfy = lattice.alfy(i);

        d_tld = alfx*dx+betx*dpx;
        sigma_x = sqrt(sigma_xbet[i]*sigma_xbet[i]+dx*dx*dp_p*dp_p);
        sigma_tmp = dp_p*sigma_xbet[i]/(gamma*sigma_x);
        q = 2*beta*gamma*sqrt(sigma_y[i]/r);

        OpticalStorage os;
        os.a = sigma_tmp*sqrt(1+alfx*alfx)/sigma_xbetp[i];
        os.b2 = sigma_tmp*sqrt(1+alfy*alfy)/sigma_yp[i];
        os.b2 *= os.b2;
        os.c2 = q*sigma_tmp;
        os.c2 *= os.c2;
        os.d2 = dp_p*dx/sigma_x;
        os.d2 *= os.d2;
        os.dtld = dp_p*d_tld/sigma_x;

        os.k1 = 1.0 / os.c2;
        os.k2 = os.a * os.a * os.k1;
        os.k3 = os.b2 * os.k1;

        storage_opt[i] = os;
    }
}

void IBSSolver_Martini::coef_f()
{
    storage_u.resize(nu_);
    storage_v.resize(nv_);

    double dv = 2*k_pi/nv_;
    double v = -0.5*dv;
    for(int i=0; i<nv_; ++i){
        v += dv;
        storage_v[i] = TrigonometryStorageV({sin(v), cos(v)});
    }

    double du = k_pi/nu_;
    double u = -0.5*du;
    for(int i=0; i<nu_; ++i){
        u += du;
        double sin_u = sin(u);
        double sin_u2 = sin_u * sin_u;
        double cos_u2 = 1 - sin_u2;
        double g3 = 1 - 3 * cos_u2;
	std::vector<TrigonometryStorageUV> uv;
	uv.resize(nv_);
        for(int j=0; j<nv_; ++j){
            const double sin_u2_cos_v2 = sin_u2 * storage_v[j].cos_v * storage_v[j].cos_v;
            const double g1 = 1 - 3 * sin_u2_cos_v2;
            const double g2_1 = 1 - 3 * sin_u2 * storage_v[j].sin_v * storage_v[j].sin_v;
            const double g2_2 = 6 * sin_u * storage_v[j].sin_v * storage_v[j].cos_v;
	    TrigonometryStorageUV tempUV = {sin_u2_cos_v2, g1, g2_1, g2_2};
	    uv[j] = tempUV;
        }
        storage_u[i] = TrigonometryStorageU({sin_u, sin_u2, cos_u2, g3, uv});
    }
}

void IBSSolver_Martini::f()
{
    const int n_element = storage_opt.size();
    f1.resize(n_element);
    f2.resize(n_element);
    f3.resize(n_element);

    if (log_c_ > 0) {
        const double duvTimes2Logc = 2*k_pi*k_pi/(nu_*nv_) * 2 * log_c_;
#ifdef _OPENMP
        // Maybe make this runtime-adjustable. SMT gives no benefit here, so choose n = physical cores
        #pragma omp parallel for num_threads(6)
#endif
        for(int ie=0; ie < n_element; ie++) {
            const OpticalStorage &os = storage_opt[ie];
            double tempf1 = 0, tempf2 = 0, tempf3 = 0;
            for(int iu=0; iu<nu_; ++iu) {
                const TrigonometryStorageU &tu = storage_u[iu];
                double sum1 = 0, sum2 = 0, sum3 = 0;
                for(int iv=0; iv<nv_; ++iv) {
                    const TrigonometryStorageV &tv = storage_v[iv];
                    const TrigonometryStorageUV &tuv = tu.uv[iv];

                    double tmp = os.a * tv.sin_v - os.dtld * tv.cos_v;
                    tmp *= tmp;
                    const double inv_int_z = (tuv.sin_u2_cos_v2 + tu.sin_u2 * tmp + os.b2 * tu.cos_u2) * os.k1;
                    sum1 += tuv.g1 / inv_int_z;
                    sum2 += (tuv.g2_1 + tuv.g2_2 * os.dtld / os.a) / inv_int_z;
                    sum3 += tu.g3 / inv_int_z;
                }
                tempf1 += tu.sin_u * sum1;
                tempf2 += tu.sin_u * sum2;
                tempf3 += tu.sin_u * sum3;
            }
            f1[ie] = tempf1 * os.k1 * duvTimes2Logc;
            f2[ie] = tempf2 * os.k2 * duvTimes2Logc;
            f3[ie] = tempf3 * os.k3 * duvTimes2Logc;
        }
    } else {
        const double duv = 2*k_pi*k_pi/(nv_*nu_);
#ifdef _OPENMP
        // Maybe make this runtime-adjustable. SMT gives no benefit here, so choose n = physical cores
        #pragma omp parallel for num_threads(6)
#endif
        for(int ie=0; ie<n_element; ++ie){
            const OpticalStorage &os = storage_opt[ie];
            double tempf1 = 0, tempf2 = 0, tempf3 = 0;
            for(int iu=0; iu<nu_; ++iu){
                const TrigonometryStorageU &tu = storage_u[iu];
                double sum1 = 0, sum2 = 0, sum3 = 0;
                for(int iv=0; iv<nv_; ++iv){
                    const TrigonometryStorageV &tv = storage_v[iv];
                    const TrigonometryStorageUV &tuv = tu.uv[iv];

                    double tmp = os.a * tv.sin_v - os.dtld * tv.cos_v;
                    tmp *= tmp;
                    const double d_uv = (tuv.sin_u2_cos_v2 + tu.sin_u2 * tmp + os.b2 * tu.cos_u2) * os.k1;
                    double int_z = 0;
                    const double dz = 20/(d_uv*nz_);
                    double z = -0.5*dz;
                    for(int iz=0; iz<nz_; ++iz){
                        z += dz;
                        int_z += exp(-d_uv*z)*log(1+z*z)*dz;
                    }
                    sum1 += int_z * tuv.g1;
                    sum2 += int_z * (tuv.g2_1 + tuv.g2_2 * os.dtld / os.a);
                    sum3 += int_z * tu.g3;
                }
                tempf1 += tu.sin_u * sum1;
                tempf2 += tu.sin_u * sum2;
                tempf3 += tu.sin_u * sum3;
            }
            f1[ie] = tempf1 * os.k1 * duv;
            f2[ie] = tempf2 * os.k2 * duv;
            f3[ie] = tempf3 * os.k3 * duv;
        }
    }
}

double IBSSolver_Martini::coef_a(const Lattice &lattice, const Beam &beam) const
{
    double lambda = beam.particle_number()/lattice.circ();
    if(beam.bunched()) lambda = beam.particle_number()/(2*sqrt(k_pi)*beam.sigma_s());

    double beta3 = beam.beta();
    beta3 = beta3*beta3*beta3;
    double gamma4 = beam.gamma();
    gamma4 *= gamma4;
    gamma4 *= gamma4;
    return k_c*beam.r()*beam.r()*lambda/(16*k_pi*sqrt(k_pi)*beam.dp_p()*beta3*gamma4)/(beam.emit_x()*beam.emit_y());
}

void IBSSolver_Martini::rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs)
{
    bunch_size(lattice, beam);
    abcdk(lattice, beam);
    if (cache_invalid) {
        coef_f();
        cache_invalid = false;
    }
    f();

    double a = coef_a(lattice, beam);

    rx = 0;
    ry = 0;
    rs = 0;
    int n=2;
    if (beam.bunched()) n=1;
    const double circ = lattice.circ();

    const int n_element = lattice.n_element();
    if(ibs_by_element) {
        std::ofstream out;
        out.open("ibs_by_element.txt");
        ibs_by_element_sddshead(out, n_element);
        out.precision(10);
        out<<std::showpos;
        out<<std::scientific;
        double s = 0;
        double inv_circ = 1/circ;
        for(int i=0; i<n_element-1; ++i){
            double l_element = lattice.l_element(i);
            const OpticalStorage &os = storage_opt[i];
            double rsi = n*a*(1-os.d2)*f1[i]*l_element*inv_circ;
            double rxi = a*(f2[i]+f1[i]*(os.d2+os.dtld*os.dtld))*l_element*inv_circ;
            double ryi = a*f3[i]*l_element*inv_circ;
            rx += rxi;
            ry += ryi;
            rs += rsi;
            out<<s<<' '<<lattice.betx(i)<<' '<<lattice.bety(i)<<' '<<lattice.alfx(i)<<' '
                <<lattice.alfy(i)<<' '<<lattice.dx(i)<<' '<<lattice.dy(i)<<' '
                <<rxi<<' '<<ryi<<' '<<rsi<<' '<<rx<<' '<<ry<<' '<<rs<<std::endl;
            s += l_element;
        }
        out.close();

    }
    else {
        for(int i=0; i<n_element-1; ++i){
            const double l_element = lattice.l_element(i);
            const OpticalStorage &os = storage_opt[i];
            rs += n*a*(1-os.d2)*f1[i]*l_element;
            rx += a*(f2[i]+f1[i]*(os.d2+os.dtld*os.dtld))*l_element;
            ry += a*f3[i]*l_element;
        }
        rx /= circ;
        ry /= circ;
        rs /= circ;
    }

    if(k_>0) ibs_coupling(rx, ry, k_, beam.emit_nx(), beam.emit_ny());
}



IBSSolver_BM::IBSSolver_BM(double log_c, double k)
    : IBSSolver(log_c, k)
{
}

void IBSSolver_BM::init_fixed_var(const Lattice& lattice, const Beam& beam) {
    int n = lattice.n_element();
    optical_strage.resize(n);

    for(int i=0; i<n; ++i) {
        OpticalStorage os;
        double dx_betax = lattice.dx(i)/lattice.betx(i);
        os.dx2 = lattice.dx(i) * lattice.dx(i);
        os.phi = lattice.dpx(i) + lattice.alfx(i)*dx_betax;
        os.dx_betax_phi_2 = dx_betax*dx_betax + os.phi*os.phi;
        os.sqrt_betay = sqrt(lattice.bety(i));
        os.gamma_phi_2 = beam.gamma() * beam.gamma() * os.phi * os.phi;
        optical_strage.at(i)= os;
    }
}

void IBSSolver_BM::calc_kernels(const Lattice& lattice, const Beam& beam) {
    auto emit_x = beam.emit_x();
    auto emit_y = beam.emit_y();
    auto sigma_p2 = beam.dp_p() * beam.dp_p();
    auto inv_sigma_p2 = 1/sigma_p2;
    auto n = lattice.n_element();
    auto gamma2 = beam.gamma()*beam.gamma();

    kernels.resize(n);

    for(int i=0; i<n; ++i) {
        Kernels knl;
        auto betx = lattice.betx(i);
        auto bety = lattice.bety(i);
        auto sigma_x = sqrt(optical_strage.at(i).dx2*sigma_p2+emit_x*betx);
        auto sigma_y = sqrt(emit_y*bety);
        knl.inv_sigma = 1/(sigma_x*sigma_y);

        auto ax = betx/emit_x;
        auto lamda_1 = bety/emit_y; //lamda_1 = ay.
        auto as = ax*optical_strage.at(i).dx_betax_phi_2 + inv_sigma_p2;
        auto a1 = gamma2*as;
        auto a2 = (ax-a1)/2;
        a1 = (ax+a1)/2;

        auto lamda_2 = sqrt(a2*a2+ax*ax*optical_strage.at(i).gamma_phi_2);
        auto lamda_3 = a1 - lamda_2;
        auto tmp1 = 3/lamda_2;
        lamda_2 = a1 + lamda_2;

        auto inv_lamda_1 = 1/lamda_1;
        auto inv_lamda_2 = 1/lamda_2;
        auto inv_lamda_3 = 1/lamda_3;

        auto r1 = rd(inv_lamda_2, inv_lamda_3, inv_lamda_1);
        auto r2 = rd(inv_lamda_3, inv_lamda_1, inv_lamda_2);
        auto r3 = 3*sqrt(lamda_1*lamda_2*lamda_3)-r1-r2;

        r1 *= inv_lamda_1*2;
        r2 *= inv_lamda_2;
        r3 *= inv_lamda_3;

        knl.psi = -r1 + r2 + r3;

        knl.sxp = tmp1*ax*optical_strage.at(i).gamma_phi_2*(r3-r2);
        tmp1 = tmp1*a2;
        auto tmp2 = 1 + tmp1;
        tmp1 = 1 - tmp1;
        knl.sp = gamma2*(r1-r2*tmp1-r3*tmp2)/2;
        knl.sx = (r1-r2*tmp2-r3*tmp1)/2;
        kernels.at(i) = knl;
    }
}

double IBSSolver_BM::coef_bm(const Lattice &lattice, const Beam &beam) const {
    double lambda = 1;
    if (beam.bunched())
        lambda /= 2*sqrt(k_pi)*beam.sigma_s();
    else
        lambda /= lattice.circ();

    double beta3 = beam.beta()*beam.beta()*beam.beta();
    double gamma5 = beam.gamma()*beam.gamma()*beam.gamma()*beam.gamma()*beam.gamma();
    return lambda*beam.particle_number()*beam.r()*beam.r()*k_c/(lattice.circ()*6*sqrt(k_pi)*beta3*gamma5);
}

void IBSSolver_BM::rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs)
{
    int n_element = lattice.n_element();

    if (cache_invalid) {
        init_fixed_var(lattice, beam);
        cache_invalid = false;
    }

    calc_kernels(lattice, beam);
    double c_bm = coef_bm(lattice, beam);
    const double lc = log_c();
    c_bm *= lc;

    rx = 0;
    ry = 0;
    rs = 0;
    int n=2;
    if (beam.bunched()) n=1;

    if(ibs_by_element) {

        std::ofstream out;
        out.open("ibs_by_element.txt");
        ibs_by_element_sddshead(out, n_element);
        out.precision(10);
        out<<std::showpos;
        out<<std::scientific;
        double s = 0;
        for(int i=0; i<n_element-1; ++i){
            double l_element = lattice.l_element(i);
            double rxi = (lattice.betx(i)*kernels.at(i).inv_sigma*(kernels.at(i).sx
                            +optical_strage.at(i).dx_betax_phi_2*kernels.at(i).sp
                            +kernels.at(i).sxp)*l_element)*c_bm/beam.emit_x();
            double ryi = (lattice.bety(i)*kernels.at(i).inv_sigma*kernels.at(i).psi*l_element)*c_bm/beam.emit_y();
            double rsi = (kernels.at(i).inv_sigma*kernels.at(i).sp*l_element)*n*c_bm/(beam.dp_p()*beam.dp_p());

            rx += rxi;
            ry += ryi;
            rs += rsi;
            out<<s<<' '<<lattice.betx(i)<<' '<<lattice.bety(i)<<' '<<lattice.alfx(i)<<' '
                <<lattice.alfy(i)<<' '<<lattice.dx(i)<<' '<<lattice.dy(i)<<' '
                <<rxi<<' '<<ryi<<' '<<rsi<<' '<<rx<<' '<<ry<<' '<<rs<<std::endl;
            s += l_element;
        }
        out.close();

    }
    else {
        for(int i=0; i<n_element-1; ++i){
            double l_element = lattice.l_element(i);
            rs += kernels.at(i).inv_sigma*kernels.at(i).sp*l_element;
            ry += lattice.bety(i)*kernels.at(i).inv_sigma*kernels.at(i).psi*l_element;
            rx += lattice.betx(i)*kernels.at(i).inv_sigma*(kernels.at(i).sx
                    +optical_strage.at(i).dx_betax_phi_2*kernels.at(i).sp
                    +kernels.at(i).sxp)*l_element;
        }

        rs *= n*c_bm;
        rx *= c_bm;
        ry *= c_bm;

        rs /= beam.dp_p()*beam.dp_p();
        rx /= beam.emit_x();
        ry /= beam.emit_y();
    }

    if(k_>0)
        ibs_coupling(rx, ry, k_, beam.emit_nx(), beam.emit_ny());
}

IBSSolver_BMZ::IBSSolver_BMZ(int nt, double log_c, double k) : IBSSolver(log_c, k), nt_(nt) {
    gw = gsl_integration_workspace_alloc(limit);
}

double IBSSolver_BMZ::core(double q, void* params){
     double a = ((P*)params)->a;
     double b = ((P*)params)->b;
     double c = ((P*)params)->c;
     double ai = ((P*)params)->ai;
     double bi = ((P*)params)->bi;

     double d = q*(q*(q+a)+b)+c;
     d *= sqrt(d);

     return sqrt(q)*(ai*q+bi)/d;

}

void IBSSolver_BMZ::init_integral(int n) {
    integral.resize(n);
    double dt = 1.0/n;
    double t = dt/2;
    for(int i=0; i<n; ++i) {
        double it = 1-t;
        itgrl o;
        o.lambda = t/it;
//        o.lambda_2 = o.lambda*o.lambda;
//        o.lambda_3 = o.lambda_2*o.lambda;
        o.lambda_sqrt = sqrt(o.lambda);
        o.ct = 1/(it*it);
        integral.at(i) = o;
        t += dt;
    }
}

void IBSSolver_BMZ::init_optical(const Lattice &lattice) {
    int n = lattice.n_element();
    optical.resize(n);
    for(int i=0; i<n; ++i) {
        optcl o;
        double alfx = lattice.alfx(i);
        double betx = lattice.betx(i);
        double dx = lattice.dx(i);
        double dpx = lattice.dpx(i);
        double alfy = lattice.alfy(i);
        double bety = lattice.bety(i);
        double dy = lattice.dy(i);
        double dpy = lattice.dpy(i);
        o.phi_x2 = dpx+alfx*dx/betx;
        o.phi_x2 *= o.phi_x2;
        o.phi_y2 = dpy+alfy*dy/bety;
        o.phi_y2 *= o.phi_y2;
        o.beta_phi_x2 = betx*betx*o.phi_x2;
        o.beta_phi_y2 = bety*bety*o.phi_y2;
        o.hx = (dx*dx + o.beta_phi_x2)/betx;
        o.hy = (dy*dy + o.beta_phi_y2)/bety;
        o.dx_2_over_beta_x = dx*dx/betx;
        o.dy_2_over_beta_y = dy*dy/bety;
//        o.beta_xy = betx*bety;
        o.beta_x_over_hx = betx/o.hx;
        o.hy_beta_x_over_hx = o.hy*o.beta_x_over_hx;
        o.hy_over_beta_y = o.hy/bety;
//        o.hx_hy_over_beta_y = o.hx*o.hy_over_beta_y;
//        o.beta_x_hy_over_beta_y = betx*o.hy/bety;
        optical.at(i) = o;
    }
}

double IBSSolver_BMZ::calc_abc(const Lattice &lattice, const Beam& beam, int i, double& a, double& b, double& c,
                               double& ax, double& bx, double& ay, double& by,double& al, double& bl) {
    double emit_x = beam.emit_x();
    double emit_y = beam.emit_y();
    double gamma = beam.gamma();
    double gamma_2 = gamma*gamma;
    double gamma_2_inv = 1/gamma_2;
    double emit_x_inv = 1/emit_x;
    double emit_y_inv = 1/emit_y;
    double emit_x2_inv = emit_x_inv*emit_x_inv;
    double emit_y2_inv = emit_y_inv*emit_y_inv;
    double dp_2_inv = 1/(beam.dp_p()*beam.dp_p());

    optcl& o = optical.at(i);
    double beta_x_over_emit_x = lattice.betx(i)*emit_x_inv;
    double beta_y_over_emit_y = lattice.bety(i)*emit_y_inv;

    double v1 = (o.hx*emit_x_inv + o.hy*emit_y_inv + dp_2_inv)*gamma_2;
    double v2 = beta_x_over_emit_x + beta_y_over_emit_y;
//    double v2 = lattice.betx(i)*emit_x_inv + lattice.bety(i)*emit_y_inv;
    a = v1 + v2;
    double v3 = (o.dx_2_over_beta_x*emit_x_inv + o.dy_2_over_beta_y*emit_y_inv + dp_2_inv) * gamma_2;

    double v4 = beta_x_over_emit_x*beta_y_over_emit_y;
//    double v4 = o.beta_xy*emit_x_inv*emit_y_inv;
    c = v3*v4;
    b = v2*v3 + v4 + v4*(o.phi_x2+o.phi_y2)*gamma_2;
    double v5 = o.hy_beta_x_over_hx*emit_y_inv;

    ax = 2*v1 - v5 + o.beta_x_over_hx*gamma_2_inv*(2*beta_x_over_emit_x - beta_y_over_emit_y - gamma_2*dp_2_inv);
    double v6 = (o.beta_phi_x2*emit_x2_inv + o.beta_phi_y2*emit_y2_inv)*gamma_2;

    double v7 = beta_x_over_emit_x - 2*beta_y_over_emit_y;
    double v8 = (2*o.beta_phi_y2*emit_y2_inv - o.beta_phi_x2*emit_x2_inv)*gamma_2;
    double v9 = v1*v2 - v6;
    bx = v9 + (beta_x_over_emit_x - 4*beta_y_over_emit_y)*beta_x_over_emit_x +
        o.beta_x_over_hx*(v7*dp_2_inv + v4*gamma_2_inv + v8*gamma_2_inv) +v7*o.hy_beta_x_over_hx*emit_y_inv;

    al = 2*v1 - v2;
    bl = v9-2*v4;

    double hy_over_emit_y = o.hy*emit_y_inv;
    double hy_over_beta_y = o.hy/lattice.bety(i);

    ay = -v1-(hy_over_emit_y+beta_x_over_emit_x*hy_over_beta_y)*gamma_2 + 2*gamma_2*hy_over_beta_y*v1 - beta_x_over_emit_x
        + 2*beta_y_over_emit_y;
    by = (beta_y_over_emit_y - 2*beta_x_over_emit_x)*v1-2*beta_x_over_emit_x*hy_over_emit_y*gamma_2 + v4
        + (2*o.beta_phi_x2*emit_x2_inv-o.beta_phi_y2*emit_y2_inv)*gamma_2 + v9*hy_over_beta_y*gamma_2;
}
//
//void IBSSolver_BMZ::calc_integral(double a, double b, double c, double ax, double bx, double ay, double by, double al,
//                                  double bl, double& ix, double& iy, double& is, int nt, std::vector<itgrl>& g) {
//    double dt = 1.0/nt;
//    ix = 0;
//    iy = 0;
//    is = 0;
//    for(int i=0; i<nt; ++i) {
//        double l = g.at(i).lambda;
//        double ls = g.at(i).lambda_sqrt;
//        double ct = g.at(i).ct;
//        double tx = ls*(ax*l+bx);
//        double ty = ls*(ay*l+by);
//        double ts = ls*(al*l+bl);
//        double d = l*(l*(l+a)+b)+c;
//        d *= ct*sqrt(d);
//        ix += tx/d;
//        iy += ty/d;
//        is += ts/d;
//    }
//    ix *= dt;
//    iy *= dt;
//    is *= dt;
//}

void IBSSolver_BMZ::calc_integral(double a, double b, double c, double ax, double bx, double ay, double by, double al,
                                  double bl, double& ix, double& iy, double& is, int nt, std::vector<itgrl>& g) {

    p.a = a;
    p.b = b;
    p.c = c;
    p.ai = ax;
    p.bi = bx;

    IBSSolver_BMZ* ptr2 = this;
    auto ptr = [=](double x)->double{return ptr2->core(x,&p);};
    gsl_function_pp<decltype(ptr)> fp(ptr);
    gsl_function *f = static_cast<gsl_function*>(&fp);

    ix = 0;
    double error = 0;
    int status =  gsl_integration_qagiu(f, 0, espabs, esprel, limit, gw, &ix, &error);
    if(status == GSL_EDIVERGE){
        status = gsl_integration_qagiu(f, 0, espabs*1e-2, esprel*1e-2, 100*limit, gw, &ix, &error);
        if(status == GSL_EDIVERGE) std::cout<<"GSL integration qagui error for transverse friction force!"<<std::endl;
    }

    p.ai = ay;
    p.bi = by;
    iy = 0;
    error = 0;
    status =  gsl_integration_qagiu(f, 0, espabs, esprel, limit, gw, &iy, &error);
    if(status == GSL_EDIVERGE){
        status = gsl_integration_qagiu(f, 0, 1e-10, 1e-10, 10*limit, gw, &iy, &error);
        if(status == GSL_EDIVERGE) std::cout<<"GSL integration qagui error for transverse friction force!"<<std::endl;
    }

    p.ai = al;
    p.bi = bl;
    is = 0;
    error = 0;
    status =  gsl_integration_qagiu(f, 0, espabs, esprel, limit, gw, &is, &error);
    if(status == GSL_EDIVERGE){
        status = gsl_integration_qagiu(f, 0, 1e-10, 1e-10, 10*limit, gw, &is, &error);
        if(status == GSL_EDIVERGE) std::cout<<"GSL integration qagui error for transverse friction force!"<<std::endl;
    }
}


double IBSSolver_BMZ::coef(const Lattice &lattice, const Beam &beam) const {
    double lambda = 2*sqrt(k_pi)*beam.particle_number()/lattice.circ();
    if(beam.bunched()) lambda = beam.particle_number()/beam.sigma_s();
     double beta3 = beam.beta();
    beta3 = beta3*beta3*beta3;
    double gamma4 = beam.gamma();
    gamma4 *= gamma4;
    gamma4 *= gamma4;
    return k_c*beam.r()*beam.r()*lambda/(8*k_pi*beam.dp_p()*beta3*gamma4*beam.emit_x()*beam.emit_y());
}

//double IBSSolver_BMZ::coef(const Lattice &lattice, const Beam &beam) const {
//    double lambda = 1;
//    if (beam.bunched())
//        lambda /= 2*sqrt(k_pi)*beam.sigma_s();
//    else
//        lambda /= lattice.circ();
//
//    double beta3 = beam.beta()*beam.beta()*beam.beta();
//    double gamma5 = beam.gamma()*beam.gamma()*beam.gamma()*beam.gamma()*beam.gamma();
//    return lambda*beam.particle_number()*beam.r()*beam.r()*k_c/(lattice.circ()*6*sqrt(k_pi)*beta3*gamma5);
//}
//

void IBSSolver_BMZ::rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs) {
    int n_element = lattice.n_element();

    if (cache_invalid) {
        init_integral(nt_);
        init_optical(lattice);
        cache_invalid = false;
    }

    double c_bmz = coef(lattice, beam);
    const double lc = log_c();
    c_bmz *= lc;

    rx = 0;
    ry = 0;
    rs = 0;
    int n=2;
    if (beam.bunched()) n=1;

    if(ibs_by_element) {
        std::ofstream out;
        out.open("ibs_by_element.txt");
        out<<"s bet_x bet_y alf_x alf_y disp_x disp_y rx_i ry_i rs_i rx_int ry_int rs_int"<<std::endl;
        out.precision(10);
        out<<std::showpos;
        out<<std::scientific;
        double s = 0;
        for(int i=0; i<n_element-1; ++i){
            double l_element = lattice.l_element(i);
            c_bmz /= lattice.circ();
            double gamma_2 = beam.gamma()*beam.gamma();
            double a, b, c, ax, bx, ay, by, as, bs;
            calc_abc(lattice, beam, i, a, b, c, ax, bx, ay, by, as, bs);
            double ix, iy, is;
            calc_integral(a, b, c, ax, bx, ay, by, as, bs, ix, iy, is, nt_, integral);

            double rxi = optical.at(i).hx*ix*l_element*c_bmz*gamma_2/beam.emit_x();
            double ryi =  lattice.bety(i)*iy*l_element*c_bmz/beam.emit_y();
            double rsi = is*l_element*n*c_bmz*gamma_2/(beam.dp_p()*beam.dp_p());

            rx += rxi;
            ry += ryi;
            rs += rsi;
            out<<s<<' '<<lattice.betx(i)<<' '<<lattice.bety(i)<<' '<<lattice.alfx(i)<<' '
                <<lattice.alfy(i)<<' '<<lattice.dx(i)<<' '<<lattice.dy(i)<<' '
                <<rxi<<' '<<ryi<<' '<<rsi<<' '<<rx<<' '<<ry<<' '<<rs<<std::endl;
            s += l_element;
        }
        out.close();
    }
    else {
        for(int i=0; i<n_element-1; ++i){
            double l_element = lattice.l_element(i);
            double a, b, c, ax, bx, ay, by, as, bs;
            calc_abc(lattice, beam, i, a, b, c, ax, bx, ay, by, as, bs);
            double ix, iy, is;
            calc_integral(a, b, c, ax, bx, ay, by, as, bs, ix, iy, is, nt_, integral);
            rx += optical.at(i).hx*ix*l_element;
            rs += is*l_element;
            ry += lattice.bety(i)*iy*l_element;
        }

        double gamma_2 = beam.gamma()*beam.gamma();
        c_bmz /= lattice.circ();
        rs *= n*c_bmz*gamma_2;
        rx *= c_bmz*gamma_2;
        ry *= c_bmz;

        rs /= beam.dp_p()*beam.dp_p();
        rx /= beam.emit_x();
        ry /= beam.emit_y();

//        std::ofstream out;
//        out.open("ibs_bmz_debug.txt");
//        out<<"a, b, c, ax, bx, ay, by, al, bl, ix, iy, is"<<std::endl;
//        out.precision(10);
//        out<<std::showpos;
//        out<<std::scientific;
//        for(int i=0; i<n_element-1; ++i) {
//            out<<my_debug.at(i).a<<' '<<my_debug.at(i).b<<' '<<my_debug.at(i).c<<' '<<my_debug.at(i).ax<<' '
//            <<my_debug.at(i).bx<<' '<<my_debug.at(i).ay<<' '<<my_debug.at(i).by<<' '<<my_debug.at(i).al<<' '
//            <<my_debug.at(i).bl<<' '<<my_debug.at(i).ix<<' '<<my_debug.at(i).iy<<' '<<my_debug.at(i).is<<std::endl;
//        }
//        out.close();
    }

    if(k_>0)
        ibs_coupling(rx, ry, k_, beam.emit_nx(), beam.emit_ny());
}



void IBSSolver_BM_Complete::init_optc(const Lattice &lattice) {
    int n = lattice.n_element();
    optc.resize(n);
    for(int i=0; i<n; ++i) {
        Optc o;
        double alfx = lattice.alfx(i);
        double betx = lattice.betx(i);
        double dx = lattice.dx(i);
        double dpx = lattice.dpx(i);
        double alfy = lattice.alfy(i);
        double bety = lattice.bety(i);
        double dy = lattice.dy(i);
        double dpy = lattice.dpy(i);
        o.phix = dpx+alfx*dx/betx;
        o.phiy = dpy+alfy*dy/bety;
        o.hx = (dx*dx + betx*betx*o.phix*o.phix)/betx;
        o.hy = (dy*dy + bety*bety*o.phiy*o.phiy)/bety;

        optc.at(i) = o;
    }
}

double IBSSolver_BM_Complete::det(std::array<std::array<double, 3>,3>& l) {
    double d = (l[0][0]*l[1][1]*l[2][2]) + (l[0][1]*l[1][2]*l[2][0]) + (l[0][2]*l[1][0]*l[2][1]) -
      (l[0][2]*l[1][1]*l[2][0]) - (l[0][0]*l[1][2]*l[2][1]) - (l[0][1]*l[1][0]*l[2][2]);
    return d;
}

double IBSSolver_BM_Complete::inv(std::array<std::array<double, 3>,3>& l, std::array<std::array<double, 3>,3>& v) {
    double l_det = det(l);

    v[0][0] = l[1][1]*l[2][2]-l[1][2]*l[2][1];
    v[0][1] = l[2][0]*l[1][2]-l[1][0]*l[2][2];
    v[0][2] = 0;
    v[1][0] = v[0][1];
    v[1][1] = l[0][0]*l[2][2]-l[0][2]*l[2][0];
    v[1][2] = l[2][0]*l[0][1]-l[0][0]*l[2][1];
    v[2][0] = 0;
    v[2][1] = v[1][2];
    v[2][2] = l[0][0]*l[1][1]-l[0][1]*l[1][0];

    for(auto& aa: v)
        for(auto& a: aa)
            a /= l_det;

    return l_det;
}

void IBSSolver_BM_Complete::calc_beam_const(const Beam& beam) {
    gamma = beam.gamma();
    gamma2 = gamma*gamma;
    inv_ex = 1/beam.emit_x();
    inv_ey = 1/beam.emit_y();
    inv_dp2 = 1/(beam.dp_p()*beam.dp_p());
}

void IBSSolver_BM_Complete::calc_l(const Lattice& lattice, int i) {
    double betx = lattice.betx(i);
    double bety = lattice.bety(i);
    double phix = optc.at(i).phix;
    double phiy = optc.at(i).phiy;
    double hx = optc.at(i).hx;
    double hy = optc.at(i).hy;

    lh[0][0] = betx*inv_ex;
    lh[0][1] = -gamma*betx*phix*inv_ex;
    lh[1][0] = lh[0][1];
    lh[1][1] = gamma2*hx*inv_ex;

    lv[1][1] = gamma2*hy*inv_ey;
    lv[1][2] = -gamma*bety*phiy/inv_ey;
    lv[2][1] = lv[1][2];
    lv[2][2] = bety*inv_ey;

    ls[1][1] = gamma2*inv_dp2;
}

double IBSSolver_BM_Complete::coef(const Lattice &lattice, const Beam &beam) const {
    double lambda;
    if(beam.bunched()) lambda = beam.particle_number()/beam.sigma_s();
    else lambda = 2*sqrt(k_pi)*beam.particle_number()/lattice.circ();
     double beta3 = beam.beta();
    beta3 = beta3*beta3*beta3;
    double gamma4 = beam.gamma();
    gamma4 *= gamma4;
    gamma4 *= gamma4;
    return k_c*beam.r()*beam.r()*lambda/(8*k_pi*beam.dp_p()*beta3*gamma4*beam.emit_x()*beam.emit_y());
}

void IBSSolver_BM_Complete::calc_itgl(int i) {
    ii = {};
    double u = 0;
    for(int j=0; j<3; ++j) {
        for(int k=0; k<3; ++k) {
            l[j][k] = lh[j][k] + lv[j][k] + ls[j][k];
            if(u<fabs(l[j][k])) u = fabs(l[j][k]);
        }
    }

    u *= factor;
    double d = u/nt_;

    for(int j=0; j<3; ++j)
            l[j][j] -= d/2;

    for(int i=0; i<nt_; ++i) {
        double lambda = (i+0.5)*d;
        for(int j=0; j<3; ++j)
            l[j][j] += d;

        double det = inv(l, ll);
        double tr= trace(ll);
        double c = sqrt(lambda)/sqrt(det);

        for(int j=0; j<3; ++j)
            for(int k=0; k<3; ++k)
                ii[j][k] += c*((j==k?1:0)*tr - 3*ll[j][k]);
    }
    for(int j=0; j<3; ++j)
        for(int k=0; k<3; ++k)
            ii[j][k] *= d;
}

void IBSSolver_BM_Complete::rate(const Lattice &lattice, const Beam &beam, double &rx, double &ry, double &rs) {
    int n_element = lattice.n_element();

    if (cache_invalid) {
        init_optc(lattice);
        cache_invalid = false;
    }

    double c_bmc = coef(lattice, beam)*log_c()/lattice.circ();

    rx = 0;
    ry = 0;
    rs = 0;
    int n=2;
    if (beam.bunched()) n=1;
    calc_beam_const(beam);

    if(ibs_by_element) {
        std::ofstream out;
        out.open("ibs_by_element.txt");
        ibs_by_element_sddshead(out, n_element);
        out.precision(10);
        out<<std::showpos;
        out<<std::scientific;
        double s = 0;


        for(int i=0; i<n_element-1; ++i){
            double l_element = lattice.l_element(i);
            calc_l(lattice, i);
            calc_itgl(i);
            double rxi{}, ryi{}, rsi{};
            for(int i=0; i<3; ++i) {
                for(int j=0; j<3; ++j) {
                    rxi += lh[i][j]*ii[i][j]*c_bmc;
                    ryi += lv[i][j]*ii[i][j]*c_bmc;
                    rsi += ls[i][j]*ii[i][j]*c_bmc*n;
                }
            }

            rx += rxi;
            ry += ryi;
            rs += rsi;
            out<<s<<' '<<lattice.betx(i)<<' '<<lattice.bety(i)<<' '<<lattice.alfx(i)<<' '
                <<lattice.alfy(i)<<' '<<lattice.dx(i)<<' '<<lattice.dy(i)<<' '
                <<rxi<<' '<<ryi<<' '<<rsi<<' '<<rx<<' '<<ry<<' '<<rs<<std::endl;
            s += l_element;
        }
        out.close();
    }
    else {
        for(int i=0; i<n_element-1; ++i){
            double l_element = lattice.l_element(i);
            calc_l(lattice, i);
            calc_itgl(i);

            for(int i=0; i<3; ++i) {
                for(int j=0; j<3; ++j) {
                    rx += lh[i][j]*ii[i][j];
                    ry += lv[i][j]*ii[i][j];
                    rs += ls[i][j]*ii[i][j];
                }
            }
        }

        rs *= n*c_bmc;
        rx *= c_bmc;
        ry *= c_bmc;
    }

    if(k_>0)
        ibs_coupling(rx, ry, k_, beam.emit_nx(), beam.emit_ny());
}

IBSSolver_BM_Complete::IBSSolver_BM_Complete(int nt, double log_c, double k) : IBSSolver(log_c, k), nt_(nt) {
}
