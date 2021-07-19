//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

//! \file gridgeod.cpp
//! This file contains several functions that sample disc photons, then send the sampled photons to
//! routines that handle single (super)photons (such as \ref tridsca_sp in \ref scatter.cpp).

#include <algorithm>
#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <memory>
#include <fstream>
#include <iterator>
#include "diskobj.h"
#include "sim5lib.h"
#include "utils.h"
#include "quadroots.h"
#include "tridgeo.h"
#include "geoinf.h"
#include "const.h"
#include "detector.h"
#include "kerr.h"
#include "scatter.h"
#include "kerr.h"
#include "photon_dist_utils.h"
                        
void mkgrid(const double rin, const double rout, const double thetamax, const long nr, const long ntheta, 
    const long nphi, std::vector<double> & radius, std::vector<double> & theta, std::vector<double> & phi) {
    
    double dlogr = std::log(rout / rin) / (double)(nr), logrmin = std::log(rin);
    std::vector<double> logradius(nr);
    
    theta.clear();
    theta.resize(ntheta);
    phi.clear();
    phi.resize(nphi);
    radius.clear();
    radius.resize(nr);
    
    long n = 0;
    std::generate(logradius.begin(), logradius.end(), [&n, dlogr, logrmin]{return double(n++) * dlogr + logrmin + dlogr / 2.;});
    std::transform(logradius.cbegin(), logradius.cend(), radius.begin(), [](double x){return std::exp(x);});
    
    double dphi = 2. * M_PI / (double)(nphi);
    double dtheta = thetamax * M_PI / (double)(ntheta);
    n = 0;
    std::generate(theta.begin(), theta.end(), [&n, dtheta]{return (double)(n++) * dtheta + dtheta / 2.;});
    n = 0;
    std::generate(phi.begin(), phi.end(), [&n, dphi]{return (double)(n++) * dphi + dphi / 2.;});
}

/*!
\brief Write to a file the parameters required to reconstruct energy and polarisation spectrum of thin disc on the equatorial plane.

Although this function has `nt` in its name, it is intended for any kind of emissivity profile.
We write the parameters into a binary parameter file. All data are in double-precision.
The file contains a *header* with 6 numbers, followed by N segments where N is the number of pixels that arrives at infinity.
The structure of the file is thus as follows:
\f{equation*}
{a, r_{\rm out}, 0.,N_r, N_\theta, N_\phi, [r, g, \mu_\infty, w, {\rm cos}\Delta\psi, {\rm sin}\Delta\psi]\times N},
\f}
where in the *header*: \f$r_{\rm out}\f$ is the disc outer edge radius, \f$N_r, N_\theta, N_\phi\f$ are number of bins (see \ref calmuobsarr).
And in each segment: \f$r\f$ is the emission radius, \f$g\equiv E_\infty/E_{\rm disc}\f$ is the redshift factor, \f$\mu_\infty\f$ is the photon
inclination at infinity, \f$w\f$ is the photon statistical weight, and \f$\Delta\psi\f$ is the rotation angle of polarisation vector from disc to
infinity. Here \f$w = {\rm cos}\theta_{\rm em} d\mathcal{S} d\Omega/
u^t\f$, where , \f$\theta_{\rm em}\f$ is the local emission angle,
\f$\mathcal{S}\f$ is the proper area, and \f$d\Omega\f$ is the solid angle,
and \f$u^t\f$ is the factor to convert local time to distant observer's time. The photon generation rate per unit distant observer's time
\f$\dot{N} = \frac{4\zeta(3) k_{\rm B}^3}{c^2h^3}\frac{f_{\rm limb} T_{\rm eff}^3}{f_{\rm col}u^t} w\f$,
where \f$\zeta()\f$ is the Riemann zeta function, \f$f_{\rm limb}\f$ is the limb-darkening factor, \f$T_{\rm eff}\f$ is the effective temperature,
and \f$f_{\rm col}\f$ is the color correction factor.

@param pol the option to turn on/off polarisation calculation. Note that the directionarity is also affected by this option. If off,
    we assume isotropic emission from the disc in the disc fluid rest frame. Otherwise the angular dependence follows Chandrasekhar 1960
    (http://adsabs.harvard.edu/abs/1960ratr.book.....C) assuming a pure scattering atmosphere with infinite optical depth, which can be 
    calculated using \ref calchandra_darkening.
@param a BH spin
@param rout disk outer disc radius, in \f$[\rm GM~c^{-2}]\f$
@param thetamax the maximum polar emission angle \f$\theta_{\rm max}\f$; see \ref calmuobsarr.
@param nr number of radial bins \f$N_r\f$
@param ntheta number of \f$\theta\f$ bins
@param nphi number of \f$\phi\f$ bins
@param paramfile the name of the output parameter file.
 */
void griddisk_nt_sp_pol_writeparam(const bool pol, const double rin, const double rout, 
    const double thetamax, const long nr, const long ntheta, const long nphi, const diskobj & disk) {
    
    long indexnow;
    double cost, dphi, dtheta, dr, one_over_dr, domega, vphi, gamma;
    double area, one_over_ut, Omega, kerrx, weightnorm, weightnow, gdisk, gdisk1, x1, omega1, knorm, muobs_re, dpsi, ang_ktodisk, ang_disktok;
    double c1, c2, cospsiinf, sinpsiinf, cospsi, sinpsi, sint;
    std::vector<double> muobsarr, radius, theta, phi, rr1arr, muplusarr, larr, qarr, kmuobssignarr, muarr(ntheta);
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::complex<double> wp;
    auto dsize = sizeof(double);

    geoinf geo;
    sim5metric met, met1;
    sim5tetrad tet;
    std::array<double, 4> umu, umu1, on_k, bl_k, kmu_bl_in, angles;
    
    //std::cout << "Calculating geodesics..." << std::endl;
    //calmuobsarr(a, disk.rms, rout, thetamax, nr, ntheta, nphi, muobsarr, radius, theta, phi, rr1arr, muplusarr, larr, qarr, kmuobssignarr);
    //std::cout << "Done!" << std::endl;
    disk.mkrgrid(nr, rin, rout, radius);
    mklingrid(ntheta, 0., thetamax * M_PI, theta);
    std::transform(theta.cbegin(), theta.cend(), muarr.begin(), [](double x){return std::abs(std::cos(x));});
    //std::cout << "theta = " << theta << std::endl;
    mklingrid(nphi, 0., 2. * M_PI, phi);
    
    dphi = phi[1] - phi[0];
    dtheta = theta[1] - theta[0];
    dr = std::sqrt(radius[1] / radius[0]);
    one_over_dr = 1. / dr;
    
    std::ofstream disc_ofile("disk_params.dat", std::ios::out | std::ofstream::binary);
    std::ofstream selfirr_ofile("selfirr_params.dat", std::ios::out | std::ofstream::binary);
    
    std::vector<double> input_params = {disk.a, rout, 0.,
        (double)(nr), (double)(ntheta), (double)(nphi)};
    disc_ofile.write((char*)&input_params[0], input_params.size() * dsize);
    selfirr_ofile.write((char*)&input_params[0], input_params.size() * dsize);
    
    on_k[0] = 1.;
    
    std::cout << "Calculating spectrum..." << std::endl;
    for (size_t ir = 0; ir < radius.size(); ++ir) {
        area = disk.calds(radius[ir] * one_over_dr, radius[ir] * dr);
        vphi = disk.vphi(radius[ir]);
        gamma = 1. / std::sqrt(1. - vphi * vphi);
        kerrx = std::sqrt(radius[ir]);
        Omega = 1. / (kerrx * kerrx * kerrx + disk.a);
        one_over_ut = std::sqrt(1. - 2. / radius[ir] + 4. * Omega * disk.a / radius[ir] - Omega * Omega * (radius[ir] * radius[ir] + 
            disk.a * disk.a + 2. * disk.a * disk.a / radius[ir]));
        weightnorm = area * gamma * one_over_ut;

        kerr_metric(disk.a, radius[ir], 0., &met);
        tetrad_azimuthal(&met, Omega, &tet);
        std::copy(std::begin(tet.e[0]), std::end(tet.e[0]), umu.begin());
        /*
        for (size_t ix = 0; ix < 4; ++ix) {
            for (size_t iy = 0; iy <4; ++iy) {
                if (iy >= ix) {
                    std::cout << "tet.e[" << ix << "] * tet.e[" << iy << "]= " << dotprod(tet.e[ix], tet.e[iy], &met) << std::endl;
                }
            }
        }
        */

        fourvelocity_azimuthal(Omega, &met, umu1.data());
        //std::cout << "umu = " << umu << std::endl;
        //std::cout << "umu1 = " << umu1 << std::endl;
        //std::cout << "================" << std::endl;
        
        //std::cout << ir + 1 << "-th out of " << radius.size() << " radial bins\r";
        for (size_t it = 0; it < theta.size(); ++it) {
            cost = std::cos(theta[it]);
            on_k[2] = cost;
            sint = std::sqrt(1. - cost * cost);
            domega = std::abs(cost) * dphi * (std::cos(theta[it] - dtheta / 2.) - std::cos(theta[it] + dtheta / 2.));
            weightnow = weightnorm * domega;
            
            for (size_t ip = 0; ip < phi.size(); ++ip) {
                on_k[1] = sint * std::cos(phi[ip]);
                on_k[3] = sint * std::sin(phi[ip]);
                on2bl(on_k.data(), bl_k.data(), &tet);
                
                // normalize wave vector such that k_t = -1;
                knorm = -1. / (met.g00 * bl_k[0] + met.g03 * bl_k[3]);
                std::transform(bl_k.cbegin(), bl_k.cend(), bl_k.begin(), std::bind2nd(std::multiplies<double>(), knorm));
                gdisk = -1. * dotprod(bl_k.data(), umu.data(), &met);
                
                //std::cout << "bl_k * bl_k = " << dotprod(bl_k.data(), bl_k.data(), &met) << std::endl;
                geo = geoinf(disk.a, radius[ir], 0., bl_k);
                muobs_re = geo.calfate(radius[ir], 0., bl_k[1], bl_k[2], disk.rms);
                indexnow = ir * theta.size() * phi.size() + it * phi.size() + ip;

                if (geo.status == 1) {
                    gdisk = -1. * dotprod(umu.data(), bl_k.data(), &met);
                    
                    if (pol) {
                        disk.polrotang_disktok(radius[ir], theta[it], phi[ip], geo.l, geo.q, c1, c2);
                        disk.polrotang_ktoinf(geo.kmuobssign, geo.l, geo.q, muobs_re, cospsiinf, sinpsiinf);
                        cospsi = c1 * cospsiinf - c2 * sinpsiinf;
                        sinpsi = c1 * sinpsiinf + c2 * cospsiinf;
                        //std::cout << cospsi * cospsi + sinpsi * sinpsi << std::endl;
                    } else {
                        cospsi = 0.;
                        sinpsi = 0.;
                    }
                    std::vector<double> params = {radius[ir], theta[it], 1. / gdisk, muobs_re, weightnow, cospsi, sinpsi,
                        geo.l, geo.q, geo.kmuobssign};
                    disc_ofile.write(reinterpret_cast<char*>(&params[0]), dsize * params.size());
                } else if (geo.status == 2) {
                    ang_disktok = disk.polrotang_disktok(radius[ir], theta[it], phi[ip], geo.l, geo.q);
                    disk.angles(muobs_re, geo.l, geo.q, geo.krobssign, angles);
                    photon_momentum(disk.a, muobs_re, 0., geo.l, geo.q, geo.krobssign, geo.kmuobssign, kmu_bl_in.data());
                    kerr_metric(disk.a, muobs_re, 0., &met1);
                    x1 = std::sqrt(muobs_re);
                    omega1 = 1. / (x1 * x1 * x1 + disk.a);
                    fourvelocity_azimuthal(omega1, &met1, umu1.data());
                    gdisk1 = -1. * dotprod(umu1.data(), kmu_bl_in.data(), &met1);
                    ang_ktodisk = -1. * disk.polrotang_disktok(muobs_re, geo.l, geo.q, angles); 
                    dpsi = ang_disktok + ang_ktodisk;
                    std::vector<double> selfirr_params = {radius[ir], theta[it], weightnow, muobs_re,
                        kmu_bl_in[0], kmu_bl_in[1], kmu_bl_in[2], kmu_bl_in[3], gdisk1 / gdisk, dpsi};
                    selfirr_ofile.write(reinterpret_cast<char*>(&selfirr_params[0]), selfirr_params.size() * dsize);
                }
                showprogress((double)(indexnow) / (double)(nr * ntheta * nphi));
            }
        }
    }
    std::cout << std::endl;
    disc_ofile.close();
    selfirr_ofile.close();
}

/*!
\brief Register blackbody photons, given blackbody temperatures, weight, inclination, and normalised Stokes parameters.

The weight written into the file is \f$w^\prime = \frac{1}{N_{\rm photon}} \frac{f_{\rm limb} {\rm cos}\theta_{\rm em} d\mathcal{S} d\Omega
 T_{\rm eff}^3}{f_{\rm col}u^t}\f$. To obtain the photon generation rate: \f$\dot{N}=\frac{4\zeta(3) k_{\rm B}^3}{c^2h^3} w^\prime\f$.

@param pol switch for polarization calculations
@param nphoton nubmer of photons 
@param tem color temperature as observed by the observer
@param weightnorm normalization of weight
@param qfac normalised Stokes parameter Q
@param ufac normalised Stokes parameter U
@param muinf \f$\mu_{\rm obs}\f$
@param gen random number generator
@param ofiles vector of output files; see \ref register_sp_nt.
*/
void register_sp_bbody(const bool pol, const long nphoton, const double tem, const double weightnorm, const double qfac, const double ufac, const
    double muinf, const double l, const double q, const double ktheta, std::mt19937_64 & gen, std::vector<std::ofstream> & ofiles) {
    
    size_t double_size = sizeof(double);
    double eninf, muinf1 = muinf, l1 = l, q1 = q, ktheta1 = ktheta;
    double weight = weightnorm / (double)(nphoton);
    double qweight = qfac * weight, uweight = ufac * weight;
    
    if (muinf > -1.5) {
        for (long i = 0; i < nphoton; ++i) {
            eninf = sample_bb(tem, gen);
            ofiles[0].write(reinterpret_cast<char*>(&eninf), double_size);
            ofiles[1].write(reinterpret_cast<char*>(&weight), double_size);
            ofiles[2].write(reinterpret_cast<char*>(&muinf1), double_size);
            ofiles[3].write(reinterpret_cast<char*>(&l1), double_size);
            ofiles[4].write(reinterpret_cast<char*>(&q1), double_size);
            ofiles[5].write(reinterpret_cast<char*>(&ktheta1), double_size);
            if (pol) {
                ofiles[6].write(reinterpret_cast<char*>(&qweight), double_size);
                ofiles[7].write(reinterpret_cast<char*>(&uweight), double_size);
            }
        }
    }
}

/*!
\brief With the geodesic information saved in the parameter file produced by \ref griddisk_nt_sp_pol_writeparam, this function
calculates the values for Novikov-Thorne disc profile.

For each pixel that arrives at infinity, this function first calcualtes the effective temperature and the statistical weight 
given \f$a, M, \dot{M}, f_{\rm col}\f$. We assume Novikov-Thorne temperature profile (Page & Thorne 1974). 
If the polarisation option is switched off, the normalised Stokes parameters are also calculated and the limb-darkening factor is for
semi-infinite pure scattering disc (Chandrasekhar 1960; see also \ref calchandra_darkening). If the option is on we assume isotropic emission.
We multiply the weight \f$w\f$ calculated by \ref griddisk_nt_sp_pol_writeparam with a factor to obtain \f$w^\prime\f$:
\f$ w^\prime = \frac{1}{N_{\rm photon}} \frac{f_{\rm limb} T_{\rm eff}^3}{f_{\rm col}} w \f$.
The temperature, weight \f$w^\prime\f$, and Stokes parameters are passed to \ref register_sp_bbody.

@param pol whether to turn on polarisation
@param parafile the input parameter file
@param m BH mass in solar mass
@param mdot mass accretion rate in Eddington rate
@param fcol color correction factor
@param nphoton number of photons
@param gen random number generator
@param ofiles vector of std::ofstream, pointing to the output files:
    - ofiles[0]: dimensionless energies of the superphotons
    - ofiles[1]: statistical weights of the superphotons
    - ofiles[2]: photon inclination angle at infinity
    - If `pol` is turned on, `ofiles` has two additional elements:
    - ofiles[3]: Stokes parameter Q
    - ofiles[4]: Stokes parameter U
    
 */
void register_sp_nt(const bool pol, const std::string parafile, const double m, const double mdot, const double fcol,
    const double rtr, const long nphoton, std::mt19937_64 & gen, std::vector<std::ofstream> & ofiles) {
    
    std::vector<double> params, theta, muarr, darkarr, poldegarr;
    double a, teff, tcol_emass, qfac, ufac, flimb, dtheta, weightnorm, ftheta, fphi, cos2psi, sin2psi, dit, poldegnow, theta1;
    rdoublevec(params, parafile);
    long ntheta, n = 0, nparam = 10, nseg, offset = 6, pos, it0, it1;
    nseg = (params.size() - offset) / nparam;
    
    // ntheta = (long)(params[4]);
    ntheta = 501;
    dtheta = 0.5 * M_PI / (double)(ntheta - 1);
    
    a = params[0];
    truncated_nt disk(a, m, mdot, rtr);
    //std::cout << "truncated disc " << std::endl;
    //std::cout << "a = " << a << std::endl;
    //std::cout << "m = " << m << std::endl;
    
    if (pol) {
        theta.resize(ntheta);
        std::generate(theta.begin(), theta.end(), [&n, dtheta]{return (double)(n++) * dtheta;});
        muarr.resize(ntheta);
        darkarr.resize(ntheta);
        std::transform(theta.cbegin(), theta.cend(), muarr.begin(), [](double x){return std::abs(std::cos(x));});
        poldegarr.resize(ntheta);
        calpoldeg(muarr, poldegarr);
        calchandra_darkening(muarr, darkarr);
        //wdoublevec(theta, "poltheta.dat");
        //wdoublevec(muarr, "polmu.dat");
        //wdoublevec(poldegarr, "polpol.dat");
        //std::cout << "theta = " << theta << std::endl;
        //std::cout << "poldeg = " << poldegarr << std::endl;
    } else {
        qfac = 0.;
        ufac = 0.;
    }
    std::cout << "nseg = " << nseg << std::endl;
    
    for (long i = 0; i < nseg; ++i) {
        pos = i * nparam + offset;
        //std::cout << "radius = " << params[pos] << std::endl;
        teff = disk.teff_kev(params[pos]);
        //std::cout << "teff = " << teff << std::endl;
        tcol_emass = fcol * teff / ME_KEV * params[pos + 2];;
        if (pol) {
            theta1 = (params[pos + 1] <= 0.5 * M_PI) ? params[pos + 1] : M_PI - params[pos + 1];
            dit = (theta1 - 0.5 * dtheta) / dtheta;
            it0 = (long)(std::floor(dit));
            it1 = it0 + 1;
            //std::cout << "theta1 = " << theta1 << std::endl;
            //std::cout << "it0 = " << it0 << ", it1 = " << it1 << std::endl;
            flimb = (darkarr[it0] * (theta1 - theta[it0]) + darkarr[it1] * (theta[it1] - theta1)) / dtheta;
            poldegnow = (poldegarr[it0] * (theta1 - theta[it0]) + poldegarr[it1] * (theta[it1] - theta1)) / dtheta;
            //std::cout << "thetanow = " << theta1 << ", poldegnow = " << poldegnow << std::endl;
            
            ftheta = -1. * params[pos + 6];
            fphi = params[pos + 5];
            sin2psi = 2. * ftheta * fphi;
            cos2psi = 2. * ftheta * ftheta - 1.;
            qfac = poldegnow * cos2psi;
            ufac = poldegnow * sin2psi;
            
        } else {
            flimb = 1.;
            qfac = 0.;
            ufac = 0.;
        }
        //std::cout << "muinf = " << params[pos + 3] << std::endl;
        weightnorm = params[pos + 4] * disk.m * disk.m  * flimb * teff * teff * teff / fcol;
        //std::cout << "params[pos + 4] = " << params[pos + 4] << ", weightnorm = " << weightnorm << std::endl;
        register_sp_bbody(pol, nphoton, tcol_emass, weightnorm, qfac, ufac, params[pos + 3], params[pos + 7], 
			params[pos + 8], params[pos + 9], gen, ofiles);
        showprogress((double)(i + 1) / (double)(nseg));
    }
    std::cout << std::endl;
}

/*!
\brief Samples disc photon with \ref calmuobsarr, and according to photon destination writes photon info into the following two files:
- `disc_params.dat`: for photons arriving at infinity without entering the corona; similar with \ref griddisk_nt_sp_pol_writeparam.
- `sca_params.dat`: for photons entering the corona: the parameter file has a *header*, followed by \f$N\f$ segments where \f$N\f$ 
    is the number of geodesic that enters the
    corona. The structure of the file: 
    \f{equation*}a, r_{\rm out}, ctype, N_r, N_\theta, N_\phi,
    [r, \theta, \mu_\infty, w, r_{\rm hit},\mu_{\rm hit}, k_{\rm hit}^t, k_{\rm hit}^r, k_{\rm hit}^\theta,
    k_{\rm hit}^\phi, g, c_1, c_2, |WP|, {\rm cos}\Delta\psi, {\rm sin}\Delta\psi]\times N,
    \f}
    here *ctype* is the corona type, \f$r\f$ is the
    emission radius, \f$\theta\f$ is the local emission angle, \f$w\f$ is the statistical weight (same as in \ref griddisk_nt_sp_pol_writeparam).
    \f$r_{\rm hit}\f$ and \f$\mu_{\rm hit}\f$ are the coordinates where the photon enters the corona, \f$\boldsymbol{k}_{\rm hit}\f$ is the
    photon wave vector in Boyer-Lindquist frame while entering the corona, \f$g\equiv E_\infty/E_{\rm disc}\f$ is the redshift factor,
    \f$c_1\f$ and \f$c_2\f$ are coefficients of the rotation matrix from polarsation vector at the disc to 
    \f$k^\prime\f$ (see \ref nt::polrotang_disktok), \f$|WP|\f$ is the modulus of the Walker-Penrose constant,
    and \f$\Delta\psi\f$ is the angle between
    \f$k^\prime\f$ to the polarisation vector at infinity (see \ref nt::polrotang_ktoinf). The file is in binary format with double precision.
- `selfirr_params.dat`: for self-irradiating photons. This file starts with a same header as `sca_params.dat`, followed by \f$N_s\f$ segments
    where \f$N_s\f$ is the number of geodesics that return to the disc. The structure of the file:
    \f{equation*}a, r_{\rm out}, ctype, N_r, N_\theta, N_\phi,
    [r_e, \theta, w, r_i, k^t, k^r, k^\theta, k^\phi, g\equiv E_i/E_e, \Delta\psi]\times N_s,
    \f}
    here the subscripts `e` and `i` denotes emission and incidence, respectively.
    
On *gtype*:
- 0: \ref zamosphere3d
- 1: \ref zamoslab3d
- 2: \ref zamosphere3d_truncated
- 3: \ref rotsphere3d
- 5: \ref kepslab3d
- 6: \ref wedge3d
- 7: \ref wedge3d_zamo_uniform
- 9: \ref conical
    
@param rin inner edge radius of the thin disc. If rin < rms then rin = rms.
@param rout outer boundary of thin disc
@param thetamax maximum local emission angle (see \ref calmuobsarr).
@param disk diskobj object
@param trid \ref tridgeo object; for different types of corona
@param nr number of radial bin on the disc plane
@param ntheta number of theta bin
@param nphi number of phi bin
*/
void grid3dcorona_sp_writeparam(const double rin, const double rout, const double thetamax, const diskobj & disk, 
    const tridgeo & trid, const long nr, const long ntheta, const long nphi, const std::string discfile, 
    const std::string selffile, const std::string scafile) {
    
    bool hit;
    int status, gtype;
    long indexnow;
    double knorm;
    double dr, one_over_dr, rnow, dphi, dtheta, gdisk, vphi, omega, area, weightnorm, weightnorm0;
    double x, sint, domega, gamma, one_over_ut;
    double rhit, muhit, wpnorm;
    double cospsi, sinpsi, c1, c2, cospsiinf, sinpsiinf;
    double d_coronatype;
    double muobs_re, dpsi, ang_disktok, ang_ktodisk, x1, omega1, gdisk1;
    std::array<double, 4> angles, kmu_bl_in, umu1, umu;
    
    if (trid.name == "zamosphere3d") {
        gtype = 0;
    } else if (trid.name == "zamoslab3d") {
        gtype = 1;
    } else if (trid.name == "zamosphere3d_truncated") {
        gtype = 2;
    } else if (trid.name == "rotsphere3d") {
        gtype = 3;
    } else if (trid.name == "kepslab3d") {
        gtype = 5;
    } else if (trid.name == "wedge3d") {
        gtype = 6;
    } else if (trid.name == "wedge3d_zamo_uniform") {
        gtype = 7;
    } else if (trid.name == "wedge3d_corotation_uniform") {
        gtype = 8;
    } else if (trid.name == "conical") {
        gtype = 9;
    } else if (trid.name == "arbitrary") {
        gtype = 10;
    } else if (trid.name == "jed") {
        gtype = 10;
	}	else if (trid.name == "zamosandwich") {
        gtype = 11;
    } else if (trid.name == "kepsandwich") {
        gtype = 12;
    } else {
        std::cerr << "Invalid corona type!" << std::endl;
        return;
    }
    d_coronatype = (double)(gtype);
    //std::cout << trid.rmax << " " << trid.rmin << " " << trid.mumin << " " << trid.mumax << std::endl;
    
    auto dsize = sizeof(double);
    
    // grid vectors
    std::vector<double> muobsarr, rr1arr, muplusarr, larr, qarr, kmuobssignarr;
    
    std::complex<double> wpnow;
    std::array<double, 4> on_k, bl_k, kmu_hit;
    std::vector<double> radius(nr), logradius(nr), theta(ntheta), phi(nphi), muarr(ntheta);
    
    std::ofstream disk_pfile(discfile, std::ios::out | std::ios::binary);
    std::ofstream sca_pfile(scafile, std::ios::out | std::ios::binary);
    std::ofstream selfirr_pfile(selffile, std::ios::out | std::ios::binary);
    
    // write params
    std::vector<double> input_params = {trid.a, rout, d_coronatype, 
        (double)(nr), (double)(ntheta), (double)(nphi)};
    disk_pfile.write((char*)&input_params[0], input_params.size() * dsize);
    sca_pfile.write((char*)&input_params[0], input_params.size() * dsize);
    selfirr_pfile.write((char*)&input_params[0], input_params.size() * dsize);
    
    /*
    std::cout << "trid.a = " << trid.a << std::endl;
    std::ofstream logfile("input.log", std::ios::out);
    logfile << "Corona type = " << trid.name << std::endl << 
        "rmin = " << trid.rmin << std::endl << "rmax = " << trid.rmax << std::endl;
    logfile << "a = " << trid.a << std::endl << "rout = " << rout << std::endl;
    logfile << "nr = " << nr << std::endl << "ntheta = " << ntheta << std::endl << "nphi = " << nphi << std::endl;
    logfile.close();
    */
    
    sim5metric met1;
    sim5tetrad tet;
    geoinf geo;
    bool print = false;
    
    //calmuobsarr(trid.a, rin, rout, thetamax, nr, ntheta, nphi, muobsarr, radius, theta, phi, rr1arr, muplusarr, larr, qarr, kmuobssignarr);
    //std::cout << "Done!" << std::endl;
    disk.mkrgrid(nr, rin, rout, radius);
    dr = std::sqrt(radius[1] / radius[0]);
    one_over_dr = 1. / dr;

    mklingrid(ntheta, 0., thetamax * M_PI, theta);
    dtheta = theta[1] - theta[0];
    std::transform(theta.cbegin(), theta.cend(), muarr.begin(), [](const double x){return std::cos(x);});
    mklingrid(nphi, 0., 2. * M_PI, phi);
    dphi = phi[1] - phi[0];
    
    if (print) {
        std::cout << "muarr = " << muarr << std::endl;
    }
    
    on_k[0] = 1.;
    
    std::cout << "Calculating geodesics to corona..." << std::endl;
    for (size_t ir = 0; ir < radius.size(); ++ir) {
        //std::cout << "ir = " << ir << std::endl;
        rnow = radius[ir];

        area = disk.calds(rnow * one_over_dr, rnow * dr);
        vphi = disk.vphi(rnow);
        gamma = 1. / std::sqrt(1. - vphi * vphi);
        
        x = std::sqrt(rnow);
        omega = 1. / (x * x * x + trid.a);
        
        one_over_ut = std::sqrt(1. - 2. / rnow + 4. * omega * trid.a / rnow - omega * omega * (rnow * rnow + 
            trid.a * trid.a + 2. * trid.a * trid.a / rnow));
        
        sim5metric met;
        kerr_metric(trid.a, radius[ir], 0., &met);
		tetrad_azimuthal(&met, omega, &tet);
        std::copy(std::begin(tet.e[0]), std::end(tet.e[0]), umu.begin());

        weightnorm0 = area * gamma * one_over_ut;

        for (size_t it = 0; it < theta.size(); ++it) {
            sint = std::sin(theta[it]);
            on_k[2] = muarr[it];
            domega = std::abs(muarr[it]) * dphi * (std::cos(theta[it] - dtheta / 2.) - std::cos(theta[it] + dtheta / 2.));
            weightnorm = weightnorm0 * domega;
            
            for (size_t ip = 0; ip < phi.size(); ++ip) {
                //print = ((it == 3 || it == 4) && ip == 5);
                indexnow = ir * theta.size() * phi.size() + it * phi.size() + ip;
                //std::cout << "ir = " << ir << ", it = " << it << ", ip = " << ip << std::endl;
                on_k[1] = sint * std::cos(phi[ip]);
                on_k[3] = sint * std::sin(phi[ip]);
                //std::cout << "on_k * on_k = " << on_k[0] * on_k[0] - on_k[1] * on_k[1] - on_k[2] * on_k[2] - on_k[3] * on_k[3] << std::endl;
                on2bl(on_k.data(), bl_k.data(), &tet);
                // normalize wave vector such that k_t = -1;
                knorm = -1. / (met.g00 * bl_k[0] + met.g03 * bl_k[3]);
                std::transform(bl_k.cbegin(), bl_k.cend(), bl_k.begin(), std::bind2nd(std::multiplies<double>(), knorm));
                geo = geoinf(trid.a, radius[ir], 0., bl_k);
                muobs_re = geo.calfate(radius[ir], 0., bl_k[1], bl_k[2], rin);
                //std::cout << "kmu * kmu = " << dotprod(bl_k.data(), bl_k.data(), &met) << std::endl;
                
                if (print) {
                    std::cout << "=================================" << std::endl;
                    std::cout << "it = " << it << "ip = " << ip << std::endl;
                    std::cout << "phinow = " << phi[ip] << std::endl;
                    std::cout << "on_k = " << on_k << std::endl;
                    std::cout << "bl_k = " << bl_k << std::endl;
                    std::cout << "geo.l = " << geo.l << ", geo.q = " << geo.q << std::endl;
                    std::cout << "geo.muplus = " << geo.muplus << ", geo.muminus = " << geo.muminus << std::endl;
                }
                
                gdisk = -1. * dotprod(umu.data(), bl_k.data(), &met);
                //std::cout << "bl_k = " << bl_k << std::endl;
                
                if (geo.l != geo.l || geo.q != geo.q) {
                    std::cout << "Inside grid3dcorona_sp; Checking geo.l:" << std::endl;
                    std::cout << "r = " << radius[ir] << ", mue = 0" << std::endl;
                    std::cout << "bl_k = " << bl_k << std::endl;
                    break;
                }
                
                disk.polrotang_disktok(radius[ir], theta[it], phi[ip], geo.l, geo.q, c1, c2);
                
                if (geo.status == 1) {
                    disk.polrotang_ktoinf(geo.kmuobssign, geo.l, geo.q, muobs_re, cospsiinf, sinpsiinf);
                    cospsi = c1 * cospsiinf - c2 * sinpsiinf;
                    sinpsi = c1 * sinpsiinf + c2 * cospsiinf;
                }
                wpnorm = std::sqrt(geo.q + (geo.l - trid.a) * (geo.l - trid.a));
                
                    
                if (gtype == 6 || gtype == 7 || gtype == 8 || (gtype == 5 && std::abs(trid.mumin) <= 1e-4)) {
                    hit = true;
                } else {
                    hit = false;
                }
                // If it is not possible for the photon to reach the corona
                if (!geo.ispossible(trid, radius[ir], 0., bl_k[1], bl_k[2])) {
                    if (print) {
                        std::cout << "not possible" << std::endl;
                    }
                    hit = false;
                } else {
                    status = geo.findxtrid(trid, radius[ir], 0., bl_k[1], bl_k[2], rin, rhit, muhit, kmu_hit);
                    //catch(const char * error) {
                    //    std::cout << error << std::endl;
                    //}
                    
                    if (status != 1) {
                        hit = false;
                        if (print) {
                            std::cout << status << std::endl;
                        }
                    } else {
                        //std::cout << "muobsarr[indexnow] = " <<  muobsarr[indexnow] << std::endl;
                        hit = true;
                        if (print) {
                            std::cout << "Hit! " << std::endl;
                        }
                        //kerr_metric(trid.a, rhit, muhit, &metnow);
                        //trid.calumu(rhit, muhit, umu, metnow);
                        //ghit = dotprod(umu.data(), kmu_hit.data(), &metnow);
                        //gfac = -1. * ghit * gdisk;
                        std::vector<double> sca_params = {radius[ir], theta[it], (geo.status == 1) ? muobs_re : -2.,
                            weightnorm, rhit, muhit, kmu_hit[0], kmu_hit[1], kmu_hit[2], kmu_hit[3], 1. / gdisk,
                            c1, c2, wpnorm, cospsi, sinpsi};
                        sca_pfile.write((char*)&sca_params[0], sca_params.size() * dsize);
                    }
                }
                
                if (!hit) {
                    if (geo.status == 1) {
                        std::vector<double> disk_params = {radius[ir], theta[it], 1. / gdisk, geo.muobs, weightnorm,
                            cospsi, sinpsi, geo.l, geo.q, geo.kmuobssign};
                        disk_pfile.write((char*)&disk_params[0], disk_params.size() * dsize);
                    } else if (geo.status == 2) {
                        ang_disktok = disk.polrotang_disktok(radius[ir], theta[it], phi[ip], geo.l, geo.q);
                        disk.angles(muobs_re, geo.l, geo.q, geo.krobssign, angles);
                        photon_momentum(trid.a, muobs_re, 0., geo.l, geo.q, geo.krobssign, geo.kmuobssign, kmu_bl_in.data());
                        kerr_metric(trid.a, muobs_re, 0., &met1);
                        x1 = std::sqrt(muobs_re);
                        omega1 = 1. / (x1 * x1 * x1 + trid.a);
                        fourvelocity_azimuthal(omega1, &met1, umu1.data());
                        gdisk1 = -1. * dotprod(umu1.data(), kmu_bl_in.data(), &met1);
                        ang_ktodisk = -1. * disk.polrotang_disktok(muobs_re, geo.l, geo.q, angles); 
                        dpsi = ang_disktok + ang_ktodisk;
                        std::vector<double> selfirr_params = {radius[ir], theta[it], weightnorm, muobs_re,
                            kmu_bl_in[0], kmu_bl_in[1], kmu_bl_in[2], kmu_bl_in[3], gdisk1 / gdisk, dpsi};
                        selfirr_pfile.write((char*)&selfirr_params[0], selfirr_params.size() * dsize);
                        std::complex<double> wp, wp1, diff;
                        wp = disk.calwp_angle(radius[ir], theta[it], phi[ip], geo.l, 30. / 180. * M_PI);
                        wp1 = disk.calwp_angle(muobs_re, geo.l, 30. / 180. * M_PI + dpsi, angles);
                        diff = wp - wp1;
                        /*
                        if (std::abs(diff) > 1e-4) {
                            std::cout << "wp = " << wp << ", wp1 = " << wp1 << std::endl;
                        }
                        */
                    }
                        
                }
                showprogress((double)(indexnow) / (double)(radius.size() * theta.size() * phi.size()));
            }
        }
    }
    std::cout << "Done!" << std::endl;
    disk_pfile.close();
    sca_pfile.close();
    selfirr_pfile.close();
}

/*!
\brief With the `sca_params.dat` produced by \ref register_sp_tridsca, propagate the photons entering the corona until they arrive at infinity
/strike at the disc/enter the BH event horizon.

@param pol whether to include polarisation. If turned off, all photons are assumed to be unpolarised when being scattered by electron
@param progress whether to print progress
@param parafile the path to `sca_params.dat`
@param m BH mass in solar mass
@param mdot mass accretion rate in Eddington rate (see \ref nt::nt(const double,const double,const double)).
@param fcol color correction factor
@param te electron temperature in \f$[\rm keV]\f$
@param mfp Thompson mean free path
@param dr step size in \f$f\f$ while propagating inside the corona
@param nphoton number of photon per geodesic
@param ofiles vector of \ref std::ofstream that point to output files; all files are in binary format and contain data in double precision.
    ofile[8] and ofiles[9] are written only if `pol` is turned on.
    - ofiles[0]: photon energy at infinity \f$E_\infty\f$
    - ofiles[1]: photon statistical weight \f$w^\prime\f$ (see \ref register_sp_nt)
    - ofiles[2]: 
        - inclination at infinity, if the photon finally arrives at infinity
        - disc radius, if the photon hits the disc;
    - ofiles[3]: \f$l\equiv L_z/E_\infty\f$, where \f$L_z\f$ is photon angular momentum about black hole spin axis
    - ofiles[4]: \f$\mathcal{Q}/E_\infty^2\f$, where \f$\mathcal{Q}\f$ is one constant of motion found by Carter 1968
    - ofiles[5]: sign of \f$k^r\f$
    - ofiles[6]: sign of \f$k^\theta\f$
    - ofiles[7]: only if `pol` is turned on; Stokes parameter Q \f$(q w^\prime)\f$, where \f$q\f$ is one normalised Stokes parametere
    - ofiles[8]: only if `pol` is turned on; Stokes parameter U \f$(u w^\prime)\f$, where \f$u\f$ is another normalised Stokes parametere
@param nscafile one \ref std::ofstream object that point to a binary file that contains photon scattering number in `int` format (note that
the actual size of `int` type depends on compiler and/or architect)
 */
void register_sp_tridsca(const bool pol, const bool KN, const bool chandra, const bool progress, 
    const std::string parafile, const double m, const double mdot, const double fcol,
    const double te, const double tau, const double dr, const long nphoton, std::vector<std::ofstream> & ofiles_inf,
    std::ofstream & nscafile_inf, std::vector<std::ofstream> & ofiles_disc, std::ofstream & nscafile_disc) {
    
    double ecut = 10., ratio = 10., ennow, weightnow, theta1, dit;
    
    std::vector<double> params, theta, muarr, darkarr, poldegarr;
    double a, teff, tcol_emass, flimb, dtheta, poldegnow, K1, K2, ftheta, fphi, cos2psi, sin2psi;
    rdoublevec(params, parafile);
    long ntheta, n = 0, nparam = 16, it0, it1, nseg, offset = 6;
    long gtype;
    nseg = (params.size() - offset) / nparam;
    std::array<double, 4> kmunow;
    std::complex<double> wpnow;
    
    double weightnorm, qfac, ufac;
    
    std::unique_ptr<tridgeo> tridptr(nullptr);
    std::vector<double> gparams;
    rdoublevec(gparams, "../gparams.dat");
    a = params[0];
    gtype = (long)(std::round(params[2]));
    
    uniform_thermal edist(te * ME_KEV);
    
    switch(gtype) {
        case 0 : {
            tridptr.reset(new zamosphere3d(a, gparams[0], gparams[1], tau));
            break;
        }
        case 1 : {
            tridptr.reset(new zamoslab3d(a, gparams[0], gparams[1], gparams[2], tau));
            break;
        }
        case 2 : {
            tridptr.reset(new zamosphere3d_truncated(a, gparams[0], gparams[1], tau));
            break;
        }
        case 3 : {
            tridptr.reset(new rotsphere3d(a, gparams[0], gparams[1], gparams[2], tau));
            break;
        }
        case 5 : {
            tridptr.reset(new kepslab3d(a, gparams[0], gparams[1], gparams[2], tau));
            break;
        }
        case 6 : {
            tridptr.reset(new wedge3d(a, gparams[0], gparams[1], gparams[2], tau));
            break;
        }
        case 7 : {
            tridptr.reset(new wedge3d_zamo_uniform(a, gparams[0], gparams[1], gparams[2], tau));
            break;
        }
        case 8 : {
            tridptr.reset(new wedge3d_corotation_uniform(a, gparams[0], gparams[1], gparams[2], tau));
            break;
        }
        case 9 : {
            tridptr.reset(new conical(a, gparams[0], gparams[1], gparams[2], gparams[3], tau));
            break;
        }
        default : {
            std::cerr << "gridgeod.cpp, register_sp_tridsca(): Invalid gtype value of " << gtype << std::endl;
            return;
        }
    }
    //zamosphere3d trid(a, gparams[0], gparams[1]);
    nt disk(a, m, mdot);
    //std::cout << "total photons = " << nphoton1 * nseg << std::endl;
    
    std::ofstream logfile("input.log", std::ios::out);
    logfile << "Corona type = " << tridptr->name << std::endl << 
        "rmin = " << tridptr->rmin << std::endl << "rmax = " << tridptr->rmax << std::endl;
    logfile << "te = " << te << std::endl << "tau = " << tau << std::endl;
    logfile << "pol = " << pol << std::endl << "a = " << tridptr->a << std::endl << "m = " << disk.m << std::endl
        << "mdot = " << disk.mdot << std::endl << "fcol = " << fcol << std::endl;
    logfile << "ntheta = " << ntheta << std::endl << "nphoton = " << nphoton << std::endl;
    logfile.close();
    
    // ntheta = (long)(params[4]);
    ntheta = 501;
    dtheta = 0.5 * M_PI / (double)(ntheta - 1);
    
    theta.resize(ntheta);
    std::generate(theta.begin(), theta.end(), [&n, dtheta]{return (double)(n++) * dtheta;});
    muarr.resize(ntheta);
    darkarr.resize(ntheta);
    std::transform(theta.cbegin(), theta.cend(), muarr.begin(), [](double x){return std::abs(std::cos(x));});
    poldegarr.resize(ntheta);
    calpoldeg(muarr, poldegarr);
    calchandra_darkening(muarr, darkarr);
    
#ifdef USEOPENMP
#pragma omp parallel
{
#pragma omp single
{
#endif
    std::mt19937_64 gen(std::random_device {}());
    for (long i = 0; i < nseg; ++i) {
        //std::cout << "i = " << i << std::endl;
        auto pos = i * nparam + offset;
        std::copy(params.cbegin() + pos + 6, params.cbegin() + pos + 10, kmunow.begin());
        
        teff = disk.teff_kev(params[pos]);
        //std::cout << "(T / Tin)^3 = " << teff / tin << std::endl;
        tcol_emass = fcol * teff / ME_KEV * params[pos + 10];;
        
        theta1 = (params[pos + 1] <= 0.5 * M_PI) ? params[pos + 1] : M_PI - params[pos + 1];
        dit = (theta1 - 0.5 * dtheta) / dtheta;
        it0 = (long)(std::floor(dit));
        it1 = it0 + 1;
        if (chandra) {
            flimb = (darkarr[it0] * (theta1 - theta[it0]) + darkarr[it1] * (theta[it1] - theta1)) / dtheta;
            poldegnow = (poldegarr[it0] * (theta1 - theta[it0]) + poldegarr[it1] * (theta[it1] - theta1)) / dtheta;
        } else {
            flimb = 1.;
            poldegnow = 0.;
        }
        //std::cout << "param[i * nparam + 1] = " << params[i * nparam + 1] << std::endl;
        //std::cout << "theta[it] = "  << theta[it] << std::endl;
        
        K1 = -1. * params[pos + 13] * params[pos + 12];
        K2 = -1. * params[pos + 13] * params[pos + 11];
        wpnow.real(K1);
        wpnow.imag(K2);
        ftheta = -1. * params[pos + 15];
        fphi = params[pos + 14];
        sin2psi = 2. * ftheta * fphi;
        cos2psi = 2. * ftheta * ftheta - 1.;
        qfac = poldegnow * cos2psi;
        ufac = poldegnow * sin2psi;
            
        weightnorm = params[pos + 3] * disk.m * disk.m  * flimb * teff * teff * teff / fcol / (double)(nphoton);
        
        for (long ip = 0; ip < nphoton; ++ip) {
            superphoton sp;
            sp.en = sample_bb(tcol_emass, gen);
            weightnow = 1.;
            //std::cout << "sp.en = " << sp.en << std::endl;
                
            sp.kmu = kmunow;
            sp.rnow = params[pos + 4];
            sp.munow = params[pos + 5];
            sp.poldeg = poldegnow;
            sp.wp = wpnow;
            sp.nscatter = 0;
            sp.weight = weightnow;
#ifdef USEOPENMP
#pragma omp task
            //std::cout << "params[pos + 2] = " << params[pos + 2] << std::endl;
            tridsca_sp_openmp(pol, KN, disk.a, disk.rms, edist, weightnorm, qfac, ufac, params[pos + 2], dr, sp,
                *tridptr, ofiles_inf, nscafile_inf, ofiles_disc, nscafile_disc);
#endif
        }
        if (progress) {
            showprogress((double)(i) / (double)(nseg));
        }
   }
    std::cout << std::endl;
#ifdef USEOPENMP
}
}
#endif
}


void plotgeo(const bool savevec, const long nphoton, const tridgeo & trid, electron_dist & edist, 
    const double xmax, const double ymin, const double ymax, std::mt19937_64 & gen, std::vector<std::ofstream> & ofiles) {
    
    std::uniform_real_distribution<> dis(0., 1.);
    double snow, hnow, rnow, munow, xnow, ynow, phinow, dennow, temnow, bfield, te_k;
    auto dsize = sizeof(double);
    std::array<double, 4> umu;
    sim5metric met;
    
    for (long i = 0; i < nphoton; ++i) {
        hnow = dis(gen) * (ymax - ymin) + ymin;
        snow = dis(gen) * xmax;
        rnow = std::sqrt(snow * snow + hnow * hnow);
        munow = hnow / rnow;
        //std::cout << "hnow = " << hnow << std::endl;
        //std::cout << "snow = " << snow << std::endl;
        
        if (trid.inside(rnow, munow)) {
            phinow = 2. * M_PI * dis(gen);
            xnow = snow * std::cos(phinow);
            ynow = snow * std::sin(phinow);
            dennow = trid.ne_sigmat(rnow, munow) / SIGMA_THOMSON;
            temnow = edist.character(rnow, munow);
            ofiles[0].write(reinterpret_cast<char *>(&xnow), dsize);
            ofiles[1].write(reinterpret_cast<char *>(&ynow), dsize);
            ofiles[2].write(reinterpret_cast<char *>(&hnow), dsize);
            ofiles[3].write(reinterpret_cast<char *>(&dennow), dsize);
            ofiles[4].write(reinterpret_cast<char *>(&temnow), dsize);
            if (trid.magnetised) {
                te_k = temnow / K_KEV;
                bfield = trid.cal_bfield(rnow, munow, te_k);
                ofiles[5].write(reinterpret_cast<char *>(&bfield), dsize);
            }
            if (savevec) {
                kerr_metric(trid.a, rnow, munow, &met);
                trid.calumu(rnow, munow, umu, met);
                ofiles[6].write(reinterpret_cast<char *>(&umu[1]), dsize);
                ofiles[7].write(reinterpret_cast<char *>(&umu[2]), dsize);
            }
        }
    }
}

void grid_ns_surface_writeparam(const bool pol, const double mns, const double rns_rg, const double pspin_ms,
    const long nincl, const long ntheta, const long nphi) {

    long indexnow;
    double cost, dphi, dtheta, domega;
    double area, weightnow, gdisk, knorm, muobs_re;
    std::vector<double> muobsarr, radius, theta, phi, rr1arr, muplusarr, larr, qarr, kmuobssignarr, muarr(ntheta), inclarr;
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::complex<double> wp;
    auto dsize = sizeof(double);
    double rns_cm = rns_rg * mns * RG;
    double a = 0.47 / pspin_ms;
    std::cout << "a = " << a << std::endl;

    double pspin_s = pspin_ms / 1000.;
    double Omega_persec = 2. * M_PI / pspin_s;
    double Omega = Omega_persec * (RG * mns / LSPEED);
    std::cout << "pspin_s = " << pspin_s << std::endl;
    std::cout << "Omega = " << Omega << std::endl;
    std::cout << "rns_rg = " << rns_rg << std::endl;
    double weightnorm = 4. * APERY_CONST / LSPEED / LSPEED / H_KEV / H_KEV / H_KEV;

    geoinf geo;
    sim5metric met;
    sim5tetrad tet;
    std::array<double, 4> umu, on_k, bl_k;

    mklingrid(nincl, 0., M_PI, inclarr);
    std::transform(theta.cbegin(), theta.cend(), muarr.begin(), [](double x){return std::abs(std::cos(x));});
    mklingrid(ntheta, 0., 0.5 * M_PI, theta);
    mklingrid(nphi, 0., 2. * M_PI, phi);

    double dincl = inclarr[1] - inclarr[0];
    dphi = phi[1] - phi[0];
    dtheta = theta[1] - theta[0];
    //std::vector<double> input_params = {trid.a, mns, rns_rg, pspin_ms, d_coronatype, inclmax,
    //    (double)(nincl), (double)(ntheta), (double)(nphi)};

    std::ofstream ofile("surface_params.dat", std::ios::out | std::ofstream::binary);
    std::vector<double> input_params = {a, mns, rns_rg, pspin_ms, -1., 1., (double)(nincl), (double)(ntheta), (double)(nphi)};
    ofile.write((char*)&input_params[0], input_params.size() * dsize);

    double cosincl;

    on_k[0] = 1.;

    std::cout << "Calculating spectrum..." << std::endl;
    for (size_t iincl = 0; iincl < inclarr.size(); ++iincl) {
        // we start from newtonian area weightnorm =
        //area = 2. * M_PI * std::sin(inclarr[iincl]) * dincl * rns_cm * rns_cm;
        area = 2. * M_PI * (std::cos(inclarr[iincl] - dincl / 2.) - std::cos(inclarr[iincl] + dincl / 2.))
            * rns_cm * rns_cm;
        cosincl = std::cos(inclarr[iincl]);
        kerr_metric(a, rns_rg, cosincl, &met);

        fourvelocity_azimuthal(Omega, &met, umu.data());
        tetrad_azimuthal(&met, Omega, &tet);

        for (size_t it = 0; it < theta.size(); ++it) {
            cost = std::cos(theta[it]);
            on_k[1] = cost;
            double sint = std::sqrt(1. - cost * cost);
            domega = dphi * (std::cos(theta[it] - dtheta / 2.) - std::cos(theta[it] + dtheta / 2.));
            weightnow = area * domega * weightnorm * cost;
            //std::cout << "area = " << area << std::endl;
            //std::cout << "domega = " << domega << std::endl;

            for (size_t ip = 0; ip < phi.size(); ++ip) {
                on_k[2] = sint * std::cos(phi[ip]);
                on_k[3] = sint * std::sin(phi[ip]);
                on2bl(on_k.data(), bl_k.data(), &tet);

                // normalize wave vector such that k_t = -1;
                knorm = -1. / (met.g00 * bl_k[0] + met.g03 * bl_k[3]);
                std::transform(bl_k.cbegin(), bl_k.cend(), bl_k.begin(), std::bind2nd(std::multiplies<double>(), knorm));
                //gdisk = -1. * dotprod(bl_k.data(), umu.data(), &met);

                //std::cout << "bl_k * bl_k = " << dotprod(bl_k.data(), bl_k.data(), &met) << std::endl;
                try {
                    geo = geoinf(a, rns_rg, cosincl, bl_k);
                    muobs_re = geo.calfate(rns_rg, cosincl, bl_k[1], bl_k[2], 1e6);
                    // here
                    indexnow = iincl * theta.size() * phi.size() + it * phi.size() + ip;

                    if (geo.status == 1) {
                        gdisk = -1. * dotprod(umu.data(), bl_k.data(), &met);
                        //std::vector<double> surface_params = {cosincl, 1. / gdisk, geo.muobs, weightnow,
                        //    geo.l, geo.q, geo.kmuobssign};

                        std::vector<double> params = {cosincl, 1. / gdisk, muobs_re, weightnow,
                            geo.l, geo.q, geo.kmuobssign};
                        ofile.write(reinterpret_cast<char*>(&params[0]), dsize * params.size());
                    }
                }
                catch (const char * error) {
                    std::cout << "incl = " << inclarr[iincl] << std::endl;
                    std::cout << "theta = " << theta[it] << std::endl;
                    std::cout << "phi = " << phi[ip] << std::endl;
                }
                showprogress((double)(indexnow) / (double)(nincl * ntheta * nphi));
            }
        }
    }
    std::cout << std::endl;
    ofile.close();
}

void register_sp_ns_surface(const std::string parafile, const long nphoton,
    const double tns_kev, std::mt19937_64 & gen, std::vector<std::ofstream> & ofiles) {

    std::vector<double> params;
    double tcol_emass;
    rdoublevec(params, parafile);
    long nparam = 7, nseg, offset = 9, pos;
    nseg = (params.size() - offset) / nparam;

    std::cout << "nseg = " << nseg << std::endl;

    //std::vector<double> surface_params = {cosincl, 1. / gdisk, geo.muobs, weightnow,
    //    geo.l, geo.q, geo.kmuobssign};

    for (long i = 0; i < nseg; ++i) {
        pos = i * nparam + offset;
        //std::cout << "radius = " << params[pos] << std::endl;
        //std::cout << "teff = " << teff << std::endl;
        tcol_emass = tns_kev / ME_KEV * params[pos + 1];;
        //std::cout << tcol_emass << std::endl;
        //std::cout << params[pos + 1] << std::endl;
        //std::cout << "muinf = " << params[pos + 3] << std::endl;
        register_sp_bbody(false, nphoton, tcol_emass, params[pos + 3] * tns_kev * tns_kev * tns_kev / (double)(nphoton),
            0., 0., params[pos + 2], params[pos + 4],
            params[pos + 5], params[pos + 6], gen, ofiles);
        showprogress((double)(i + 1) / (double)(nseg));
    }
    std::cout << std::endl;
}

void grid3dcorona_sp_writeparam_ns(const double mns, const double rns_rg, const double pspin_ms, const double rin,
    const double inclmax, const tridgeo & trid, const long nincl, 
    const long ntheta, const long nphi, const std::string scafile, const std::string discfile) {
    
    bool hit;
    int status, gtype;
    long indexnow;
    double knorm, dphi, dtheta, gdisk, area, sint, weightnow, domega;
    double rhit, muhit, d_coronatype, muobs_re;
    std::array<double, 4> umu;
    double rns_cm = RG * mns * rns_rg;

    if (trid.name == "zamosphere3d") {
        gtype = 0;
    } else if (trid.name == "zamoslab3d") {
        gtype = 1;
    } else if (trid.name == "zamosphere3d_truncated") {
        gtype = 2;
    } else if (trid.name == "rotsphere3d") {
        gtype = 3;
    } else if (trid.name == "kepslab3d") {
        gtype = 5;
    } else if (trid.name == "wedge3d") {
        gtype = 6;
    } else if (trid.name == "wedge3d_zamo_uniform") {
        gtype = 7;
    } else if (trid.name == "wedge3d_corotation_uniform") {
        gtype = 8;
    } else if (trid.name == "conical") {
        gtype = 9;
    } else if (trid.name == "arbitrary") {
        gtype = 10;
    } else if (trid.name == "jed") {
        gtype = 10;
	}	else if (trid.name == "zamosandwich") {
        gtype = 11;
    } else if (trid.name == "kepsandwich") {
        gtype = 12;
    } else {
        std::cerr << "Invalid corona type!" << std::endl;
        return;
    }
    d_coronatype = (double)(gtype);
    
    //std::cout << trid.rmax << " " << trid.rmin << " " << trid.mumin << " " << trid.mumax << std::endl;
    double weightnorm = 4. * APERY_CONST / LSPEED / LSPEED / H_KEV / H_KEV / H_KEV;
    double dntotal = (double)(nincl * ntheta * nphi);
    
    auto dsize = sizeof(double);
    
    // grid vectors
    std::vector<double> muobsarr, rr1arr, muplusarr, larr, qarr, kmuobssignarr, inclarr;
    
    std::array<double, 4> on_k, bl_k, kmu_hit;
    std::vector<double> incl(nincl), theta(ntheta), phi(nphi), muarr(ntheta);
    
    std::ofstream disk_pfile(discfile, std::ios::out | std::ios::binary);
    std::ofstream sca_pfile(scafile, std::ios::out | std::ios::binary);
    
    // write params
    std::vector<double> input_params = {trid.a, mns, rns_rg, pspin_ms, d_coronatype, inclmax,
        (double)(nincl), (double)(ntheta), (double)(nphi)};
    disk_pfile.write((char*)&input_params[0], input_params.size() * dsize);
    sca_pfile.write((char*)&input_params[0], input_params.size() * dsize);
    
    /*
    std::cout << "trid.a = " << trid.a << std::endl;
    std::ofstream logfile("input.log", std::ios::out);
    logfile << "Corona type = " << trid.name << std::endl << 
        "rmin = " << trid.rmin << std::endl << "rmax = " << trid.rmax << std::endl;
    logfile << "a = " << trid.a << std::endl << "rout = " << rout << std::endl;
    logfile << "nr = " << nr << std::endl << "ntheta = " << ntheta << std::endl << "nphi = " << nphi << std::endl;
    logfile.close();
    */
    
    sim5tetrad tet;
    geoinf geo;
    bool print = false;

    mklingrid(nincl, 0., inclmax * M_PI, inclarr);
    double dincl = inclarr[1] - inclarr[0];
    mklingrid(ntheta, 0., 0.5 * M_PI, theta);
    dtheta = theta[1] - theta[0];
    std::transform(theta.cbegin(), theta.cend(), muarr.begin(), [](const double x){return std::cos(x);});
    mklingrid(nphi, 0., 2. * M_PI, phi);
    dphi = phi[1] - phi[0];
    
    if (print) {
        std::cout << "muarr = " << muarr << std::endl;
    }
    
    on_k[0] = 1.;

    double pspin_s = pspin_ms / 1000.;
    double Omega_persec = 2. * M_PI / pspin_s;
    double Omega = Omega_persec * (RG * mns / LSPEED);
    
    std::cout << "Calculating geodesics to corona..." << std::endl;
    for (size_t iincl = 0; iincl < inclarr.size(); ++iincl) {
        area = 2. * M_PI * (std::cos(inclarr[iincl] - dincl / 2.) - std::cos(inclarr[iincl] + dincl / 2.))
            * rns_cm * rns_cm;
        double cosincl = std::cos(inclarr[iincl]);

        sim5metric met;
        kerr_metric(trid.a, rns_rg, cosincl, &met);
        tetrad_azimuthal(&met, Omega, &tet);
        fourvelocity_azimuthal(Omega, &met, umu.data());

        for (size_t it = 0; it < theta.size(); ++it) {
            sint = std::sin(theta[it]);
            on_k[1] = muarr[it];
            domega = std::abs(muarr[it]) * dphi * (std::cos(theta[it] - dtheta / 2.) - std::cos(theta[it] + dtheta / 2.));
            weightnow = area * domega * weightnorm;
            
            for (size_t ip = 0; ip < phi.size(); ++ip) {
                indexnow = iincl * theta.size() * phi.size() + it * phi.size() + ip;
                on_k[2] = sint * std::cos(phi[ip]);
                on_k[3] = sint * std::sin(phi[ip]);
                on2bl(on_k.data(), bl_k.data(), &tet);
                // normalize wave vector such that k_t = -1;
                knorm = -1. / (met.g00 * bl_k[0] + met.g03 * bl_k[3]);
                std::transform(bl_k.cbegin(), bl_k.cend(), bl_k.begin(), std::bind2nd(std::multiplies<double>(), knorm));
                //std::cout << "bl_k = " << bl_k << std::endl;
                geo = geoinf(trid.a, rns_rg, cosincl, bl_k);
                muobs_re = geo.calfate(rns_rg, cosincl, bl_k[1], bl_k[2], rin);
                //std::cout << "kmu * kmu = " << dotprod(bl_k.data(), bl_k.data(), &met) << std::endl;
                
                if (print) {
                    std::cout << "=================================" << std::endl;
                    std::cout << "it = " << it << "ip = " << ip << std::endl;
                    std::cout << "phinow = " << phi[ip] << std::endl;
                    std::cout << "on_k = " << on_k << std::endl;
                    std::cout << "bl_k = " << bl_k << std::endl;
                    std::cout << "geo.l = " << geo.l << ", geo.q = " << geo.q << std::endl;
                    std::cout << "geo.muplus = " << geo.muplus << ", geo.muminus = " << geo.muminus << std::endl;
                }
                
                gdisk = -1. * dotprod(umu.data(), bl_k.data(), &met);
                //std::cout << "bl_k = " << bl_k << std::endl;
                
                if (geo.l != geo.l || geo.q != geo.q) {
                    std::cout << "Inside grid3dcorona_sp; Checking geo.l:" << std::endl;
                    std::cout << "r = " << rns_rg<< ", mue = " << cosincl << std::endl;
                    std::cout << "bl_k = " << bl_k << std::endl;
                    break;
                }
                
                hit = false;
                // If it is not possible for the photon to reach the corona
                if (!geo.ispossible(trid, rns_rg, cosincl, bl_k[1], bl_k[2])) {
                    if (print) {
                        std::cout << "not possible" << std::endl;
                    }
                    hit = false;
                } else {
                    status = geo.findxtrid(trid, rns_rg, cosincl, bl_k[1], bl_k[2], rin, rhit, muhit, kmu_hit);
                    /*
                    std::array<double, 4> pos = {0., blr, blmu, 0.}, kmu1(bl_k);
                    int status1 = geo.findxtrid_sim5(rin, trid, pos, kmu1);
                    if ((status == 1 || status1 == 1) && status != status1) {
                        std::cout << pos << ", " << kmu1 << std::endl;
                        std::cout << trid.inside(pos[1], pos[2]) << std::endl;
                        std::cout << "status = " << status << ", status1 = " << status1 << std::endl;
                    }
                    */
                    //catch(const char * error) {
                    //    std::cout << error << std::endl;
                    //}
                    
                    if (status != 1) {
                        hit = false;
                        if (print) {
                            std::cout << status << std::endl;
                        }
                    } else {
                        //std::cout << "muobsarr[indexnow] = " <<  muobsarr[indexnow] << std::endl;
                        hit = true;
                        if (print) {
                            std::cout << "Hit! " << std::endl;
                        }
                        //kerr_metric(trid.a, rhit, muhit, &metnow);
                        //trid.calumu(rhit, muhit, umu, metnow);
                        //ghit = dotprod(umu.data(), kmu_hit.data(), &metnow);
                        //gfac = -1. * ghit * gdisk;
                        std::vector<double> sca_params = {cosincl, (geo.status == 1) ? muobs_re : -2.,
                            weightnow, rhit, muhit, kmu_hit[0], kmu_hit[1], kmu_hit[2], kmu_hit[3], 1. / gdisk};
                        sca_pfile.write((char*)&sca_params[0], sca_params.size() * dsize);
                    }
                }

                if (!hit && (geo.status == 1)) {
                    std::vector<double> surface_params = {cosincl, 1. / gdisk, geo.muobs, weightnow,
                        geo.l, geo.q, geo.kmuobssign};
                    disk_pfile.write((char*)&surface_params[0], surface_params.size() * dsize);
                }
                showprogress((double)(indexnow) / dntotal);
            }
        }
    }
    std::cout << "Done!" << std::endl;
    disk_pfile.close();
    sca_pfile.close();
}
