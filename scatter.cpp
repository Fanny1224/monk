//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

//! \file scatter.cpp
//! Contains various functions dealing with compton scattering.
#include <vector>
#include <array>
#include <iostream>
#include <random>
#include <algorithm>
#include <complex>
//#include <omp.h>
#include "superphoton.h"
#include "const.h"
#include "scatter.h"
#include "rootsearch.h"
#include "utils.h"
#include "geoinf.h"
#include "tridgeo.h"
#include "diskobj.h"
#include "detector.h"
#include "kerr.h"
#include "photon_dist.h"
#include "sim5lib.h"
#include "electron_population.h"

const double RMAX = 5e4;

//! Differential Klein-Nishina cross section in electron rest frame
//! @param a dimensionless photon energy before scattering
//! @param ap dimensionless photon energy after scattering
double klein_nishina(double a, double ap) {
    double ch, kn;

    ch = 1. + 1. / a - 1. / ap;
    kn = (a / ap + ap / a - 1. + ch * ch) / (a * a);

    return (kn);
}

/*! \brief Probability density distribution of `xp`, dimensionless energy of the scattered photon
@param x dimensionless energy of the incoming photon
@param xp dimensionless energy of the scattered photon
 */
double KN_pdfx(const double x, const double xp) {
    double mu = 1. + 1. / x - 1. / xp;
    return 0.375 / (x * x) * (x / xp + xp / x - 1. + mu * mu) / KN_xsection(x);
}

/*!
\brief Sample differential KN cross section (\ref KN_pdfx),
following Kahn 1964 (https://www.rand.org/pubs/research_memoranda/RM1237.html).
@param x dimensionless energy of the incoming photon
@param gen random number generator
*/
double sample_KN_xp(const double x, std::mt19937_64 & gen) {
    std::uniform_real_distribution<> dis(0., 1.);
    double mu, z, reci;
    double r1, r2, r3;
    bool found = false;
    
    while (!found) {
        r1 = dis(gen);
        r2 = dis(gen);
        r3 = dis(gen);
        if (r1 <= (2. * x + 1) / (2. * x + 9)) {
            z = 1. + 2. * x * r2;
            reci = 1. / z;
            if (r3 <= 4. * (reci - reci * reci)) {
                found = true;
            }
        } else {
            z = (2. * x + 1) / (1. + 2. * x * r2);
            mu = 1. + (1. - z) / x;
            if (r3 <= 0.5 * (mu * mu + 1. / z)) {
                found = true;
            }
        }
    }
    return x / z;
}

/*! 
\brief Evaluate the function \f$cdf(\mu) - \lambda\f$, where \f$cdf(\mu)\f$ is the cdf. of
\f$\mu\f$, and \f$\lambda\f$ is a random number drawn from an uniform distribution between 0 and 1.

Reference: Schnittman & Krolik 2013 (http://adsabs.harvard.edu/abs/2013ApJ...777...11S).
@param mu cosine of polar angle
@param lambda a uniform random number between 0 and 1
*/
double xsecmucdf(double mu, const double & lambda) {
    return mu * mu * mu + 3. * mu + 4. - 8. * lambda;
}

/*!
\brief Calculate cdf. of azimuthal angle for Thompson scattering cross section, given
cosine of polar angle and polarization degree

Reference: Schnittman & Krolik 2013.
@param phi azimuthal phi
@param params struct that contains mu and polarization degree
*/
double xsecphicdf(double phi, const xsectphicdfparams & params) {
    return phi / 2. / M_PI + params.poldeg * (params.mu * params.mu - 1.) * std::sin(2. * phi) /
        4. / M_PI / (params.mu * params.mu + 1.) - params.result;
}

/*!
\brief Total Klein-Nishina cross section.

Reference: PSS83. Return \f$\sigma(x) / \sigma_{\rm T}\f$, where
\f$\sigma_{\rm T}=\frac{8}{3}\pi r_0^2\f$ is the Thompson scattering cross section, and \f$r_0\f$ is the classical electron radius.
@param x dimensionless photon energy, as measured in the electron rest frame
*/
double KN_xsection(const double x) {
    double sigma_KN;
    double temp = 1. + 2. * x;
    double reci = 1. / x, reci2 = reci * reci;

    if (x > 1e-4) {
        sigma_KN = 0.375 * reci * ((1. - 2. * reci - 2. * reci2) * std::log(temp) + 0.5 + 4. * reci - 0.5 / temp / temp);
    } else {
        sigma_KN = 1. - 2. * x;
    }

    return sigma_KN;
}

/*
\brief abc
@param gen random number generator
void sample_thompsonx_fixmu(std::mt19937_64 & gen, const double poldeg, const double mu, double & phi) {
    double lambda;
    std::uniform_real_distribution<> dis(0., 1.);

    lambda = dis(gen);
    xsectphicdfparams params;
    params.mu = mu;
    params.poldeg = poldeg;
    params.result = lambda;
    rtbis<xsectphicdfparams>(0., 2. * M_PI, 1e-6, &xsecphicdf, params, phi);
}
*/

/*! \brief 
Sample angular distribution from thompson scattering cross section, following
Equations C11-C16 of Schnittman & Krolik 2013.
@param gen Mersenne-Twister pseudo-random generator
@param poldeg polarization degree
@param mu cosine of scattering polar angle
@param phi scattering azimuthal angle
*/
void sample_thompsonx(std::mt19937_64 & gen, const double poldeg, double & mu, double & phi) {
    double lambda;
    std::uniform_real_distribution<> dis(0., 1.);

    // sampling theta
    lambda = dis(gen);
    rtbis<double>(-1., 1., 1e-6, &xsecmucdf, lambda, mu);
    lambda = dis(gen);

    xsectphicdfparams params;
    params.mu = mu;
    params.poldeg = poldeg;
    params.result = lambda;
    rtbis<xsectphicdfparams>(0., 2. * M_PI, 1e-6, &xsecphicdf, params, phi);
}

/*!
\brief An ancillary function for \ref sample_phi_pol.

Returns \f$cdf(\phi) - \lambda\f$, where \f$cdf(\phi)\f$ is cdf. of \f$\phi\f$, and \f$\lambda\f$ is a random number uniformly distribued 
between 0 and 1.
@param phi scattering azimuthal angle
@param params `pol_phi_params` objects
*/
double pol_phi_cdf(const double phi, const pol_phi_params & params) {
    return 0.5 * phi / M_PI - params.par1 * std::sin(2. * phi) - params.rndnum;
}

/*!
\brief Sample scattering azimuthal angle
@param gen random number generator
@param par1 a function of polarization degree and scattering polar angle; For both Thompson and Klein-Nishina cross sections,
 the cdf. of \f$\phi\f$ can be written into such a form: \f$\rm cdf(\phi) = \phi/2\pi- par1 * {\rm sin}2\phi\f$
*/
double sample_phi_pol(std::mt19937_64 & gen, const double par1) {
    double phi;
    std::uniform_real_distribution<> dis(0., 1.);

    pol_phi_params params;
    params.par1 = par1;
    params.rndnum = dis(gen);
    
    rtbis<pol_phi_params>(0., 2. * M_PI, 1e-6, &pol_phi_cdf, params, phi);
    return phi;
}

/*! \brief 
Sample scattering azimuthal angle for polarized beam, with Thompson cross section
@param gen random number generator
@param poldeg polarization degree
@param mu cosine of scattering polar angle
*/
double sample_thompson_phi(std::mt19937_64 & gen, const double poldeg, const double mu) {
    double par1, musqr = mu * mu, phi;
    par1 = 0.25 * poldeg * (1. - musqr) / (1. + musqr) / M_PI;
    phi = sample_phi_pol(gen, par1);
    return phi;
}

/*! \brief Sample scattering azimuthal angle for polarized beam, with Klein-Nishina cross section
@param gen random number generator
@param poldeg polarization degree
@param mu cosine of scattering polar angle
@param x1_over_x energy ratio of scattered photon to incoming photon
*/
double sample_kn_phi(std::mt19937_64 & gen, const double poldeg, const double mu, const double x1_over_x) {
    double par1, sintsqr = 1. - mu * mu, phi;
    par1 = 0.25 * poldeg * sintsqr / M_PI / (x1_over_x + 1. / x1_over_x - sintsqr);
    phi = sample_phi_pol(gen, par1);
    return phi;
}
    


/*!
\brief Given dimensionless electron velocity \f$\beta\f$, sample \f$\mu\f$, cosine of the angle made by the electron
and the incoming photon.

Given \f$\beta\f$, the pdf. of \f$\mu\f$ is:
\f{equation}
P(\mu) \propto (1 - \mu \beta)\sigma(x^\prime),
\f}  
where \f$\sigma(x^\prime)\f$ the cross section and \f$x^\prime\f$ is the dimensionless photon energy in electron rest frame.
In this function we just sample the pdf.
\f{equation}
P^\prime(\mu) = (1 - \mu \beta)/2.
\f} 
The part of the pdf which is dependent on cross section will be handled in \ref sample_electron_KN or \ref sample_electron_Thompson
with a rejection method.

@param beta electron velocity
@param gen Mersenne-Twister pseudo-random number generator
*/
double sample_scatter_mu(const double beta, std::mt19937_64 & gen) {
    double y, mu;
    // generate random number
    std::uniform_real_distribution<> dis(0., 1.);
    y = dis(gen);
    // solve mu
    mu = (1. - std::sqrt(1. - beta * (4. * y - 2. - beta))) / beta;
    return mu;
}

/*!
\brief Sample momentum of the scattering electron assuming Thompson scattering cross section.
@param te dimensionless electron temperature
@param gen MT random number generator
@param gamma output; electron Lorentz factor
@param mu output; cosine of angle made by photon and electron
@param phi output; electron azimuthal angle with respect to photon polarisation vector
*/
void sample_electron_Thompson(const double te, std::mt19937_64 & gen, double & gamma, double & mu, double & phi) {
    double beta, rndnum;
    std::uniform_real_distribution<> dis(0., 1.);
    rndnum = dis(gen);
    phi = 2. * M_PI * rndnum;
    gamma = sample_maxjutt1(te, gen);
    beta = std::sqrt(1. - 1. / gamma / gamma);
        // then sample mu
    mu = sample_scatter_mu(beta, gen);
}

/*!
\brief Sample momentum of the scattering electron assuming Klein-Nishina scattering cross section.

The probability density of electron momentum
\f{equation}
P(\boldsymbol{p}_e)\propto \frac{dN_e}{d^3p} (1-\mu_e \beta_e)\sigma_{\rm KN}.
\f}
The procedure follows PSS83.
@param te dimensionless eletron temperature
@param x dimensionless photon energy
@param gen Mersenne-Twister pseudo-random number generator
@param gamma output; electron Lorentz factor
@param mu output; cosine of angle made by photon and electron
@param phi output; electron azimuthal angle with respect to photon polarisation vector
*/
void sample_electron_KN(const double te, const double x, std::mt19937_64 & gen, double & gamma, double & mu, double & phi) {
    double beta, xe, total_xsect, rndnum;
    std::uniform_real_distribution<> dis(0., 1.);
    rndnum = dis(gen);
    phi = 2. * M_PI * rndnum;

    while(true) {
        // first sample gamma
        gamma = sample_maxjutt1(te, gen);
        beta = std::sqrt(1. - 1. / gamma / gamma);
        // then sample mu
        mu = sample_scatter_mu(beta, gen);
        xe = x * gamma * (1. - beta * mu);
        //total_xsect = 0.5 * hc_klein_nishina(xe) * (1. - beta * mu);
        total_xsect = KN_xsection(xe);
        rndnum = dis(gen);
        if (rndnum <= total_xsect)
            break;
    }
}

/*!
\brief Given results by \ref sample_electron_KN or \ref sample_electron_Thompson, reconstructs the electron momentum vector.

Note that the angles by the two functions are in a coordinate system with photon wave and polarisation vectors pointing to
\f$z-\f$ and \f$x-\f$ directions, respectively. Hence at the end we need to rotate the electron photon momentum vector back
to fluid frame.

@param emu cosine of angle made by photon and electron
@param ephi electron azimuthal angle with respect to photon polarisation vector
@param kmu photon wave vector
@param pe output; electron four momentum in fluid frame
*/
void calpe(const double emu, const double ephi, const double * kmu, std::array<double, 4> & pe) {
    double sinetheta;
    std::array<double, 3> xvec, yvec, pe_photon;

    dummy_vec(kmu, xvec.data());

    arrcrossprod(kmu, xvec.data(), yvec.data());
    sinetheta = std::sqrt(1. - emu * emu);

    pe_photon[0] = sinetheta * std::cos(ephi);
    pe_photon[1] = sinetheta * std::sin(ephi);
    pe_photon[2] = emu;

    backward_rotation(xvec.data(), yvec.data(), kmu, pe_photon.data(), pe.data() + 1);
}

/*!
\brief Compton scattering in electron rest frame.

Given the photon wave and polarisation vectors in the electron rest frame, polarisation degree, 
and scattering angles, this function calculates the wave and polarisation vectors, and polarisation degree
of the scattered photon. If KN option is on, the polarisation degree and vector are calculated following
Connors+1980 (http://adsabs.harvard.edu/abs/1980ApJ...235..224C); otherwise we use Rayleigh phase matrix.

@param KN whether to assume Klein-Nishia cross section. In this function this option only affects the calculation of
    Stokes parameters of the scattered (super)photon.
@param pk pointer to the spatial part of incoming photon wave vector
@param pf pointer to the spatial part of incoming photon polarization vector
@param coskt cosine of the scattering polar angle
@param kphi cosine of the scattering azimuthal angle
@param poldeg polarization degree of the incoming radiation
@param poldeg_p output; polarzation degree of the outgoing radiation
@param z \f$=x^\prime/x\f$, the ratio of energies of scattered photon to incident photon. 
@param pk_e_p output; pointer to the spatial part of outgoing photon wave vector;
    note that this vector has been already rotated back to electron frame
@param pf_e_p output; pointer to the spatial part of outgoing photon polarization vector;
    note that this vector has been already rotated back to electron frame
*/
void scatter_restframe(const bool KN, const double * pk, const double * pf, const double coskt,
    const double kphi, const double poldeg, double & poldeg_p, const double z, double * pk_e_p, double * pf_e_p) {

    double sinkt, sinkp, coskp, polflux, cosktsqr, sinktsqr, sinkpsqr, coskpsqr;
    double iprime, qprime, uprime, polang_p, cosphipol_p, sinphipol_p;
    std::array<double, 3> yaxis;

    arrcrossprod(pk, pf, yaxis.data());

    cosktsqr = coskt * coskt;
    sinktsqr = 1. - cosktsqr;
    sinkt = std::sqrt(sinktsqr);
    
    sinkp = std::sin(kphi);
    coskp = std::cos(kphi);
    
    sinkpsqr = sinkp * sinkp;
    coskpsqr = coskp * coskp;
    
    //double iprime1, qprime1, uprime1;
    //double cosphipol2 = coskp * coskp, sinphipol = sinkp, cosphi = coskp, sinphipol2 = sinkp * sinkp, cosphipol = coskp;
    // cospsi = coskp; sinpsi = sinkp;

    // spatial components of photon wave vector after scattering.
    std::array<double, 3> kmu_p = {sinkt * coskp, sinkt * sinkp, coskt};
    // epsilon_para; see Figure 5 of Schnittman & Krolik 2013;
    //std::array<double, 3> epsilon_para = {coskp, sinkp, 0.};
    // right-hand coordinate system; epsilon_perp = kmu x epsilon_para (the spatial components only)
    std::array<double, 3> epsilon_perp = {sinkp, -1. * coskp, 0.};
    // epsilon_para_prime, unit vector in the ki-kf plane, and is parallel with ki; note the handedness
    std::array<double, 3> epsilon_para_p = {coskt * coskp, coskt * sinkp, -1. * sinkt};

    //iprime = 0.5 * (1. - poldeg) * (1. + coskt * coskt) + poldeg * (1. - sinkt * sinkt * cosphipol2);
    //qprime = poldeg * (coskt * coskt * cosphipol2 - sinphipol2) - (1. - poldeg) * sinkt * sinkt / 2.;
    //uprime = coskt * poldeg * 2. * sinphipol * cosphipol;
    
    //iprime = 0.5 * (1. + cosktsqr - poldeg * cos2phi * sinktsqr);
    //qprime = 0.5 * (poldeg * (1. + cosktsqr) - sinktsqr);
    //uprime = coskt * poldeg * sin2phi;
    
    if (!KN) {
        // Rayleigh phase matrix
        double sin2phi = 2. * sinkp * coskp;
        iprime = 0.5 * (1. - poldeg) * (1. + cosktsqr) + poldeg * (sinkpsqr + coskpsqr * cosktsqr);
        qprime = -0.5 * (1. - poldeg) * sinktsqr + poldeg * (coskpsqr * cosktsqr - sinkpsqr);
        uprime = poldeg * coskt * sin2phi;
    } else {
        // The following lines for evaluting Q' and U' follow Connors+1980, Eq. II.(4)
        // As z -> 1, i(kn) -> 2i, q(kn) -> -2q (!!!), u(kn) -> 2u
        double cos2phi = 2. * coskpsqr - 1., sin2phi = 2. * sinkp * coskp;
        iprime = z + 1. / z - sinktsqr - sinktsqr * poldeg * cos2phi;
        qprime = (1. + cosktsqr) * poldeg * cos2phi - sinktsqr;
        uprime = 2. * coskt * poldeg * sin2phi;
    }
    
//     std::cout << "iprime = " << iprime << ", iprime1/2 = " << iprime1 / 2. << std::endl;
//     std::cout << "qprime = " << qprime << ", qprime1/2 = " << qprime1 / 2. << std::endl;
//     std::cout << "uprime = " << uprime << ", uprime1/2 = " << uprime1 / 2. << std::endl;
    
    polflux = std::sqrt(qprime * qprime + uprime * uprime);
    poldeg_p = polflux / iprime;
    polang_p = 0.5 * std::acos(qprime / polflux);
    if (uprime < 0.) {
        polang_p = M_PI - polang_p;
    }
    cosphipol_p = std::cos(polang_p);
    sinphipol_p = std::sin(polang_p);

    // photon polarization vector after scattering
    std::array<double, 3> fmu_p =
        {cosphipol_p * epsilon_para_p[0] + sinphipol_p * epsilon_perp[0],
        cosphipol_p * epsilon_para_p[1] + sinphipol_p * epsilon_perp[1],
        cosphipol_p * epsilon_para_p[2] + sinphipol_p * epsilon_perp[2]};
    //std::cout << "Inside scatter_restframe: f * f = " << arrdotprod(3, fmu_p.data(), fmu_p.data()) << std::endl;

    // we here assume thompson scattering, such that k[0] = k^\prime[0]
    // Then we transform photon wave and polarization vectors all the way back to fluid frame, step-by-step
    backward_rotation(pf, yaxis.data(), pk, kmu_p.data(), pk_e_p);
    backward_rotation(pf, yaxis.data(), pk, fmu_p.data(), pf_e_p);
    
    //std::cout << "Inside scatter_restframe: f * f = " << arrdotprod(3, pf_e_p, pf_e_p) << std::endl;
}

/*!
\brief Sample scattering angle assuming Thompson scattering cross section,
and calculates the photon energy, photon momentum, polarization vector of the scattered photon.
In this function we assume that the electron are *stationary*.

Suffix of vector names: _p: scattered photon.
@param pol if this option is off, we assume that the incoming photon is unpolarised.
@param poldeg polarization degree of the incoming photon
@param en dimensionless photon energy of the incoming photon 
@param kmu incoming photon wave-vector in the fluid frame
@param fmu photon polarization vector in the fluid frame
@param en_p dimensionless photon energy after scattering
@param poldeg_p polarization degree after scatter (note that polarization degree is invariant)
@param kmu_p photon wave vector after scattering
@param fmu_p photon polarization vector after scattering
@param gen Mersenne-Twister pseudo-random number generator
*/
void thompson_scatter(const bool pol, const double en, const double poldeg,
    std::array<double, 4> & kmu, std::array<double, 4> & fmu, double & en_p, double & poldeg_p,
    std::array<double, 4> & kmu_p, std::array<double, 4> & fmu_p, std::mt19937_64 & gen) {
    
    double kphi, coskt, kfdot, ffdot;

    // dummy vector for x-axis
    std::array<double, 3> xvec;

    std::uniform_real_distribution<> dis(0., 1.);

    if (pol) {
        sample_thompsonx(gen, poldeg, coskt, kphi);
    } else {
        sample_thompsonx(gen, 0., coskt, kphi);
    }
        
    en_p = en / (1. + en * (1. - coskt));
    
    if (!pol) {
        dummy_vec(kmu.data() + 1, xvec.data());
        scatter_restframe(false, kmu.data() + 1, xvec.data(), coskt, kphi, 0., poldeg_p, 1., kmu_p.data() + 1, fmu_p.data() + 1);
        poldeg_p = 0.;
    } else {
        scatter_restframe(false, kmu.data() + 1, fmu.data() + 1, coskt, kphi, poldeg, poldeg_p, 1., kmu_p.data() + 1, fmu_p.data() + 1);
    }
    
    // now we correct fmu the keep fmu and kmu are perpendicular
    kfdot = arrdotprod(3, kmu_p.data() + 1, fmu_p.data() + 1);
    ffdot = std::sqrt(arrdotprod(3, fmu_p.data() + 1, fmu_p.data() + 1));
    //std::cout << "k * f = " << kfdot << std::endl;
    //std::cout << "f * f = " << ffdot << std::endl;
    for (int i = 1; i <=3; ++i)
        fmu_p[i] -= kfdot * kmu_p[i];
    // normalize fmu
    ffdot = std::sqrt(arrdotprod(3, fmu_p.data() + 1, fmu_p.data() + 1));
    for (int i = 1; i <=3; ++i)
        fmu_p[i] /= ffdot;

    // To keep kmu as a null vector, set the time component to be 1;
    kmu_p[0] = 1.;
    fmu_p[0] = 0.;
}

/*!
\brief Sample scattering angle and calculates the photon energy, photon momentum, polarization vector of one photon after scatteirng. 
The electron momentum is given.

All the input and output are measured in the rest frame of the fluid.
Suffix of vector names inside this function: 
    - _p: after scattering
    - _f: fluid frame
    - _e: electron frame
    - _s: scattering frame
    
@param recoil whether to include modification of photon energy due to recoil effect
@param pol if off, we assume that the incoming photon is unpolarised
@param KN true - Klein-Nishina X-section; else Thompson. This option will override `recoil`: if we assume
    Klein-Nishina cross section, we always take into account recoil effect
@param en dimensionless photon energy before scattering, measured in the fluid frame
@param poldeg polarization angle before scattering
@param pe pe[0] - electron Lorentz factor; pe[1:3] - unit vector along electron 3-momentum
@param kmu_f photon wave-vector in the fluid frame
@param fmu_f photon polarization vector in the fluid frame
@param en_p dimensionless photon energy after scattering, measured in the fluid frame
@param poldeg_p polarization degree after scatter (note that polarization degree is invariant)
@param kmu_f_p photon wave vector after scatter, in the fluid frame
@param fmu_f_p photon polarization vector after scattering, in the fluid frame
@param gen Mersenne-Twister pseudo-random number generator
*/
void scatter(const bool pol, const bool recoil, const bool KN, const double en, const double poldeg, const std::array<double, 4> & pe,
    std::array<double, 4> & kmu_f, std::array<double, 4> & fmu_f, double & en_p, double & poldeg_p,
    std::array<double, 4> & kmu_f_p, std::array<double, 4> & fmu_f_p, std::mt19937_64 & gen) {
    
    if (!recoil && KN) {
        std::cerr << "no-recoil option not valid in KN" << std::endl;
        return;
    }
    //beta = v_e / c;
    double beta = std::sqrt(1. - 1. / pe[0] / pe[0]);
    //std::cout << "beta = " << beta << std::endl;
    // direction of scattered light; should be drawn from scattering cross-section
    double kphi, coskt, en_e, kfdot, ffdot;
    // photon wave and polarization vectors
    double enboost, knorm;
    std::array<double, 4> kmu_e, kmu_e_norm, fmu_e;
    double fk0ratio;
    
    // post scattering wave and polarization vector
    std::array<double, 4> kmu_e_p, fmu_e_p;
    // dummy vector for x-axis
    std::array<double, 3> xvec;

    std::uniform_real_distribution<> dis(0., 1.);

    enboost = loboost(pe[1], pe[2], pe[3], beta, pe[0], kmu_f, kmu_e);
    // Lorentz boost
    en_e = en * enboost;
    knorm = 1. / kmu_e[0];
    // get kmu direction vector; note that kmu is a null vector
    std::transform(kmu_e.cbegin(), kmu_e.cend(), kmu_e_norm.begin(), std::bind2nd(std::multiplies<double>(), knorm));

    //std::cout << "Inside scatter(), f_f * f_f = " << arrdotprod(3, fmu_f.data() + 1, fmu_f.data() + 1) << std::endl;
    if (pol) {
        loboost(pe[1], pe[2], pe[3], beta, pe[0], fmu_f, fmu_e);
        // adjust fmu that that fmu[0] = 0; remember that we are left with freedom to add fmu with factors of kmu
        fk0ratio = fmu_e[0] / kmu_e[0];
        for (int i = 0; i < 4; ++i)
            fmu_e[i] -= (fk0ratio * kmu_e[i]);
    }
    //std::cout << "Inside scatter(), f_e * f_e = " << arrdotprod(3, fmu_e.data() + 1, fmu_e.data() + 1) << std::endl;
    //std::cout << "Inside scatter(), f_e * kmu = " << arrdotprod(3, kmu_e.data() + 1, fmu_e.data() + 1) << std::endl;

    if (KN) {
        // sample KN X-section
        en_p = sample_KN_xp(en_e, gen);
        coskt = 1. + 1. / en_e - 1. / en_p;
        if (pol) {
            kphi = sample_kn_phi(gen, poldeg, coskt, en_p / en_e);
        } else {
            kphi = dis(gen) * 2. * M_PI;
        }
    } else {
        if (pol) {
            sample_thompsonx(gen, poldeg, coskt, kphi);
        } else {
            sample_thompsonx(gen, 0., coskt, kphi);
        }
        // recoil
        en_p = recoil ? en_e / (1. + en_e * (1. - coskt)) : en_e;
    }
    
    if (!pol) {
        dummy_vec(kmu_e_norm.data() + 1, xvec.data());
        scatter_restframe(KN, kmu_e_norm.data() + 1, xvec.data(), coskt, kphi, 0., 
            poldeg_p, en_p / en_e, kmu_e_p.data() + 1, fmu_e_p.data() + 1);
    } else {
        scatter_restframe(KN, kmu_e_norm.data() + 1, fmu_e.data() + 1, coskt, kphi, poldeg,
            poldeg_p, en_p / en_e, kmu_e_p.data() + 1, fmu_e_p.data() + 1);
    }
    
    if (pol) {
        // now we correct fmu the keep fmu and kmu are perpendicular
        kfdot = arrdotprod(3, kmu_e_p.data() + 1, fmu_e_p.data() + 1);
        for (int i = 1; i <=3; ++i)
            fmu_e_p[i] -= kfdot * kmu_e_p[i];
        // normalize fmu
        ffdot = std::sqrt(arrdotprod(3, fmu_e_p.data() + 1, fmu_e_p.data() + 1));
        for (int i = 1; i <=3; ++i)
            fmu_e_p[i] /= ffdot;
    }

    // To keep kmu as a null vector, set the time component to be 1;
    kmu_e_p[0] = 1.;
    fmu_e_p[0] = 0.;

    enboost = loboost(-1. * pe[1], -1. * pe[2], -1. * pe[3], beta, pe[0], kmu_e_p, kmu_f_p);
    // get kmu direction vector; note that kmu is a null vector
    knorm = 1. / kmu_f_p[0];
    std::transform(kmu_f_p.cbegin(), kmu_f_p.cend(), kmu_f_p.begin(), std::bind2nd(std::multiplies<double>(), knorm));

    en_p *= enboost;
    if (pol) {
        loboost(-1. * pe[1], -1. * pe[2], -1. * pe[3], beta, pe[0], fmu_e_p, fmu_f_p);

        fk0ratio = fmu_f_p[0] / kmu_f_p[0];
        for (int i = 0; i < 4; ++i)
            fmu_f_p[i] -= (fk0ratio * kmu_f_p[i]);
        //std::cout << "Inside scatter(), f * f = " << arrdotprod(4, fmu_f_p.data(), fmu_f_p.data()) << std::endl;
    }
}

/*!
\brief Propagate a superphoton inside a corona.
@param brem_self_abs whether to include bremsstrahlung self-absorption
@param pol whether to calculate polarization
@param a BH spin
@param rin inner edge radius of the disc on the equatorial plane
@param rg \f$GM/c^2\f$ in \f$[cm]\f$
@param mfp scattering mean free path
@param te dimensionless electron temperature
@param weightnorm photon weight
@param qfac0 normalised Stokes parameter q at infinity if the photon finally escape the corona without any interaction
@param ufac0 normalised Stokes parameter u at infinity if the photon finally escape the corona without any interaction
@param muinf_n0scatter photon inclination at infinity if the photon finally escape the corona without any interaction
@param dr maximum step size in \f$r\f$
@param sp \ref superphoton object; superphoton
@param trid \ref tridgeo object; corona object
@param ofiles see \ref register_sp_tridsca
@param nscafile see \ref register_sp_tridsca
*/
void tridsca_sp(const bool pol, const bool KN, const double a, 
    const double rin, electron_dist & edist, const double weightnorm, 
    const double qfac0, const double ufac0, const double muinf_n0scatter, const double dr_ratio, superphoton sp, 
    const tridgeo & trid, std::vector<std::vector<double>> & vecarr_inf, std::vector<int> & nsca_inf, 
    std::vector<std::vector<double>> & vecarr_disc, std::vector<int> & nsca_disc, std::mt19937_64 & gen,
	const bool ns, const double rns_rg) {
//  std::vector<std::ofstream> & ofiles, std::ofstream & nscafile) {
    
    //std::cout << "len(ofiles) = " << ofiles.size() << std::endl;
    //std::random_device rd;
    //std::mt19937_64 gen(rd());

    double inte, dinte, re, mue, rprev, muprev, kre, kmue, taue,
        xie, krsignnow, kmusignnow, gfacnow, gfacprev,
        lambdaprev, lambdanow, qweight, uweight, qfac, ufac, ktheta, dr;
    double bias, maxtau, weight1, en1, ephi, mfp;
    bool usebias;
    double dlambda_r, dlambda_mu, dlambda, dlambda_mue, dlambda_re;
    double rndnum, dtau, interact_prob;
    double knorm, poldeg_p, emu, ennow;
    std::array<double, 4> kmunow, kmu_on, umunow, fmu_p, kmu_p, pe, fmu, kmu_temp, fmu_on;
    int nrr2, nrr3, nmuplus, nmuminus;
    double temp, muinf, a2, rfunc, rfunc0, rfunc2;
	double alpha_now, alpha_prev;
    superphoton sp1;
    sim5complex sim5wp;
    sim5complex * sim5wp_pointer;
    //auto dsize = sizeof(double);
    //auto isize = sizeof(int);
    //thread_local std::mt19937_64 gen(std::random_device {}());
    a2 = a * a;
    //std::cout << "muinf_n0scatter = " << muinf_n0scatter << std::endl;

    // ktarr: array of -k_(t) = U^\mu dot k^\mu, which is the time component of
    // k^\mu measured by local tetrad
    // loboostarr: arr of Lorentz boost to energy
    // inloboostarr: arr of inverse Lorentz boost to energy
    // scattermuarr: array of cos(catter polar angle)
    
    if (sp.weight < MINWEIGHT) {
        //std::cout << "sp.weight = " << sp.weight << std::endl;
        return;
    }

    //mfp = trid.mean_free_path(tau, rnow, munow);
    //maxtau = trid.length() / mfp;
    maxtau = trid.length();
    bias = 1. / maxtau;
    if (1. / maxtau <= 1.) {
        usebias = false;
        bias = 1.;
    } else {
        usebias = true;
    }
    //std::cout << "usebias = " << usebias << std::endl;
    //dinte = 1e-4;
    // dinte will be chosen adaptively

    sim5tetrad tet;
    sim5metric metnow;

    krsignnow = (sp.kmu[1] >= 0.) ? 1. : -1.;
    kmusignnow = (sp.kmu[2] >= 0.) ? 1. : -1.;

    re = sp.rnow;
    mue = sp.munow;
    inte = 0.;
    kre = sp.kmu[1];
    kmue = sp.kmu[2];
    kmunow = sp.kmu;
    rprev = 0.;
    muprev = 0.;

    kerr_metric(a, re, mue, &metnow);
    trid.calumu(sp.rnow, sp.munow, umunow, metnow);
    gfacnow = -1. * dotprod(kmunow.data(), umunow.data(), &metnow);
    //std::cout << "kmu * kmu = " << dotprod(sp.kmu.data(), sp.kmu.data(), &metnow) << std::endl;
    
    //prod = dotprod(sp.kmu.data(), sp.kmu.data(), &metnow);
    /*
    if (std::abs(sp.kmu[0]) + std::abs(sp.kmu[1]) + std::abs(sp.kmu[2]) + std::abs(sp.kmu[3]) < 1e-8 || std::abs(prod) > 1e-8) {
        std::cout << "k * k = " << prod << std::endl;
        std::cout << "Inside tridsca_sp();" << std::endl;
        std::cout << "rnow = " << sp.rnow << ", munow = " << sp.munow << std::endl;
        std::cout << "sp.kmu = " << sp.kmu << std::endl;
        std::cout << "sp.nscatter = " << sp.nscatter << std::endl;
        return;
    }
    */

    lambdanow = 0.;

    if (trid.inside(sp.rnow, sp.munow)) {
        mfp = trid.mean_free_path(sp.rnow, sp.munow);
        alpha_now = 1. / mfp;
        if (KN) {
            alpha_now *= edist.hotcross(sp.rnow, sp.munow, sp.en * gfacnow);
        }

        dr = dr_ratio * mfp;
    } else {
          dr = sp.rnow * dr_ratio;
          alpha_now = 0.;
    }

    geoinf geo;
    try {
        geo = geoinf(a, re, mue, kmunow);
    } 
    catch (const char * error) {
        std::cerr << "Inside scatter.cpp, tridsca_sp()" << std::endl;
        std::cerr << "geo = geoinf(a, re, mue, kmunow)" << std::endl;
        std::cerr << error << std::endl;
        return;
    }
    if (geo.l != geo.l || geo.q != geo.q) {
        std::cout << "Inside tridsca_sp(); Checking geo.l:" << std::endl;
        std::cout << "re = " << re << ", mue = " << mue << std::endl;
        std::cout << "kmunow = " << kmunow << std::endl;
        std::cout << "sp.nscatter = " << sp.nscatter << std::endl;
        return;
    }
    //std::cout << "starting caltaue():" << std::endl;
    try {
        taue = geo.caltaue(re);
        xie = geo.calximu(mue);
        //std::cout << "starting caldlambda_re():" << std::endl;
        dlambda_re = geo.caldlambda_re(re);
        dlambda_mue = geo.caldlambda_muplus(mue);
    } 
    catch (const char * error) {
        std::cerr << error << std::endl;
        return;
    }
    
    std::array<double, 5> rfunc_coeff = {1., 0., a2 - geo.l * geo.l - geo.q, 2. * (geo.q + (geo.l -a) * (geo.l - a)), -1. * a2 * geo.q};
    
    std::uniform_real_distribution<> dis(0., 1.);

	double realrh;
    if (ns) {
        realrh = rns_rg;
    } else {
        realrh = geo.rh;
    }
    
    while(true) {
        if (trid.inside(sp.rnow, sp.munow)) {
            mfp = 1. / alpha_now;
            dr = dr_ratio * mfp;
        } else {
            dr = sp.rnow * dr_ratio;
        }
        //std::cout << "dr = " << dr << std::endl;
        
        rfunc = polyval(rfunc_coeff, sp.rnow);
        rfunc0 = polyval(rfunc_coeff, sp.rnow + dr);
        rfunc2 = polyval(rfunc_coeff, sp.rnow - dr);
        if (rfunc0 > rfunc) rfunc = rfunc0;
        if (rfunc2 > rfunc) rfunc = rfunc2;
        
        dinte = dr / std::sqrt(rfunc);
        
        //if (dinte > 1e-3) dinte = 1e-3;
        
        inte += dinte;
        muprev = sp.munow;
        rprev = sp.rnow;
        lambdaprev = lambdanow;
        gfacprev = gfacnow;
        alpha_prev = alpha_now;

        //std::cout << "starting solver_inte():" << std::endl;
        try {
            sp.rnow = geo.solver_inte(inte, re, kre, taue, krsignnow, nrr2, nrr3);
        }
        catch(const char * error) {
            std::cout << error << std::endl;
            std::cout << "re = " << re << ", kre = " << kre << std::endl;
            std::cout << "rnow = " << sp.rnow << ", munow = " << sp.munow << std::endl;
            std::cout << "l = " << geo.l << ", q = " << geo.q << std::endl;
            std::cout << "nrr = " << geo.nrr << std::endl;
            std::cout << "rr1 = " << geo.rr1 << ", rr2 = " << geo.rr2 << ", rr3 = " << geo.rr3 << ", rr4 = " << geo.rr4 << std::endl;
            std::cout << "dinte = " << dinte << std::endl;
            std::cout << "rfunc = " << rfunc << std::endl;
            std::cout << "rfunc(r - 0.01) = " << polyval(rfunc_coeff, sp.rnow - dr) << std::endl;
        }
        
        //std::cout << "rnow - rprev = " << sp.rnow - rprev << std::endl;

        //std::cout << "starting solvemu_inte():" << std::endl;
        try {
            sp.munow = geo.solvemu_inte(inte, mue, kmue, xie, kmusignnow, nmuminus, nmuplus);
        }
        catch(const char * error) {
            std::cout << error << std::endl;
            std::cout << "rnow = " << sp.rnow << ", munow = " << sp.munow << std::endl;
            std::cout << "l = " << geo.l << ", q = " << geo.q << std::endl;
            std::cout << "muplus = " << geo.muplus << ", muminus = " << geo.muminus << std::endl;
            std::cout << "taumu0 = " << geo.taumu0 << std::endl;
            std::cout << "dinte = " << dinte << std::endl;
        }
        
        if (sp.rnow <= 1.01 * realrh) {//|| sp.munow <= 0. || (geo.q > 0. && nmuminus > 0)) {
            break;
        }
        if (sp.munow * muprev < 0. && sp.rnow >= rin) {
            nsca_disc.push_back(sp.nscatter);
            vecarr_disc[0].push_back(sp.en);
            vecarr_disc[1].push_back(sp.weight * weightnorm);
            vecarr_disc[2].push_back(0.5 * (sp.rnow + rprev));
            vecarr_disc[3].push_back(geo.l);
            vecarr_disc[4].push_back(geo.q);
            vecarr_disc[5].push_back(krsignnow);
            vecarr_disc[6].push_back(kmusignnow);
            
            if (pol) {
                vecarr_disc[7].push_back(sp.poldeg);
                vecarr_disc[8].push_back(sp.wp.real());
                vecarr_disc[9].push_back(sp.wp.imag());
            }
            break;
        }
        
        if (sp.rnow != sp.rnow || sp.munow != sp.munow) {
            std::cout << "rprev = " << rprev << ", muprev = " << muprev << std::endl;
            std::cout << "l = " << geo.l << ", q = " << geo.q << std::endl;
            std::cout << "muplus = " << geo.muplus << ", muminus = " << geo.muminus << std::endl;
            std::cout << "taumu0 = " << geo.taumu0 << std::endl;
            std::cout << "dinte = " << dinte << std::endl;
            std::cout << "sp.rnow or sp.munow NaN!" << std::endl;
            break;
        }

        //std::cout << "sp.rnow = " << sp.rnow << ", munow = " << sp.munow << std::endl;
        if (!trid.inside(sp.rnow, sp.munow) && !geo.ispossible(trid, sp.rnow, sp.munow, krsignnow, kmusignnow)) {
            // find out muinf
            temp = sp.weight * weightnorm;
            
            muinf = geo.calfate(sp.rnow, sp.munow, krsignnow, kmusignnow, rin);
            //std::cout << "kmuobssign = " << geo.kmuobssign << std::endl;
            if (pol && geo.status == 1) {
                cal_pol_inf(muinf, geo.kmuobssign, geo.a, geo.l, geo.q, sp.poldeg, sp.wp, qfac, ufac);
                if (qfac != qfac || ufac != ufac) {
                    std::cout << "Inside scatter()" << std::endl;
                    std::cout << "sp.wp = " << sp.wp << std::endl;
                    std::cout << "muinf = " << muinf << std::endl;
                    std::cout << "ktheta = " << ktheta << std::endl;
                    std::cout << "geo.l = " << geo.l << std::endl;
                    std::cout << "geo.q = " << geo.q << std::endl;
                    std::cout << "sp.poldeg = " << sp.poldeg << std::endl;
                    std::cout << "qfac or ufac is NaN" << std::endl;
                }
                qweight = qfac * temp;
                uweight = ufac * temp;
            } else {
                qweight = 0.;
                uweight = 0.;
            }
            
            if (geo.status == 1) {
                nsca_inf.push_back(sp.nscatter);
                vecarr_inf[0].push_back(sp.en);
                vecarr_inf[1].push_back(temp);
                vecarr_inf[2].push_back(muinf);
                vecarr_inf[3].push_back(geo.l);
                vecarr_inf[4].push_back(geo.q);
                vecarr_inf[5].push_back(geo.kmuobssign);
                
                if (pol) {
                    vecarr_inf[6].push_back(qweight);
                    vecarr_inf[7].push_back(uweight);
                }
            } else if (geo.status == 2) {
                nsca_disc.push_back(sp.nscatter);
                vecarr_disc[0].push_back(sp.en);
                vecarr_disc[1].push_back(temp);
                vecarr_disc[2].push_back(muinf);
                vecarr_disc[3].push_back(geo.l);
                vecarr_disc[4].push_back(geo.q);
                vecarr_disc[5].push_back(geo.krobssign);
                vecarr_disc[6].push_back(geo.kmuobssign);
                
                if (pol) {
                    vecarr_disc[7].push_back(sp.poldeg);
                    vecarr_disc[8].push_back(sp.wp.real());
                    vecarr_disc[9].push_back(sp.wp.imag());
                }
            }
            return;
        }
                
        if (trid.inside(sp.rnow, sp.munow) && trid.inside(rprev, muprev)) {
            try {
                dlambda_r = geo.caldlambda_rall(re, sp.rnow, kre, nrr2, nrr3, dlambda_re, inte);
                dlambda_mu = geo.caldlambda_muall(mue, sp.munow, kmue, nmuminus, nmuplus, dlambda_mue, inte);
                lambdanow = dlambda_r + dlambda_mu;
            }
            catch (const char * error) {
                std::cerr << error << std::endl;
                break;
            }

            dlambda = lambdanow - lambdaprev;
            //std::cout << "dlambda = " << dlambda << std::endl;

            kerr_metric(a, sp.rnow, sp.munow, &metnow);
            photon_momentum(a, sp.rnow, sp.munow, geo.l, geo.q, krsignnow, kmusignnow, kmunow.data());
            if (kmunow != kmunow) return;

            trid.calumu(sp.rnow, sp.munow, umunow, metnow);
            gfacnow = -1. * dotprod(kmunow.data(), umunow.data(), &metnow);
            
            /*
            std::cout << "mu = " << sp.munow << ", umu = " << umunow << std::endl;
            std::array<double, 4> tempumu;
            sim5tetrad temptetrad;
            tetrad_zamo(&metnow, &temptetrad);
            bl2on(umunow.data(), tempumu.data(), &temptetrad);
            std::cout << "umuon = " << tempumu << std::endl;
            */
            mfp = trid.mean_free_path(sp.rnow, sp.munow);
            alpha_now = 1. / mfp;
            if (KN) {
                alpha_now *= edist.hotcross(sp.rnow, sp.munow, sp.en * gfacnow);
            }
            
            dtau = 0.5 * (gfacprev * alpha_prev + gfacnow * alpha_now) * dlambda;

            // scattering optical depth in the fluid frame
            interact_prob = bias * dtau;
            //std::cout << "scatterprob = " << scatterprob << std::endl;
            rndnum = dis(gen);

            if (rndnum <= interact_prob) {
                if (usebias) {
                    weight1 = sp.weight / bias;
                    sp.weight *= (1. - 1. / bias);
                } else {
                    weight1 = sp.weight;
                }

                trid.caltet(sp.rnow, sp.munow, tet);
                
                if (pol) {
                    sim5wp_pointer = reinterpret_cast<sim5complex *>(&sp.wp);
                    polarization_vector(kmunow.data(), *sim5wp_pointer, &metnow, fmu.data());
                }

                bl2on(kmunow.data(), kmu_on.data(), &tet);
                knorm = 1. / kmu_on[0];
                std::transform(kmu_on.begin(), kmu_on.end(), kmu_on.begin(),
                    std::bind2nd(std::multiplies<double>(), knorm));

                // gravitational redshift
                ennow = sp.en * gfacnow;
                
                // sample electron
                //sample_electron_KN(te, ennow, gen, pe[0], emu, ephi);
                edist.sample_mu(sp.rnow, sp.munow, ennow, pe[0], emu, gen);
                ephi = 2. * M_PI * dis(gen);
                calpe(emu, ephi, kmu_on.data() + 1, pe);
                
                if (pol) {
                    bl2on(fmu.data(), fmu_on.data(), &tet);
                }

                scatter(pol, true, KN, ennow, sp.poldeg, pe, kmu_on, fmu_on, en1, poldeg_p, kmu_p, fmu_p, gen);

                on2bl(kmu_p.data(), kmu_temp.data(), &tet);
                knorm = -1. / (metnow.g00 * kmu_temp[0] + metnow.g03 * kmu_temp[3]);
                if (knorm < 0.) {
                    return;
                }
                std::transform(kmu_temp.begin(), kmu_temp.end(), kmu_temp.begin(),
                    std::bind2nd(std::multiplies<double>(), knorm));

                double gfactemp = -1. * dotprod(kmu_temp.data(), umunow.data(), &metnow);
                
                if (pol) {
                    on2bl(fmu_p.data(), fmu.data(), &tet);
                    sim5wp = photon_wp_const(kmu_temp.data(), fmu.data(), &metnow);
                    sp1.poldeg = poldeg_p;
                    sp1.wp.real(creal(sim5wp));
                    sp1.wp.imag(cimag(sim5wp));
                }

                sp1.rnow = sp.rnow;
                sp1.munow = sp.munow;
                sp1.kmu = kmu_temp;
                sp1.weight = weight1;
                sp1.nscatter = sp.nscatter + 1;
                sp1.en = en1 / gfactemp;
                //std::cout << "After scattering; sp.poldeg = " << sp1.poldeg << std::endl;
                tridsca_sp(pol, KN, a, rin, edist, weightnorm, qfac0, ufac0, muinf_n0scatter, 
                    dr_ratio, sp1, trid, vecarr_inf, nsca_inf, vecarr_disc, nsca_disc, gen, ns, rns_rg);
                if (!usebias)
                    return;
            }
        }
    }
}
