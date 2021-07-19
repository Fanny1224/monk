//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

//! \file tridgeo.cpp
#include <cmath>
#include <array>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "tridgeo.h"
#include "sim5lib.h"
#include "utils.h"
#include "kerr.h"
#include "const.h"

//const double EulerConstant = std::exp(1.0);
//const double efac = 1. - 1. / EulerConstant;

/*!
\brief Constructor.
@param aa BH spin
@param ss corona radius
@param hminn corona minimum height along the rotation axis
@param hmaxx corona maximum height along the rotation axis
 */
zamoslab3d::zamoslab3d(const double aa, const double ss, const double hminn, const double hmaxx, const double tauu) {
    name = "zamoslab3d";
    a = aa;
    s = ss;
    tau = tauu;
    
    hmin = hminn;
    hmax = hmaxx;
    
    rmax = std::sqrt(s * s + hmax * hmax);
    rmin = hmin;
    
    mumax = 1.;
    mumin = hmin / std::sqrt(s * s + hmin * hmin);
    magnetised = false;
}

/*!
\brief Tell if the photon is inside the corona.
*/
bool zamoslab3d::inside(double rnow, double munow) const {
    double hnow = rnow * munow, snow = std::sqrt(rnow * rnow - hnow * hnow);
    return (hnow >= hmin && hnow <= hmax && snow <= s);
}

/*!
\brief Returns tetrad.
*/
void zamoslab3d::caltet(double r, double mu, sim5tetrad & tet) const {
    sim5metric met;
    kerr_metric(a, r, mu, &met);
    tetrad_zamo(&met, &tet);
}

/*!
\brief Calculates four velocity.
*/
void zamoslab3d::calumu(const double r, const double mu, std::array<double, 4> & umu, sim5metric & met) const {
    fourvelocity_zamo(&met, umu.data());
}

void zamoslab3d::genpos(std::mt19937_64 & gen, double & rnow, double & munow) const {
    double hnow, snow;
    std::uniform_real_distribution<> dis(0., 1.);
    hnow = dis(gen) * (hmax - hmin) + hmin;
    snow = s * std::sqrt(dis(gen));
    rnow = std::sqrt(hnow * hnow + snow * snow);
    munow = hnow / rnow;
}

double zamoslab3d::mean_free_path(const double r, const double mu) const{
    return 0.5 * (hmax - hmin) / tau;
}

double zamoslab3d::ne_sigmat(const double r, const double mu) const{
    return 2. * tau / (hmax - hmin);
}

zamosandwich::zamosandwich(const double aa, const double ss, const double sminn, const double hminn, const double hmaxx, const double tauu) {
    name = "zamosandwich";
    a = aa;
    tau = tauu;
    hmin = hminn;
    hmax = hmaxx;
    smin = sminn;
    s = ss;

    rmax = std::sqrt(s * s + hmax * hmax);
    rmin = std::sqrt(smin * smin + hmax * hmax);

    mumax = hmax / std::sqrt(smin * smin + hmax * hmax);
    mumin = hmin / std::sqrt(s * s + hmin * hmin);
    magnetised = false;

}

bool zamosandwich::inside(const double rnow, const double munow) const {
    double hnow = rnow * munow, snow = std::sqrt(rnow * rnow - hnow * hnow);
    return (hnow >= hmin && hnow <= hmax && snow <= s && snow >= smin);
}

void zamosandwich::calumu(const double r, const double mu, std::array<double, 4> & umu, sim5metric & met) const {
    fourvelocity_zamo(&met, umu.data());
}

void zamosandwich::caltet(double r, double mu, sim5tetrad & tet) const {
    sim5metric met;
    kerr_metric(a, r, mu, &met);
    tetrad_zamo(&met, &tet);
}

kepsandwich::kepsandwich(const double aa, const double ss, const double sminn, const double hminn, const double hmaxx, const double tauu) {
    name = "kepsandwich";
    a = aa;
    tau = tauu;
    hmin = hminn;
    hmax = hmaxx;
    smin = sminn;
    s = ss;

    rmax = std::sqrt(s * s + hmax * hmax);
    rmin = std::sqrt(smin * smin + hmax * hmax);

    mumax = hmax / std::sqrt(smin * smin + hmax * hmax);
    mumin = hmin / std::sqrt(s * s + hmin * hmin);
    magnetised = false;

}

bool kepsandwich::inside(const double rnow, const double munow) const {
    double hnow = rnow * munow, snow = std::sqrt(rnow * rnow - hnow * hnow);
    return (hnow >= hmin && hnow <= hmax && snow <= s && snow >= smin);
}

void kepsandwich::calumu(const double r, const double mu, std::array<double, 4> & umu, sim5metric & met) const {
    double x, Omega, snow;
    snow = r * std::sqrt(1. - mu * mu);
    x = std::sqrt(snow);
    Omega = 1. / (x * x * x + a);
    fourvelocity_azimuthal(Omega, &met, umu.data());
}

void kepsandwich::caltet(double r, double mu, sim5tetrad & tet) const {
    double x, Omega, snow;
    sim5metric met;
    snow = r * std::sqrt(1. - mu * mu);
    x = std::sqrt(snow);
    Omega = 1. / (x * x * x + a);
    kerr_metric(a, r, mu, &met);
    tetrad_azimuthal(&met, Omega, &tet);
}

    
/*!
\brief Constructor.
@param aa BH spin
@param ss corona radius
@param hminn corona minimum height along the rotation axis
@param hmaxx corona maximum height along the rotation axis
*/
kepslab3d::kepslab3d(const double aa, const double ss, const double hminn, const double hmaxx, const double tauu) {
    name = "kepslab3d";
    a = aa;
    s = ss;
    hmin = hminn;
    hmax = hmaxx;
    tau = tauu;
    
    rmax = std::sqrt(s * s + hmax * hmax);
    rmin = hmin;
    
    mumax = 1.;
    mumin = hmin / std::sqrt(s * s + hmin * hmin);
    magnetised = false;
}

double kepslab3d::mean_free_path(const double, const double) const
{
    return 0.5 * (hmax - hmin) / tau;
}
double kepslab3d::ne_sigmat(const double r, const double mu) const{
    return 2. * tau / (hmax - hmin);
}

bool kepslab3d::inside(double rnow, double munow) const {
    double hnow = rnow * munow, snow = std::sqrt(rnow * rnow - hnow * hnow);
    return (hnow >= hmin && hnow <= hmax && snow <= s);
}

void kepslab3d::caltet(double r, double mu, sim5tetrad & tet) const {
    double x, Omega, snow;
    sim5metric met;
    snow = r * std::sqrt(1. - mu * mu);
    x = std::sqrt(snow);
    Omega = 1. / (x * x * x + a);
    kerr_metric(a, r, mu, &met);
    tetrad_azimuthal(&met, Omega, &tet);
}

void kepslab3d::calumu(const double r, const double mu, std::array<double, 4> & umu, sim5metric & met) const {
    double x, Omega, snow;
    snow = r * std::sqrt(1. - mu * mu);
    x = std::sqrt(snow);
    Omega = 1. / (x * x * x + a);
    fourvelocity_azimuthal(Omega, &met, umu.data());
}

void kepslab3d::genpos(std::mt19937_64 & gen, double & rnow, double & munow) const {
    double hnow, snow;
    std::uniform_real_distribution<> dis(0., 1.);
    hnow = dis(gen) * (hmax - hmin) + hmin;
    snow = s * std::sqrt(dis(gen));
    rnow = std::sqrt(hnow * hnow + snow * snow);
    munow = hnow / rnow;
}

/*!
\brief Constructor.
@param aa BH spin
@param hh height
@param RR radius
*/
zamosphere3d::zamosphere3d(const double aa, const double hh, const double RR, const double tauu) {
    name = "zamosphere3d";
    a = aa;
    h = hh;
    R = RR;
    tau = tauu;
    
    rmax = h + R;
    mumax = 1.;
    if (h > 1e-8) {
        rmin = h - R;
        rmax = h + R;
        mumax = 1.;
        mumin = R / h;
    } else {
        rmin = 0.;
        rmax = R;
        mumin = -1.;
    }
    magnetised = false;
}

bool zamosphere3d::inside(double rnow, double munow) const {
    double dissqr = rnow * rnow + h * h - 2. * rnow * h * munow;
    //dissqr = rnow * rnow * (1. - munow * munow) + (rnow * munow - h) * (rnow * munow - h);
    return (dissqr <= R * R);
}

double zamosphere3d::mean_free_path(const double r, const double mu) const {
    return R / tau;
}

double zamosphere3d::ne_sigmat(const double r, const double mu) const {
    return tau / R;
}

void zamosphere3d::caltet(double r, double mu, sim5tetrad & tet) const {
    sim5metric met;
    kerr_metric(a, r, mu, &met);
    tetrad_zamo(&met, &tet);
}

void zamosphere3d::calumu(const double r, const double mu, std::array<double, 4> & umu, sim5metric & met) const {
    if (mu < 1.) {
        fourvelocity_zamo(&met, umu.data());
    } else {
        umu[0] = std::sqrt(-1. / met.g00);
        umu[1] = 0.;
        umu[2] = 0.;
        umu[3] = 0.;
    }
}

void zamosphere3d::genpos(std::mt19937_64 & gen, double & rnow, double & munow) const {
    double r0, mu0, hnow, snow;
    std::uniform_real_distribution<> dis(0., 1.);
    
    r0 = R * std::cbrt(dis(gen));
    mu0 = 2. * dis(gen) - 1.;
    hnow = h + r0 * mu0;
    snow = r0 * std::sqrt(1. - mu0 * mu0);
    rnow = std::sqrt(hnow * hnow + snow * snow);
    munow = hnow / rnow;
}

void zamosphere3d::write_brem_sp(const size_t ndim1, const size_t ndim2, std::ofstream & ofile) const {
}

/*! \brief Constructor.
@param aa BH spin
@param rr outer radius of the shell
@param rinn inner radius of the shell
 */
zamosphere3d_truncated::zamosphere3d_truncated(const double aa, const double rr, const double rinn, const double tauu) {
    name = "zamosphere3d_truncated";
    a = aa;
    R = rr;
    rin = rinn;
    tau = tauu;
    
    rmax = R;
    rmin = rin;
    
    mumax = 1.;
    mumin = -1.;
    magnetised = false;
}

bool zamosphere3d_truncated::inside(double rnow, double munow) const {
    //dissqr = rnow * rnow * (1. - munow * munow) + (rnow * munow - h) * (rnow * munow - h);
    return (rnow <= R && rnow >= rin);
}

void zamosphere3d_truncated::caltet(double r, double mu, sim5tetrad & tet) const {
    sim5metric met;
    kerr_metric(a, r, mu, &met);
    tetrad_zamo(&met, &tet);
}

void zamosphere3d_truncated::calumu(const double r, const double mu, std::array<double, 4> & umu, sim5metric & met) const {
    if (mu < 1.) {
        fourvelocity_zamo(&met, umu.data());
    } else {
        umu[0] = std::sqrt(-1. / met.g00);
        umu[1] = 0.;
        umu[2] = 0.;
        umu[3] = 0.;
    }
}

void zamosphere3d_truncated::genpos(std::mt19937_64 & gen, double & rnow, double & munow) const {
    bool found = false;
    double r0, mu0, hnow, snow;
    std::uniform_real_distribution<> dis(0., 1.);
    
    while (!found) {
        r0 = R * std::cbrt(dis(gen));
        if (r0 >= rin) 
            found = true;
    }
    mu0 = 2. * dis(gen) - 1.;
    hnow = r0 * mu0;
    snow = r0 * std::sqrt(1. - mu0 * mu0);
    rnow = std::sqrt(hnow * hnow + snow * snow);
    munow = hnow / rnow;
}

double zamosphere3d_truncated::mean_free_path(const double r, const double mu) const {
    return (R - rin) / tau;
}

double zamosphere3d_truncated::ne_sigmat(const double r, const double mu) const {
    return tau / (R - rin);
}

bool rotsphere3d::inside(double rnow, double munow) const {
    double dissqr = rnow * rnow + h * h - 2. * rnow * h * munow;
    //dissqr = rnow * rnow * (1. - munow * munow) + (rnow * munow - h) * (rnow * munow - h);
    return (dissqr <= R * R);
}

double rotsphere3d::mean_free_path(const double r, const double mu) const {
    return R / tau;
}

double rotsphere3d::ne_sigmat(const double r, const double mu) const {
    return tau / R;
}

/*!
\brief Object for corona rotating with given linear velocity
@param aa BH spin
@param hh height
@param RR radius
@param vv linear circular velocity 
*/
rotsphere3d::rotsphere3d(const double aa, const double hh, const double RR, const double vv, const double tauu) {
    name = "rotsphere3d";
    a = aa;
    h = hh;
    R = RR;
    tau = tauu;
    
    velocity = vv;
    gamma = 1. / std::sqrt(1. - velocity * velocity);
    rmax = h + R;
    mumax = 1.;
    if (h > 1e-8) {
        rmin = h - R;
        rmax = h + R;
        mumax = 1.;
        mumin = R / h;
    } else {
        rmin = 0.;
        rmax = R;
        mumin = -1.;
    }
    magnetised = false;
}

void rotsphere3d::caltet(double r, double mu, sim5tetrad & tet) const {
    double omega;
    sim5metric met;
    sim5tetrad zamotet;
    kerr_metric(a, r, mu, &met);
    std::array<double, 4> umu_on = {gamma, 0., 0., gamma * velocity}, umu_bl;
    if (mu < 1.) {
        tetrad_zamo(&met, &zamotet);
        on2bl(umu_on.data(), umu_bl.data(), &zamotet);
        omega = umu_bl[3] / umu_bl[0];
        tetrad_azimuthal(&met, omega, &tet);
    } else {
        tetrad_zamo(&met, &tet);
    }
}

void rotsphere3d::calumu(const double r, const double mu, std::array<double, 4> & umu, sim5metric & met) const {
    sim5tetrad zamotet;
    std::array<double, 4> umu_on = {gamma, 0., 0., gamma * velocity};
    if (mu < 1.) {
        tetrad_zamo(&met, &zamotet);
        on2bl(umu_on.data(), umu.data(), &zamotet);
    } else {
        fourvelocity_zamo(&met, umu.data());
    }
}

void rotsphere3d::genpos(std::mt19937_64 & gen, double & rnow, double & munow) const {
    double r0, mu0, hnow, snow;
    std::uniform_real_distribution<> dis(0., 1.);
    
    r0 = R * std::cbrt(dis(gen));
    mu0 = 2. * dis(gen) - 1.;
    hnow = h + r0 * mu0;
    snow = r0 * std::sqrt(1. - mu0 * mu0);
    rnow = std::sqrt(hnow * hnow + snow * snow);
    munow = hnow / rnow;
}

bool sphere3d_sr::inside(const std::array<double, 4> & xmu) const {
    double rnowsqr;
    rnowsqr = xmu[3] * xmu[3] + xmu[1] * xmu[1] + xmu[2] * xmu[2];
    return (rnowsqr <= radius * radius);
}


bool slab3d_sr::inside(const std::array<double, 4> & xmu) const {
    double snowsqr = xmu[2] * xmu[2] + xmu[1] * xmu[1];
    return (std::abs(xmu[3]) <= h && snowsqr <= s * s);
}

wedge3d::wedge3d(const double aa, const double sminn, const double smaxx, const double mumaxx, const double tauu) {
    double cos_opening_angle; // cosine of the opening angle; the opening angle is \f$\pi / 2 - \theta\f$, where \f$\theta\$ is one of the Boyer-Lindquist coordinates
    name = "wedge3d";
    a = aa;
    smin = sminn;
    smax = smaxx;
    tau = tauu;
    
    mumin = 0.;
    mumax = mumaxx;
    cos_opening_angle = std::sqrt(1. - mumax * mumax);
    tan_opening_angle = mumax / cos_opening_angle;
    
    rmin = sminn;
    rmax = smax / cos_opening_angle;
    magnetised = false;
}

bool wedge3d::inside(const double r, const double mu) const {
    bool insidenow = false;
    double s = r * std::sqrt(1. - mu * mu);
    insidenow = (s >= smin && s <= smax && mu >= 0. && mu <= mumax);
    return insidenow;
}

void wedge3d::caltet(const double r, const double mu, sim5tetrad & tet) const {
    double x, Omega, snow;
    sim5metric met;
    snow = r * std::sqrt(1. - mu * mu);
    x = std::sqrt(snow);
    Omega = 1. / (x * x * x + a);
    kerr_metric(a, r, mu, &met);
    tetrad_azimuthal(&met, Omega, &tet);
}

void wedge3d::calumu(const double r, const double mu, std::array<double, 4> & umu, sim5metric & met) const {
    double x, Omega, snow;
    snow = r * std::sqrt(1. - mu * mu);
    x = std::sqrt(snow);
    Omega = 1. / (x * x * x + a);
    fourvelocity_azimuthal(Omega, &met, umu.data());
}

//! This is a 
double wedge3d::mean_free_path(const double r, const double mu) const
{
    double snow = r * std::sqrt(1. - mu * mu);
    double Hnow = snow * tan_opening_angle;
    return (Hnow * std::exp(r * mu / Hnow) / tau);
}

double wedge3d::ne_sigmat(const double r, const double mu) const {
    double snow = r * std::sqrt(1. - mu * mu);
    double Hnow = snow * tan_opening_angle;
    return (tau / (Hnow * std::exp(r * mu / Hnow)));
}

void wedge3d::genpos(std::mt19937_64 & gen, double & rnow, double & munow) const {
    rnow = 0.;
    munow = 0.;
}

wedge3d_zamo_uniform::wedge3d_zamo_uniform(const double aa, const double sminn, const double smaxx, const double mumaxx, const double tauu) {
    double cos_opening_angle; // cosine of the opening angle; the opening angle is \f$\pi / 2 - \theta\f$, where \f$\theta\$ is one of the Boyer-Lindquist coordinates
    name = "wedge3d_zamo_uniform";
    a = aa;
    smin = sminn;
    smax = smaxx;
    tau = tauu;
    
    mumin = 0.;
    mumax = mumaxx;
    cos_opening_angle = std::sqrt(1. - mumax * mumax);
    tan_opening_angle = mumax / cos_opening_angle;
    
    rmin = sminn;
    rmax = smax / cos_opening_angle;
    magnetised = false;
}

bool wedge3d_zamo_uniform::inside(const double r, const double mu) const {
    double s = r * std::sqrt(1. - mu * mu);
    return (s >= smin && s <= smax && mu >= 0. && mu <= mumax);
}

void wedge3d_zamo_uniform::caltet(const double r, const double mu, sim5tetrad & tet) const {
    sim5metric met;
    kerr_metric(a, r, mu, &met);
    tetrad_zamo(&met, &tet);
}

void wedge3d_zamo_uniform::calumu(const double r, const double mu, std::array<double, 4> & umu, sim5metric & met) const {
    fourvelocity_zamo(&met, umu.data());
}

//! This is a 
double wedge3d_zamo_uniform::mean_free_path(const double r, const double mu) const {
    double snow = r * std::sqrt(1. - mu * mu);
    double Hnow = snow * tan_opening_angle;
    return (Hnow / tau);
}

double wedge3d_zamo_uniform::ne_sigmat(const double r, const double mu) const {
    double snow = r * std::sqrt(1. - mu * mu);
    double Hnow = snow * tan_opening_angle;
    return (tau / Hnow);
}

void wedge3d_zamo_uniform::genpos(std::mt19937_64 & gen, double & rnow, double & munow) const {
    rnow = 0.;
    munow = 0.;
}

wedge3d_corotation_uniform::wedge3d_corotation_uniform(const double aa, const double sminn, const double smaxx, const double mumaxx, const double tauu) {
    double cos_opening_angle; // cosine of the opening angle; the opening angle is \f$\pi / 2 - \theta\f$, where \f$\theta\$ is one of the Boyer-Lindquist coordinates
    name = "wedge3d_corotation_uniform";
    a = aa;
    smin = sminn;
    smax = smaxx;
    tau = tauu;
    
    mumin = 0.;
    mumax = mumaxx;
    cos_opening_angle = std::sqrt(1. - mumax * mumax);
    tan_opening_angle = mumax / cos_opening_angle;
    
    rmin = sminn;
    rmax = smax / cos_opening_angle;
    magnetised = false;
}

bool wedge3d_corotation_uniform::inside(const double r, const double mu) const {
    double s = r * std::sqrt(1. - mu * mu);
    return (s >= smin && s <= smax && mu >= 0. && mu <= mumax);
}

void wedge3d_corotation_uniform::caltet(const double r, const double mu, sim5tetrad & tet) const {
    double snow = std::sqrt(1. - mu * mu) * r, xnow, omega;
    xnow = std::sqrt(snow);
    omega = 1. / (xnow * xnow * xnow + a);
    
    sim5metric met;
    kerr_metric(a, r, mu, &met);
    tetrad_azimuthal(&met, omega, &tet);
}

void wedge3d_corotation_uniform::calumu(const double r, const double mu, std::array<double, 4> & umu, sim5metric & met) const {
    double snow = std::sqrt(1. - mu * mu) * r, xnow, omega;
    xnow = std::sqrt(snow);
    omega = 1. / (xnow * xnow * xnow + a);
    fourvelocity_azimuthal(omega, &met, umu.data());
}

//! This is a 
double wedge3d_corotation_uniform::mean_free_path(const double r, const double mu) const {
    double snow = r * std::sqrt(1. - mu * mu);
    double Hnow = snow * tan_opening_angle;
    return (Hnow / tau);
}

double wedge3d_corotation_uniform::ne_sigmat(const double r, const double mu) const {
    double snow = r * std::sqrt(1. - mu * mu);
    double Hnow = snow * tan_opening_angle;
    return (tau / Hnow);
}

void wedge3d_corotation_uniform::genpos(std::mt19937_64 & gen, double & rnow, double & munow) const {
    rnow = 0.;
    munow = 0.;
}

conical::conical(const double aa, const double hminn, const double thicknesss, const double angle, const double vv, const double tauu) {
    double snow, hmax;
    
    name = "conical";
    a = aa;
    hmin = hminn;
    thickness = thicknesss;
    opening_angle = angle;
    tan_opening_angle = std::tan(opening_angle / 360. * M_PI);
    v = vv;
    gamma = 1. / std::sqrt(1. - v * v);
    tau = tauu;
    
    rmin = hmin;
    hmax = hmin + thickness;
    snow = hmax * tan_opening_angle;
    rmax = std::sqrt(snow * snow + hmax * hmax);
    mumin = hmax / rmax;
    mumax = 1.;
    diameter_base = 2. * hmin * tan_opening_angle;
    magnetised = false;
}

bool conical::inside(const double r, const double mu) const {
    double hnow;
    hnow = r * mu;
    bool isinside = ((hnow <= (hmin + thickness)) && (mu >= mumin));
    
    return isinside;
}

double conical::mean_free_path(const double r, const double mu) const {
    return diameter_base/ tau;
}

double conical::ne_sigmat(const double r, const double mu) const {
    return tau / diameter_base;
}

void conical::caltet(const double r, const double mu, sim5tetrad & tet) const {
    sim5tetrad tet_zamo;
    sim5metric met, met_up;
    kerr_metric(a, r, mu, &met);
    kerr_metric_contravariant(a, r, mu, &met_up);
    tetrad_zamo(&met, &tet_zamo);
    
    //double sinmu = std::sqrt(1. - mu * mu);
    //if (sinmu != sinmu) {
    //    sinmu = 0.;
    //}
    
    /*
    double tan2, sin2, vx, vy, mu2 = mu * mu;
    sin2 = 1. - mu2;
    tan2 = sin2 / mu2;
    vx = v / std::sqrt(1. + tan2);
    vy = vx * std::sqrt(tan2);
    */
    
    std::array<double, 4> umu_zamo = {gamma, gamma * v, 0., 0.}, umu_bl;
    on2bl(umu_zamo.data(), umu_bl.data(), &tet_zamo);
    
    caltetrad(umu_bl, &met, &met_up, &tet);
}

void conical::calumu(const double r, const double mu, std::array<double, 4> & umu, sim5metric & met) const {
    sim5tetrad tet;
    kerr_metric(a, r, mu, &met);
    tetrad_zamo(&met, &tet);
    //double sinmu = std::sqrt(1. - mu * mu);
    //if (sinmu != sinmu) {
    //    sinmu = 0.;
    //}
    
    /*
    double tan2, sin2, vx, vy, mu2 = mu * mu;
    sin2 = 1. - mu2;
    tan2 = sin2 / mu2;
    vx = v / std::sqrt(1. + tan2);
    vy = vx * std::sqrt(tan2);
    //std::array<double, 4> umu_zamo = {gamma, gamma * vx, gamma * vy, 0.};
    */
    
    std::array<double, 4> umu_zamo = {gamma, gamma * v, 0., 0.};
    on2bl(umu_zamo.data(), umu.data(), &tet);
    //std::cout << "gamma = " << gamma << ", v = " << v << std::endl;
}
