//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

//! \file scatter.h
#ifndef _SLABSCATTER_H
#define _SLABSCATTER_H
#include <vector>
#include <array>
#include <random>
#include "detector.h"
#include "tridgeo.h"
#include "diskobj.h"
#include "superphoton.h"
#include "photon_dist.h"
#include "electron_population.h"

class xsectphicdfparams {
public:
    double mu;
    double poldeg;
    double result;
};

class pol_phi_params {
public:
    double par1;
    double rndnum;
};

double KN_pdfx(const double x, const double xp);
void sample_omega(std::mt19937_64 & gen, double & theta, double & ephi);
double sample_maxjutt(const double te, std::mt19937_64 & gen);
double sample_maxjutt1(const double te, std::mt19937_64 & gen);
void sample_thompsonx(std::mt19937_64 & gen, const double poldeg, double & mu, double & phi);
//void sample_thompsonx_fixmu(std::mt19937_64 & gen, const double poldeg, const double mu, double & phi);
double sample_scatter_mu(const double beta, std::mt19937_64 & gen);
double sample_kn_phi(std::mt19937_64 & gen, const double poldeg, const double mu, const double x1_over_x);
double sample_thompson_phi(std::mt19937_64 & gen, const double poldeg, const double mu);
double loboost(const double nx, const double ny, const double nz, const double beta, const double gamma,
    const std::array<double, 4> & kmu_f, std::array<double, 4> & kmu_e);
double sample_KN_xp(const double x, std::mt19937_64 & gen);

void scatter(const bool pol, const bool recoil, const bool KN, const double en, const double poldeg, const std::array<double, 4> & pe,
    std::array<double, 4> & kmu_f, std::array<double, 4> & fmu_f, double & en_p, double & poldeg_p,
    std::array<double, 4> & kmu_f_p, std::array<double, 4> & fmu_f_p, std::mt19937_64 & gen);

void scatter_restframe(const bool KN, const double * pk, const double * pf, const double coskt,
    const double kphi, const double poldeg, double & poldeg_p, const double z, double * pk_e_p, double * pf_e_p);
void tridsca_sp(const bool pol, const bool KN, const double a, const double rin,
    electron_dist & edist, const double weightnorm, 
    const double qfac0, const double ufac0, const double muinf_n0scatter, const double dr_ratio, superphoton sp, 
    const tridgeo & trid, std::vector<std::vector<double>> & vecarr_inf, std::vector<int> & nsca_inf, 
    std::vector<std::vector<double>> & vecarr_disc, std::vector<int> & nsca_disc, std::mt19937_64 & gen,
	const bool ns=false, const double rns_rg=5.);
double KN_xsection(const double x);
double klein_nishina(double a, double ap);

void sample_electron_KN(const double te, const double x, std::mt19937_64 & gen, double & gamma, double & mu, double & phi);
void calpe(const double emu, const double ephi, const double * kmu, std::array<double, 4> & pe);

#endif
