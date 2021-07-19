//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

/*! \file electron_population.cpp
 */

#include <string>
#include <vector>
#include <experimental/filesystem>
#include "electron_population.h"
#include "utils.h"
#include "scatter.h"
#include "electron_population_utils.h"

namespace fs = std::experimental::filesystem;

thermal_electron_data isotropic_thermal::edata;

void thermal_electron_data::readdata() {
	std::cout << "hotdir = " << hotdir << std::endl;
    std::string tpath = hotdir + "logthetat.dat", xpath = hotdir + "logx.dat", hotpath = hotdir + "hotx.dat";
    if (!fs::exists(tpath)) {
        throw(tpath + " does not exist!");
        return;
    }
    if (!fs::exists(xpath)) {
        throw(xpath + " does not exist!");
        return;
    }
    if (!fs::exists(hotpath)) {
        throw(hotpath + " does not exist!");
        return;
    }
    std::cout << "hotpath = " << hotpath << std::endl;
    
    rdoublevec(logthetatarr, tpath);
    rdoublevec(logxarr, xpath);
    rdoublevec(hotxarr, hotpath);
    
    nx = logxarr.size();
    nthetat = logthetatarr.size();
    logxmin = logxarr[0];
    logxmax = logxarr[nx - 1];
    logthetatmin = logthetatarr[0];
    logthetatmax = logthetatarr[nthetat - 1];
    dlogx = logxarr[1] - logxarr[0];
    dlogthetat = logthetatarr[1] - logthetatarr[0];
    xmin = std::pow(10., logxmin);
    xmax = std::pow(10., logxmax);
    thetatmin = std::pow(10., logthetatmin);
    thetatmax = std::pow(10., logthetatmax);
    
std::cout << "Hot cross section read!" << std::endl;
}

double isotropic_thermal::cross_section(const double x) {
    double logx, logthetat, hotx, hotx_t0, hotx_t1;
    long it0, it1, ix0, ix1;
    long index00, index01, index10, index11;
    logthetat = std::log10(thetat);
    
    if (thetat * x < 1.e-6)
        return 1.;
    
    if (logthetat < -4.)
        return KN_xsection(x);
    
    if (x > edata.xmin && x < edata.xmax && thetat > edata.thetatmin && thetat < edata.thetatmax) {
        logx = std::log10(x);
        ix0 = (long)(std::floor((logx - edata.logxmin) / edata.dlogx));
        ix1 = ix0 + 1;
        it0 = (long)(std::floor((logthetat - edata.logthetatmin) / edata.dlogthetat));
        it1 = it0 + 1;
        /*
        std::cout << "logx = " << logx << " " << ix0 << " " << ix1 << " " << it0 << " " << it1 << std::endl;
        std::cout << "logxmin = " << logxmin << std::endl;
        std::cout << "logxarr[ix0] = " << logxarr[ix0] << ", logxarr[ix1] = " << logxarr[ix1] << std::endl;
        */
        
        index00 = ix0 * edata.nthetat + it0;
        index01 = ix0 * edata.nthetat + it1;
        index10 = ix1 * edata.nthetat + it0;
        index11 = ix1 * edata.nthetat + it1;
        
        hotx_t0 = (edata.logxarr[ix1] - logx) * edata.hotxarr[index00] + (logx - edata.logxarr[ix0]) * edata.hotxarr[index10];
        hotx_t1 = (edata.logxarr[ix1] - logx) * edata.hotxarr[index01] + (logx - edata.logxarr[ix0]) * edata.hotxarr[index11];
        hotx = ((edata.logthetatarr[it1] - logthetat) * hotx_t0 +
             (logthetat - edata.logthetatarr[it0]) * hotx_t1) / edata.dlogthetat / edata.dlogx;
             /*
        std::cout << "index00 = " << index00 << std::endl;
        std::cout << "index10 = " << index10 << std::endl;
        std::cout << "arrsize = " << edata.hotxarr.size() << std::endl;
        std::cout << "hotx_t0 = " << hotx_t0 << ", edata.hotxarr[index00] = " << edata.hotxarr[index00] << std::endl;
        std::cout << "hotx_t1 = " << hotx_t1 << std::endl;
        std::cout << "hotx = " << hotx << std::endl;
        */
        
    } else {
        hotx = isotropic_thermal_calbetainte(thetat, x);
    }
    
    return hotx;
}

void isotropic_thermal::update(const double tte) {
    te = tte;
    thetat = tte / ME_KEV;
}

void isotropic_thermal::sample_mu(const double x, double & gamma, double & emu, std::mt19937_64 & gen) const {
    double beta, xe, total_xsect, rndnum;
    std::uniform_real_distribution<> dis(0., 1.);
    rndnum = dis(gen);

    while(true) {
        // first sample gamma
        gamma = sample_maxjutt1(thetat, gen);
        beta = std::sqrt(1. - 1. / gamma / gamma);
        // then sample mu
        emu = sample_scatter_mu(beta, gen);
        xe = x * gamma * (1. - beta * emu);
        total_xsect = KN_xsection(xe);
        rndnum = dis(gen);
        if (rndnum <= total_xsect)
            break;
    }
}

void isotropic_nonthermal::calc() {
    double logxmin = std::log10(xmin), logxmax = std::log10(xmax);
    nx = std::floor((logxmax - logxmin) / dlogx);
    logxarr.resize(nx);
    xsecarr.resize(nx);
    std::cout << "starting calculating cross section:" << std::endl;
    std::cout << "betamin = " << betamin << ", betamax = " << betamax << std::endl;
    
    for (long i = 0; i < nx; ++i) {
        //std::cout << "i = " << i << " / " << nx << std::endl;
        logxarr[i] = logxmin + (double)(i) * dlogx;
        xsecarr[i] = isotropic_nonthermal_calbetainte(betamin, betamax, index, std::pow(10., logxarr[i]));
        //std::cout << "xsecarr[i] = " << xsecarr[i] << std::endl;
    }
    wdoublevec(logxarr, "nonthermal_logx.dat");
    wdoublevec(xsecarr, "nonthermal_xsec.dat");
    iscalc = true;
    std::cout << "Cross section calculated!" << std::endl;
}

double isotropic_nonthermal::cross_section(const double x) {
    double logx = std::log10(x);
    if (!iscalc) {
        calc();
    }
    if (x >= xmin && x <= xmax) {
        return linear_interpolate(logx, logxarr, xsecarr);
    } else {
        return isotropic_nonthermal_calbetainte(betamin, betamax, index, std::pow(10., logx));
    }
}

void isotropic_nonthermal::sample_mu(const double x, double & gamma, double & emu, std::mt19937_64 & gen) const {
    double beta, xe, total_xsect, rndnum;
    std::uniform_real_distribution<> dis(0., 1.);
    rndnum = dis(gen);

    while(true) {
        // first sample gamma
        gamma = sample_pl(-1. * index, gammamin, gammamax, gen);
        beta = std::sqrt(1. - 1. / gamma / gamma);
        // then sample mu
        emu = sample_scatter_mu(beta, gen);
        xe = x * gamma * (1. - beta * emu);
        total_xsect = KN_xsection(xe);
        rndnum = dis(gen);
        std::cout << "gamma = " << gamma << std::endl;
        /*
        std::cout << "xe = " << xe << std::endl;
        std::cout << "beta = " << beta << std::endl;
        std::cout << "total_xsect = " << total_xsect << std::endl;
        std::cout << "emu = " << emu << std::endl;
        std::cout << "rndnum = " << rndnum << std::endl;
        std::cout << "================ " << std::endl;
        */
        if (rndnum <= total_xsect)
            break;
    }
}

