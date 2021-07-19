//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

/*! \file electron_population.h

 */
#ifndef ELECTRON_POPULATION_H
#define ELECTRON_POPULATION_H

#include <random>
#include <string>
#include "const.h"
#include "photon_dist.h"

//! The class that saves the data of hot Klein-Nishina cross section
//! Now we save the data in this normal class and define a static `thermal_electron_data` variable
//! in classes that use the data, therefore we only need to read the data once
class thermal_electron_data {
public:
    const std::string hotdir = "/home/wdzhang/data/sim5corona/src/data/";
    //const std::string hotdir = std::string(std::getenv("DATADIR"));
    std::vector<double> logxarr, logthetatarr, hotxarr;
    long nthetat, nx;
    double dlogx, dlogthetat, xmin, xmax, thetatmin, thetatmax, logxmin, logxmax, logthetatmin, logthetatmax;
    thermal_electron_data() {readdata();};
    ~thermal_electron_data(){};
    void readdata();
};

//! virtual electron velocity distribution class
class electron_population {
public:
    //!
    virtual ~electron_population(){};
    //! Given photon energy, return the mean scattering cross section with respect to a electron population.
    virtual double cross_section(const double x) = 0;
    //! Given photon energy, sample the Lorentz factor of the scattering electron and the angle between the electron and the photon.
    virtual void sample_mu(const double x, double & gamma, double & emu, std::mt19937_64 & gen) const = 0;
};

class isotropic_thermal : public electron_population {
public:
    double te, thetat;
    ~isotropic_thermal(){};
    isotropic_thermal(const double tte) : te(tte), thetat(tte / ME_KEV) {};
    double cross_section(const double x);
    void sample_mu(const double x, double & gamma, double & emu, std::mt19937_64 & gen) const;
    void update(const double tte);
    static thermal_electron_data edata;
};

/*!
   Nonthermal electron population. For powerlaw distribution with number density \f$n_e\f$, index \f$p\f$ and lower and upper boundaries
    \f$\gamma_1\f$ and \f$\gamma_2\f$, the distribution of \f$\gamma\f$, the Lorentz factor, is
    \f{equation*}
    \frac{dN_e}{d\gamma} = \frac{n_e (1-p)}{\gamma_2^{1-p} - \gamma_1^{1-p}} \gamma^{-p}.
    \f}
    As \f$\gamma \equiv (1-\beta^2)^{-1/2}\f$, we have 
    \f{equation*}
    \frac{d\gamma}{d\beta} = \gamma^3 \beta.
    \f}
    Therefore in terms of \f$\beta\f$,
    \f{equation*}
    \frac{dN_e}{d\beta} = \frac{n_e (1-p)}{\gamma_2^{1-p} - \gamma_1^{1-p}} \gamma^{3-p} \beta.
    \f}
    If the electron is isotropic, the mean cross section is therefore
    \f{equation*}
    \sigma = \frac{1}{2n_e} \int_0^1 d\beta \int_{-1}^1 (1-\mu \beta) \sigma_{\rm KN}(x^\prime) \frac{dN_e}{d\beta} d\mu = \frac{1}{2}
        \frac{(1-p)}{\gamma_2^{1-p} - \gamma_1^{1-p}} \int_0^1 d\beta \int_{-1}^1 (1-\mu \beta) \gamma^{3-p}\beta \sigma_{KN}(x^\prime) d\mu,
    \f}
    where \f$x^\prime \equiv \gamma x (1-\mu \beta) \f$ is the dimensionless photon energy in the electron rest frame.
*/
class isotropic_nonthermal : public electron_population {
public:
    double gammamin, gammamax, betamin, betamax, index;
    ~isotropic_nonthermal(){};
    isotropic_nonthermal(const double gammaminn, const double gammamaxx, const double indexx) : gammamin(gammaminn),
        gammamax(gammamaxx), betamin(std::sqrt(1. - 1. / gammaminn / gammaminn)), betamax(std::sqrt(1. - 1. / gammamaxx / gammamaxx)),
        index(indexx) {};
    double cross_section(const double x);
    void sample_mu(const double x, double& gamma, double& emu, std::mt19937_64& gen) const;
private:
    bool iscalc = false;
    void calc();
    std::vector<double> logxarr, xsecarr;
    long nx;
    const double dlogx = 0.1, xmin = 1e-12, xmax = 1e6;
};

//! Now we offer the scattering functions the temperature profile of the corona with the `electron_dist` class virtually defined below, 
//! instead of the previously used `epop` class. `electron_dist` is `epop` plus info of electron population spatial distribution.
class electron_dist {
public:
    bool magnetised;
    electron_dist(const bool mag) : magnetised(mag) {};
    //! Return the hot Klein-Nishina cross-section in \f$[\sigma_T]\f$, given the photon energy in the fluid rest frame, and photon position
    virtual double hotcross(const double r, const double mu, const double x) = 0;
    ///*! Sample the four-velocity of the scattering electron,
    //given the photon energy in the fluid rest frame, and photon position */
    virtual void sample_mu(const double, const double, const double x, double & gamma, double & emu, std::mt19937_64 & gen) = 0;
    //! Return the characteristic physical value of the electron distribution; for thermal distribution, the electron temperature
    virtual double character(const double, const double) = 0;
    //! Virtual destructor
    ~electron_dist(){};
};

//! Thermal electron with uniform temperature
class uniform_thermal : public electron_dist {
public:
    // BH mass, in \f$M_\odot\f$
    isotropic_thermal epop;
    uniform_thermal(const double tee) :electron_dist(false), epop(tee) {};
    ~uniform_thermal(){};
    double hotcross(const double r, const double mu, const double x) { return epop.cross_section(x);};
    void sample_mu(const double r, const double mu, const double x, double & gamma, double & emu, std::mt19937_64 & gen)
        { epop.sample_mu(x, gamma, emu, gen);};
    double character(const double r, const double mu) { return epop.te;};
};

class uniform_nonthermal : public electron_dist {
    isotropic_nonthermal epop;
    uniform_nonthermal(const double gammaminn, const double gammamaxx, const double indexx) :
        electron_dist(false), epop(gammaminn, gammamaxx, indexx) {};
    ~uniform_nonthermal(){};
    double hotcross(const double r, const double mu, const double x) { return epop.cross_section(x);};
    void sample_mu(const double r, const double mu, const double x, double& gamma, double& emu, std::mt19937_64& gen)
        { epop.sample_mu(x, gamma, emu, gen);};
    double character(const double r, const double mu) { return epop.index;};
    // to do list
};

#endif
