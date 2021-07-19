//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

//! \file photon_dist.h
#ifndef PHOTON_DIST_H
#define PHOTON_DIST_H
#include <random>
#include <vector>
#include "photon_dist_utils.h"
#include "utils.h"
#include "const.h"

class photon_dist {
public:
    //! Virtual destructor
    virtual ~photon_dist(){};
    //! Constructor
    photon_dist(){};
    //! Sampling method
    virtual double genen(std::mt19937_64 & gen) const = 0;
    virtual double cal_weight() const = 0;
};

//! Planckian distribution
class bb : public photon_dist {
public:
     //! blackbody temperature in keV 
    double tbb;
    //! Destructor
    ~bb(){};
    //! Constructor
    bb(const double tbbb) : tbb(tbbb) {};
    double genen(std::mt19937_64 & gen) const {return sample_bb(tbb, gen);};
    double cal_weight() const {return 0.;};
};

//! \brief Powerlaw distribution.
//! The probability density distribution \f$p(E)\propto E^\gamma\f$ where \f$E\in[E_{\rm min}, E_{\rm max}]\f$.
//! Note that \f$\gamma\f$ is the opposite of photon index.
class pl : public photon_dist {
public:
    //! Opposite of photon index
    double gamma;
    //! Lower bound energy
    double emin;
    //! Upper bound energy
    double emax;
    //! Destructor
    ~pl(){};
    //! Constructor
    pl(const double gammaa, const double eminn, const double emaxx) : gamma(gammaa), emin(eminn), emax(emaxx) {};
    double genen(std::mt19937_64 & gen) const {return sample_pl(gamma, emin, emax, gen);};
    double cal_weight() const {return 0.;};
};

class cutoffpl : public photon_dist {
public:
    double gamma;
    //! Lower bound energy
    double emin;
    //! Upper bound energy
    double ecut;
    //! Destructor
    ~cutoffpl(){};
    //! Constructor
    cutoffpl(const double gammaa, const double eminn, const double ecutt) : gamma(gammaa), emin(eminn), ecut(ecutt) {};
    double genen(std::mt19937_64 & gen) const {return sample_cutoffpl(gamma, emin, ecut, gen);};
    double cal_weight() const {return 0.;};
};

class monoenergetic : public photon_dist {
    double en;
    ~monoenergetic(){};
    monoenergetic(const double enn) : en(enn) {};
    double genen(std::mt19937_64 & gen) const {return en;};
    double cal_weight() const {return 0.;};
};

#endif
