//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

//! \file tridgeo.h
//! Here we define an ABC `tridgeo` as a base class for classes of all kinds of 3D geometries.
#ifndef _TRIDGEO_H
#define _TRIDGEO_H
#include <string>
#include <random>
#include <array>
#include "sim5lib.h"

/*!
An abstract base class for all kinds of 3D corona. We expect that a general corona, regardless of its shape,
would be able to tell us the following: 1) if a photon is inside the corona given position; 2) the four-velocity of
corona fluid anywhere in the corona; 3) the local tetrad attached to corona fluid; 4) the size. Therefore we define the corresponding
virtual functions that will be implemented by the derived classes. 
The corona should have finite size, such that it has finite lower and upper bounds in \f$r\f$ and \f$\mu\f$.
These can be used to tell if the photon can reach the corona given current position and wave vector.
*/
class tridgeo {
public:
    //! geometry name
    std::string name;
    //! BH spin
    double a;
    //! Minimum \f$r\f$ of corona
    double rmin;
    //! Maximum \f$r\f$ of corona
    double rmax;
    //! Minimum \f$\mu\f$
    double mumin;
    //! Maximum \f$\mu\f$
    double mumax;
    //! whether if the corona is magnetised
    bool magnetised;
    //! Tell if the photon is inside
    virtual bool inside(const double r, const double mu) const = 0;
    //! Calculate the tetrad attached to the corona fluid
    virtual void caltet(const double, const double, sim5tetrad &) const = 0;
    //! Calculate the four velocity of the corona fluid
    virtual void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const = 0;
    //! Virtual destructor
    virtual ~tridgeo(){};
    //! Return the length scale of the corona
    virtual double length() const = 0;
    virtual void genpos(std::mt19937_64 & gen, double &, double &) const = 0;
    virtual double mean_free_path(const double, const double) const = 0;
    //! returns \f$n_e * \sigma\f$, given r, mu
    virtual double ne_sigmat(const double, const double) const = 0;
    virtual void write_brem_sp(const size_t ndim1, const size_t ndim2, std::ofstream & ofile) const = 0;
    virtual double cal_bfield(const double r, const double mu, const double te_K) const = 0;
};

//! ZAMO slab, uniform temperature and density
class zamoslab3d : public tridgeo {
public:
    //! Corona radius
    double s;
    //! Minimum height
    double hmin;
    //! Maximum height
    double hmax;
    //! \f$\tau \equiv (h_{\rm max} - h_{\rm min}) n_e /2\f$, the optical depth
    double tau;
    zamoslab3d(const double, const double, const double, const double, const double);
    bool inside(const double r, const double mu) const;
    void caltet(const double, const double, sim5tetrad &) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    ~zamoslab3d(){};
    double length() const { return s * tau;};
    void genpos(std::mt19937_64 & gen, double &, double &) const;
    double mean_free_path(const double, const double) const;
    double ne_sigmat(const double, const double) const;
    void write_brem_sp(const size_t ndim1, const size_t ndim2, std::ofstream & ofile) const {};
    double cal_bfield(const double, const double, const double) const {return 0.;};
};

//! ZAMO sandwich, uniform temperature and density
class zamosandwich : public tridgeo {
public:
    //! Corona radius
    double s;
    //! inner edge of the corona
    double smin;
    //! Minimum height
    double hmin;
    //! Maximum height
    double hmax;
    //! \f$\tau \equiv (h_{\rm max} - h_{\rm min}) n_e /2\f$, the optical depth
    double tau;
    zamosandwich(const double, const double, const double, const double, const double, const double);
    bool inside(const double r, const double mu) const;
    void caltet(const double, const double, sim5tetrad &) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    ~zamosandwich(){};
    double length() const { return s * tau;};
    void genpos(std::mt19937_64 & gen, double &, double &) const {};
    double mean_free_path(const double, const double) const {return 0.5 * (hmax - hmin) / tau;};
    double ne_sigmat(const double, const double) const {return 2. * tau / (hmax - hmin);};
    void write_brem_sp(const size_t ndim1, const size_t ndim2, std::ofstream & ofile) const {};
    double cal_bfield(const double, const double, const double) const {return 0.;};
};

//! ZAMO sandwich, uniform temperature and density
class kepsandwich : public tridgeo {
public:
    //! Corona radius
    double s;
    //! inner edge of the corona
    double smin;
    //! Minimum height
    double hmin;
    //! Maximum height
    double hmax;
    //! \f$\tau \equiv (h_{\rm max} - h_{\rm min}) n_e /2\f$, the optical depth
    double tau;
    kepsandwich(const double, const double, const double, const double, const double, const double);
    bool inside(const double r, const double mu) const;
    void caltet(const double, const double, sim5tetrad &) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    ~kepsandwich(){};
    double length() const { return s * tau;};
    void genpos(std::mt19937_64 & gen, double &, double &) const {};
    double mean_free_path(const double, const double) const {return 0.5 * (hmax - hmin) / tau;};
    double ne_sigmat(const double, const double) const {return 2. * tau / (hmax - hmin);};
    void write_brem_sp(const size_t ndim1, const size_t ndim2, std::ofstream & ofile) const {};
    double cal_bfield(const double, const double, const double) const {return 0.;};
};

//! Slab that is co-rotating with the underlying Keplerian disc.
class kepslab3d : public tridgeo {
public:
    //! Corona radius
    double s;
    //! Minimum height
    double hmin;
    //! Maximum height
    double hmax;
    double tau;
    kepslab3d(const double, const double, const double, const double, const double);
    bool inside(const double r, const double mu) const;
    void caltet(const double, const double, sim5tetrad &) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    double length() const { return s * tau;};
    ~kepslab3d() {};
    void genpos(std::mt19937_64 & gen, double &, double &) const;
    double mean_free_path(const double, const double) const;
    double ne_sigmat(const double, const double) const;
    void write_brem_sp(const size_t ndim1, const size_t ndim2, std::ofstream & ofile) const {};
    double cal_bfield(const double, const double, const double) const {return 0.;};
};
    
//! ZAMO sphere
class zamosphere3d : public tridgeo {
public:
    //! Height
    double h;
    //! Radius
    double R;
    double tau;
    zamosphere3d(const double, const double, const double, const double);
    bool inside(const double r, const double mu) const;
    void caltet(const double, const double, sim5tetrad &) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    ~zamosphere3d() {};
    double length() const { return tau;};
    void genpos(std::mt19937_64 & gen, double &, double &) const;
    double mean_free_path(const double, const double) const;
    double ne_sigmat(const double, const double) const;
    void write_brem_sp(const size_t ndim1, const size_t ndim2, std::ofstream & ofile) const;
    double cal_bfield(const double, const double, const double) const {return 0.;};
};

//! ZAMO spherical shell.
class zamosphere3d_truncated : public tridgeo {
public:
    //! Radius of outer boundary
    double R;
    //! Radius of inner boundary
    double rin;
    double tau;
    zamosphere3d_truncated(const double, const double, const double, const double);
    bool inside(const double r, const double mu) const;
    void caltet(const double, const double, sim5tetrad &) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    ~zamosphere3d_truncated() {};
    double length() const { return tau;};
    void genpos(std::mt19937_64 & gen, double &, double &) const;
    double mean_free_path(const double, const double) const;
    double ne_sigmat(const double, const double) const;
    void write_brem_sp(const size_t ndim1, const size_t ndim2, std::ofstream & ofile) const {};
    double cal_bfield(const double, const double, const double) const {return 0.;};
};

//! Spherical corona that is rotating with uniform linear circular velocity
class rotsphere3d: public tridgeo {
public:
    //! Height
    double h;
    //! Radius
    double R;
    //! Linear circular velocity
    double velocity;
    //! Lorentz factor, \f$\gamma=1/\sqrt{1 - v^2}\f$.
    double gamma;
    double tau;
    rotsphere3d(const double, const double, const double, const double, const double);
    ~rotsphere3d() {};
    bool inside(const double r, const double mu) const;
    void caltet(const double, const double, sim5tetrad &) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    double length() const { return tau;};
    void genpos(std::mt19937_64 & gen, double &, double &) const;
    double mean_free_path(const double, const double) const;
    double ne_sigmat(const double, const double) const;
    void write_brem_sp(const size_t ndim1, const size_t ndim2, std::ofstream & ofile) const {};
    double cal_bfield(const double, const double, const double) const {return 0.;};
};

//! Abstract base class for corona in flat spacetime. No tetrad/\f$U^\mu\f$ calculation compared with \ref tridgeo
class tridgeo_sr {
public:
    //! Tell if the photon is inside. Accept std:array<double,3> argument, in accordance with \ref superphoton_sr.
    virtual bool inside(const std::array<double, 4> & xmu) const = 0;
    //! Return the length scale of the corona
    virtual double length() const = 0;
    //! Virtual constructor
    virtual ~tridgeo_sr(){};
};

//! Spherical corona in flat spacetime
class sphere3d_sr : public tridgeo_sr {
public:
    //! Radius
    double radius;
    //! Destructor
    ~sphere3d_sr(){};
    //! Constructor
    sphere3d_sr(const double rr) : radius(rr) {};
    //! Constructor without argument
    sphere3d_sr(){};
    bool inside(const std::array<double, 4> & xmu) const;
    double length() const { return radius;};
};

//! Slab corona in flat spacetime
class slab3d_sr : public tridgeo_sr {
public:
    //! Thickness
    double h;
    //! Radius
    double s;
    //! Destructor
    ~slab3d_sr(){};
    //! Constructor without argument
    slab3d_sr(const double hh, const double ss) : h(hh), s(ss) {};
    //! Constructor
    slab3d_sr(){};
    bool inside(const std::array<double, 4> & xmu) const;
    double length() const { return std::sqrt(s * h);};
};

class infislab3d_sr : public tridgeo_sr {
public:
    double h;
    ~infislab3d_sr(){};
    //! Constructor without argument
    infislab3d_sr(const double hh) : h(hh) {};
    //! Constructor
    infislab3d_sr(){};
    bool inside(const std::array<double, 4> & xmu) const {return (std::abs(xmu[3]) < h);};
    double length() const { return h * 40.;};
};
    
/*
// co-rotating slab
class coroslab3d: public tridgeo {
public:
    //! Corona radius
    double s;
    //! Minimum height
    double hmin;
    //! Maximum height
    double hmax;
    coroslab3d(const double, const double, const double, const double);
    ~coroslab3d(){};
    
    bool inside(const double r, const double mu) const;
    sim5tetrad caltet(const double, const double) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    double length() const { return std::sqrt((hmax - hmin) * s);};
    double mean_free_path(const double, const double, const double) const;
    double ne_sigmat(const double, const double, const double) const;
};
*/

//! Wedge geometry; see Schnittman & Krolik 2010. The corona is co-rotating with the underlying disc
class wedge3d : public tridgeo {
public:
    // outer boundary in the horizontal direction
    double smax;
    // inner boundary in the horizontal direction
    double smin;
    double tan_opening_angle;
    double tau;
    wedge3d(const double, const double, const double, const double, const double);
    ~wedge3d(){};
    bool inside(const double r, const double mu) const;
    void caltet(const double, const double, sim5tetrad &) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    double length() const { return (smax - smin) * tau;};
    double mean_free_path(const double r, const double mu) const;
    double ne_sigmat(const double, const double) const;
    void genpos(std::mt19937_64 & gen, double &, double &) const;
    void write_brem_sp(const size_t ndim1, const size_t ndim2, std::ofstream & ofile) const {};
    double cal_bfield(const double, const double, const double) const {return 0.;};
};

//! 3D wedge, and the electrons are uniformly distributed, unlike the `wedge3d` class
class wedge3d_zamo_uniform : public tridgeo {
public:
    // outer boundary in the horizontal direction
    double smax;
    // inner boundary in the horizontal direction
    double smin;
    double tan_opening_angle;
    double tau;
    wedge3d_zamo_uniform(const double, const double, const double, const double, const double);
    ~wedge3d_zamo_uniform(){};
    bool inside(const double r, const double mu) const;
    void caltet(const double, const double, sim5tetrad &) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    double length() const { return (smax - smin) * tau;};
    double mean_free_path(const double, const double) const;
    double ne_sigmat(const double, const double) const;
    void genpos(std::mt19937_64 & gen, double &, double &) const;
    void write_brem_sp(const size_t ndim1, const size_t ndim2, std::ofstream & ofile) const {};
    double cal_bfield(const double, const double, const double) const {return 0.;};
};

//! 3D wedge, and the electrons are uniformly distributed, unlike the `wedge3d` class
class wedge3d_corotation_uniform : public tridgeo {
public:
    // outer boundary in the horizontal direction
    double smax;
    // inner boundary in the horizontal direction
    double smin;
    double tan_opening_angle;
    double tau;
    wedge3d_corotation_uniform(const double aa, const double sminn, const double smaxx, const double mumaxx, const double);
    ~wedge3d_corotation_uniform(){};
    bool inside(const double r, const double mu) const;
    void caltet(const double, const double, sim5tetrad &) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    double length() const { return (smax - smin) * tau;};
    double mean_free_path(const double, const double) const;
    double ne_sigmat(const double, const double) const;
    void genpos(std::mt19937_64 & gen, double &, double &) const;
    void write_brem_sp(const size_t ndim1, const size_t ndim2, std::ofstream & ofile) const {};
    double cal_bfield(const double, const double, const double) const {return 0.;};
};

//! Conical 
class conical : public tridgeo {
public:
    // height
    double hmin;
    // thickness
    double thickness;
    // *half* opening angle *in degree*
    double opening_angle;
    // tan of the opening angle
    double tan_opening_angle;
    // velocity of the fluid along the +z direction, as measured by a zamo
    double v;
    // Lorentz factor
    double gamma;
    // jet base diameter
    double diameter_base;
    double tau;
    conical(const double aa, const double hminn, const double thicknesss, const double angle, const double vv, const double);
    ~conical(){};
    bool inside(const double r, const double) const;
    void caltet(const double, const double, sim5tetrad &) const;
    void calumu(const double, const double, std::array<double, 4> &, sim5metric &) const;
    double mean_free_path(const double, const double) const;
    double ne_sigmat(const double, const double) const;
    void genpos(std::mt19937_64 & gen, double &, double &) const {};
    double length() const {return thickness * tau;};
    void write_brem_sp(const size_t ndim1, const size_t ndim2, std::ofstream & ofile) const {};
    void write_syn_sp(const size_t ndim1, const size_t ndim2, const double bfield, std::ofstream & ofile) const {};
    double cal_bfield(const double, const double, const double) const {return 0.;};
};
#endif
