//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

//! \file gridgeod.h
#ifndef _GRIDGEOD_H
#define _GRIDGEOD_H
#include "diskobj.h"
#include "sim5lib.h"
#include "tridgeo.h"
#include "detector.h"
#include "electron_population.h"

void griddisk_nt_sp_pol_writeparam(const bool pol, const double rin, const double rout, 
    const double thetamax, const long nr, const long ntheta, const long nphi, const diskobj & disk);
//void grid3dcorona_sp_writeparam(const double rin, const double rout, const double thetamax, const tridgeo & trid,
//    const long nr, const long ntheta, const long nphi);
void grid3dcorona_sp_writeparam(const double rin, const double rout, const double thetamax, const diskobj & disk, 
    const tridgeo & trid, const long nr, const long ntheta, const long nphi, const std::string discfile, 
    const std::string selffile, const std::string scafile);
void register_sp_nt(const bool pol, const std::string parafile, const double m, const double mdot, const double fcol,
    const double rtr, const long nphoton, std::mt19937_64 & gen, std::vector<std::ofstream> & ofiles);
void register_sp_tridsca(const bool pol, const bool KN, const bool chandra, const bool progress, const std::string parafile,
    const double m, const double mdot, const double fcol,
    const double te, const double mfp, const double dr, const long nphoton, std::vector<std::ofstream> & ofiles_inf,
    std::ofstream & nscafile_inf, std::vector<std::ofstream> & ofiles_disc, std::ofstream & nscafile_disc);
void plotgeo(const bool savevec, const long nphoton, const tridgeo & trid, electron_dist & edist, 
    const double xmax, const double ymin, const double ymax, std::mt19937_64 & gen, std::vector<std::ofstream> & ofiles);
void mkgrid(const double rin, const double rout, const double thetamax, const long nr, const long ntheta, 
    const long nphi, std::vector<double> & radius, std::vector<double> & theta, std::vector<double> & phi);
void grid_ns_surface_writeparam(const bool pol, const double mns, const double rns_rg, const double pspin_ms,
	const long nincl, const long ntheta, const long nphi);
void register_sp_ns_surface(const std::string parafile, const long nphoton,
	const double tns_kev, std::mt19937_64 & gen, std::vector<std::ofstream> & ofiles);
void grid3dcorona_sp_writeparam_ns(const double mns, const double rns_rg, const double pspin_ms, const double rin,
	const double inclmax, const tridgeo & trid, const long nincl,
	const long ntheta, const long nphi, const std::string scafile, const std::string discfile);

#endif
