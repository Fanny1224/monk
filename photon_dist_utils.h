//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

#ifndef PHOTON_DIST_UTILS_H
#define PHOTON_DIST_UTILS_H
//! Sample energy of photons that follow Planckian distribution. 
double sample_bb(const double tbb, std::mt19937_64 & gen);
//! Weighted sample of photon energy that follow Planckian distribution. 
void sample_bb1(const double tbb, std::mt19937_64 & gen, const double ecut, const double ratio, double & en, double & weight);
#endif
