//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

#ifndef ELECTRON_POPULATON_UTILS_H
#define ELECTRON_POPULATON_UTILS_H
class Ggammaparam {
public:
    double thetat, result;
};
double isotropic_thermal_calbetainte(const double thetat, const double x);
double isotropic_thermal_calbetainte_lowtem(const double thetat, const double x);
double isotropic_nonthermal_calbetainte(const double beta1, const double beta2, const double index, const double x);
double maxwell_jutt(const double gamma, const double thetat);
#endif
