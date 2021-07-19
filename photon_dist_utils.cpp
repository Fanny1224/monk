//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

//! \file photon_dist_utils.cpp
#include <cmath>
#include <random>
#include <iostream>
#include "electron_population_utils.h"
#include "const.h"
#include "utils.h"

/*!
Following Pozdnyakov+1983 (http://adsabs.harvard.edu/abs/1983ASPRv...2..189P).
@param tbb photon temperature.
@param gen 64-bit MT random number generator
 */
double sample_bb(const double tbb, std::mt19937_64 & gen) {
    std::uniform_real_distribution<> dis(0., 1.);
    double xi1, xi2, xi3, xi4, alpha, m = 2., en;
    double series0, series1;
    
    series0 = 1.;
    series1 = 9. / 8.;
    
    xi1 = dis(gen);
    xi2 = dis(gen);
    xi3 = dis(gen);
    xi4 = dis(gen);
    
    if (xi1 < 1. / APERY_CONST) {
        alpha = 1.;
    } else {
        while(true) {
            if (xi1 >= series0 / APERY_CONST && xi1 < series1 / APERY_CONST) {
                alpha = m;
                break;
            } else {
                ++m;
                series0 = series1;
                series1 += 1. / m / m / m;
            }
        }
    }
    
    en = -1. * tbb * std::log(xi2 * xi3 * xi4) / alpha;
    
    return en;
}

/*!
Following Pozdnyakov+1983 (http://adsabs.harvard.edu/abs/1983ASPRv...2..189P).
To increase the statistics at high energy, we increase the probability of photons above `ecut` and decrease their weight correspondingly.
We use a highly inefficient rejection method (1/ratio of \ref sample_bb).
@param tbb photon temperature.
@param gen 64-bit MT random number generator
@param ecut the critical energy above which the probability will be increased by a factor of `ratio`, and correspondingly weight
    dereased by the same factor
@param ratio the factor by which the probability increased
@param en energy of the sampled photon
@param weight statistical weight of the sampled photon
*/
void sample_bb1(const double tbb, std::mt19937_64 & gen, const double ecut, const double ratio, double & en, double & weight) {
    std::uniform_real_distribution<> dis(0, 1);
    bool found = false;
    double r1;
    while (!found) {
        en = sample_bb(tbb, gen);
        if (en >= ecut) {
            found = true;
            weight = 1. / ratio;
        } else {
            r1 = dis(gen);
            if (r1 < 1. / ratio) {
                found = true;
                weight = 1.;
            }
        }
    }
}

