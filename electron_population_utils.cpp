//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

//! \file electron_population_utils.cpp

#include <cmath>
#include <array>
#include <random>
#include "utils.h"
#include "const.h"
#include "scatter.h"
#include "electron_population_utils.h"
#include "rootsearch.h"

/*! \brief \f$\mu\f$ integrand: \f$(1-\mu\beta)\sigma_{\rm KN}(x^\prime)\f$. Note that this integrand is independent of the distribution
once the electron velocity distribution is isotropic

@param mu cosine of angle made by photon and electron
@param beta the dimensionless electron velocity
@param x dimensionless photon energy
*/
double KN_integrand_mu(const double mu, const double beta, const double x) {
    //std::cout << "beta = " << beta << std::endl;
    double xp, gamma;
    gamma = 1. / std::sqrt(1. - beta * beta);
    xp = gamma * (1. - mu * beta) * x;
    return (1. - mu * beta) * KN_xsection(xp);
}

/*! \brief the same with \ref KN_integrand_mu, but the arguments are passed with a void pointer,
to be compatible with \ref int_trape1.
@param mu cosine of angle made by photon and electron
@param params a reinterpreted pointer to std:array<double, 2> which contains \f$[\beta, x]\f$
 */
double KN_integrand_mu1(const double mu, void * params) {
    std::array<double, 2> * paramarr = reinterpret_cast<std::array<double, 2> *> (params);
    return KN_integrand_mu(mu, (*paramarr)[0], (*paramarr)[1]);
}

/*! \brief The part of the integrand which is function of \f$\beta\f$ only; i.e., \f$\gamma^5 \beta^2 e^{-\gamma/\theta_T}\f$
@param thetat dimensionless electron temperature
@param beta the dimensionless electron velocity
*/
double isotropic_thermal_integrand_beta(const double thetat, const double beta) {
    double gamma = 1. / std::sqrt(1. - beta * beta);
    double gamma2 = gamma * gamma;
    return gamma2 * gamma2 * gamma * beta * beta * std::exp(-1. * gamma / thetat);
}

/*! \brief The \f$\beta\f$ integrand; \f$\gamma^5 \beta^2 e^{-\gamma/\theta_T}\int_{-1}^1 d\mu (1-\mu\beta) \sigma_{\rm KN}(x^\prime) \f$
@param beta the dimensionless electron velocity
@param params a reinterpreted pointer to std:array<double, 2> which contains \f$[\theta_T, x]\f$
 */
double isotropic_thermal_integrand_beta1(const double beta, void * params) {
    std::array<double, 2> * paramarr = reinterpret_cast<std::array<double, 2> *> (params);
    std::array<double, 2> muparams = {beta, (*paramarr)[1]};
    void * intparams = reinterpret_cast<void *>(&muparams);
    double muint = int_trape1(&KN_integrand_mu1, -1., 1., 1e-3, intparams);
    return isotropic_thermal_integrand_beta((*paramarr)[0], beta) * muint;
}

/*! \brief Evaluating \f$\beta\f$ integration. We first evaluate \f$\mu\f$ integral given \f$beta\f$.
@param thetat dimensionless electron temperature
@param x dimensionless photon energy
 */
double isotropic_thermal_calbetainte(const double thetat, const double x) {
    //double abissca[NGAUSS1], weights[NGAUSS1];
    //double integauss = 0.;
    //gauleg(0., 1., abissca, weights, NGAUSS1);
    //
    //for (int i = 0; i < NGAUSS1; ++i) {
    //    integauss += weights[i] * hotcross_integrand_beta(thetat, abissca[i]) * hotcross_calmuinte(abissca[i], x);
    //}
    //return 0.5 * integauss / thetat / bessel_k2(1. / thetat);
    std::array<double, 2> paramarr = {thetat, x};
    void * param = reinterpret_cast<void *>(&paramarr);
    double inte = int_trape1(&isotropic_thermal_integrand_beta1, 0., 1., 1e-3, param);

    //std::cout << "bessel_k2 = " << bessel_k2(1. / thetat) << std::endl;
    
    return 0.5 * inte / thetat / bessel_k2(1. / thetat);
}

double isotropic_thermal_integrand_beta1_lowtem(const double beta, void * params) {
    std::array<double, 2> * paramarr = reinterpret_cast<std::array<double, 2> *> (params);
    std::array<double, 2> muparams = {beta, (*paramarr)[1]};
    void * intparams = reinterpret_cast<void *>(&muparams);
    double muint = int_trape1(&KN_integrand_mu1, -1., 1., 1e-3, intparams);
    return beta * beta * std::exp(-0.5 * beta * beta / (*paramarr)[0]) * muint;
}

/*! \brief Evaluating \f$\beta\f$ integration. We first evaluate \f$\mu\f$ integral given \f$beta\f$.
@param thetat dimensionless electron temperature
@param x dimensionless photon energy
 */
double isotropic_thermal_calbetainte_lowtem(const double thetat, const double x) {
    std::array<double, 2> paramarr = {thetat, x};
    void * param = reinterpret_cast<void *>(&paramarr);
    double inte = int_trape1(&isotropic_thermal_integrand_beta1_lowtem, 0., 1., 1e-3, param);

    return std::sqrt(0.5 / (M_PI * thetat * thetat * thetat)) * inte;
}


/* \brief The part of the integrand which is function of \f$\beta\f$ only; i.e., \f$\gamma^5 \beta^2 e^{-\gamma/\theta_T}\f$
@param beta the dimensionless electron velocity
double isotropic_nonthermal_integrand_beta(const double beta, const double index) {
    double gamma = 1. / std::sqrt(1. - beta * beta);
    return beta * std::pow(gamma, 3. - index);
}
*/


/*! \brief The \f$\beta\f$ integrand; \f$\gamma^{3-p} \beta \int_{-1}^1 d\mu (1-\mu\beta) \sigma_{\rm KN}(x^\prime) \f$
@param beta the dimensionless electron velocity
@param params a reinterpreted pointer to std:array<double, 2> which contains \f$[p, x]\f$, where \f$p\f$ is the index of the electron
    velocity distribution
 */
double isotropic_nonthermal_integrand_beta(const double beta, void * params) {
    double gamma = 1. / std::sqrt(1. - beta * beta);
    std::array<double, 2> * paramarr = reinterpret_cast<std::array<double, 2> *> (params);
    std::array<double, 2> muparams = {beta, (*paramarr)[1]};
    void * intparams = reinterpret_cast<void *>(&muparams);
    double muint = int_trape1(&KN_integrand_mu1, -1., 1., 1e-3, intparams);
    return muint * beta * std::pow(gamma, 3. - (*paramarr)[0]);
}

/*! \brief Evaluating \f$\beta\f$ integration. We first evaluate \f$\mu\f$ integral given \f$beta\f$.
@param thetat dimensionless electron temperature
@param x dimensionless photon energy
 */
double isotropic_nonthermal_calbetainte(const double beta1, const double beta2, const double index, const double x) {
    //double abissca[NGAUSS1], weights[NGAUSS1];
    //double integauss = 0.;
    //gauleg(0., 1., abissca, weights, NGAUSS1);
    //
    //for (int i = 0; i < NGAUSS1; ++i) {
    //    integauss += weights[i] * hotcross_integrand_beta(thetat, abissca[i]) * hotcross_calmuinte(abissca[i], x);
    //}
    //return 0.5 * integauss / thetat / bessel_k2(1. / thetat);
    double gamma1 = 1. / std::sqrt(1. - beta1 * beta1);
    double gamma2 = 1. / std::sqrt(1. - beta2 * beta2);
    std::array<double, 2> paramarr = {index, x};
    void * param = reinterpret_cast<void *>(&paramarr);
    double inte = int_trape1(&isotropic_nonthermal_integrand_beta, beta1, beta2, 1e-3, param);
    
    if ((std::abs(index - 1.)) > 1e-6) {
        return 0.5 * (1. - index) / (std::pow(gamma2, 1. - index) - std::pow(gamma1, 1. - index)) * inte;
    } else {
        return 0.5 / std::log(gamma2 / gamma1) * inte;
    }
}

/*! \brief Return Maxwell Jutter distribution.
@param gamma electron Lorentz factor
@param thetat plasma dimensionless temperature
 */
double maxwell_jutt(const double gamma, const double thetat) {
    double beta = std::sqrt(1. - 1. / gamma / gamma);
    return gamma * gamma * beta * std::exp(-1. * gamma / thetat) / thetat / bessel_k2(1. / thetat);
}

/*! 
\brief Ancillary function for \ref sample_maxjutt.
@param gamma Lorentz factor
@param param structer that contains parameters.
*/
double Ggamma(const double gamma, const Ggammaparam & param) {
    double nom, denom;
    nom = std::exp(-1. * gamma / param.thetat) * (2. * param.thetat * param.thetat +
        2. * param.thetat * gamma + gamma * gamma);
    denom = std::exp(-1. / param.thetat) * (2. * param.thetat * param.thetat +
        2. * param.thetat + 1.);

    return 1. - nom / denom - param.result;
}

/*!
\brief Sampling Maxwell-Juttner distribution, following Schnittman & Krolik 2013. This function is empolyed 
when the electron temperature is greater than 150 keV, see \ref sample_maxjutt1.
//! @param te dimensionless electron temperature
//! @param gen Mersenne-Twister pseudo-random number generator
*/
double sample_maxjutt(const double te, std::mt19937_64 & gen) {
    double fx, gx, gammamax, temp, beta;
    double x0, lambda0, lambda1;
    Ggammaparam params;
    params.thetat = te;

    std::uniform_real_distribution<> dis(0., 1.);

    while(true) {
        lambda0 = dis(gen);
        lambda1 = dis(gen);
        // find upper limit of gamma, for root finding
        gammamax = 2.;
        params.result = lambda0;
        while(Ggamma(gammamax, params) < 0.)
            gammamax *= 2.;
        rtbis<Ggammaparam>(1., gammamax, 1e-6, &Ggamma, params, x0);
        temp = std::exp(-1 * x0 / te);
        beta = std::sqrt(1. - 1. / x0 / x0);
        gx = x0 * x0 * temp;
        fx = beta * gx;
        if (lambda1 < fx / gx)
            break;
    }
    return x0;
}

/*!
\brief Sampling Maxwell-Juttner distribution, with two different methods for electron temperature below and above 150 keV.

If the plasma is hotter than 150 keV, we use \ref sample_maxjutt to sample the distribution. Other we follow
the low temperature method that remains efficient as temperature \f$\rightarrow\f$ 0.
Reference: Pozdnyakov, Sobol, Sunyaev 1983 (PSS83; http://adsabs.harvard.edu/abs/1983ASPRv...2..189P).

@param te electron temperature in electron rest mass
@param gen Mersenne-Twister pseudo-random number generator
*/
double sample_maxjutt1(const double te, std::mt19937_64 & gen) {
    double x0, xi1, xi2, zeta, temp;
    std::uniform_real_distribution<> dis(0., 1.);
    if (te >= 150. / ME_KEV) {
        x0 = sample_maxjutt(te, gen);
    } else {
        while(true) {
            xi1 = dis(gen);
            xi2 = dis(gen);
            zeta = -1.5 * std::log(xi1);
            temp = te * zeta;
            
            if (xi2 * xi2 < 0.151 * (1. + temp) * (1. + temp) * zeta * xi1 * (2. + temp)) {
                x0 = std::sqrt(1. + temp * (2. + temp));
                break;
            }
        }
                
    }
    return x0;
}
