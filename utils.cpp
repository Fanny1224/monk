//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

//! \file utils.cpp

#include <cmath>
#include <random>
#include <complex>
#include <limits>
#include "utils.h"
#include "const.h"

const int NMAX = 23;

/*! \brief Trapezoid integral of a given order. Taken from `sim5` package; re-written to be able to pass additional arguments with a void pointer.
@param pt2func pointer to function to be integrated
@param a lower boundary of integral
@param b upper boundary of integral
@param n order
@param params void pointer to pass additional parameters
@param res result
 */
void int_trapezoid(double(*pt2func)(const double, void *), const double a, const double b, const int n, void * params, double & res) {
    double x, tnm, sum, del;
    int it, j;

    if(n==1){
       res = 0.5 * (b-a) * (pt2func(a, params) + pt2func(a, params));
    } else {
        for(it=1, j=1; j < n-1; j++)
            it <<= 1;
        tnm = (double) it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for(sum=0.0, j=1; j<=it; j++, x+=del) {
            sum += pt2func(x, params); 
        }
        res = 0.5 * (res + del * sum);
    }
}

/*! \brief Trapezoid integral. The function will return the result if 1) the given accurancy is reached; 2) the maximum
order is reached. Taken from `sim5` package; re-written to be able to pass additional arguments with a void pointer.
@param pt2func pointer to function to be integrated
@param a lower boundary of integral
@param b upper boundary of integral
@param acc the desired accurancy
@param params void pointer to pass additional parameters
 */
double int_trape1(double (*pt2func)(const double, void *), const double a, const double b, const double acc, void * params) {
    int n;
    double s = 0.0;
    double olds = -1000.;
    
    for (n = 0; n <= NMAX; ++n) {
        int_trapezoid(pt2func, a, b, n, params, s);
        
        if (n > 3) {
            if (std::abs(s - olds) < acc * std::abs(olds) || ((s == 0.) && (olds == 0.))) 
                break;
        }
        olds = s;
    }
    
    //std::cout << "n = " << n << std::endl;
    return s;
}

/*! \brief Returns the modified Bessel function \f$K_0(x)\f$ for positive real x.

Reference: Press+2007, Numerical Recipes 3rd Edition, Chapter 6.5.1
 */
double bessel_k0(const double x) {
    double z, term;
    if (x <= 1.0) {
        z = x * x;
        term = polyval1(k0pi, z) * std::log(x) / polyval1(k0qi, 1. - z);
        return polyval1(k0p, z) / polyval1(k0q, 1. - z) - term;
    } else {
        z = 1.0 / x;
        return std::exp(-1. * x) * polyval1(k0pp, z) / (polyval1(k0qq, z) * std::sqrt(x));
    }
}

/*! \brief Returns the modified Bessel function \f$K_1(x)\f$ for positive real x.

Reference: Press+2007, Numerical Recipes 3rd Edition, Chapter 6.5.1
 */
double bessel_k1(const double x) {
    double z, term;
    if (x <= 1.0) {
        z = x * x;
        term = polyval1(k1pi, z) * std::log(x) / polyval1(k1qi, 1. - z);
        return x * (polyval1(k1p, z) / polyval1(k1q, 1.- z) + term) + 1. / x;
    } else {
        z = 1.0 / x;
        return std::exp(-x) * polyval1(k1pp, z) / (polyval1(k1qq, z) * std::sqrt(x));
    }
}

/*! \brief Returns the modified Bessel function \f$K_2(x)\f$ for positive real x.

Reference: Press+2007, Numerical Recipes 3rd Edition, Chapter 6.5.1
 */
double bessel_k2(const double x) {
    return (2. / x) * bessel_k1(x) + bessel_k0(x);
}

/*! \brief simple linear interpolation
@param x input
@param xarr x array; (we presume that xarr is evenly spaced and do not check this) we just assume that 
@param yarr y array
*/
double linear_interpolate(const double x, const std::vector<double> & xarr, const std::vector<double> & yarr) {
    double ynow = 0.;

    if (std::abs(x - xarr[0]) <= 1e-8) {
        return yarr[0];
    }
    if (std::abs(x - xarr[xarr.size() - 1]) <= 1e-8) {
        return yarr[xarr.size() - 1];
    }

    
    if ((x - xarr[0]) * (x - xarr[xarr.size() -1]) > 0.) {
        throw("utils::linear_interpolate: x out of range!");
    }
    
    for (unsigned i = 0; i < xarr.size() - 1; ++i) {
        if ((x - xarr[i]) * (x - xarr[i + 1]) <= 0.) {
            ynow = (xarr[i + 1] - x) * yarr[i] + (x - xarr[i]) * yarr[i + 1];
            ynow /= (xarr[i + 1] - xarr[i]);
            break;
        }
    }
    return ynow;
}

double bilinear_interpolate_1d_base (size_t ix0, size_t iy0, const double x, const double y,
    const std::vector<double> & xarr, const std::vector<double> & yarr, const std::vector<double> & zarr) {

    size_t ix1, iy1, index00, index01, index10, index11, ny = yarr.size();
    double z0, z1, znow;

    ix1 = ix0 + 1;
    iy1 = iy0 + 1;
    index00 = ix0 * ny + iy0;
    index01 = ix0 * ny + iy1;
    index10 = ix1 * ny + iy0;
    index11 = ix1 * ny + iy1;
            
    z0 = (xarr[ix1] - x) * zarr[index00] + (x - xarr[ix0]) * zarr[index10];
    z1 = (xarr[ix1] - x) * zarr[index01] + (x - xarr[ix0]) * zarr[index11];
    znow = ((yarr[iy1] - y) * z0 + (y - yarr[iy0]) * z1) / ((xarr[ix1] - xarr[ix0]) * (yarr[iy1] - yarr[iy0]));
    
    return znow;
}

/*!
Bi-linear interpolation. In this function there are a few assumptions (which we do not check):
* Both xarr and yarr are increasing (this function does not require xarr/yarr to be evenly-spaced)
* The zarr is a one-dimensional array (as suggested by `_1d` in the function name), and continuous with y
*/
double bilinear_interpolate_1d_uneven(const double x, const double y, const std::vector<double> & xarr, const std::vector<double> & yarr,
    const std::vector<double> & zarr) {

    size_t ix0, iy0;

    if (std::abs(xarr[0] - x) < 1e-8) {
        ix0 = 0;
    } else if (std::abs(xarr[xarr.size() - 1] - x) < 1e-8) {
        ix0 = xarr.size() - 1;
    } else {
        if (xarr[0] > x || xarr[xarr.size() - 1] < x) {
            std::cout << "x value out of range!" << std::endl;
            throw("x value out of range!");
        }
        ix0 = std::lower_bound(xarr.cbegin(), xarr.cend(), x) - xarr.cbegin();
    }

    if (std::abs(yarr[0] - y) < 1e-8) {
        iy0 = 0;
    } else if (std::abs(yarr[yarr.size() - 1] - y) < 1e-8) {
        iy0 = yarr.size() - 1;
    } else {
        if (yarr[0] > y || yarr[yarr.size() - 1] < y) {
            std::cerr << "Inside bilinear_interpolate_1d_uneven" << std::endl;
            std::cout << "y = " << y << "; value out of range!" << std::endl;
            throw("y value out of range!");
        }
        iy0 = std::lower_bound(yarr.cbegin(), yarr.cend(), y) - yarr.cbegin();
    }

    //std::cout << ix0 << ", " << iy0 << std::endl;
    if (ix0 > 0) ix0 -= 1;
    if (iy0 > 0) iy0 -= 1;

    return bilinear_interpolate_1d_base(ix0, iy0, x, y, xarr, yarr, zarr);
}

/*!
Bi-linear interpolation. In this function there are a few assumptions (which we do not check):
* Both xarr and yarr have regular intervals
* Both xarr and yarr are in increasing order
* The zarr is a one-dimensional array (as suggested by `_1d` in the function name), and continuous with y
*/
double bilinear_interpolate_1d(const double x, const double y, const std::vector<double> & xarr, const std::vector<double> & yarr,
    const std::vector<double> & zarr) {
    
    double xmin = xarr[0], xmax = xarr[xarr.size() - 1], ymin = yarr[0], ymax = yarr[yarr.size() - 1];
    double dx = xarr[1] - xarr[0], dy = yarr[1] - yarr[0];
    long ix0, iy0;
    //, ix1, iy1, index00, index01, index10, index11, ny = yarr.size(); 
    //double z0, z1, znow;
    
    if (x > xmax || x < xmin) {
        throw("utils::bilinear_interpolate(): x out of range!");
    }
    
    if (y > ymax || y < ymin) {
        throw("utils::bilinear_interpolate(): y out of range!");
    }
    
    ix0 = (long)(std::floor((x - xmin) / dx));
    iy0 = (long)(std::floor((y - ymin) / dy));

    return bilinear_interpolate_1d_base(ix0, iy0, x, y, xarr, yarr, zarr);
    
}

/*!
The probability density distribution
\f$p(x) \propto x^\gamma\f$ between \f$x_0\f$ and \f$x_1\f$.
@param gamma powlaw index \f$\gamma\f$
@param x0 lower limit
@param x1 upper limit
@param gen 64-bit MT19937 random number generator
*/
double sample_pl(const double gamma, const double x0, const double x1, std::mt19937_64 & gen) {
    double x, y, temp;
    std::uniform_real_distribution<> dis(0, 1);
    y = dis(gen);
    
    if (std::abs(gamma + 1.) >= 1e-6) {
        temp = (std::pow(x1, gamma + 1.) - std::pow(x0, gamma + 1.)) * y + std::pow(x0, gamma + 1.);
        x = std::pow(temp, 1. / (gamma + 1.));
    } else {
        x = std::exp(y * std::log(x1 / x0)) * x0;
    }
    return x;
}

double sample_cutoffpl(const double alpha, const double emin, const double ecut, std::mt19937_64 & gen) {
    double emax = 10. * ecut;
    std::uniform_real_distribution<> dis(0., 1.);
    double x, xi;
    bool find = false;
    while (!find) {
        x = sample_pl(alpha, emin, emax, gen);
        xi = dis(gen);
        if (xi < std::exp(-1. * x / ecut))
            find = true;
    }
    return x;
}

double sample_expo(const double lambda, std::mt19937_64 & gen) {
    std::uniform_real_distribution<> dis(0., 1.);
    double xi = dis(gen);
    return -1. * lambda * std::log(1. - xi);
}

/*!
\brief Given an unit vector z, calculate one dummy, unit vector x that is perpendicular to z
@param itz the pointer to the vector z
@param itx the pointer to the vector x
*/
void dummy_vec(const double * itz, double * itx) {
    double prod;
    *itx = 1.;
    *(itx + 1) = 0.;
    *(itx + 2) = 0.;

    for (int i = 0; i < 3; ++i)
        *(itx + i) -= *(itz + i) * *(itz);
    prod = *itx  * *itx + *(itx + 1) * *(itx + 1) + *(itx + 2) * *(itx + 2);
    if (std::abs(prod) < DEPS) {
        *itx = 0.;
        *(itx + 1) = 1.;
        *(itx + 2) = 0.;
        for (int i = 0; i < 3; ++i)
            *(itx + i) -= *(itz + i) * *(itz);
        prod = *itx  * *itx + *(itx + 1) * *(itx + 1) + *(itx + 2) * *(itx + 2);
    }
    prod = std::sqrt(prod);
    for (int i = 0; i < 3; ++i)
        *(itx + i) /= prod;
}

void random_dummy_vec(const double * itz, double * itx, std::mt19937_64 & gen) {
    double prod, phi, cosp, sinp;
    std::uniform_real_distribution<> dis(0., 1.);
    std::array<double, 3> xvec = {*(itz + 1) - *(itz + 2), *(itz + 2) - *itz,  *itz - *(itz + 1)}, yvec;
    
    prod = 1. / std::sqrt(xvec[0] * xvec[0] + xvec[1] * xvec[1] + xvec[2] * xvec[2]);
    std::transform(xvec.cbegin(), xvec.cend(), xvec.begin(), std::bind2nd(std::multiplies<double>(), prod));
    
    arrcrossprod(itz, xvec.data(), yvec.data());
    phi = dis(gen) * 2. * M_PI;
    cosp = std::cos(phi);
    sinp = std::sin(phi);
    
    *itx = cosp * xvec[0] + sinp * yvec[0];
    *(itx + 1) = cosp * xvec[1] + sinp * yvec[1];
    *(itx + 2) = cosp * xvec[2] + sinp * yvec[2];
}

/*!
\brief Inverse transformation of \ref forward_rotation.
@param itx pointer to x'-axis unit vector
@param ity pointer to y'-axis unit vector
@param itz pointer to z'-axis unit vector
@param it pointer to vector to be rotated
@param it1 pointer to vector after rotation
*/
void backward_rotation(const double * itx, const double * ity, const double * itz, const double * it, double * it1) {
    *it1 = *(itx) * (*it) + *(ity) * *(it + 1) + *(itz) * *(it + 2);
    *(it1 + 1) = *(itx + 1) * *(it) + *(ity + 1) * *(it + 1) + *(itz + 1) * *(it + 2);
    *(it1 + 2) = *(itx + 2) * *(it) + *(ity + 2) * *(it + 1) + *(itz + 2) * *(it + 2);
}

/*!
\brief Transformation that rotates vector `x`, `y`, `z` to [1, 0, 0], [0, 1, 0], and [0, 0, 1], respectively.

Note that this routine does not check if the three vectors are orthogonal unit vectors; nor handedness.
@param itx pointer to x'-axis unit vector
@param ity pointer to y'-axis unit vector
@param itz pointer to z'-axis unit vector
@param it pointer to vector to be rotated
@param it1 pointer to vector after rotation
*/
void forward_rotation(const double * itx, const double * ity, const double * itz, const double * it, double * it1) {
    *(it1) = *(itx) * *(it) + *(itx + 1) * *(it + 1) + *(itx + 2) * *(it + 2);
    *(it1 + 1) = *(ity) * *(it) + *(ity + 1) * *(it + 1) + *(ity + 2) * *(it + 2);
    *(it1 + 2) = *(itz) * *(it) + *(itz + 1) * *(it + 1) + *(itz + 2) * *(it + 2);
}

/*!
\brief Sample polar and azimuthal angles for isotropic emission.
@param gen Mersenne Twister pseudo-random generator
@param theta output; polar angle
@param ephi output; azimuthal angle
*/
void sample_omega(std::mt19937_64 & gen, double & theta, double & ephi) {
    std::uniform_real_distribution<> dis(0., 1.);
    double mu;
    ephi = 2. * M_PI * dis(gen);
    mu = 2. * dis(gen) - 1.;
    theta = std::acos(mu);
}

/*!
\brief Vector dot product in Minkovski space; return \f$\bm{v}_1 \cdot \bm{v}_2\f$
@param v1 \f$\bm{v}_1\f$
@param v2 \f$\bm{v}_2\f$
*/
double min_dotprod(const std::array<double, 4> & v1, const std::array<double, 4> & v2) {
    return v1[0] * v2[0] - v1[1] * v2[1] - v1[2] * v2[2] - v1[3] * v2[3];
}

/*! 
\brief Generic Lorentz boost
@param nx x-component of electron velocity unit vector
@param ny y-commponent
@param nz z-component
@param beta \f$\beta\equiv v/c\f$
@param gamma Lorentz factor
@param kmu_f four vector in the fluid frame
@param kmu_e output; four vector in the electron frame
*/
double loboost(const double nx, const double ny, const double nz,
    const double beta, const double gamma, const std::array<double, 4> & kmu_f,
    std::array<double, 4> & kmu_e) {

    double kmudotne, klotemp, enboost;
    kmudotne = nx * kmu_f[1] + ny * kmu_f[2] + nz * kmu_f[3];
    klotemp = -1. * beta * gamma * kmu_f[0] + (gamma - 1.) * kmudotne;
    kmu_e[0] = gamma * (kmu_f[0] - beta * kmudotne);
    kmu_e[1] = klotemp * nx + kmu_f[1];
    kmu_e[2] = klotemp * ny + kmu_f[2];
    kmu_e[3] = klotemp * nz + kmu_f[3];
    enboost = kmu_e[0] / kmu_f[0];
    return enboost;
}

    
/*! A very simple read column function. Read a ASCII data file and save the columns in a std::vector<std::vector<double>>.
@param filename path to the ASCII data file
@param vec the 2D vector 
*/
void readcol(std::string filename, std::vector<std::vector<double>> & vec) {
    std::string buffer;
    std::ifstream infile(filename);
    std::vector<std::string> strvec;
    long count = 0;
    
    while(std::getline(infile, buffer)) {
        strvec = str_split(buffer, ' ');
        if (count == 0) {
            vec.resize(strvec.size());
        }
        for (unsigned i = 0; i < vec.size(); ++i) {
            vec[i].push_back(std::stod(strvec[i])); 
        }
           
        ++count;
    }
}

/*! The same with readcol(), but the data are saved in a 1D vector instead of a 2D vector
@param filename path to the ASCII data file
@param vec the 1D vector 
*/
void readcol_1d(std::string filename, std::vector<double> & vec) {
    std::string buffer;
    std::ifstream infile(filename);
    std::vector<std::string> strvec;
    
    while(std::getline(infile, buffer)) {
        strvec = str_split(buffer, ' ');
        for (unsigned i = 0; i < strvec.size(); ++i) {
            vec.push_back(std::stod(strvec[i])); 
        }
    }
}

//! Polynomial evaluation, for std::vector
double polyval(const std::vector<double> coeffs, const double x) {
    size_t n = coeffs.size();
    double res = coeffs[0];
    for (size_t i = 1; i < n; ++i) {
        res = (res * x + coeffs[i]);
    }
    return res;
}

//! Polynomial evaluation, for std::vector
double polyval1(const std::vector<double> coeffs, const double x) {
    size_t n = coeffs.size();
    double res = coeffs[n - 1];
    for (int j = n - 2; j >= 0; --j) {
        res = (res * x + coeffs[j]);
    }
    return res;
}

void mklingrid(const long ntheta, const double thetamin, const double thetamax, std::vector<double> & theta) {
    theta.clear();
    theta.resize(ntheta);
    double dtheta = (thetamax - thetamin) / (double)(ntheta);
    double arr0 = 0.5 * dtheta + thetamin;
    long n = 0;
    std::generate(theta.begin(), theta.end(), [&n, dtheta, arr0]{return (double)(n++) * dtheta + arr0;});
}

void mkloggrid(const long nr, const double rin, const double rout, std::vector<double> & radius) {
    radius.clear();
    radius.resize(nr);
    std::vector<double> logradius(nr);
    double dlogr = std::log(rout / rin) / (double)(nr), logrmin = std::log(rin);

    long n = 0;
    std::generate(logradius.begin(), logradius.end(), [&n, dlogr, logrmin]{return double(n++) * dlogr + logrmin + 0.5 * dlogr;});
    std::transform(logradius.cbegin(), logradius.cend(), radius.begin(), [](double x){return std::exp(x);});
}
