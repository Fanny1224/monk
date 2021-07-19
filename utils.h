//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

//! \file utils.h
//! Utilities
#ifndef _UTILS_H
#define _UTILS_H
#include <iostream>
#include <iterator>
#include <array>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <experimental/filesystem>
#include <random>

//! Print all elements in a std::vector
template <typename T>
std::ostream & operator << (std::ostream & out, const std::vector<T> & vec) {
  out << "[";
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(out, ", "));
  out << "]";
  return out;
};

//! Print all elements in a std::array
template <typename T, size_t n>
std::ostream & operator << (std::ostream & out, const std::array<T, n> & vec) {
  out << "[";
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(out, ", "));
  out << "]";
  return out;
}

//! Multiply std::arrays by a scalar
template <typename T, size_t n>
void mul(const double x, const std::array<T, n> & arr1, std::array<T, n> & arr2) {
    std::transform(arr1.cbegin(), arr1.cend(), arr2.begin(), std::bind2nd(std::multiplies<double>(), x));
}

//! Multiply one std::array by another
template <typename T, size_t n>
void mul(const std::array<T, n> & arr1, const std::array<T, n> & arr2, std::array<T, n> & arr3) {
    std::transform(arr1.cbegin(), arr1.cend(), arr2.cbegin(), arr3.begin(), std::multiplies<double>());
}

//! Array addition
template <typename T, size_t n>
void add(const std::array<T, n> & arr1, const std::array<T, n> & arr2, std::array<T, n> & arr3) {
    std::transform(arr1.cbegin(), arr1.cend(), arr2.cbegin(), arr3.begin(), std::plus<double>());
}

//! Array dot product
template <typename T, size_t n>
double arrdotprod(const std::array<T, n> & arr1, const std::array<T, n> & arr2) {
    double prod = 0.;
    for (size_t i = 0; i < n; ++i)
        prod += arr1[i] * arr2[i];
    return prod;
}

//! Dot product with arrays passed by pointer to fit all kinds of container
template <typename T>
double arrdotprod(size_t n, const T * arr1, const T * arr2) {
    double prod = 0.;
    for (size_t i = 0; i < n; ++i)
        prod += *(arr1 + i) * *(arr2 + i);
    return prod;
}

//! Cross product. Here we limit ourself in cross product in 3 dimension.
template <typename T>
void arrcrossprod(const std::array<T, 3> & u, const std::array<T, 3> & v, std::array<T, 3> & prod) {
    prod[0] = u[1] * v[2] - u[2] * v[1];
    prod[1] = u[2] * v[0] - u[0] * v[2];
    prod[2] = u[0] * v[1] - u[1] * v[0];
}

//! Overloaded `arrcrossprod()` with vectors being passed by pointer
template <typename T>
void arrcrossprod(const T * pu, const T * pv, T * pp) {
    *(pp) = *(pu + 1) * *(pv + 2) - *(pu + 2) * *(pv + 1);
    *(pp + 1) = *(pu + 2) * *(pv) - *(pu) * *(pv + 2);
    *(pp + 2) = *(pu) * *(pv + 1) - *(pu + 1) * *(pv);
}
    
//! Write std::vector to binary file
template <typename T>
void wdoublevec(std::vector<T> & vec, const std::string & path) {
    std::ofstream ofile(path, std::ios::out | std::ofstream::binary);
    vec.push_back(0.);
    std::copy(reinterpret_cast<const char*>(&vec.front()), reinterpret_cast<const char*>(&vec.back()),
        std::ostreambuf_iterator<char>(ofile));
    vec.pop_back();
    ofile.flush();
    ofile.close();
}

//! Append std::vector to binary file
template <typename T>
void wdoublevec_app(std::vector<T> & vec, const std::string & path) {
    std::ofstream ofile(path, std::ios::out | std::ofstream::binary | std::ofstream::app);
    vec.push_back(0.);
    std::copy(reinterpret_cast<const char*>(&vec.front()), reinterpret_cast<const char*>(&vec.back()),
        std::ostreambuf_iterator<char>(ofile));
    vec.pop_back();
    ofile.flush();
    ofile.close();
}

//! Read a whole binary file and save data in a std::vector
template <typename T>
void rdoublevec(std::vector<T> & vec, const std::string & path) {
    std::ifstream ifile(path, std::ios::in | std::ios::binary | std::ios::ate);
    
    ifile.seekg(0, std::ios::end);
    auto fileSize = ifile.tellg();
    ifile.seekg(0, std::ios::beg);
    auto buffer = new char[fileSize];
    ifile.read(buffer, fileSize);
    ifile.close();
    
    T * double_values = (T*)buffer;
    vec =  std::vector<T>(double_values, double_values + (fileSize / sizeof(T)));
    delete[] buffer;
}

/*!
\brief Read a segment of a binary file and save data in a std::vector. This is necessary while reading a large file
@param vec the vector to contain the data
@param path path to binary file
@param rstart the offset, in number of elements (not in bytes)
@param ndata number of elements to be read
*/
template <typename T>
void rdoublevec(std::vector<T> & vec, const std::string & path, size_t rstart, size_t ndata) {
    if (!std::experimental::filesystem::exists(path)) {
        std::cout << path + " does not exist!" << std::endl;
        return;
    }
    size_t ndata_char = (sizeof(T)) / (sizeof(char)) * ndata;
    size_t pos_start = (sizeof(T)) / (sizeof(char)) * rstart;
    size_t pos_end = pos_start + ndata_char;
    std::ifstream ifile(path, std::ios::in | std::ios::binary | std::ios::ate);
    
    auto fileSize = ifile.tellg();
    if (pos_start > fileSize) {
        std::cerr << "Starting point out of range!" << std::endl;
        return;
    }
    
    if (pos_end > fileSize) {
        ndata_char = (size_t)(fileSize) - pos_start;
    }
    
    ifile.seekg(pos_start, std::ios::beg);
    auto buffer = new char[ndata_char];
    ifile.read(buffer, ndata_char);
    ifile.close();
    T * double_values = (T*)buffer;
    vec =  std::vector<T>(double_values, double_values + (ndata_char / sizeof(T)));
    delete[] buffer;
}

/*!
\brief Print progress on screen. Use this function **ONLY** when working in an interactive shell (other than submitting a job).
@param progress fraction of work having been done.
*/
inline void showprogress(double progress) {
    int barWidth = 70, pos;
    std::cout << "[";
    pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

/*!
\brief Find out file size in byte, with `std::string` argument.
*/
inline std::ifstream::pos_type get_filesize(const std::string filename) {
    std::ifstream ifile(filename, std::ios::in | std::ios::binary | std::ios::ate);
    auto fileSize = ifile.tellg();
    ifile.close();
    return fileSize;
}

/*!
\brief Find out file size in byte, with `const char *` argument.
*/
inline std::ifstream::pos_type get_filesize(const char * filename) {
    std::ifstream ifile(filename, std::ios::in | std::ios::binary | std::ios::ate);
    auto fileSize = ifile.tellg();
    ifile.close();
    return fileSize;
}

/*!
\brief Polynomial evaluation.
*/
template <size_t n>
double polyval(const std::array<double, n> & coeffs, const double x) {
    double res = coeffs[0];
    for (size_t i = 1; i < n; ++i) {
        res = (res * x + coeffs[i]);
    }
    return res;
}

/*!
\brief Polynomial evaluation.
*/
template <size_t n>
double polyval1(const std::array<double, n> & coeffs, const double x) {
    double res = coeffs[n - 1];
    for (int j = n - 2; j >= 0; --j) {
        res = (res * x + coeffs[j]);
    }
    return res;
}

double polyval(const std::vector<double> coeffs, const double x);
double polyval1(const std::vector<double> coeffs, const double x);

/*!
\brief Return sign of `x` in double
*/
inline double dsign(const double x) {
    double y;
    y = (x >= 0.) ? 1 : -1;
    if (x == 0.) y = 0.;
    return y;
}

inline std::vector<std::string> str_split(const std::string & s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

void readcol(std::string filename, std::vector<std::vector<double>> & vec);
void readcol_1d(std::string filename, std::vector<double> & vec);

double int_trape1(double (*pt2func)(const double, void *), const double a, const double b, const double acc, void * params);
void int_trapezoid(double(*pt2func)(const double, void *), const double a, const double b, const int n, void * params, double & res);
double bessel_k0(const double x);
double bessel_k1(const double x);
double bessel_k2(const double x);
double linear_interpolate(const double x, const std::vector<double> & xarr, const std::vector<double> & yarr);
double bilinear_interpolate_1d_uneven(const double x, const double y, const std::vector<double> & xarr, const std::vector<double> & yarr,
    const std::vector<double> & zarr);
double bilinear_interpolate_1d(const double x, const double y, const std::vector<double> & xarr, const std::vector<double> & yarr,
    const std::vector<double> & zarr);
//! Sample powerlaw distribution.
double sample_pl(const double gamma, const double x0, const double x1, std::mt19937_64 & gen);
double sample_cutoffpl(const double alpha, const double emin, const double ecut, std::mt19937_64 & gen);
double sample_expo(const double lambda, std::mt19937_64 & gen);

void dummy_vec(const double * itz, double * itx);
void random_dummy_vec(const double * itz, double * itx, std::mt19937_64 & gen);
double loboost(const double nx, const double ny, const double nz,
    const double beta, const double gamma, const std::array<double, 4> & kmu_f,
    std::array<double, 4> & kmu_e);
void forward_rotation(const double * itx, const double * ity, const double * itz, const double * it, double * it1);
void backward_rotation(const double * itx, const double * ity, const double * itz, const double * it, double * it1);
void myslowdft(const std::vector<double> & xt, std::vector<double> & res);
void mysimplefft(std::vector<double> & data);
double expint(const double x, const int n);
void mklingrid(const long ntheta, const double thetamin, const double thetamax, std::vector<double> & theta);
void mkloggrid(const long ntheta, const double thetamin, const double thetamax, std::vector<double> & theta);

const std::array<double, 5> k0pi = {1.0, 2.346487949187396e-1, 1.187082088663404e-2, 2.150707366040937e-4, 1.425433617130587e-6};
const std::array<double, 3> k0qi = {9.847324170755358e-1, 1.518396076767770e-2, 8.362215678646257e-5};
const std::array<double, 5> k0p = {1.159315156584126e-1, 2.770731240515333e-1, 2.066458134619875e-2, 4.574734709978264e-4, 3.454715527986737e-6};
const std::array<double, 3> k0q = {9.836249671709183e-1, 1.627693622304549e-2, 9.809660603621949e-5};
const std::array<double, 8> k0pp = {1.253314137315499, 1.475731032429900e1,  6.123767403223466e1, 1.121012633939949e2,
    9.285288485892228e1, 3.198289277679660e1, 3.595376024148513, 6.160228690102976e-2};
const std::array<double, 8> k0qq = {1.0, 1.189963006673403e1, 5.027773590829784e1, 9.496513373427093e1,
    8.318077493230258e1, 3.181399777449301e1, 4.443672926432041, 1.408295601966600e-1};
const std::array<double, 5> k1pi = {0.5, 5.598072040178741e-2, 1.818666382168295e-3, 2.397509908859959e-5, 1.239567816344855e-7};
const std::array<double, 3> k1qi = {9.870202601341150e-1, 1.292092053534579e-2, 5.881933053917096e-5};
const std::array<double,5> k1p = {-3.079657578292062e-1, -8.109417631822442e-2, -3.477550948593604e-3, -5.385594871975406e-5, -3.110372465429008e-7};
const std::array<double, 3> k1q = {9.861813171751389e-1, 1.375094061153160e-2, 6.774221332947002e-5};
const std::array<double, 8> k1pp = {1.253314137315502, 1.457171340220454e1, 6.063161173098803e1, 1.147386690867892e2,
    1.040442011439181e2,  4.356596656837691e1, 7.265230396353690, 3.144418558991021e-1};
const std::array<double, 8> k1qq = {1.0, 1.125154514806458e1, 4.427488496597630e1, 7.616113213117645e1, 5.863377227890893e1,
    1.850303673841586e1, 1.857244676566022, 2.538540887654872e-2};
    
#endif
