

#ifndef _SIM5CONFIG_H
#define _SIM5CONFIG_H




#ifdef __CUDACC__
    #define CUDA

    #define DEVICEFUNC   __device__
    #define HOSTFUNC     __host__
    #define INLINE
#else
    #define DEVICEFUNC
    #define HOSTFUNC
    #define INLINE       
#endif


#endif
#ifndef _SIM5INCLUDE_H
#define _SIM5INCLUDE_H

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stddef.h>
#include <string.h>
#include <float.h>
#include <time.h>
#ifndef CUDA
#include <complex.h>
#endif

#endif
/* -------------------------------------------------------------- */

/*
Mersenne Twister pseudorandom number generator
--
NOTE: This unit redistributes code for 64bit Mersenne Twister
pseudorandom number generator.by Nishimura and Matsumoto
No changes has been made to the original code except renaming
the routines (adding mt19937 prefix) and introducing thread_local
specifier to globally declared static variables (requires compiler
to implement C11 standard)
*/

#ifndef __CUDACC__

/*
   The original copyright notice:

   A C-program for MT19937-64 (2004/9/29 version).
   Coded by Takuji Nishimura and Makoto Matsumoto.

   This is a 64-bit version of Mersenne Twister pseudorandom number
   generator.

   Before using, initialize the state by using init_genrand64(seed)
   or init_by_array64(init_key, key_length).

   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   References:
   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
     ACM Transactions on Modeling and
     Computer Simulation 10. (2000) 348--357.
   M. Matsumoto and T. Nishimura,
     ``Mersenne Twister: a 623-dimensionally equidistributed
       uniform pseudorandom number generator''
     ACM Transactions on Modeling and
     Computer Simulation 8. (Jan. 1998) 3--30.

   Any feedback is very welcome.
   http:
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
*/


/* initializes mt[NN] with a seed */
void mt19937_init(unsigned long long seed);

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void mt19937_init_by_array64(unsigned long long init_key[],
		     unsigned long long key_length);

/* generates a random number on [0, 2^64-1]-interval */
unsigned long long mt19937_int64(void);


/* generates a random number on [0, 2^63-1]-interval */
long long mt19937_int63(void);

/* generates a random number on [0,1]-real-interval */
double mt19937_real1(void);

/* generates a random number on [0,1)-real-interval */
double mt19937_real2(void);

/* generates a random number on (0,1)-real-interval */
double mt19937_real3(void);

#endif 

/* -------------------------------------------------------------- */



#ifndef _SIM5CONST_H
#define _SIM5CONST_H

#define TINY                1e-40                

#define TRUE                1                    
#define FALSE               0                    


#define grav_radius         1.476716e+05         
#define speed_of_light      2.997925e+10         
#define speed_of_light2     8.987554e+20         
#define boltzmann_k         1.380650e-16         
#define sb_sigma            5.670400e-05         
#define sigma_thomson       6.652458e-25         
#define parsec              3.085680e+18         
#define mass_proton         1.672622e-24         
#define mass_electron       9.109382e-28         
#define solar_mass          1.988920e+33         
#define grav_const          6.673000e-08         
#define planck_h            6.626069e-27         
#define atomic_mass_unit    1.660539e-24         
#define avogadro_number     6.022141e+23         


#define Mdot_Edd        2.225475942e+18          
#define L_Edd           1.257142540e+38          




#define si_grav_radius      1.476716e+03         
#define si_speed_of_light   2.997925e+08         
#define si_speed_of_light2  8.987554e+16         
#define si_boltzmann_k      1.380650e-23         
#define si_sb_sigma         5.670400e-08         
#define si_electronvolt     1.602177e-19         
#define si_parsec           3.085680e+16         
#define si_grav_const       6.673000e-11         
#define si_erg              1.000000e-07         
#define si_solar_mass       1.988920e+30         
#define si_angstrom         1.000000e-10         
#define si_sigma_T          6.652459e-29         
#define si_mass_proton      1.672622e-27         
#define si_mass_electron    9.109382e-31         
#define si_planck_h         6.626069e-34         


#define gu_sb_sigma         1.562539e-60         

#define erg2kev             6.241507e+08         
#define kev2erg             1.602177e-09         
#define joule2kev           6.241507e+15         
#define joule2erg           1.000000e+07         
#define erg2joule           1.000000e-07         
#define kev2joule           1.602177e-16         
#define freq2kev            4.135667e-18         
#define freq2erg            6.626069e-27         
#define kev2freq            2.417990e+17         
#define msq2cmsq            1.000000e+04         
#define cmsq2msq            1.000000e-04         
#define kelvin2kev          8.617342e-08         
#define kev2kelvin          1.160451e+07         
#define m2cm                1.000000e+02         
#define cm2m                1.000000e-02         
#define kev2ev              1.000000e+03         
#define ev2kev              1.000000e-03         


#endif




#ifndef _SIM5MATH_H
#define _SIM5MATH_H

#ifndef CUDA
#include <complex.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif




#define PI      3.14159265359                                   
#define PI2     6.28318530718                                   
#define PI4     12.5663706144                                   
#define PI_half 1.57079632679                                   

#define sqr(a)   ((a) * (a))                                    
#define sqr2(a)  ((a) * (a))                                    
#define sqr3(a)  ((a) * (a) * (a))                              
#define sqr4(a)  ((a) * (a) * (a) * (a))                        
#define sqrt3(a) cbrt(a)                                        
#define sqrt4(a) pow(a,0.25)                                    
#define odd(a) ((a%2==1)?1:0)                                   
#define sign(a) ((a) >= 0.0 ? (+1.0) : (-1.0))                  
#define deg2rad(a) ((a)/180.0*M_PI)                             
#define rad2deg(a) ((a)*180.0/M_PI)                             
#define EE(a) pow(10.0,a)                                       
#define ave(a, b, w) ((1.0-(w))*(a) + (w)*(b))                  
#define logave(a, b, w) (exp((1.0-(w))*log(a) + (w)*log(b)))    
#define inrange(a, min, max) (((a)>=(min))&&((a)<=(max)))       



#ifdef CUDA
    #define ComplexI    makeComplex(0.0,1.0)                    
    typedef double2     sim5complex;                            
#else
    #define ComplexI _Complex_I                                 
    typedef double _Complex sim5complex;                         
#endif


DEVICEFUNC INLINE long sim5round(double num);

DEVICEFUNC INLINE long int factorial(long int n);

DEVICEFUNC INLINE double reduce_angle_pi(double phi);
DEVICEFUNC INLINE double reduce_angle_2pi(double phi);

DEVICEFUNC int ensure_range(double *val, double min, double max, double acc);


DEVICEFUNC void cartesian2spherical1(double x, double y, double z, double Vx, double Vy, double Vz, double* Vr, double* Vh, double* Vf);
DEVICEFUNC void cartesian2spherical2(double cos_h, double sin_f, double cos_f, double Vx, double Vy, double Vz, double* Vr, double* Vh, double* Vf);


DEVICEFUNC INLINE void sim5seed();
DEVICEFUNC INLINE unsigned long long sim5rand();
DEVICEFUNC INLINE double sim5urand();


DEVICEFUNC INLINE sim5complex makeComplex(double r, double i);
DEVICEFUNC INLINE sim5complex nullComplex();
DEVICEFUNC INLINE double sim5creal (sim5complex a);
DEVICEFUNC INLINE double sim5cimag (sim5complex a);


#ifdef __cplusplus
}
#endif


#endif
#ifndef _SIM5UTILS_H_
#define _SIM5UTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

void gprintf(FILE* file, const char *templatex, ...);

void error(const char *templatex, ...);
void warning(const char *templatex, ...);

void sort_array(double *array, int N);
void sort_array_f(float *array, int N);

#ifndef CUDA
void* sim5_array_alloc(size_t capacity, size_t element_size);
void sim5_array_free(void* array);
void* sim5_array_realloc(void* array, size_t new_capacity);
long sim5_array_count(void* arrry);
long sim5_array_capa(void* array);
size_t sim5_array_esize(void* array);
void sim5_array_push(void** array_ptr, const void* data);
void sim5_array_push_int(void** array_ptr, const int data);
void sim5_array_push_long(void** array_ptr, const long data);
void sim5_array_push_double(void** array_ptr, const double data);
void sim5_array_push_float(void** array_ptr, const float data);
int  sim5_array_exists(void* array, const void* data);
void sim5_array_push_if_not_exists(void** array_ptr, const void* data);
void sim5_array_reverse(void* array);
#endif

char* key_value_get(const char *string, const char *key);

void backtrace();

#ifdef __cplusplus
}
#endif


#endif



#ifndef _SIM5INTEGRATION_H
#define _SIM5INTEGRATION_H

#ifdef __cplusplus
extern "C" {
#endif

DEVICEFUNC double integrate_trapezoid(double(*f)(double), double a, double b, double acc);
DEVICEFUNC double integrate_simpson(double (*f)(double), double a, double b, double acc);

DEVICEFUNC void gauleg(double x1, double x2, double x[], double w[], int n);

#ifdef __cplusplus
}
#endif

#endif



#ifndef _SIM5INTERPOLATION_H
#define _SIM5INTERPOLATION_H

#ifndef CUDA

#ifdef __cplusplus
extern "C" {
#endif

#define INTERP_DATA_REF                 0       
#define INTERP_DATA_COPY                1       
#define INTERP_DATA_BUILD               2       

#define INTERP_OPT_ACCEL                1       
#define INTERP_OPT_CAN_EXTRAPOLATE      2       

#define INTERP_TYPE_LINLIN              0       
#define INTERP_TYPE_LINLOG              1       
#define INTERP_TYPE_LOGLIN              2       
#define INTERP_TYPE_LOGLOG              3       
#define INTERP_TYPE_SPLINE              4       


typedef struct sim5interp {
    long    N;                  
    long    capa;               
    double* X;                  
    double* Y;                  
    double* d2Y;                
    int     datamodel;          
    int     type;               
    int     options;            
    double  xmin;               
    double  xmax;               

    long last_index;            
} sim5interp;


DEVICEFUNC sim5interp* sim5_interp_alloc();
DEVICEFUNC void sim5_interp_init(sim5interp* interp, double xa[], double ya[], long N, int data_model, int interp_type, int interp_options);
DEVICEFUNC void sim5_interp_data_push(sim5interp* interp, double x, double y);
DEVICEFUNC double sim5_interp_eval(sim5interp* interp, double x);
DEVICEFUNC void sim5_interp_done(sim5interp* interp);
DEVICEFUNC void sim5_interp_free(sim5interp* interp);

DEVICEFUNC INLINE long sim5_interp_search(const double x_array[], double x, long index_lo, long index_hi);

#ifdef __cplusplus
}
#endif

#endif 

#endif




#ifndef _SIM5DISTRIBUTIONS_H
#define _SIM5DISTRIBUTIONS_H

#ifndef CUDA

#ifdef __cplusplus
extern "C" {
#endif

typedef struct sim5distrib {
    double x_min;
    double x_max;
    double norm;
    sim5interp pdf;                     
    sim5interp cdf;                     
    sim5interp icd;                     
} sim5distrib;



DEVICEFUNC void distrib_init(sim5distrib* d, double(*pdf)(double), double x_min, double x_max, int N);
DEVICEFUNC void distrib_done(sim5distrib* d);
DEVICEFUNC INLINE double distrib_hit(sim5distrib* d);

#ifdef __cplusplus
}
#endif

#endif 

#endif

#ifndef _SIM5ROOTFINDING_H
#define _SIM5ROOTFINDING_H

#ifdef __cplusplus
extern "C" {
#endif

long rtbis(double x1, double x2, double xacc, double (*fx)(double), double* result);

#ifdef __cplusplus
}
#endif

#endif

#ifndef _SIM5ELLIPTIC_H
#define _SIM5ELLIPTIC_H


#ifdef __cplusplus
extern "C" {
#endif

DEVICEFUNC double rf(double x, double y, double z);
DEVICEFUNC double rd(double x, double y, double z);
DEVICEFUNC double rc(double x, double y);
DEVICEFUNC double rj(double x, double y, double z, double p);
DEVICEFUNC INLINE double elliptic_k(double m);
DEVICEFUNC double elliptic_e_sin(double sin_phi, double m);
DEVICEFUNC double elliptic_e_cos(double cos_phi, double m);
DEVICEFUNC double elliptic_f(double phi, double m);
DEVICEFUNC double elliptic_f_cos(double cos_phi, double m);
DEVICEFUNC double elliptic_f_sin(double sin_phi, double m);
DEVICEFUNC double elliptic_pi_complete(double n, double m);
DEVICEFUNC double elliptic_pi_cos(double cos_phi, double n, double m);
DEVICEFUNC double elliptic_pi_sin(double sin_phi, double n, double m);
DEVICEFUNC sim5complex elliptic_pi(double phi, double n, double m);
DEVICEFUNC INLINE double jacobi_isn(double z, double m);
DEVICEFUNC INLINE double jacobi_icn(double z, double m);
DEVICEFUNC INLINE double jacobi_itn(double z, double m);
DEVICEFUNC void  jacobi_sncndn(double u, double m, double *sn, double *cn, double *dn);
DEVICEFUNC INLINE double jacobi_sn(double u, double m);
DEVICEFUNC INLINE double jacobi_cn(double u, double m);
DEVICEFUNC INLINE double jacobi_dn(double u, double m);

DEVICEFUNC double integral_R_r0_re(double a, double b, double c, double d, double X);
DEVICEFUNC double integral_R_r0_cc(double a, double b, sim5complex c, double X);
DEVICEFUNC double integral_R_r0_re_inf(double a, double b, double c, double d);
DEVICEFUNC double integral_R_r0_cc_inf(double a, double b, sim5complex c);
DEVICEFUNC double integral_R_r1_re(double a, double b, double c, double d, double X);
DEVICEFUNC double integral_R_r1_cc(double a, double b, sim5complex c, double X1, double X2);
DEVICEFUNC double integral_R_r2_re(double a, double b, double c, double d, double X);
DEVICEFUNC double integral_R_r2_cc(double a, double b, sim5complex c, double X1, double X2);
DEVICEFUNC double integral_R_rp_re(double a, double b, double c, double d, double p, double X);
DEVICEFUNC double integral_R_rp_cc2(double a, double b, sim5complex c, double p, double X1, double X2);
DEVICEFUNC double integral_R_rp_re_inf(double a, double b, double c, double d, double p);
DEVICEFUNC double integral_R_rp_cc2_inf(double a, double b, sim5complex c, double p, double X1);
DEVICEFUNC double integral_T_m0(double a2, double b2, double X);
DEVICEFUNC double integral_T_m2(double a2, double b2, double X);
DEVICEFUNC double integral_T_mp(double a2, double b2, double p, double X);

#ifdef __cplusplus
}
#endif

#endif

#ifndef _SIM5POLYROOTS_H
#define _SIM5POLYROOTS_H

#ifdef __cplusplus
extern "C" {
#endif

DEVICEFUNC int quadratic_eq(double pr, double pi, double qr, double qi, double zr[2], double zi[2]);
DEVICEFUNC int cubic_eq(double p, double q, double r, double zr[3], double zi[3]);
DEVICEFUNC int quartic_eq(double a3, double a2, double a1, double a0, double zr[4], double zi[4]);
DEVICEFUNC void sort_roots_re(double *r1, double *r2, double *r3, double *r4);
DEVICEFUNC void sort_mix(double *r1, double *r2, int *s);
DEVICEFUNC void sort_mix2(double *r1, double *r2, int *s);
DEVICEFUNC void sort_roots(int *s, sim5complex *z1, sim5complex *z2, sim5complex *z3, sim5complex *z4);

DEVICEFUNC void quartic_eq_c(double a3, double a2, double a1, double a0, int *nr, sim5complex *z1, sim5complex *z2, sim5complex *z3, sim5complex *z4);

#ifdef __cplusplus
}
#endif

#endif


#ifndef _SIM5_RAYTRACE_H
#define _SIM5_RAYTRACE_H

#ifdef __cplusplus
extern "C" {
#endif


#define RTOPT_NONE              0         
#define RTOPT_FLAT              1         
#define RTOPT_POLARIZATION      2         


typedef struct raytrace_data {
    int opt_gr;             
    int opt_pol;            
    double step_epsilon;    

    double bh_spin;         
    double E;               
    double Q;               
    sim5complex WP;         

    int pass;               
    int refines;            
    double dk[4];           
    double df[4];           
    double kt;              
    float error;            
} raytrace_data;


DEVICEFUNC
void raytrace_prepare(double bh_spin, double x[4], double k[4], double presision_factor, int options, raytrace_data* rtd);

DEVICEFUNC
void raytrace(double x[4], double k[4], double *step, raytrace_data* rtd);

DEVICEFUNC
double raytrace_error(double x[4], double k[4], raytrace_data* rtd);

#ifdef __cplusplus
}
#endif


#endif



#ifndef _SIM5KERR_H
#define _SIM5KERR_H

#ifdef __cplusplus
extern "C" {
#endif

struct sim5metric {
    double a,r,m;
    double g00;
    double g11;
    double g22;
    double g33;
    double g03;
};
typedef struct sim5metric sim5metric;

struct sim5tetrad {
    double e[4][4];
    sim5metric metric;
};
typedef struct sim5tetrad sim5tetrad;



DEVICEFUNC
void flat_metric(double r, double m, sim5metric *metric);

DEVICEFUNC
void flat_metric_contravariant(double r, double m, sim5metric *metric);

DEVICEFUNC
void kerr_metric(double a, double r, double m, sim5metric *metric);

DEVICEFUNC
void kerr_metric_contravariant(double a, double r, double m, sim5metric *metric);

DEVICEFUNC
void flat_connection(double r, double m, double G[4][4][4]);

DEVICEFUNC
void kerr_connection(double a, double r, double m, double G[4][4][4]);

DEVICEFUNC INLINE
void Gamma(double G[4][4][4], double U[4], double V[4], double result[4]);

DEVICEFUNC INLINE
void vector_set(double x[4], double x0, double x1, double x2, double x3);

DEVICEFUNC INLINE
void vector_copy(double src[4], double dst[4]);

DEVICEFUNC INLINE
void vector_covariant(double V1[4], double V2[4], sim5metric* m);

DEVICEFUNC INLINE
double vector_norm(double V[4], sim5metric* m);

DEVICEFUNC INLINE
double vector_3norm(double V[4]);

DEVICEFUNC INLINE
void vector_norm_to(double V[4], double norm, sim5metric* m);

DEVICEFUNC
void vector_norm_to_null(double V[4], double V0, sim5metric* m);

DEVICEFUNC INLINE
void vector_multiply(double V[4], double factor);

DEVICEFUNC INLINE
double dotprod(double V1[4], double V2[4], sim5metric* m);

DEVICEFUNC
void tetrad_general(sim5metric *m, double U[], sim5tetrad *t);

DEVICEFUNC
void tetrad_zamo(sim5metric *m, sim5tetrad *t);

DEVICEFUNC
void tetrad_radial(sim5metric *m, double v_r, sim5tetrad *t);

DEVICEFUNC
void tetrad_azimuthal(sim5metric *m, double Omega, sim5tetrad *t);

DEVICEFUNC
void tetrad_surface(sim5metric *m, double Omega, double vr, double dhdr, sim5tetrad *t);

DEVICEFUNC
void bl2on(double Vin[4], double Vout[4], sim5tetrad* t);

DEVICEFUNC
void on2bl(double Vin[4], double Vout[4], sim5tetrad* t);




DEVICEFUNC INLINE
double r_bh(double a);

DEVICEFUNC INLINE
double r_ms(double a);

DEVICEFUNC INLINE
double r_mb(double a);

DEVICEFUNC INLINE
double r_ph(double a);

DEVICEFUNC INLINE
double OmegaK(double r, double a);

DEVICEFUNC INLINE
double ellK(double r, double a);

DEVICEFUNC INLINE
double omega_r(double r, double a);

DEVICEFUNC INLINE
double omega_z(double r, double a);

DEVICEFUNC INLINE
double Omega_from_ell(double ell, sim5metric *m);

DEVICEFUNC INLINE
double ell_from_Omega(double Omega, sim5metric *m);

DEVICEFUNC INLINE
double gfactorK(double r, double a, double l);

DEVICEFUNC
void photon_momentum(double a, double r, double m, double l, double q, double r_sign, double m_sign, double k[4]);

DEVICEFUNC
void photon_motion_constants(double a, double r, double m, double k[4], double* L, double* Q);

DEVICEFUNC
double photon_carter_const(double k[4], sim5metric *metric);

DEVICEFUNC INLINE
void fourvelocity_zamo(sim5metric *m, double U[4]);

DEVICEFUNC INLINE
void fourvelocity_azimuthal(double Omega, sim5metric *m, double U[4]);

DEVICEFUNC INLINE
void fourvelocity_radial(double vr, sim5metric *m, double U[4]);

DEVICEFUNC INLINE
double fourvelocity_norm(double U1, double U2, double U3, sim5metric *m);

DEVICEFUNC
void fourvelocity(double U1, double U2, double U3, sim5metric *m, double U[]);

#ifdef __cplusplus
}
#endif

#endif


#ifndef _SIM5KERR_GEOD_H
#define _SIM5KERR_GEOD_H

#ifdef __cplusplus
extern "C" {
#endif

#define GEOD_TYPE_RR               40       
#define GEOD_TYPE_RR_DBL           41       
#define GEOD_TYPE_RR_BH            42       
#define GEOD_TYPE_RC                2       
#define GEOD_TYPE_CC                0       


#define GD_OK                           0
#define GD_ERROR_Q_ZERO                 1
#define GD_ERROR_BOUND_GEODESIC         2
#define GD_ERROR_UNKNOWN_SOLUTION       3
#define GD_ERROR_TYPE_RR_DOUBLE         4
#define GD_ERROR_TYPE_CC                5
#define GD_ERROR_Q_RANGE                7
#define GD_ERROR_MUPLUS_RANGE           8
#define GD_ERROR_MU0_RANGE              9
#define GD_ERROR_MM_RANGE              10
#define GD_ERROR_INCL_RANGE            11
#define GD_ERROR_SPIN_RANGE            12



typedef struct geodesic {
	double a;                                
	double alpha;                            
	double beta;                             
	double incl;                            
	double cos_i;                            

	double l;                                
	double q;                                
	sim5complex r1,r2,r3,r4;                 
	int    nrr;                              
	int    type;                             
	double m2p,m2m,mm,mK;                    

	double rp;                               
	double dmdp_inf;                         

	double Rpc;                              
	double Tpp;                              
	double Tip;                              

	double k[4];                             
	double p;                                
} geodesic;


DEVICEFUNC int geodesic_init_inf(double i, double a, double alpha, double beta, geodesic *g, int *error);
DEVICEFUNC int geodesic_init_src(double a, double r, double m, double k[4], int bpa, geodesic *g, int *error);
DEVICEFUNC double geodesic_P_int(geodesic *g, double r, int bpa);
DEVICEFUNC void geodesic_position(geodesic *g, double P, double x[4]);
DEVICEFUNC double geodesic_position_rad(geodesic *g, double P);
DEVICEFUNC double geodesic_position_pol(geodesic *g, double P);
DEVICEFUNC double geodesic_position_pol_sign_k_theta(geodesic *g, double P);
DEVICEFUNC double geodesic_position_azm(geodesic *g, double r, double m, double P);
DEVICEFUNC double geodesic_P_int(geodesic *g, double r, int bpa);
DEVICEFUNC double geodesic_dm_sign(geodesic *g, double P);
DEVICEFUNC void geodesic_momentum(geodesic *g, double P, double r, double m, double k[]);
DEVICEFUNC double geodesic_find_midplane_crossing(geodesic *g, int order);
DEVICEFUNC void geodesic_follow(geodesic *g, double step, double *P, double *r, double *m, int *status);
DEVICEFUNC double geodesic_timedelay(geodesic *g, double P1, double r1, double m1, double P2, double r2, double m2);

#ifdef __cplusplus
}
#endif

#endif


#ifndef _SIM5DISKNT_H
#define _SIM5DISKNT_H

#ifndef CUDA

#define DISK_NT_OPTION_LUMINOSITY     1

#ifdef __cplusplus
extern "C" {
#endif

DEVICEFUNC int disk_nt_setup(double M, double a, double mdot_or_L, double alpha, int options);
DEVICEFUNC void disk_nt_done();
DEVICEFUNC double disk_nt_r_min();
DEVICEFUNC double disk_nt_flux(double r);
DEVICEFUNC double disk_nt_lumi();
DEVICEFUNC double disk_nt_mdot();
DEVICEFUNC double disk_nt_sigma(double r);
DEVICEFUNC double disk_nt_ell(double r);
DEVICEFUNC double disk_nt_vr(double r);
DEVICEFUNC double disk_nt_h(double r);
DEVICEFUNC double disk_nt_dhdr(double r);
DEVICEFUNC void disk_nt_dump(char* filename);

#ifdef __cplusplus
}
#endif

#endif 

#endif





#ifndef _SIM5_POLARIZATION_H
#define _SIM5_POLARIZATION_H

#ifdef __cplusplus
extern "C" {
#endif

DEVICEFUNC void polarization_vector(double k[4], sim5complex wp, sim5metric *metric, double f[4]);
DEVICEFUNC sim5complex polarization_constant(double k[4], double f[4], sim5metric *metric);
DEVICEFUNC sim5complex polarization_constant_infinity(double a, double alpha, double beta, double incl);

DEVICEFUNC double polarization_angle_rotation(double a, double inc, double alpha, double beta, sim5complex kappa);

#ifdef __cplusplus
}
#endif


#endif


#ifndef _SIM5RADIATION_H
#define _SIM5RADIATION_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct stokes_params {
  double i;
  double q;
  double u;
  double v;
  double tau;         
} stokes_params;


static const stokes_params stokes_null = {0.0, 0.0, 0.0, 0.0, 0.0};



DEVICEFUNC double blackbody_Iv(double T, double hardf, double cos_mu, double E);
DEVICEFUNC void blackbody(double T, double hardf, double cos_mu, double E[], double Iv[], int en_bins);
DEVICEFUNC INLINE double blackbody_photons(double T, double hardf, double cos_mu, double E);
DEVICEFUNC double blackbody_photons_total(double T, double hardf);
DEVICEFUNC double blackbody_photon_energy_random(double T);

#ifdef __cplusplus
}
#endif

#endif