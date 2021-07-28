#include "sim5lib.h"
/* -------------------------------------------------------------- */

/*
Mersenne Twister pseudorandom number generator
--
NOTE: This unit redistributes code for 64bit Mersenne Twister
pseudorandom number generator.by Nishimura and Matsumoto.
No changes has been made to the original code except renaming
the routines (adding mt19937 prefix) and introducing thread_local
specifier to globally declared static variables (requires compiler
to implement C11 standard).
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



#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */


/* The array for the state vector */
static unsigned long long mt[NN];
/* mti==NN+1 means mt[NN] is not initialized */
static int mti=NN+1;

/* initializes mt[NN] with a seed */
void mt19937_init(unsigned long long seed)
{
 mt[0] = seed;
 for (mti=1; mti<NN; mti++)
     mt[mti] =  (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void mt19937_init_by_array64(unsigned long long init_key[],
		  unsigned long long key_length)
{
 unsigned long long i, j, k;
 mt19937_init(19650218ULL);
 i=1; j=0;
 k = (NN>key_length ? NN : key_length);
 for (; k; k--) {
     mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 3935559000370003845ULL))
       + init_key[j] + j; /* non linear */
     i++; j++;
     if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
     if (j>=key_length) j=0;
 }
 for (k=NN-1; k; k--) {
     mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 2862933555777941757ULL))
       - i; /* non linear */
     i++;
     if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
 }

 mt[0] = 1ULL << 63; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0, 2^64-1]-interval */
unsigned long long mt19937_int64(void)
{
 int i;
 unsigned long long x;
 static unsigned long long mag01[2]={0ULL, MATRIX_A};

 if (mti >= NN) { /* generate NN words at one time */

     /* if mt19937_init() has not been called, */
     /* a default initial seed is used     */
     if (mti == NN+1)
         mt19937_init(5489ULL);

     for (i=0;i<NN-MM;i++) {
         x = (mt[i]&UM)|(mt[i+1]&LM);
         mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
     }
     for (;i<NN-1;i++) {
         x = (mt[i]&UM)|(mt[i+1]&LM);
         mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
     }
     x = (mt[NN-1]&UM)|(mt[0]&LM);
     mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];

     mti = 0;
 }

 x = mt[mti++];

 x ^= (x >> 29) & 0x5555555555555555ULL;
 x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
 x ^= (x << 37) & 0xFFF7EEE000000000ULL;
 x ^= (x >> 43);

 return x;
}

/* generates a random number on [0, 2^63-1]-interval */
long long mt19937_int63(void)
{
 return (long long)(mt19937_int64() >> 1);
}

/* generates a random number on [0,1]-real-interval */
double mt19937_real1(void)
{
 return (mt19937_int64() >> 11) * (1.0/9007199254740991.0);
}

/* generates a random number on [0,1)-real-interval */
double mt19937_real2(void)
{
 return (mt19937_int64() >> 11) * (1.0/9007199254740992.0);
}

/* generates a random number on (0,1)-real-interval */
double mt19937_real3(void)
{
 return ((mt19937_int64() >> 12) + 0.5) * (1.0/4503599627370496.0);
}


#endif 

/* -------------------------------------------------------------- */



DEVICEFUNC INLINE
long sim5round(double num) {
 return (long)(num+0.5);
}


DEVICEFUNC INLINE
long int factorial(long int n)
{
	if (n<=1) return(1);	else n=n*factorial(n-1);
	return(n);
}


DEVICEFUNC INLINE
double reduce_angle_pi(double phi)
{
 while (phi < 0.0)  phi += 2.*PI;
 while (phi > PI) phi -= PI;
 return phi;
}


DEVICEFUNC INLINE
double reduce_angle_2pi(double phi)
{
 while (phi >= +PI2) phi -= PI2;
 while (phi <     0.0) phi += PI2;
 return phi;
}


DEVICEFUNC
int ensure_range(double *val, double min, double max, double acc)
{
	if (*val<min-acc) return 0;
	if (*val>max+acc) return 0;

	if (*val<min) *val = min;
	if (*val>max) *val = max;
	return 1;
}



DEVICEFUNC INLINE
void sim5seed()
{
 #ifndef CUDA
 mt19937_init(time(NULL));
 #else
 
 #endif
}

DEVICEFUNC INLINE
unsigned long long sim5rand()
{
 #ifndef CUDA
 return mt19937_int64();
 #else

 static long rndseed = 0x4d544750;
 #if LONG_MAX > (16807*2147483647)
 int const a    = 16807;      
 int const m    = 2147483647; 
	rndseed = ((long)(rndseed * a))%m;
	return rndseed;
 #else
 double const a    = 16807;      
 double const m    = 2147483647; 
	double temp = rndseed * a;
	rndseed = (int) (temp - m * floor ( temp / m ));
	return rndseed;
 #endif

 #endif
}


DEVICEFUNC INLINE
double sim5urand()
{
 #ifndef CUDA
 return mt19937_real1();
 #else

 static long rndseed = 0x4d544750;
 #if LONG_MAX > (16807*2147483647)
 int const a    = 16807;      
 int const m    = 2147483647; 
	rndseed = ((long)(rndseed * a))%m;
	return (double)(rndseed)/(double)(m);
 #else
 double const a    = 16807;      
 double const m    = 2147483647; 
	double temp = rndseed * a;
	rndseed = (int) (temp - m * floor ( temp / m ));
	return (double)(rndseed)/m;
 #endif

 #endif
}



/*
#ifdef CUDA
 __device__
#endif
long rndseed = 0x4d544750;
DEVICEFUNC
double urand()
{
 #if LONG_MAX > (16807*2147483647)
 int const a    = 16807;      
 int const m    = 2147483647; 
	rndseed = ((long)(rndseed * a))%m;
	return (double)(rndseed)/(double)(m);
 #else
 double const a    = 16807;      
 double const m    = 2147483647; 
	double temp = rndseed * a;
	rndseed = (int) (temp - m * floor ( temp / m ));
	return (double)(rndseed)/m;
 #endif
}
*/


DEVICEFUNC
void cartesian2spherical1(double x, double y, double z, double Vx, double Vy, double Vz, double* Vr, double* Vh, double* Vf)
{
 double r = sqrt(x*x + y*y + z*z);
 double cos_h = z/r;
 double sin_h = sqrt(1.-sqr(cos_h));
 double cos_f = x/r/sin_h;
 double sin_f = y/r/sin_h;
 (*Vr) = sin_h*cos_f*Vx + sin_h*sin_f*Vy + cos_h*Vz;
 (*Vh) = cos_h*cos_f*Vx + cos_h*sin_f*Vy - sin_h*Vz;
 (*Vf) =      -sin_f*Vx +       cos_f*Vy;
}

DEVICEFUNC
void cartesian2spherical2(double cos_h, double sin_f, double cos_f, double Vx, double Vy, double Vz, double* Vr, double* Vh, double* Vf)
{
 double sin_h = sqrt(1.-sqr(cos_h));
 (*Vr) = sin_h*cos_f*Vx + sin_h*sin_f*Vy + cos_h*Vz;
 (*Vh) = cos_h*cos_f*Vx + cos_h*sin_f*Vy - sin_h*Vz;
 (*Vf) =      -sin_f*Vx +       cos_f*Vy;
}







DEVICEFUNC INLINE
sim5complex makeComplex(double r, double i)
{
#ifdef CUDA
 sim5complex res;
 res.x = r;
 res.y = i;
 return res;
#else
 sim5complex res;
 res = r + ComplexI*i;
 return res;
#endif
}


DEVICEFUNC INLINE
sim5complex nullComplex()
{
#ifdef CUDA
 return makeComplex(0.0,0.0);
#else
 return 0.0;
#endif
}

DEVICEFUNC INLINE
double sim5creal (sim5complex a)
{
 return creal(a);
}

DEVICEFUNC INLINE
double sim5cimag (sim5complex a)
{
 return cimag(a);
}


#ifdef CUDA

DEVICEFUNC INLINE
double creal (sim5complex a)
{
 return a.x;
}


DEVICEFUNC INLINE
double cimag (sim5complex a)
{
 return a.y;
}


DEVICEFUNC INLINE
sim5complex cconj (sim5complex a)
{
 return makeComplex(creal(a), -cimag(a));
}


DEVICEFUNC INLINE
double cabs(sim5complex a)
{
 

 double x = creal(a);
 double y = cimag(a);
 double v, w, t;
 x = (double)fabs(x);
 y = (double)fabs(y);
 if (x > y) {
     v = x;
     w = y;
 } else {
     v = y;
     w = x;
 }
 t = w / v;
 t = 1.0 + t * t;
 t = v * (double)sqrt(t);
 if ((v == 0.0) || (v > 3.402823466e38) || (w > 3.402823466e38)) {
     t = v + w;
 }
 return t;
}


DEVICEFUNC INLINE
sim5complex cadd(sim5complex a, sim5complex b)
{
 return makeComplex(creal(a)+creal(b), cimag(a)+cimag(b));
}


DEVICEFUNC INLINE
sim5complex csub(sim5complex a, sim5complex b)
{
 return makeComplex(creal(a)-creal(b), cimag(a)-cimag(b));
}


DEVICEFUNC INLINE
sim5complex cmul(sim5complex a, sim5complex b)
{
 sim5complex prod;
 prod = makeComplex(
     creal(a)*creal(b) - cimag(a)*cimag(b),
     creal(a)*cimag(b) + cimag(a)*creal(b)
 );
 return prod;
}


DEVICEFUNC INLINE
sim5complex cdiv(sim5complex a, sim5complex b)
{
 

 sim5complex quot;
 double s = ((double)fabs((double)creal(b))) +
            ((double)fabs((double)cimag(b)));
 double oos = 1.0 / s;
 double ars = creal(a) * oos;
 double ais = cimag(a) * oos;
 double brs = creal(b) * oos;
 double bis = cimag(b) * oos;
 s = (brs * brs) + (bis * bis);
 oos = 1.0 / s;
 quot = makeComplex(
     (ars*brs + ais*bis) * oos,
     (ais*brs - ars*bis) * oos
 );
 return quot;
}


DEVICEFUNC INLINE
sim5complex csqrt(sim5complex a)
{
	sim5complex c;
	double x,y,w,r;
	if ((a.x == 0.0) && (a.y == 0.0)) {
	 return nullComplex();
	} else {
		x=fabs(a.x);
		y=fabs(a.y);
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (a.x >= 0.0) {
			c.x=w;
			c.y=a.y/(2.0*w);
		} else {
			c.y=(a.y >= 0) ? w : -w;
			c.x=a.y/(2.0*c.y);
		}
		return c;
	}
}


DEVICEFUNC INLINE
sim5complex catan(sim5complex a)
{
/*
 double re = creal(a), im = cimag(a);
 sim5complex z;

 if (im == 0) return makeComplex(atan(re), 0.0);

 
 
 
 double r = hypot(re,im);
 double imag;
 double u = 2.*im / (1. + r*r);

 
 if (fabs(u) < 0.1) {
     imag = 0.25 * (log1p(u) - log1p(-u));
 }
 else {
     double A = hypot(re, im+1.);
     double B = hypot(re, im-1.);
     imag = 0.5*log(A/B);
 }

 if (re == 0.0) {
     if (I > 1) {
         z = makeComplex(PI_half, im);
     }
     else if (I < -1.0) {
         z = makeComplex(-PI_half, im);
     }
     else {
         z = makeComplex(0.0, im);
     };
 } else {
     z = makeComplex(0.5*atan2(2.*R, (1.+r)*(1.-r)), im);
 }

  return z;
*/
 return nullComplex();
}





inline DEVICEFUNC sim5complex operator+(sim5complex a, sim5complex b)
{
 return cadd(a,b);
}


inline DEVICEFUNC sim5complex operator+(sim5complex a, double b)
{
 return makeComplex(creal(a)+b, cimag(a));
}


inline DEVICEFUNC sim5complex operator+(double a, sim5complex b)
{
 return makeComplex(a+creal(b), cimag(b));
}


inline DEVICEFUNC void operator+=(sim5complex a, double b)
{
 a.x += b;
}


inline DEVICEFUNC void operator+=(sim5complex a, sim5complex b)
{
 a.x += b.x;
 a.y += b.y;
}


inline DEVICEFUNC sim5complex operator-(sim5complex a, sim5complex b)
{
 return csub(a,b);
}


inline DEVICEFUNC sim5complex operator-(sim5complex a, double b)
{
 return makeComplex(creal(a)-b, cimag(a));
}


inline DEVICEFUNC sim5complex operator-(double a, sim5complex b)
{
 return makeComplex(a-creal(b), -cimag(b));
}


inline DEVICEFUNC void operator-=(sim5complex a, double b)
{
 a.x -= b;
}


inline DEVICEFUNC void operator-=(sim5complex a, sim5complex b)
{
 a.x -= b.x;
 a.y -= b.y;
}

inline DEVICEFUNC sim5complex operator-(sim5complex &a)
{
 return makeComplex(-creal(a), -cimag(a));
}


inline DEVICEFUNC sim5complex operator*(sim5complex a, sim5complex b)
{
 return cmul(a,b);
}


inline DEVICEFUNC sim5complex operator*(double a, sim5complex b)
{
 return makeComplex(a*creal(b), a*cimag(b));
}


inline DEVICEFUNC sim5complex operator*(sim5complex a, double b)
{
 return makeComplex(b*creal(a), b*cimag(a));
}


inline DEVICEFUNC sim5complex operator/(sim5complex a, sim5complex b)
{
 return cdiv(a,b);
}


inline DEVICEFUNC sim5complex operator/(sim5complex a, double b)
{
 return makeComplex(creal(a)/b, cimag(a)/b);
}


#endif



void gprintf(FILE* file, const char *templatex, ...)
{
 #ifndef CUDA
 va_list ap;

 
 va_start (ap, templatex);
 vfprintf (stderr, templatex, ap);
 va_end (ap);

 
 va_start (ap, templatex);
 vfprintf (file, templatex, ap);
 va_end (ap);
 #endif
}


void warning(const char *templatex, ...)
{
 #ifndef CUDA
 va_list ap;
 fprintf (stderr, "WRN: ");
 va_start (ap, templatex);
 vfprintf (stderr, templatex, ap);
 va_end (ap);
	fprintf(stderr, "\n");
 #endif
}


void error(const char *templatex, ...)
{
 #ifndef CUDA
 va_list ap;
 fprintf (stderr, "ERROR: ");
 va_start (ap, templatex);
 vfprintf (stderr, templatex, ap);
 va_end (ap);
 fprintf(stderr, "\n");
 
 #endif
}




int sort_array_d_compare_func(const void *x, const void *y) {
 return (*(double*)x - *(double*)y);
}

void sort_array(double *array, int N)
{
 qsort(array, N, sizeof(double), sort_array_d_compare_func);
}



int sort_array_f_compare_func(const void *x, const void *y) {
 return (*(float*)x - *(float*)y);
}

void sort_array_f(float *array, int N)
{
 qsort(array, N, sizeof(float), sort_array_f_compare_func);
}


#ifndef CUDA
void* sim5_array_alloc(size_t capacity, size_t element_size)
{
 size_t size_total = element_size*capacity    
                   + sizeof(size_t)           
                   + sizeof(size_t)           
                   + sizeof(size_t);          

 void* ptr = malloc(size_total);
 memset(ptr, '\0', size_total);

 
 *(size_t*)(ptr) = element_size;
 ptr += sizeof(size_t);

 
 *(size_t*)(ptr) = capacity;
 ptr += sizeof(size_t);

 
 *(size_t*)(ptr) = 0;
 ptr += sizeof(size_t);

 return ptr;
}


void sim5_array_free(void* ptr)
{
 ptr -= 3*sizeof(size_t);
 free(ptr);
}


void* sim5_array_realloc(void* array, size_t new_capacity)
{
 size_t* count = (size_t*)(array-1*sizeof(size_t));
 size_t* capa  = (size_t*)(array-2*sizeof(size_t));
 size_t* esize = (size_t*)(array-3*sizeof(size_t));

 size_t size_total = new_capacity * (*esize)    
                   + sizeof(size_t)             
                   + sizeof(size_t)             
                   + sizeof(size_t);            

 (*capa)  = new_capacity;
 (*count) = fmin(*count, new_capacity);

 void* src_ptr = (array-3*sizeof(size_t));
 src_ptr = realloc(src_ptr, size_total);
 return src_ptr+3*sizeof(size_t);
}


inline
long sim5_array_count(void* array)
{
 array -= 1*sizeof(size_t);
 return *(size_t*)(array);
}


inline
long sim5_array_capa(void* array)
{
 array -= 2*sizeof(size_t);
 return *(size_t*)(array);
}


inline
size_t sim5_array_esize(void* array)
{
 array -= 3*sizeof(size_t);
 return *(size_t*)(array);
}


void sim5_array_push(void** array_ptr, const void* data)
{
 void* array = *array_ptr;

 size_t* count = (size_t*)(array-1*sizeof(size_t));
 size_t* capa  = (size_t*)(array-2*sizeof(size_t));
 size_t* esize = (size_t*)(array-3*sizeof(size_t));

 if (*count+1 > *capa) {
     array = *array_ptr = sim5_array_realloc(array, 2*(*capa));
     count = (size_t*)(array-1*sizeof(size_t));
     capa  = (size_t*)(array-2*sizeof(size_t));
     esize = (size_t*)(array-3*sizeof(size_t));
 }

 void* dptr = array + (*esize)*(*count);
 memcpy(dptr, data, *esize);
 (*count)++;
}


inline
void sim5_array_push_int(void** array_ptr, const int data)
{
 sim5_array_push(array_ptr, &data);
}


inline
void sim5_array_push_long(void** array_ptr, const long data)
{
 sim5_array_push(array_ptr, &data);
}


inline
void sim5_array_push_float(void** array_ptr, const float data)
{
 sim5_array_push(array_ptr, &data);
}


inline
void sim5_array_push_double(void** array_ptr, const double data)
{
 sim5_array_push(array_ptr, &data);
}


int  sim5_array_exists(void* array, const void* data)
{
 size_t count = *(size_t*)(array-1*sizeof(size_t));
 size_t esize = *(size_t*)(array-3*sizeof(size_t));

 size_t i;
 for (i=0; i<count; i++) {
     if (memcmp(array, data, esize) == 0) return 1;
     array += esize;
 }

 return 0;
}

inline
void sim5_array_push_if_not_exists(void** array_ptr, const void* data)
{
 if (!sim5_array_exists(*array_ptr, data)) sim5_array_push(array_ptr, data);
}


void sim5_array_reverse(void* array)
{
 long count = *(size_t*)(array-1*sizeof(size_t));
 long esize = *(size_t*)(array-3*sizeof(size_t));

 void* frst_ptr = array;
 void* last_ptr = array + esize*(count-1);
 void* tmpdata1 = malloc(esize);
 void* tmpdata2 = malloc(esize);

 while(frst_ptr < last_ptr) {
     memcpy(tmpdata1, frst_ptr, esize);
     memcpy(tmpdata2, last_ptr, esize);
     memcpy(frst_ptr, tmpdata2, esize);
     memcpy(last_ptr, tmpdata1, esize);
     frst_ptr += esize;
     last_ptr -= esize;
 }

 free(tmpdata1);
 free(tmpdata2);
}
#endif


char* key_value_get(const char *string, const char *key)
{
 static char value[256];
 char* p;

 if ((!string) || (!key)) return NULL;

 
 char* string_copy = (char*)malloc(strlen(string)+1);
 strcpy (string_copy, string);

 for (p = strtok(string_copy,","); p != NULL; p = strtok(NULL, ",")) {
     if ((strncmp(p, key, strlen(key))==0) && (p[strlen(key)]=='=')) {
         p += strlen(key) + 1;
         strncpy(value, p, sizeof(value)-1);
         free(string_copy);
         return value;
     }
 }

 
 free(string_copy);
 return NULL;
}




#define MAXTOKENS    1024
#define MAXLINE      8096     
#define MINLEN       3        

char **split(char *string, char *delim) {
 char **tokens = NULL;
 char *working = NULL;
 char *token = NULL;
 int idx = 0;

 tokens  = (char**)malloc(sizeof(char *) * MAXTOKENS);
 if(tokens == NULL) return NULL;
 working = (char*)malloc(sizeof(char) * strlen(string) + 1);
 if(working == NULL) return NULL;

 
 strcpy(working, string);
 for(idx = 0; idx < MAXTOKENS; idx++) tokens[idx] = NULL;

 token = strtok(working, delim);
 idx = 0;

 
 while((idx < (MAXTOKENS - 1)) && (token != NULL)) {
     tokens[idx] = (char*)malloc(sizeof(char) * strlen(token) + 1);
     if(tokens[idx] != NULL) {
         strcpy(tokens[idx], token);
         idx++;
         token = strtok(NULL, delim);
     }
 }

 free(working);
 return tokens;
}



long getlinecount(FILE* f) {
 long fpos = ftell(f);
 char line[8192];
 long result = 0;
 fseek(f,0,SEEK_SET);
 while (fgets(line, 8192, f) != NULL) result++;
 fseek(f,fpos,SEEK_SET);
 return result;
}






void backtrace() {
/*
 unw_cursor_t cursor;
 unw_context_t context;

 unw_getcontext(&context);
 unw_init_local(&cursor, &context);

 int n=0;
 while ( unw_step(&cursor) ) {
     unw_word_t ip, sp, off;

     unw_get_reg(&cursor, UNW_REG_IP, &ip);
     unw_get_reg(&cursor, UNW_REG_SP, &sp);

     char symbol[256] = {"<unknown>"};
     char *name = symbol;

     off = 0;
     
     
     
     
     

     printf("#%-2d 0x%016" PRIxPTR " sp=0x%016" PRIxPTR " %s + 0x%" PRIxPTR "\n", ++n, ip, sp, name, off);

     if (name != symbol) free(name);
 }
*/
}






#define NMAX_TRAPEZOID 23
#define NMAX_SIMPSON 23


DEVICEFUNC
static void integrate_trapezoid_rule(double(*f)(double), double a, double b, int n, double *s)
{
 double x, tnm, sum, del;
 int it, j;

 if(n==1){
     *s = 0.5 * (b-a) * (f(b) + f(a));
 } else {
     for(it=1, j=1; j < n-1; j++)  it <<= 1;
     tnm = (double) it;
     del = (b-a) / tnm;
     x = a + 0.5 * del;
     for(sum=0.0, j=1; j<=it; j++, x+=del) { sum += f(x); }
     *s = 0.5 * (*s + del * sum);
 }
}


DEVICEFUNC
double integrate_trapezoid(double(*f)(double), double a, double b, double acc)
{
 int n;
 double s = 0.0;
 double olds = DBL_MIN;  

 for (n=1; n<=NMAX_TRAPEZOID; n++) {
     integrate_trapezoid_rule(f, a, b, n, &s);
     if (n > 3) {        
         if ((fabs(s-olds) < acc*fabs(olds)) || ((s==0.) && (olds==0.))) return s;
     }
     olds = s;
 }

 #ifndef CUDA
 warning("too many steps in integrate_trapezoid()");
 #endif

 return s;
}



DEVICEFUNC
double integrate_simpson(double (*f)(double), double a, double b, double acc)
{
 int n;
 double s, st=0.0, ost, os;

 ost = os = -1.e50;

 for (n=1; n<=NMAX_SIMPSON; n++){
     integrate_trapezoid_rule(f, a, b, n, &st);
     s = (4.*st - ost) / 3.;
     if (n > 3) {        
         if ((fabs(s-os) < acc*fabs(os)) || ((s==0.) && (os==0.))) return s;
     }
     os = s;
     ost = st;
 }

 #ifndef CUDA
 warning("too many steps in integrate_simpson()");
 #endif

 return s;
}



double gammln(double xx)
{
  double x,tmp,ser;
  static double cof[6]={76.18009173,-86.50532033,24.01409822,
			-1.231739516,0.120858003e-2,-0.536382e-5};
  int j;

  x=xx-1.0;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<=5;j++) {
 x += 1.0;
 ser += cof[j]/x;
  }
  return -tmp+log(2.50662827465*ser);
}



DEVICEFUNC
void gauleg(double x1, double x2, double x[], double w[], int n)
{
 const double EPS = 3.0e-11;
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;

	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=1;i<=m;i++) {
		z=cos(3.141592654*(i-0.25)/(n+0.5));
		do {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1) > EPS);
		x[i-1]=xm-xl*z;
		x[n+1-i-1]=xm+xl*z;
		w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n+1-i-1]=w[i-1];
	}
}



float qgaus(float (*func)(float), float a, float b)
{
 int j;
 float xr,xm,dx,s;
 static float x[]={0.0,0.1488743389,0.4333953941,0.6794095682,0.8650633666,0.9739065285}; 
 static float w[]={0.0,0.2955242247,0.2692667193,0.2190863625,0.1494513491,0.0666713443};
 xm=0.5*(b+a);
 xr=0.5*(b-a);
 s=0;  
 for (j=1;j<=5;j++) {
     dx=xr*x[j];  
     s += w[j]*((*func)(xm+dx)+(*func)(xm-dx));
 }
 return s *= xr;  
}



double qgaus_general(double w[], double y[], int N, double a, double b)
{
 int j;
 double s  = 0.0;
 for (j=0; j<N; j++) s += w[j]*y[j];
 return s*(b-a);  
}


#undef NMAX_TRAPEZOID
#undef NMAX_SIMPSON


#ifndef CUDA




DEVICEFUNC INLINE
long sim5_interp_search(const double x_array[], double x, long index_lo, long index_hi)
{
 long ilo = index_lo;
 long ihi = index_hi;
 while(ihi > ilo + 1) {
     long i = (ihi + ilo)/2;
     if (x < x_array[i])
         ihi = i;
     else
         ilo = i;
 }
 return ilo;
}




DEVICEFUNC INLINE
long sim5_interp_search_accel(sim5interp* interp, double x)
{
 long x_index = interp->last_index;
 if(x < interp->X[x_index]) {
     
     interp->last_index = sim5_interp_search(interp->X, x, 0, x_index);
 } else
 if(x >= interp->X[x_index + 1]) {
     
     interp->last_index = sim5_interp_search(interp->X, x, x_index, interp->N-1);
 } 
     
 
 return interp->last_index;
}



DEVICEFUNC
static void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
 int i,k;
 double p, qn, sig, un, *u;

 
 u = (double*)malloc((n-1)*sizeof(double));
 
 if (u == NULL) exit(EXIT_FAILURE);
 
 
 
 

 if(yp1 > 0.99e30)
     y2[0] = u[0] = 0.0;
 else{
     y2[0] = -0.5;
     u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
 }

 for(i = 1; i < n-1; i++){
     sig = (x[i] - x[i-1])/(x[i+1] - x[i-1]);
     p = sig*y2[i-1] + 2.0;
     y2[i] = (sig - 1.0)/p;
     u[i] = (y[i+1] - y[i])/(x[i+1] - x[i]) - (y[i] - y[i-1])/(x[i] - x[i-1]);
     u[i] = (6.0*u[i]/(x[i+1] - x[i-1]) - sig*u[i-1])/p;
 }

 if(ypn > 0.99e30)
     qn = un = 0.0;
 else{
     qn = 0.5;
     un = (3.0/(x[n-1] - x[n-2]))*(ypn - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]));
 }

 y2[n-1] = (un - qn*u[n-2])/(qn*y2[n-2] + 1.0);

 for(k = n-2; k >= 0; k--){
     y2[k] = y2[k]*y2[k+1] + u[k];
 }

 free(u);
}



DEVICEFUNC
static double splint(double xa[], double ya[], double y2a[], int n, double x)
{
 int  klo,khi,k;
 double h,b,a;
 static int pklo=0, pkhi=1;

 #pragma omp threadprivate(pklo,pkhi)

 
 
 
 
 if(xa[pklo] <= x && xa[pkhi] > x){
     klo = pklo;
     khi = pkhi;
 }
 else{
     klo = 0;
     khi = n-1;
     while (khi-klo > 1){
         k = (khi + klo) >> 1;
         if(xa[k] > x) khi = k; else klo = k;
     }
     pklo = klo;
     pkhi = khi;
 }

 h = xa[khi] - xa[klo];
 
 
 
 
 
 a = (xa[khi] - x)/h;
 b = (x - xa[klo])/h;
 return  a*ya[klo] + b*ya[khi] + ((a*a*a - a)*y2a[klo] + (b*b*b - b)*y2a[khi])*(h*h)/6.0;
}



DEVICEFUNC
void sim5_interp_init(sim5interp* interp, double xa[], double ya[], long N, int data_model, int interp_type, int interp_options)

{
 if ((interp_type==INTERP_TYPE_SPLINE) && (interp_options & INTERP_OPT_CAN_EXTRAPOLATE)) {
     
     fprintf(stderr, "ERR (sim5_interp_init): spline interpolation cannot be used with extrapolation option\n");
     
     return;
 }


 interp->datamodel = data_model;
 interp->type      = interp_type;
 interp->options   = interp_options;
 interp->d2Y       = NULL;

 
 if ((interp->datamodel==INTERP_DATA_REF) || (interp->datamodel==INTERP_DATA_COPY)) {
     long i;
     for (i=0; i<N-1; i++) {
         if (xa[i] >= xa[i+1]) {
             
             fprintf(stderr, "ERR (sim5_interp_init): unordered X grid (x[%ld]=%.4e, x[%ld]=%.4e, N=%ld, opt=%d)\n", i, xa[i], i+1, xa[i+1], N, interp_options);
             backtrace();
             

             interp->N = 0;
             interp->X = NULL;
             interp->Y = NULL;
             
             exit(-1);
             
             
             
         }
     }
 }

 switch (interp->datamodel) {
     case INTERP_DATA_REF:
         
         interp->N = N;
         interp->capa = 0;
         interp->X = xa;
         interp->Y = ya;
         interp->xmin = interp->X[0];
         interp->xmax = interp->X[N-1];
         interp->last_index = (N-1)/2;
         break;

     case INTERP_DATA_COPY:
         
         interp->N = N;
         interp->capa = N;
         interp->X = (double*)malloc(N*sizeof(double));
         interp->Y = (double*)malloc(N*sizeof(double));
         memcpy (interp->X, xa, N*sizeof(double));
         memcpy (interp->Y, ya, N*sizeof(double));
         interp->xmin = interp->X[0];
         interp->xmax = interp->X[N-1];
         interp->last_index = (N-1)/2;
         break;

     case INTERP_DATA_BUILD:
         
         interp->N = 0;
         interp->capa = N>0 ? N : 100;
         interp->X = (double*)calloc(interp->capa, sizeof(double));
         interp->Y = (double*)calloc(interp->capa, sizeof(double));
         interp->xmin = 0.0;
         interp->xmax = 0.0;
         interp->last_index = 0;
         break;

     default:
         
         fprintf(stderr, "ERR (sim5_interp_init): unimplemented data model (%d)\n", interp->datamodel);
         exit(-1);
         
         
         
 }

}



DEVICEFUNC
void sim5_interp_data_push(sim5interp* interp, double x, double y)
{
 if (interp->datamodel != INTERP_DATA_BUILD) {
     
     fprintf(stderr, "ERR (sim5_interp_data_push): you can only push in data with INTERP_DATA_BUILD data model\n");
     
     return;
 }

 long i = interp->N;

 if ((i>0) && (interp->X[i-1] >= x)) {
     
     fprintf(stderr, "ERR (sim5_interp_data_push): unordered X grid (x[%ld]=%.4e, x[%ld]=%.4e)\n", i-1, interp->X[i-1], i, x);
     exit(-1);
     
     
     
 }

 interp->X[i] = x;
 interp->Y[i] = y;
 interp->N++;

 if (interp->N >= interp->capa) {
     interp->capa *= 2;
     interp->X = (double*)realloc(interp->X, interp->capa*sizeof(double));
     interp->Y = (double*)realloc(interp->Y, interp->capa*sizeof(double));
 }

 interp->xmin = interp->X[0];
 interp->xmax = interp->X[i];
 interp->last_index = i/2;
}



DEVICEFUNC
double sim5_interp_eval(sim5interp* interp, double x)
{
 double x_lo, x_hi;
 double y_lo, y_hi;
 long index;

 
 if (interp->type == INTERP_TYPE_SPLINE) {
     
     if (!interp->d2Y) {
         interp->d2Y = (double*) malloc(interp->N*sizeof(double));
         spline(interp->X, interp->Y, interp->N, 1e50, 1e50, interp->d2Y);
     }
     return splint(interp->X, interp->Y, interp->d2Y, interp->N, x);
 }



 if ((!(interp->options & INTERP_OPT_CAN_EXTRAPOLATE)) && ((x < interp->xmin) || (x > interp->xmax))) {
     
     fprintf(stderr, "WRN (sim5_interp_eval): unwarranted extrapolation (x=%.4e, xmin=%.4e, xmax=%.4e)\n", x, interp->xmin, interp->xmax);
     
 }

 if (interp->options & INTERP_OPT_ACCEL) {
     
     index = sim5_interp_search_accel(interp, x);
 } else {
     
     index = sim5_interp_search(interp->X, x, 0, interp->N-1);
 }

 x_lo = interp->X[index];
 x_hi = interp->X[index + 1];
 y_lo = interp->Y[index];
 y_hi = interp->Y[index + 1];

 
 
 
 
 

 switch (interp->type) {
     case INTERP_TYPE_LINLIN:
         return y_lo + (x-x_lo)/(x_hi-x_lo) * (y_hi-y_lo);

     case INTERP_TYPE_LOGLOG:
         return exp(log(y_lo) + (log(x)-log(x_lo)) / (log(x_hi) - log(x_lo)) * (log(y_hi)-log(y_lo)));
         

     case INTERP_TYPE_LOGLIN:
         return y_lo + log(x/x_lo)/log(x_hi/x_lo) * (y_hi-y_lo);
         

     default:
         
         fprintf(stderr, "ERR (sim5_interp_eval): unimplemented interpolation type (%d)\n", interp->type);
         
         return NAN;
 }
}


/*
double sim5_interp_integral(sim5interp* interp, double a, double b)
{
 int i, N = 500;
 double result = 0.0;
 for (i=0; i<N; i++) result += sim5_interp_eval(interp, a+(i+0.5)*(b-a)/(N));
 return result*(b-a)/N;
}
*/



DEVICEFUNC
void sim5_interp_done(sim5interp* interp)
{
 if ((interp->datamodel==INTERP_DATA_COPY) || (interp->datamodel==INTERP_DATA_BUILD)){
     free(interp->X);
     free(interp->Y);
 }

 if (interp->d2Y) free(interp->d2Y);

 interp->N = 0;
 interp->capa = 0;
 interp->X = NULL;
 interp->Y = NULL;
}




DEVICEFUNC
sim5interp* sim5_interp_alloc()
{
 return (sim5interp*)calloc(1, sizeof(sim5interp));
}




DEVICEFUNC
void sim5_interp_free(sim5interp* interp)
{
 sim5_interp_done(interp);
 free(interp);
}





#ifdef SIM5FILEIO_TESTING
int main() {
 double X[5] = {1.,2.,3.,4.,5.};
 double Y[5] = {2.,4.,6.,8.,10.};

 sim5interp interp;
 sim5_interp_init(&interp, X, Y, 5, INTERP_TYPE_LINLIN, INTERP_OPT_ALLOW_EXTRAPOLATION+INTERP_OPT_ACCEL);

 double x;
 for (x=0.; x<10.; x+=0.1) printf("%e %e\n", x, sim5_interp_eval(&interp, x));
 sim5_interp_free(&interp);
 return 0;
}
#endif

#endif 





#ifndef CUDA

DEVICEFUNC
void distrib_init(sim5distrib* d, double(*pdf)(double), double x_min, double x_max, int N)
{
 int i;
 double *tmp_x = (double*)calloc(N, sizeof(double));
 double *tmp_w = (double*)calloc(N, sizeof(double));
 double *tmp_pdf = (double*)calloc(N, sizeof(double));
 double *tmp_cdf = (double*)calloc(N, sizeof(double));

 d->x_min = x_min;
 d->x_max = x_max;
 d->norm  = integrate_simpson(pdf, x_min, x_max, 1e-5);

 
 tmp_x[0]   = x_min;
 tmp_pdf[0] = pdf(x_min);
 tmp_cdf[0] = 0.0;

 
 
 
 gauleg(x_min, x_max, &tmp_x[1], tmp_w, N-2);

 
 for (i=1; i<N-1; i++) {
     
     tmp_pdf[i] = pdf(tmp_x[i])/d->norm;
     tmp_cdf[i] = integrate_simpson(pdf, x_min, tmp_x[i], 1e-6)/d->norm;
     if ((i>1) && (tmp_cdf[i]<tmp_cdf[i-1])) fprintf(stderr, "CDF problem at %d/%d (%e  %e)\n", i,N,tmp_cdf[i],tmp_cdf[i-1]);
 }

 
 tmp_x[N-1]   = x_max;
 tmp_pdf[N-1] = pdf(x_max);
 tmp_cdf[N-1] = 1.0;

 sim5_interp_init(&d->pdf, tmp_x, tmp_pdf, N, INTERP_DATA_COPY, INTERP_TYPE_LINLIN, 0);
 sim5_interp_init(&d->cdf, tmp_x, tmp_cdf, N, INTERP_DATA_COPY, INTERP_TYPE_LINLIN, 0);
 sim5_interp_init(&d->icd, tmp_cdf, tmp_x, N, INTERP_DATA_COPY, INTERP_TYPE_LINLIN, 0);

 free(tmp_x);
 free(tmp_w);
 free(tmp_pdf);
 free(tmp_cdf);
}



DEVICEFUNC
void distrib_done(sim5distrib* d)
{
 sim5_interp_done(&d->pdf);
 sim5_interp_done(&d->cdf);
 sim5_interp_done(&d->icd);
}



DEVICEFUNC INLINE
double distrib_hit(sim5distrib* d)
{
 return sim5_interp_eval(&d->icd, sim5urand());
}

#endif






#define MAX_STEPS  500

long rtbis(double x1, double x2, double xacc, double (*fx)(double), double* result)
{
	double dx, f, fmid, xmid, rtb;
	long j=0;

	fmid = (*fx)(x2);
	f = (*fx)(x1);
	if ((f*fmid) >= 0.0) return(0);

	if (f < 0.0) {
		rtb = x1;
		dx  = x2-x1;
	}
	else {
		rtb = x2;
		dx  = x1-x2;
	}

	for (j=0; j<MAX_STEPS; j++) {
		dx = dx*0.5;
		xmid = rtb+dx;
		fmid = (*fx)(xmid);
		if (fmid <= 0.0) rtb = xmid;
		if ((fabs(dx) < xacc) || (fmid == 0.0)) break;
	}
	if (j >= MAX_STEPS) {
		error("rtbis: too many steps");
		return(0);
	}

	*result = rtb;
	return(1);
}


#undef MAX_STEPS
/*
#include "sim5config.h"
#ifndef CUDA
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "sim5math.h"
#include "sim5utils.h"
#include "sim5elliptic.h"
#endif
*/

DEVICEFUNC
double rf(double x, double y, double z)
{
	const double ERRTOL=0.0003, rfTINY=3.*DBL_MIN, rfBIG=DBL_MAX/3., THIRD=1.0/3.0;
	const double C1=1.0/24.0, C2=0.1, C3=3.0/44.0, C4=1.0/14.0;
	double alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt;

	if ((fmin(fmin(x,y),z) < 0.0) || (fmin(fmin(x+y,x+z),y+z) < rfTINY) || (fmax(fmax(x,y),z) > rfBIG)) {
     #ifndef CUDA
		warning("%e/%e/%e\n", x,y,z);
		error("invalid arguments in rf %e/%e/%e", x,y,z);
		backtrace();
     #endif
	}

	xt=x;
	yt=y;
	zt=z;
	do {
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		ave=THIRD*(xt+yt+zt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
	} while (fmax(fmax(fabs(delx),fabs(dely)),fabs(delz)) > ERRTOL);
	e2=delx*dely-delz*delz;
	e3=delx*dely*delz;
	return (1.0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave);
}


DEVICEFUNC
double rd(double x, double y, double z) {
	const double ERRTOL=0.0003, rdTINY=3.*DBL_MIN, rdBIG=DBL_MAX/3.;
	const double C1=3.0/14.0, C2=1.0/6.0, C3=9.0/22.0, C4=3.0/26.0, C5=0.25*C3, C6=1.5*C4;
	double alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt;

	if ((fmin(x,y) < 0.0) || (fmin(x+y,z) < rdTINY) || (fmax(fmax(x,y),z) > rdBIG)) {
     #ifndef CUDA
		fprintf(stderr, "%e/%e/%e\n", x,y,z);
		error("invalid arguments in rd (%e/%e/%e)\n", x,y,z);
     #endif
	}

	xt=x;
	yt=y;
	zt=z;
	sum=0.0;
	fac=1.0;
	do {
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		sum += fac/(sqrtz*(zt+alamb));
		fac=0.25*fac;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		ave=0.2*(xt+yt+3.0*zt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
	} while (fmax(fmax(fabs(delx),fabs(dely)),fabs(delz)) > ERRTOL);
	ea=delx*dely;
	eb=delz*delz;
	ec=ea-eb;
	ed=ea-6.0*eb;
	ee=ed+ec+ec;
	return 3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)
		+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave));
}


DEVICEFUNC
double rc(double x, double y) {
	const double ERRTOL=0.0003, rcTINY=3.*DBL_MIN, rcBIG=DBL_MAX/3.;
	const double THIRD=1.0/3.0, C1=0.3, C2=1.0/7.0, C3=0.375, C4=9.0/22.0;
	const double COMP1=2.236/sqrt(rcTINY),COMP2=sqr(rcTINY*rcBIG)/25.0;

	double alamb,ave,s,w,xt,yt;

	if ((x < 0.0) || (y == 0.0) || ((x+fabs(y)) < rcTINY) || ((x+fabs(y)) > rcBIG) ||
	  ((y<-COMP1) && (x > 0.0) && (x < COMP2))) {
     #ifndef CUDA
		fprintf(stderr, "%e/%e\n", x,y);
		error("invalid arguments in rc (%e/%e)\n", x, y);
     #endif
	}

	if (y > 0.0) {
		xt=x;
		yt=y;
		w=1.0;
	} else {
		xt=x-y;
		yt= -y;
		w=sqrt(x)/sqrt(xt);
	}
	do {
		alamb=2.0*sqrt(xt)*sqrt(yt)+yt;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		ave=THIRD*(xt+yt+yt);
		s=(yt-ave)/ave;
	} while (fabs(s) > ERRTOL);
	return w*(1.0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave);
}


DEVICEFUNC
double rj(double x, double y, double z, double p) {
	const double ERRTOL=0.0003, rjTINY=pow(5.0*DBL_MIN,1./3.), rjBIG=0.3*pow(0.1*DBL_MAX,1./3.);
	const double C1=3.0/14.0, C2=1.0/3.0, C3=3.0/22.0, C4=3.0/26.0,
		C5=0.75*C3, C6=1.5*C4, C7=0.5*C2, C8=C3+C3;
	double a,alamb,alpha,ans,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,ed,ee,
		fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt;
	if ((fmin(fmin(x,y),z) < 0.0) || (fmin(fmin(x+y,x+z),fmin(y+z,fabs(p))) < rjTINY) ||
	  (fmax(fmax(x,y),fmax(z,fabs(p))) > rjBIG)) {
     #ifndef CUDA
		fprintf(stderr, "WRN: invalid arguments in rj (%e/%e/%e/%e)\n", x,y,z,p);
     #endif
     return 0.0;
	}

 a=b=rcx=0.0;
	sum=0.0;
	fac=1.0;
	if (p > 0.0) {
		xt=x;
		yt=y;
		zt=z;
		pt=p;
	} else {
		xt=fmin(fmin(x,y),z);
		zt=fmax(fmax(x,y),z);
		yt=x+y+z-xt-zt;
		a=1.0/(yt-p);
		b=a*(zt-yt)*(yt-xt);
		pt=yt+b;
		rho=xt*zt/yt;
		tau=p*pt/yt;
		rcx=rc(rho,tau);
	}
	do {
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		alpha=sqr(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz);
		beta=pt*sqr(pt+alamb);
		sum += fac*rc(alpha,beta);
		fac=0.25*fac;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		pt=0.25*(pt+alamb);
		ave=0.2*(xt+yt+zt+pt+pt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
		delp=(ave-pt)/ave;
	} while (fmax(fmax(fabs(delx),fabs(dely)),fmax(fabs(delz),fabs(delp))) > ERRTOL);
	ea=delx*(dely+delz)+dely*delz;
	eb=delx*dely*delz;
	ec=delp*delp;
	ed=ea-3.0*ec;
	ee=eb+2.0*delp*(ea-ec);
	ans=3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4))
		+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave));
	if (p <= 0.0) ans=a*(b*ans+3.0*(rcx-rf(xt,yt,zt)));
	return ans;
}



DEVICEFUNC INLINE
double elliptic_k(double m)
{
  if (m==1.0) m = 1.0-1e-8;
  #ifndef CUDA
  if (m>1.0) error("Error in routine elliptic_k(m): m >= 1.0 (%e)", m);
  #endif
  return rf(0,1.0-m,1.0);
}



DEVICEFUNC
double elliptic_f(double phi, double m)
{
 #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_f(m): m >= 1.0");
 #endif
	if (m==1.0) m = 0.99999999;
	if (phi==0.0) return 0.0;

 int k = 0;
 while (fabs(phi) > M_PI/2.) { (phi>0)?k++:k--; phi += (phi>0)?-M_PI:+M_PI;  }

	double s2 = pow(sin(phi), 2);
	double ell_f = (phi>0?+1:-1) * sqrt(s2)*rf(1-s2, 1.0-s2*m, 1.0);
 if (k!=0) ell_f += 2.*k*elliptic_k(m);
 return ell_f;
}

DEVICEFUNC
double elliptic_f_cos(double cos_phi, double m)
{
 #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_f_cos(m): m >= 1.0");
 #endif
	if (m==1.0) m = 0.99999999;
	if (cos_phi==1.0) return 0.0;

	double X = 0.0;
	if (cos_phi<0.0) {
		cos_phi = - cos_phi;
		X = 2.0 * rf(0.0, 1.0-m, 1.0);
	}

	double s2 = 1.0-sqr(cos_phi);
	return X + ((X==0.0)?(+1):(-1))*sqrt(s2)*rf(1.0-s2, 1.0-s2*m, 1.0);
}

DEVICEFUNC
double elliptic_f_sin(double sin_phi, double m)
{
 #ifndef CUDA
	if ((sin_phi<0.0)||(sin_phi>1.0)) error("elliptic_f_sin: invalid input (sin_phi<0)");
	if (m>1.0) error("Error in routine elliptic_f_sin(m): m >= 1.0");
 #endif
	if (m==1.0) m = 0.99999999;
	if (sin_phi==0.0) return 0.0;
	double s2 = sqr(sin_phi);
	return sin_phi*rf(1.-s2, 1.0-s2*m, 1.0);
}





DEVICEFUNC
double elliptic_e(double phi, double m)
{
 #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_e(m): m >= 1.0");
 #endif
	if (m==1.0) m = 0.99999999;
	if (phi==0.0) return 0.0;

 #ifndef CUDA
	if ((phi<0.0)||(phi>M_PI)) error("elliptic_e: invalid input");
 #endif
	double X = 0.0;
	if (phi>M_PI/2.) {
		phi = M_PI - phi;
		X = 2.0 * ( rf(0.0,1.0-m,1.0) - m*rd(0.0,1.0-m,1.0)/3.0 );
	}

	double s  = sin(phi);
	double c2 = sqr(cos(phi));
	double q  = 1.0 - s*s*m;
	return X + ((X==0.0)?(+1):(-1))*s*( rf(c2,q,1.0) - sqr(s*sqrt(m))*rd(c2,q,1.0)/3.0 );
}

DEVICEFUNC
double elliptic_e_cos(double cos_phi, double m)
{
 #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_e_cos: m >= 1.0");
 #endif
	if (m==1.0) m = 0.99999999;
	if (cos_phi==1.0) return 0.0;

	double X = 0.0;
	if (cos_phi<0.0) {
		cos_phi = -cos_phi;
		X = 2.0 * ( rf(0.0,1.0-m,1.0) - m*rd(0.0,1.0-m,1.0)/3.0 );
	}

	double c2 = sqr(cos_phi);
	double s  = sqrt(1.0-c2);
	double q  = 1.0 - m + c2*m;
	return X + ((X==0.0)?(+1):(-1))*s*( rf(c2,q,1.0) - sqr(s*sqrt(m))*rd(c2,q,1.0)/3.0);
}

DEVICEFUNC
double elliptic_e_sin(double sin_phi, double m)
{
 #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_e_sin: m >= 1.0");
 #endif
	if (m==1.0) m = 0.99999999;
	if (sin_phi==0.0) return 0.0;

 #ifndef CUDA
	if ((sin_phi<0.0)||(sin_phi>1.0)) error("elliptic_e_sin: invalid input (sin_phi<0)");
 #endif
	double s2 = sin_phi*sin_phi;
	double c2 = 1.0 - s2;
	double q  = 1.0 - s2*m;
	return sin_phi*( rf(c2,q,1.0) - sqr(sin_phi*sqrt(m))*rd(c2,q,1.0)/3.0);
}



DEVICEFUNC
double elliptic_pi_complete(double n, double m)
{
 #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_pi_complete: m>1.0");
	if (n>1.0) error("Error in routine elliptic_pi_complete: n>1.0");
 #endif
 if (isinf(n)) return 0.0;
	if (m==1.0) m = 0.99999999;
	if (n==1.0) n = 0.99999999;

	double q  = 1.0-m;
	return rf(0.0,q,1.0) + n*rj(0.0,q,1.0,1.0-n)/3.0;
}


DEVICEFUNC
sim5complex elliptic_pi(double phi, double n, double m)
{
 #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_pi: m>1.0");
 #endif
 if (isinf(n)) return makeComplex(0.0,0.0);
	if (m==1.0) m = 0.99999999;
	if (phi==0.0) return makeComplex(0.0,0.0);

	int hyperbolicIII = (n > 1.0);
 double p=0,dn=0,la=0;
 if (hyperbolicIII) {
     p  = sqrt((n-1.)*(1.-m/n));
     dn = sqrt(1. - m*sin(phi)*sin(phi));
     la = (dn+p*tan(phi)) / (dn-p*tan(phi));
     n = m/n;
 }

 int k = 0;
 while (fabs(phi) > M_PI/2.) { (phi>0)?k++:k--; phi += (phi>0)?-M_PI:+M_PI;  }

 double s = sin(phi);
	double c2  = 1.0 - s*s;
	double q   = 1.0 - s*s*m;
	double ns2 = -n*s*s;

 #ifndef CUDA
 if (rj(c2,q,1.0,1.0+ns2)==0.0) error("Error in routine elliptic_pi: rj==0 (%e/%e/%e)", phi,n,m);
 #endif

 
 sim5complex ell_pi = makeComplex(s*(rf(c2,q,1.0) - ns2*rj(c2,q,1.0,1.0+ns2)/3.0), 0.0);
 if ((k!=0) && (!hyperbolicIII)) ell_pi += 2.*k*elliptic_pi_complete(n,m);

 
 if (hyperbolicIII) ell_pi = -ell_pi + elliptic_f(phi,m) + log(fabs(la))/(2*p) + ((la<0)?ComplexI*M_PI/(2*p):nullComplex());

 return ell_pi;
}


DEVICEFUNC
double elliptic_pi_cos(double cos_phi, double n, double m)
{
	

 #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_pi_cos: m>1.0");
 #endif
 if (isinf(n)) return 0.0;
 if (cos_phi==1.0) return 0.0;
 if (cos_phi==0.0) return elliptic_pi_complete(n, m);
	if (m==1.0) m = 0.99999999;
 

	double X = 0.0;
	if (cos_phi<0.0) {
		cos_phi = -cos_phi;
		X = 2.0 * ((rf(0.0,1.0-m,1.0) + n*rj(0.0,1.0-m,1.0,1.0-n)/3.0));
	}

	double c2  = sqr(cos_phi);
	double s   = sqrt(1.0-c2);
	double ns2 = -n*(1.0-c2);
	double q  = 1.0 - (1.0-c2)*m;
	return X + ((X==0.0)?(+1):(-1))*s*(rf(c2,q,1.0) - ns2*rj(c2,q,1.0,1.0+ns2)/3.0);
}


DEVICEFUNC
double elliptic_pi_sin(double sin_phi, double n, double m)
{
	double s2 = sin_phi*sin_phi;

 #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_pi_sin: m>1.0");
 #endif
 if (isinf(n)) return 0.0;
	if (m==1.0) m = 0.99999999;
	if (sin_phi==0.0) return 0.0;
 if (sin_phi==1.0) return elliptic_pi_complete(n, m);
	

 #ifndef CUDA
	if ((sin_phi<0.0)||(sin_phi>1.0)) error("elliptic_pi_sin: invalid input (sin_phi<0)");
 #endif
	double c2  = 1.0-s2;
	double ns2 = -n*s2;
	double q  = 1.0 - s2*m;
	return sin_phi*(rf(c2,q,1.0) - ns2*rj(c2,q,1.0,1.0+ns2)/3.0);
}


DEVICEFUNC INLINE
double jacobi_isn(double z, double m)
{
 if (fabs(m-0.0)<1e-8) return asin(z);                       
 if (fabs(m-1.0)<1e-8) return log(sqrt((1.+z)/(1.-z)));      
 return z*rf(1.0-z*z,1.0-m*z*z,1.0);
}


DEVICEFUNC INLINE
double jacobi_icn(double z, double m)
{
 if ((z>+1.0) && (z<+1.0+1e-8)) z = +1.0;
 if ((z<-1.0) && (z>-1.0-1e-8)) z = -1.0;
 if ((m>+1.0) && (m<+1.0+1e-8)) m =  1.0;
 if ((m< 0.0) && (m> 0.0-1e-8)) m =  0.0;

 if (z==0.0) return elliptic_k(m);
 if (z==1.0) return 0.0;
 if (m==0.0) return acos(z);                     
 if (m==1.0) return log((1.+sqrt(1.-z))/z);      

 #ifndef CUDA
 if (fabs(z)>1.0) error("jacobi_icn: u>1 (%.10e)", z);
 if (m>1.0) error("jacobi_cn: m>1.0 (%.10e)", m);
 #endif

 double icn1 = sqrt(1.-z*z) * rf(z*z, 1.0 - m*(1.-z*z), 1.0);

 
 return (z>0.0) ? icn1 : 2./sqrt(1.-m)*elliptic_f_sin(-z, m/(m-1.))+icn1;
}


DEVICEFUNC INLINE
double jacobi_itn(double z, double m)
{
 if (m==0.0) return atan(z);              
 if (m==1.0) return log(z+sqrt(1.+z*z));  
 return jacobi_isn(sqrt(z*z/(1.+z*z)), m);
}



DEVICEFUNC
void jacobi_sncndn(double u, double m, double *sn, double *cn, double *dn)
{
 #ifndef CUDA
	if (m>1.0) error("Error in routine sncndn: m1>1.0");
	if (u>2.*elliptic_k(m)) error("sncndn: invalid input (u>2K(m))");
 #endif
	if (m==1.0) m = 0.999999999;

	const double CA=1.0e-8;
	int bo;
	int i,ii,l;
	double a,b,c,d,emc;
	double em[13],en[13];

 d=1.0;
	emc=1.0-m; 
	if (emc != 0.0) {
		bo=(emc < 0.0);
		if (bo) {
			d=1.0-emc;
			emc /= -1.0/d;
			u *= (d=sqrt(d));
		}
		a=1.0;
		*dn=1.0;
		for (i=0;i<13;i++) {
			l=i;
			em[i]=a;
			en[i]=(emc=sqrt(emc));
			c=0.5*(a+emc);
			if (fabs(a-emc) <= CA*a) break;
			emc *= a;
			a=c;
		}
		u *= c;
		*sn=sin(u);
		*cn=cos(u);
		if (*sn != 0.0) {
			a=(*cn)/(*sn);
			c *= a;
			for (ii=l;ii>=0;ii--) {
				b=em[ii];
				a *= c;
				c *= *dn;
				*dn=(en[ii]+a)/(b+a);
				a=c/b;
			}
			a=1.0/sqrt(c*c+1.0);
			*sn=((*sn) >= 0.0 ? a : -a);
			*cn=c*(*sn);
		}
		if (bo) {
			a=*dn;
			*dn=*cn;
			*cn=a;
			*sn /= d;
		}
	} else {
		*cn=1.0/cosh(u);
		*dn=*cn;
		*sn=tanh(u);
	}
}

DEVICEFUNC INLINE
double jacobi_sn(double u, double m)
{
  double snx, cnx, dnx;
  jacobi_sncndn(u,m,&snx,&cnx,&dnx);
  return snx;
}

DEVICEFUNC INLINE
double jacobi_cn(double u, double m)
{
  double snx, cnx, dnx;
  jacobi_sncndn(u,m,&snx,&cnx,&dnx);
  return cnx;
}


DEVICEFUNC INLINE
double jacobi_dn(double u, double m)
{
  double sn, cn, dn;
  jacobi_sncndn(u,m,&sn,&cn,&dn);
  return dn;
}




DEVICEFUNC INLINE
double integral_C0(double u, double m)
{
	return u;
}


DEVICEFUNC INLINE
double integral_C1(double u, double m)
{
	double sn, cn, dn;
	jacobi_sncndn(u,m,&sn,&cn,&dn);
	return acos(dn)/sqrt(m);
}


DEVICEFUNC INLINE
double integral_C2(double u, double m)
{
	double sn, cn, dn;
	jacobi_sncndn(u,m,&sn,&cn,&dn);
	return 1./m*(elliptic_e_cos(cn,m) - (1.-m)*u);
}


DEVICEFUNC INLINE
double integral_C2_cos(double cn_u, double m)
{
	return 1./m*(elliptic_e_cos(cn_u,m) - (1.-m)*elliptic_f_cos(cn_u,m));
}


DEVICEFUNC
double integral_Z1(double a, double b, double u, double m)
{
 #ifndef CUDA
	if (m>1.0) error("Error in routine integral_Z1 (m>1.0)");
	if (u>2.*elliptic_k(m)) error("integral_Z1: invalid input (u>2K(m))");
 #endif
	double sn, cn, dn;
	jacobi_sncndn(u, m, &sn, &cn, &dn);
	return 1./a*((a-b)*elliptic_pi_cos(cn,a,m) + b*u);
}


DEVICEFUNC
double integral_Z2(double a, double b, double u, double m)
{
 #ifndef CUDA
	if (m>1.0) error("Error in routine integral_Z2: m>1.0");
	if (u>2.*elliptic_k(m)) error("integral_Z2: invalid input (u>2K(m))");
 #endif
	double sn, cn, dn;
	jacobi_sncndn(u, m, &sn, &cn, &dn);
	
	double V1 = elliptic_pi_cos(cn,a,m);
	double V2 = 0.5/((a-1.)*(m-a)) * (
	           a*elliptic_e_cos(cn,m) + (m-a)*u +
	           (2.*a*m+2.*a-a*a-3.*m)*V1 -
	           (a*a*sn*cn*dn)/(1.-a*sn*sn)
            );
     double ab = a-b;
	return 1./sqr(a)*( sqr(b)*u + 2.*b*ab*V1 + ab*ab*V2 );
}


DEVICEFUNC INLINE
double integral_Rm1(double a, double u, double m)
{
 #ifndef CUDA
	if (m>1.0) error("Error in routine integral_Rm1: m>1.0");
	
 #endif

	return u + a/sqrt(m)*acos(jacobi_dn(u,m));
}


DEVICEFUNC
double integral_Rm2(double a, double u, double m)
{
 #ifndef CUDA
	if (m>1.0) error("Error in routine integral_Rm2: m>1.0");
	
 #endif
 double a2 = sqr(a);
	double sn, cn, dn;
	jacobi_sncndn(u, m, &sn, &cn, &dn);
	return 1/m*( (m-a2*(1.-m))*u + a2*elliptic_e_cos(cn,m) + 2*a*sqrt(m)*acos(dn) );
}


DEVICEFUNC INLINE
double integral_R0(double u, double m)
{
	return u;
}


DEVICEFUNC
double integral_R1(double a, double u, double m)
{
 #ifndef CUDA
	if (m>1.0) error("Error in routine integral_R1: m>1.0");
	
 #endif

	double a2 = sqr(a);
	double n = a2/(a2-1.);
	double sn, cn, dn;
	jacobi_sncndn(u, m, &sn, &cn, &dn);
	double mma = (m + (1.-m)*a2) / (1.-a2);
	sim5complex f1 = (fabs(mma)>1e-5) ? csqrt(makeComplex(1./mma,0.0))*catan(csqrt(makeComplex(mma,0.0))*sn/dn) : makeComplex(sn/dn,0.0);
	sim5complex ellpi;
	#ifdef CUDA
	ellpi.x = elliptic_pi_cos(cn, n, m); ellpi.y = 0.0;
	#else
	ellpi = elliptic_pi_cos(cn, n, m);
	#endif
	sim5complex res   = 1./(1.-a2) * (ellpi + a*f1);
	/*
 fprintf(stderr, "R1 - a   = %f\n",a);
	fprintf(stderr, "R1 - u   = %f\n",u);
	fprintf(stderr, "R1 - m   = %f\n",m);
 fprintf(stderr, "R1 - f1  = %e + %e I\n", creal(f1), cimag(f1));
 fprintf(stderr, "R1 - af1  = %e + %e I\n", a*creal(f1), a*cimag(f1));
 fprintf(stderr, "R1 - epi = %e + %e I\n", creal(ellpi), cimag(ellpi));
 fprintf(stderr, "R1 - res = %e + %e I\n", creal(res), cimag(res));
 fprintf(stderr, "R1-tot = %e + %e I\n", -1.114784e-01-creal(res), cimag(res));
  */
	return creal(res);
}


DEVICEFUNC
double integral_R2(double a, double u, double m)
{
 #ifndef CUDA
	if (m>1.0) error("Error in routine integral_R2: m>1.0");
	
 #endif

 double a2  = sqr(a);
	double mma = (m + (1.-m)*a2);
	double sn, cn, dn;
	jacobi_sncndn(u, m, &sn, &cn, &dn);

	return 1/(a2-1.)/mma * (
	 (a2*(2.*m-1.)-2.*m)*integral_R1(a,u,m) +
	 2.*m*integral_Rm1(a,u,m) -
	 m*integral_Rm2(a,u,m) +
	 a*a2*sn*dn/(1.+a*cn)
 );
}







DEVICEFUNC
double integral_R_r0_re(double a, double b, double c, double d, double X)
{
 #ifndef CUDA
	if ((a<=b)||(b<=c)||(c<=d)) error("integral_R_r0_re: invalid input");
 #endif
	double m4 = ((b-c)*(a-d))/((a-c)*(b-d));
	double sn = sqrt(((b-d)*(X-a))/((a-d)*(X-b)));
	
	return 2.0/sqrt((a-c)*(b-d)) * jacobi_isn(sn, m4);    
}


DEVICEFUNC
double integral_R_r0_re_inf(double a, double b, double c, double d)
{
 #ifndef CUDA
	if ((a<=b)||(b<=c)||(c<=d)) error("integral_R_r0_re: invalid input");
 #endif
	double m4 = ((b-c)*(a-d))/((a-c)*(b-d));
	double sn = sqrt((b-d)/(a-d));
	
	return 2.0/sqrt((a-c)*(b-d)) * jacobi_isn(sn, m4);
}


DEVICEFUNC
double integral_R_r0_cc(double a, double b, sim5complex c, double X)
{
	double u  = creal(c);
	double v2 = sqr(cimag(c));
	double A  = sqrt(sqr(a-u) + v2);
	double B  = sqrt(sqr(b-u) + v2);
	double m2 = (sqr(A+B) - sqr(a-b))/(4.*A*B);
	double cn = (X*(A-B)+a*B-b*A) / (X*(A+B)-a*B-b*A);
	
	return 1./sqrt(A*B) * jacobi_icn(cn, m2);

}


DEVICEFUNC
double integral_R_r0_cc_inf(double a, double b, sim5complex c)
{
	double u  = creal(c);
	double v2 = sqr(cimag(c));
	double A  = sqrt(sqr(a-u) + v2);
	double B  = sqrt(sqr(b-u) + v2);
	double m2 = (sqr(A+B) - sqr(a-b))/(4.*A*B);
	double cn = (A-B)/(A+B);
	
	return 1./sqrt(A*B) * jacobi_icn(cn, m2);
}


DEVICEFUNC
double integral_R_r1_re(double a, double b, double c, double d, double X)
{
	double m2 = ((b-c)*(a-d))/((a-c)*(b-d));
	double sn = sqrt(((b-d)*(X-a))/((a-d)*(X-b)));
	double u  = jacobi_isn(sn, m2);
	double a2 = (a-d)/(b-d);
	double b2 = ((a-d)*b)/(a*(b-d));
	double Z  = integral_Z1(a2,b2,u,m2) - integral_Z1(a2,b2,0,m2);
	return a*2.0/sqrt((a-c)*(b-d)) * Z;
}


DEVICEFUNC
double integral_R_r1_cc(double a, double b, sim5complex c, double X1, double X2)
{
	double u  = creal(c);
	double v2 = sqr(cimag(c));
	double A  = sqrt(sqr(a-u) + v2);
	double B  = sqrt(sqr(b-u) + v2);
	double m  = (sqr(A+B) - sqr(a-b))/(4.*A*B);
	double g  = 1./sqrt(A*B);
	double alpha1 = (B*a+b*A)/(B*a-b*A);
	double alpha2 = (B+A)/(B-A);

	double u1 = elliptic_f_cos((X1*(A-B)+a*B-b*A)/(X1*(A+B)-a*B-b*A),m);
	double u2 = elliptic_f_cos((X2*(A-B)+a*B-b*A)/(X2*(A+B)-a*B-b*A),m);

	double t0 = alpha1 * (integral_R0(u2,m) - integral_R0(u1,m));
	double t1 = (alpha2-alpha1) * (integral_R1(alpha2,u2,m) - integral_R1(alpha2,u1,m));

	return (B*a-b*A)/(B+A)*g * (t0+t1);

/*
 double cr = creal(c);
	double ci = cimag(c);

	double func(double x) {
		return x/sqrt((x-a)*(x-b)*(x*x-2.*x*cr+cr*cr+ci*ci));
	}

	if (X1<a) error("integral_R_r2_cc: X1<a");
	if (X1<b) error("integral_R_r2_cc: X1<b");
	if (X1>X2) error("integral_R_r2_cc: X1>X2");

	double qromb(double (*func)(double), double a, double b);
	double x = qromb(&func, X1, X2);
	
	return x;
*/

}


DEVICEFUNC
double integral_R_r2_re(double a, double b, double c, double d, double X)
{
	double m2 = ((b-c)*(a-d))/((a-c)*(b-d));
	double sn = sqrt(((b-d)*(X-a))/((a-d)*(X-b)));
	double u  = jacobi_isn(sn, m2);
	double a2 = (a-d)/(b-d);
	double b2 = ((a-d)*b)/(a*(b-d));
	double Z  = integral_Z2(a2,b2,u,m2) - integral_Z2(a2,b2,0,m2);
	return sqr(a)*2.0/sqrt((a-c)*(b-d)) * Z;
}


DEVICEFUNC
double integral_R_r2_cc(double a, double b, sim5complex c, double X1, double X2)
{
	double u  = creal(c);
	double v2 = sqr(cimag(c));
	double A  = sqrt(sqr(a-u) + v2);
	double B  = sqrt(sqr(b-u) + v2);
	double m  = (sqr(A+B) - sqr(a-b))/(4.*A*B);
	double g  = 1./sqrt(A*B);
	double alpha1 = (B*a+b*A)/(B*a-b*A);
	double alpha2 = (B+A)/(B-A);

	double u1 = elliptic_f_cos((X1*(A-B)+a*B-b*A)/(X1*(A+B)-a*B-b*A),m);
	double u2 = elliptic_f_cos((X2*(A-B)+a*B-b*A)/(X2*(A+B)-a*B-b*A),m);

	double t0 = pow(alpha1,2.)         * (integral_R0(u2,m)        - integral_R0(u1,m));
	double t1 = 2.*alpha1*(alpha2-alpha1) * (integral_R1(alpha2,u2,m) - integral_R1(alpha2,u1,m));
	double t2 = pow(alpha2-alpha1,2.)  * (integral_R2(alpha2,u2,m) - integral_R2(alpha2,u1,m));

	return pow((B*a-b*A)/(B+A),2.)*g * (t0+t1+t2);
/*
	double cr = creal(c);
	double ci = cimag(c);

	double func(double x) {
		return x*x/sqrt((x-a)*(x-b)*(x*x-2.*x*cr+cr*cr+ci*ci));
	}

	if (X1<a) error("integral_R_r2_cc: X1<a");
	if (X1<b) error("integral_R_r2_cc: X1<b");
	if (X1>X2) error("integral_R_r2_cc: X1>X2");

	double qromb(double (*func)(double), double a, double b);
	double x = qromb(&func, X1, X2);
	
	return x;
*/


}


DEVICEFUNC
double integral_R_rp_re(double a, double b, double c, double d, double p, double X)
{
	double m2 = ((b-c)*(a-d))/((a-c)*(b-d));
	double sn = sqrt(((b-d)*(X-a))/((a-d)*(X-b)));
	double u1 = jacobi_isn(sn, m2);
	double a2 = (a-d)/(b-d);
	double c2 = ((p-b)*(a-d))/((p-a)*(b-d));
	return -2.0/sqrt((a-c)*(b-d))/(p-a) * (integral_Z1(c2,a2,u1,m2) - integral_Z1(c2,a2,0,m2));
}


DEVICEFUNC
double integral_R_rp_re_inf(double a, double b, double c, double d, double p)
{
	double m2 = ((b-c)*(a-d))/((a-c)*(b-d));
	double sn = sqrt((b-d)/(a-d));
	double u1 = jacobi_isn(sn, m2);
	double a2 = (a-d)/(b-d);
	double c2 = ((p-b)*(a-d))/((p-a)*(b-d));
	return -2.0/sqrt((a-c)*(b-d))/(p-a) * (integral_Z1(c2,a2,u1,m2) - integral_Z1(c2,a2,0,m2));
}


DEVICEFUNC
double integral_R_rp_cc2(double a, double b, sim5complex c, double p, double X1, double X2)
{
 #ifndef CUDA
 if (X1<a) error("integral_R_rp_cc2: X1<a (X1=%e, a=%e)", X1, a);
 if (X1<b) error("integral_R_rp_cc2: X1<b (X1=%e, b=%e)", X1, b);
 if (X1<p) error("integral_R_rp_cc2: X1<p (X1=%e, p=%e)", X1, p);
 #endif

	double u  = creal(c);
	double v2 = sqr(cimag(c));
	double A  = sqrt(sqr(a-u) + v2);
	double B  = sqrt(sqr(b-u) + v2);
	double m  = (sqr(A+B) - sqr(a-b))/(4.*A*B);
	double g  = 1./sqrt(A*B);
	double alpha1 = (B*a+b*A-p*A-p*B)/(B*a-b*A+p*A-p*B);
	double alpha2 = (B+A)/(B-A);


	double u1 = elliptic_f_cos((X1*(A-B)+a*B-b*A)/(X1*(A+B)-a*B-b*A),m);
	double u2 = elliptic_f_cos((X2*(A-B)+a*B-b*A)/(X2*(A+B)-a*B-b*A),m);

	double t0 = alpha2 * (integral_R0(u2,m) - integral_R0(u1,m));
	double t1 = (alpha1-alpha2) * (integral_R1(alpha1,u2,m) - integral_R1(alpha1,u1,m));

	return (B-A)*g/(B*a+b*A-p*A-p*B) * (t0+t1);
}


DEVICEFUNC
double integral_R_rp_cc2_inf(double a, double b, sim5complex c, double p, double X1)
{
 #ifndef CUDA
 if (X1<a) error("integral_R_rp_cc2_inf: X1<a (X1=%e, a=%e)", X1, a);
 if (X1<b) error("integral_R_rp_cc2_inf: X1<b (X1=%e, b=%e)", X1, b);
 if (X1<p) error("integral_R_rp_cc2_inf: X1<p (X1=%e, p=%e)", X1, p);
 #endif

	double u  = creal(c);
	double v2 = sqr(cimag(c));
	double A  = sqrt(sqr(a-u) + v2);
	double B  = sqrt(sqr(b-u) + v2);
	double m  = (sqr(A+B) - sqr(a-b))/(4.*A*B);
	double g  = 1./sqrt(A*B);
	double alpha1 = (B*a+b*A-p*A-p*B)/(B*a-b*A+p*A-p*B);
	double alpha2 = (B+A)/(B-A);


	double u1 = elliptic_f_cos((X1*(A-B)+a*B-b*A)/(X1*(A+B)-a*B-b*A),m);
	double u2 = elliptic_f_cos((A-B)/(A+B),m);

	double t0 = alpha2 * (integral_R0(u2,m) - integral_R0(u1,m));
	double t1 = (alpha1-alpha2) * (integral_R1(alpha1,u2,m) - integral_R1(alpha1,u1,m));

	return (B-A)*g/(B*a+b*A-p*A-p*B) * (t0+t1);
}





DEVICEFUNC
double integral_T_m0(double a2, double b2, double X)
{
	double m = b2/(a2+b2);
	return 1./sqrt(a2+b2) * jacobi_icn(X/sqrt(b2),m);
}


DEVICEFUNC
double integral_T_m2(double a2, double b2, double X)
{
	double m  = b2/(a2+b2);
	double cn = X/sqrt(b2);
 return b2/sqrt(a2+b2) * (integral_C2_cos(cn,m) - integral_C2(0,m));
}


DEVICEFUNC
double integral_T_mp(double a2, double b2, double p, double X)
{
 #ifndef CUDA
	if (fabs(X) > sqrt(b2)) error("integral_T_mp: invalid input ((X<0)||(X>b))");
	if (p==b2) error("integral_T_mp: invalid input (p==b2)");
 #endif
	double m = b2/(a2+b2);
 double n = b2/(b2-p);

 if (X >= 0.0)
     return 1./sqrt(a2+b2)/(p-b2) * elliptic_pi_cos(X/sqrt(b2), n, m);
 else
     return 1./sqrt(a2+b2)/(p-b2) * (2.*elliptic_pi_complete(n, m) - elliptic_pi_cos(-X/sqrt(b2), n, m));
}


DEVICEFUNC
int quadratic_eq(double pr, double pi, double qr, double qi, double zr[2], double zi[2])
{
/*
  complex bb, bb2, cc, del, qq, z1, z2;
  double tst;
  int i, nu_re=0;

  bb  = cmake(pr, +pi);
  bb2 = ccc(bb);
  cc  = cmake(qr, qi);

  del = csub(csqr(bb), cmulr(cc,4.0));
  tst = cre(cmulc(bb2,csqrt(del)));

  if (tst>=0.0)
 qq=cmulr(cadd(bb,csqrt(del)),-0.5);
  else
 qq=cmulr(csub(bb,csqrt(del)),-0.5);

  z1 = qq;
  z2 = cdiv(cc,qq);

  zr[0]=z1.re;
  zi[0]=z1.im;
  zr[1]=z2.re;
  zi[1]=z2.im;

  for (i=0;i<2;i++) {
	 if (zi[i]==0.0) ++nu_re;
  }

  return nu_re;
*/
  
  
  sim5complex bb, bb2, cc, del, qq, z1, z2;
  double tst;
  int i, nu_re=0;

  
  
  
  bb  = makeComplex(pr, pi);
  bb2 = makeComplex(pr, -pi);
  cc  = makeComplex(qr, qi);

  del=bb*bb-4.*cc;
  tst=creal(bb2*csqrt(del));

  if (tst>=0.0)
 qq=-0.5*(bb+csqrt(del));
  else
 qq=-0.5*(bb-csqrt(del));

  z1=qq;
  z2=cc/qq;

  zr[0]=creal(z1);
  zi[0]=cimag(z1);
  zr[1]=creal(z2);
  zi[1]=cimag(z2);

  for (i=0;i<2;i++)
 {
	 if (zi[i]==0.0)
	   ++nu_re;
 }
  return nu_re;
}


DEVICEFUNC
int cubic_eq(double p, double q, double r, double zr[3], double zi[3])
{
  double x1, x2, x3, y1, y2, y3;
  double theta, aa, bb, qq, rr;
  int sgnrr, i, nu_re=0;

  qq=(p*p-3.*q)/9.;
  rr=(2*p*p*p-9.*p*q+27*r)/54.;

  if (rr>=0.0)
 sgnrr=1;
  else
 sgnrr=-1;

  if ((rr*rr)<(qq*qq*qq))
 {
	 theta=acos(rr/sqrt(qq*qq*qq));

	 x1=-2*sqrt(qq)*cos(theta/3.)-p/3.;
	 x2=-2*sqrt(qq)*cos((theta+2*M_PI)/3.)-p/3.;
	 x3=-2*sqrt(qq)*cos((theta-2*M_PI)/3.)-p/3.;

	 y1=0.0;
	 y2=0.0;
	 y3=0.0;
 }
  else
 {
	 aa=-sgnrr*pow(fabs(rr)+sqrt(rr*rr-qq*qq*qq),1/3.);
	 if (aa!=0.0)
	   bb=qq/aa;
	 else
	   bb=0.0;

	 x1=(aa+bb)-p/3.;
	 x2=-0.5*(aa+bb)-p/3;
	 x3=x2;

	 y1=0;
	 y2=0.5*sqrt(3.)*(aa-bb);
	 y3=-y2;
 }

  zr[0]=x1;
  zi[0]=y1;
  zr[1]=x2;
  zi[1]=y2;
  zr[2]=x3;
  zi[2]=y3;

  for (i=0;i<3;i++)
 {
	 if (zi[i]==0.0)
	   ++nu_re;
 }
  return nu_re;
}


DEVICEFUNC
void sort_roots_re(double *r1, double *r2, double *r3, double *r4)
{
	double t;
	if (*r2 > *r1) { t=*r1; *r1=*r2; *r2=t; }
	if (*r3 > *r1) { t=*r1; *r1=*r3; *r3=t; }
	if (*r4 > *r1) { t=*r1; *r1=*r4; *r4=t; }
	if (*r3 > *r2) { t=*r2; *r2=*r3; *r3=t; }
	if (*r4 > *r2) { t=*r2; *r2=*r4; *r4=t; }
	if (*r4 > *r3) { t=*r3; *r3=*r4; *r4=t; }
}



DEVICEFUNC
void sort_mix(double *r1, double *r2, int *s)
{
  const int N=4;
  double rr1[N], rr2[N], tmp;
  int i, j, k=0;

  *s=0;

  for (i=0; i<N; i++) {
 if (r2[i] != 0.) {
     rr1[*s]=r1[i];
     rr2[*s]=r2[i];
     *s+=1;
	}
   else
	{
     rr1[N-1-k]=r1[i];
     rr2[N-1-k]=0.0;
     k+=1;
	}
  }

  for (i=0; i<N; i++){
   r1[i]=rr1[i];
   r2[i]=rr2[i];
  }

  for (i=*s; i<N; i++) {
   for (j=0; j<(N-i); j++) {
     if (r1[i+j]>r1[i]) {
         tmp=r1[i+j];
         r1[i+j]=r1[i];
         r1[i]=tmp;
	 }
   }
 }
}


DEVICEFUNC
void sort_mix2(double *r1, double *r2, int *s)
{
  const int N=4;
  double rr1[N], rr2[N], tmp;
  int i, j, k=0;

  *s=0;

  for (i=0; i<N; i++) {
 if (r2[i] == 0.) {
     rr1[*s]=r1[i];
     rr2[*s]=0.0;
     *s+=1;
	}
   else
	{
     rr1[N-1-k]=r1[i];
     rr2[N-1-k]=r2[i];
     k+=1;
	}
  }

  for (i=0; i<N; i++){
   r1[i]=rr1[i];
   r2[i]=rr2[i];
  }

  for (i=0; i<*s; i++) {
   for (j=0; j<(*s-i); j++) {
     if (r1[i+j]>r1[i]) {
         tmp=r1[i+j];
         r1[i+j]=r1[i];
         r1[i]=tmp;
	 }
   }
 }
}



DEVICEFUNC
void sort_roots(int *s, sim5complex *z1, sim5complex *z2, sim5complex *z3, sim5complex *z4)
{
  const int N=4;
  
  sim5complex rr1[N], rr2[N], tmp;
  int i, j, k=0;

  *s=0;

  rr1[0] = *z1;
  rr1[1] = *z2;
  rr1[2] = *z3;
  rr1[3] = *z4;

  
  for (i=0; i<N; i++) {
   if (cimag(rr1[i])==0.) {
       rr2[*s] = rr1[i];
       *s += 1;
   }
  }

  
  k=*s;
  for (i=0; i<N; i++) {
   if (cimag(rr1[i])!=0.) {
       rr2[k] = rr1[i];
       k += 1;
   }
  }

  for (i=0; i<*s; i++) {
   for (j=0; j<(*s-i); j++) {
     if (creal(rr2[i+j])>creal(rr2[i])) {
         tmp=rr2[i+j];
         rr2[i+j]=rr2[i];
         rr2[i]=tmp;
	 }
   }
 }

 *z1 = rr2[0];
 *z2 = rr2[1];
 *z3 = rr2[2];
 *z4 = rr2[3];
}



DEVICEFUNC
int quartic_eq(double a3, double a2, double a1, double a0, double zr[4], double zi[4])
{
  double u1, pp, qq, rr, sup, del;
  double x1, x2, x3, x4, y1, y2, y3, y4;
  double zr2[2], zi2[2], zr3[3], zi3[3];
  int i, nu_re;

  pp=a2-3.*a3*a3/8.;
  qq=a1-0.5*a2*a3+0.125*pow(a3,3);
  rr=a0-0.25*a1*a3+a2*a3*a3/16.-3.*pow(a3,4)/256.;

  if (qq!=0)                  
 {
   cubic_eq(-pp,-4.*rr,4.*rr*pp-qq*qq,zr3,zi3);
   u1=zr3[0];                 

   if (u1-pp>0.0)
	{
	  sup=sqrt(u1-pp);

	  quadratic_eq(sup,0,0.5*(u1-qq/sup),0,zr2,zi2);
	  x1=zr2[0];
	  x2=zr2[1];
	  y1=zi2[0];
	  y2=zi2[1];

	  quadratic_eq(-sup,0,0.5*(u1+qq/sup),0,zr2,zi2);
	  x3=zr2[0];
	  x4=zr2[1];
	  y3=zi2[0];
	  y4=zi2[1];
	}

   else
	{
	  sup=sqrt(pp-u1);

	  quadratic_eq(0,sup,0.5*u1,0.5*qq/sup,zr2,zi2);
	  x1=zr2[0];
	  x2=zr2[1];
	  y1=zi2[0];
	  y2=zi2[1];

	  quadratic_eq(0,-sup,0.5*u1,-0.5*qq/sup,zr2,zi2);
	  x3=zr2[0];
	  x4=zr2[1];
	  y3=zi2[0];
	  y4=zi2[1];
	}
 }
  else                        
 {
   del=pp*pp-4.*rr;

   if (del>=0)
	{
	  quadratic_eq(0,0,0.5*(pp-sqrt(del)),0,zr2,zi2);
	  x1=zr2[0];
	  x2=zr2[1];
	  y1=zi2[0];
	  y2=zi2[1];

	  quadratic_eq(0,0,0.5*(pp+sqrt(del)),0,zr2,zi2);
	  x3=zr2[0];
	  x4=zr2[1];
	  y3=zi2[0];
	  y4=zi2[1];
	}
   else
	{
	  quadratic_eq(0,0,0.5*pp,-0.5*sqrt(-del),zr2,zi2);
	  x1=zr2[0];
	  x2=zr2[1];
	  y1=zi2[0];
	  y2=zi2[1];

	  quadratic_eq(0,0,0.5*pp,0.5*sqrt(-del),zr2,zi2);
	  x3=zr2[0];
	  x4=zr2[1];
	  y3=zi2[0];
	  y4=zi2[1];
	}
 }

  zr[0]=x1-0.25*a3;
  zi[0]=y1;
  zr[1]=x2-0.25*a3;
  zi[1]=y2;
  zr[2]=x3-0.25*a3;
  zi[2]=y3;
  zr[3]=x4-0.25*a3;
  zi[3]=y4;

  nu_re=0;
  for (i=0;i<4;i++)
 {
	 if (zi[i]==0.0)
	   ++nu_re;
 }

  if (nu_re == 4) sort_roots_re(&zr[0], &zr[1], &zr[2], &zr[3]);
  else
  if (nu_re == 2) sort_mix2(zr, zi, &i);

  return nu_re;
}


DEVICEFUNC
void quartic_eq_c(
	double a3, double a2, double a1, double a0,
	int *nr,
	sim5complex *z1, sim5complex *z2, sim5complex *z3, sim5complex *z4)
{
	double zr[4];
	double zi[4];
	*nr = quartic_eq(a3,a2,a1,a0,zr,zi);
	*z1 = makeComplex(zr[0], zi[0]);
	*z2 = makeComplex(zr[1], zi[1]);
	*z3 = makeComplex(zr[2], zi[2]);
	*z4 = makeComplex(zr[3], zi[3]);
}








/*
#include "sim5config.h"
#ifndef CUDA
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sim5config.h"
#include "sim5math.h"
#include "sim5kerr.h"
#include "sim5polyroots.h"
#include "sim5photontrace.h"
#endif
*/


#define frac_error(a,b) (fabs((b)-(a))/(fabs(b)+1e-40))
#define vect_copy(v1,v2) {v2[0]=v1[0]; v2[1]=v1[1]; v2[2]=v1[2]; v2[3]=v1[3];}

#define raytrace_max_error 1e-2


DEVICEFUNC
void raytrace_prepare(double bh_spin, double x[4], double k[4], double presision_factor, int options, raytrace_data* rtd)
{
 
 rtd->opt_gr  = !((options & RTOPT_FLAT) == RTOPT_FLAT);
 rtd->step_epsilon = sqrt(presision_factor)/10.;   

 
 sim5metric m;
 double G[4][4][4];
 rtd->opt_gr ? kerr_metric(bh_spin, x[1], x[2], &m) : flat_metric(x[1], x[2], &m);
 rtd->opt_gr ? kerr_connection(bh_spin, x[1], x[2], G) : flat_connection(x[1], x[2], G);

 
 double kk = dotprod(k, k, &m);
 #ifndef CUDA
 if (fabs(kk) > 1e-10) fprintf(stderr,"ERR (kerr_raytrace_prepare): k is not a null vector (k.k=%.3e)\n", kk);
 #endif

 
 rtd->bh_spin = bh_spin;
 rtd->E = k[0]*m.g00 + k[3]*m.g03;
 rtd->Q  = photon_carter_const(k, &m);

 
 rtd->pass = 0;
 rtd->refines = 0;
 rtd->kt = rtd->E;
 rtd->error = 0.0;

 
 Gamma(G, k, k, rtd->dk);
}

#ifdef CUDA
DEVICEFUNC
inline double k_deriv(int j, double _k[4], double G[4][4][4]) {
 int a,b;
 double _dki = 0.0;
 for (a=0;a<4;a++) for (b=a;b<4;b++) _dki -= G[j][a][b]*_k[a]*_k[b];
 return _dki;
}
#endif


DEVICEFUNC
void raytrace(double x[4], double k[4], double *step, raytrace_data* rtd)
{
 int i;
 
 sim5metric m;
 double G[4][4][4];
 double* dk = rtd->dk;
 double x_orig[4], k_orig[4];
 double xp[4], kp[4], kp_prev[4];
 double kk=0.0, kt = rtd->kt;
 float k_frac_error;

 vect_copy(x, x_orig);
 vect_copy(k, k_orig);

 #ifndef CUDA
 
 
 
 
 
 inline double k_deriv(int j, double _k[4]) {
     int a,b;
     double _dki = 0.0;
     for (a=0;a<4;a++) for (b=a;b<4;b++) _dki -= G[j][a][b]*_k[a]*_k[b];
     return _dki;
 }
 #endif

 
 
 
 
 
 double stepsize = rtd->step_epsilon/(fabs(dk[0])/(fabs(k[0])+TINY) + fabs(dk[1])/(fabs(k[1])+TINY) + fabs(dk[2])/(fabs(k[2])+TINY) + fabs(dk[3])/(fabs(k[3])+TINY) + TINY);
 double dl = fmin(*step, stepsize);
 if (dl < 1e-3) dl = 1e-3;

 rtd->pass++;

 
 
 
 double half_dl  = 0.5 * dl;
 double half_dl2 = 0.5 * dl * dl;
 xp[0] = x[0] + k[0]*dl + dk[0]*half_dl2;
 xp[1] = x[1] + k[1]*dl + dk[1]*half_dl2;
 xp[2] = cos(acos(x[2]) +(k[2]*dl + dk[2]*half_dl2));
 xp[3] = x[3] + k[3]*dl + dk[3]*half_dl2;

 
 for (i=0;i<4;i++) k[i] += dk[i]*half_dl;

 
 if (rtd->opt_gr) {
     kerr_metric(rtd->bh_spin, xp[1], xp[2], &m);
     kerr_connection(rtd->bh_spin, xp[1], xp[2], G);
 } else {
     flat_metric(xp[1], xp[2], &m);
     flat_connection(xp[1], xp[2], G);
 };

 
 for (i=0;i<4;i++) kp[i] = k[i] + dk[i]*half_dl;


 
 int k_iter = 0;
 do {
     k_frac_error = 0.0;

     vect_copy(kp, kp_prev);
     for (i=0;i<4;i++) {
         
         #ifdef CUDA
         kp[i] = k[i] + k_deriv(i,kp_prev,G)*half_dl;
         #else
         kp[i] = k[i] + k_deriv(i,kp_prev)*half_dl;
         #endif
         k_frac_error += frac_error(kp[i],kp_prev[i]);
     }

     k_iter++;
	} while (k_frac_error>raytrace_max_error*1e-3 && k_iter<3);


 
 kt = kp[0]*m.g00 + kp[3]*m.g03;
 kk = fabs(dotprod(kp, kp, &m));
 rtd->error = fmax(frac_error(kt,rtd->kt), kk);
 if ((k_frac_error>raytrace_max_error*1e-2) || (rtd->error>raytrace_max_error*1e-2)) {
     vect_copy(x_orig, x);
     vect_copy(k_orig, k);
     DEVICEFUNC void raytrace_rk4(double x[4], double k[4], double dl, raytrace_data* rtd);
     raytrace_rk4(x, k, dl, rtd);
     *step = dl;
     return;
 }

 
 for (i=0;i<4;i++) {
     x[i]  = xp[i];
     k[i]  = kp[i];
     #ifdef CUDA
     dk[i] = k_deriv(i,kp,G);
     #else
     dk[i] = k_deriv(i,kp);
     #endif
 }

 
 rtd->kt = kt;

 
 *step = dl;
}



DEVICEFUNC
void raytrace_rk4(double x[4], double k[4], double dl, raytrace_data* rtd)
{
 int i;
 sim5metric m;
 double G[4][4][4];
 double xp[4];
 double k1[4], dk1[4], k2[4], dk2[4], k3[4], dk3[4], k4[4], dk4[4];
	double dl_2 = 0.5 * dl;

 
 
 
 


 double kt0 = rtd->kt;

 
 x[2] = acos(x[2]);

	for (i=0; i<4; i++) xp[i] = x[i];
 rtd->opt_gr ? kerr_connection(rtd->bh_spin, xp[1], cos(xp[2]), G) : flat_connection(xp[1], cos(xp[2]), G);
	for (i=0; i<4; i++) k1[i] = k[i];
 Gamma(G, k1, k1, dk1);

	for (i=0; i<4; i++) xp[i] = x[i] + k1[i]*dl_2;
 rtd->opt_gr ? kerr_connection(rtd->bh_spin, xp[1], cos(xp[2]), G) : flat_connection(xp[1], cos(xp[2]), G);
	for (i=0; i<4; i++) k2[i] = k[i] + dk1[i]*dl_2;
 Gamma(G, k2, k2, dk2);

	for (i=0; i<4; i++) xp[i] = x[i] + k2[i]*dl_2;
 rtd->opt_gr ? kerr_connection(rtd->bh_spin, xp[1], cos(xp[2]), G) : flat_connection(xp[1], cos(xp[2]), G);
	for (i=0; i<4; i++) k3[i] = k[i] + dk2[i]*dl_2;
 Gamma(G, k3, k3, dk3);

	for (i=0; i<4; i++) xp[i] = x[i] + k3[i]*dl;
 rtd->opt_gr ? kerr_connection(rtd->bh_spin, xp[1], cos(xp[2]), G) : flat_connection(xp[1], cos(xp[2]), G);
	for (i=0; i<4; i++) k4[i] = k[i] + dk3[i]*dl;
 Gamma(G, k4, k4, dk4);

 
	for (i=0; i<4; i++) {
 	x[i] += dl/6. * ( k1[i] + 2.* k2[i] + 2.* k3[i] +  k4[i]);
 	k[i] += dl/6. * (dk1[i] + 2.*dk2[i] + 2.*dk3[i] + dk4[i]);
 }

 
	x[2] = cos(x[2]);

 
 rtd->opt_gr ? kerr_connection(rtd->bh_spin, x[1], x[2], G) : flat_connection(x[1], x[2], G);
 Gamma(G, k, k, rtd->dk);


 kerr_metric(rtd->bh_spin, x[1], x[2], &m);
 double kt1 = k[0]*m.g00 + k[3]*m.g03;

 rtd->error = frac_error(kt1,kt0);
 /*
	if (rtd->error>1e-4) {
	 fprintf(stderr,"WRN: RK4 rtd->error=%.3e (%d/%d) dl=%.2e ke=%.2e kp[1]=%.5e kp[2]=%.5e\n", rtd->error, rtd->pass, rtd->refines, dl, 0.0, k[1], k[2]);
	 if (dl>1e-8) {
         vect_copy(x_orig, x);
         vect_copy(k_orig, k);
         
         rtd->refines++;
	     for (i=0; i<10; i++) raytrace_rk4(x, k, f, dl/10., rtd);
	 } else {
	 fprintf(stderr,"WRN: RK4 step too small rtd->error=%.3e (%d/%d) dl=%.2e ke=%.2e kp[1]=%.5e kp[2]=%.5e\n", rtd->error, rtd->pass, rtd->refines, dl, 0.0, k[1], k[2]);
	 }
	}
 */
}


DEVICEFUNC
double raytrace_error(double x[4], double k[4], raytrace_data* rtd)
{
 sim5metric m;
 rtd->opt_gr ? kerr_metric(rtd->bh_spin, x[1], x[2], &m) : flat_metric(x[1], x[2], &m);
 return frac_error(rtd->Q, photon_carter_const(k,&m));
}


#undef frac_error
#undef vect_copy
#undef raytrace_max_error










DEVICEFUNC
void flat_metric(double r, double m, sim5metric *metric)
{
 metric->a = 0.0;
 metric->r = r;
 metric->m = m;
 metric->g00 = -1.0;
 metric->g11 = +1.0;
 metric->g22 = +r*r;
 metric->g33 = +r*r*(1.-m*m);
 metric->g03 = 0.0;
}



DEVICEFUNC
void flat_metric_contravariant(double r, double m, sim5metric *metric)
{
 metric->a = 0.0;
 metric->r = r;
 metric->m = m;
 metric->g00 = -1.0;
 metric->g11 = +1.0;
 metric->g22 = +1./(r*r);
 metric->g33 = +1./(r*r)/(1.-m*m);
 metric->g03 = 0.0;
}



DEVICEFUNC
void kerr_metric(double a, double r, double m, sim5metric *metric)
{
 double r2  = sqr(r);
 double a2  = sqr(a);
 double m2  = sqr(m);
 double S   = r2 + a2*m2;
 double s2S = (1.0-m2)/S;
 
 
 metric->a = a;
 metric->r = r;
 metric->m = m;
 metric->g00 = -1. + 2.0*r/S;
 metric->g11 = S/(r2-2.*r+a2); 
 metric->g22 = S;
 metric->g33 = (r2*r2 + a2*a2*m2 + a2*r2*(1.+m2) + 2.*r*a2*s2S*S)*s2S;
 metric->g03 = -2.*a*r*s2S;
}



DEVICEFUNC
void kerr_metric_contravariant(double a, double r, double m, sim5metric *metric)
{
 double r2  = sqr(r);
 double a2  = sqr(a);
 double m2  = sqr(m);
 double S   = r2 + a2*m2;
 double SD  = S*(r2 - 2.*r + a2);  
 
 
 metric->a = a;
 metric->r = r;
 metric->m = m;
 metric->g00 = -sqr(r2+a2)/SD + a2*(1.-m2)/S;  
 metric->g11 = (r2-2.*r+a2)/S; 
 metric->g22 = 1./S;
 metric->g33 = 1./S/(1.-m2) - a2/SD;
 metric->g03 = -2.*a*r/SD;
}



DEVICEFUNC
void kerr_connection(double a, double r, double m, double G[4][4][4])
{
 
 
 

 double s  = sqrt(1.-m*m);
 double cs = s*m;
 double c2 = m*m;
 double s2 = s*s;
 double cc = c2-s2; 
 double CC = 8.*c2*c2-8.*c2+1.;   
 double a2 = a*a;
 double a4 = a2*a2;
 double a2cc = a2*cc;
 double a2c2 = a2*c2;
 double a2cs = a2*cs;
 double a4CC = a4*CC;
 double two_r = 2.0*r;
 double r2 = r*r;
 double r3 = r2*r;
 double r4 = r2*r2;
 double a2r2 = a2*r2;
 double R  = pow(a2 + 2.*r2 + a2cc, 2.);
 double D  = r2 - 2.*r + a2;
 double S  = r2 + a2c2;
 double S_1 = 1./S;
 double S_3 = 1./(S*S*S);
 double D_1 = 1./D;
 double R_1 = 1./R;
 double m_s = m/s;
 double DR_1 = D_1*R_1;
 double DS_1 = D_1*S_1;
 double dbl_r2 = 2.*r2;

 memset(G, 0, 4*4*4*sizeof(double));

 G[0][0][1] =  -4.0 * (a2 + r2)*(a2 - dbl_r2 + a2cc)*DR_1;
 G[0][0][2] = -16.0 * a2cs*r*R_1;
 G[0][1][3] =   4.0 * a*s2*(a4 - 3.*a2r2 - 6.*r4 + a2cc*(a2 - r2))*DR_1;
 G[0][2][3] = -G[0][0][2]*s2*a;

 G[1][0][0] = -D*(a2c2-r2)*S_3;
 G[1][0][3] = -2.0 * a * G[1][0][0] * s2;
 G[1][1][1] = (r*(a2 - r) + a2*(1. - r)*c2)*DS_1;
 G[1][1][2] = -2.0*a2cs*S_1;
 G[1][2][2] = -r*D*S_1;
 G[1][3][3] = -D*s2*(2.*a2c2*r3 + r2*r3 + a2*a2c2*s2 + a2c2*a2c2*r - a2r2*s2)*S_3;

 G[2][0][0] = -two_r*a2cs*S_3;
 G[2][0][3] = 4.0 * a*r*cs*(a2 + r2)*S_3;
 G[2][1][1] = a2cs*DS_1;
 G[2][1][2] = two_r*S_1;
 G[2][2][2] = -G[2][1][1]*D;
 G[2][3][3] = -0.125*cs*(3.*a2*a4 + 10.*a4*r + 11.*a4*r2 + 16.*a2*r3 +
         16.*a2*r4 + 8.*r3*r3 + 4.*(a2 + dbl_r2)*D*a2cc +
         a4CC*D)*S_3;

 G[3][0][1] =   2.0 * a*(r2 - a2c2)*DS_1*S_1;
 G[3][0][2] = -16.0 * a*r*m_s*R_1;
 G[3][1][3] = (a4 + 3.*a4*r - 12.*a2r2 + 8.*a2*r3 -
         16.*r4 + 8.*r2*r3 + 4.*r*(dbl_r2 -r + a2)*a2cc -
         a4CC*(1. - r))*DR_1;
 G[3][2][3] = ((3.*a4 + 8.*a2*r + 8.*a2r2 + 8.*r4 +
         4.*(dbl_r2 -2.*r + a2)*a2cc + a4CC)*m_s)*R_1;

 
 
}


DEVICEFUNC
void flat_connection(double r, double m, double G[4][4][4])
{
 double s = sqrt(1.-m*m);
 memset(G, 0, 4*4*4*sizeof(double));

 G[1][2][2] = -r;
 G[1][3][3] = -r*s*s;

 G[2][1][2] = 2.0 * 1./r;
 G[2][3][3] = -m*s;

 G[3][1][3] = 2.0 * 1./r;
 G[3][2][3] = 2.0 * m/s;
}

/*
DEVICEFUNC INLINE
void Gamma(double G[4][4][4], double V[4], double result[4])
{
 int i,j,k;
 for (i=0;i<4;i++) {
     result[i] = 0.0;
     for (j=0;j<4;j++) for (k=j;k<4;k++) result[i] -= G[i][j][k]*V[j]*V[k];
 }
}
*/


DEVICEFUNC INLINE
void Gamma(double G[4][4][4], double U[4], double V[4], double result[4])
{
 int i,j,k;
 for (i=0;i<4;i++) {
     result[i] = 0.0;
     for (j=0;j<4;j++) for (k=j;k<4;k++) result[i] -= 0.5*G[i][j][k]*(U[j]*V[k] + U[k]*V[j]);
 }
}


DEVICEFUNC INLINE
void vector_set(double x[4], double x0, double x1, double x2, double x3)
{
 x[0] = x0;
 x[1] = x1;
 x[2] = x2;
 x[3] = x3;
}


DEVICEFUNC INLINE
void vector_copy(double src[4], double dst[4])
{
 dst[0] = src[0];
 dst[1] = src[1];
 dst[2] = src[2];
 dst[3] = src[3];
}


DEVICEFUNC INLINE
void vector_covariant(double V1[4], double V2[4], sim5metric* m)
{
 if (m) {
     V2[0] = V1[0]*m->g00 + V1[3]*m->g03;
     V2[1] = V1[1]*m->g11;
     V2[2] = V1[2]*m->g22;
     V2[3] = V1[3]*m->g33 + V1[0]*m->g03;
 } else {
     V2[0] = -V1[0];
     V2[1] = +V1[1];
     V2[2] = +V1[2];
     V2[3] = +V1[3];
 }
}


DEVICEFUNC INLINE
double vector_norm(double V[4], sim5metric* m)
{
 return sqrt(dotprod(V, V, m));
}



DEVICEFUNC INLINE
double vector_3norm(double V[4])
{
 return sqrt(V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
}



DEVICEFUNC INLINE
void vector_multiply(double V[4], double factor)
{
 V[0] *= factor;
 V[1] *= factor;
 V[2] *= factor;
 V[3] *= factor;
}



DEVICEFUNC INLINE
void vector_norm_to(double V[4], double norm, sim5metric* m)
{
 double N = dotprod(V, V, m);
 V[0] *= sqrt(norm/N);
 V[1] *= sqrt(norm/N);
 V[2] *= sqrt(norm/N);
 V[3] *= sqrt(norm/N);
}


DEVICEFUNC
void vector_norm_to_null(double V[4], double V0, sim5metric* m)
{
 double a,b,c,alpha;
 if (m) {
     a = V[1]*V[1]*m->g11 + V[2]*V[2]*m->g22 + V[3]*V[3]*m->g33;
     b = V0*V[3]*m->g03;  
     c = V0*V0*m->g00;
     alpha = fmax(-b/a+sqrt(b*b-a*c)/a, -b/a-sqrt(b*b-a*c)/a);
 } else {
     a = V[1]*V[1] + V[2]*V[2] + V[3]*V[3];
     c = -V0*V0;
     alpha = sqrt(-c/a);
 }
 V[0]  = V0;
 V[1] *= alpha;
 V[2] *= alpha;
 V[3] *= alpha;
}



DEVICEFUNC INLINE
double dotprod(double V1[4], double V2[4], sim5metric* m)
{
 if (m)
     return V1[0]*V2[0]*m->g00 + V1[1]*V2[1]*m->g11 + V1[2]*V2[2]*m->g22 +
            V1[3]*V2[3]*m->g33 + V1[0]*V2[3]*m->g03 + V1[3]*V2[0]*m->g03;
 else
     return -V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2] + V1[3]*V2[3];
}



DEVICEFUNC
void tetrad_general(sim5metric *m, double U[], sim5tetrad *t)
{
 
 double u[4];
 double D = sqr(m->r) - 2*(m->r) + sqr(m->a);
 vector_covariant(U, u, m);  

 
 double N1 = sqrt( -m->g11*(u[0]*U[0]+u[3]*U[3])*(1.+u[2]*U[2]) );
 double N2 = sqrt( +m->g22*(1.+u[2]*U[2]) );
 double N3 = sqrt( -(u[0]*U[0]+u[3]*U[3])*D*(1.-sqr(m->m)) );

 t->e[0][0] = U[0];
 t->e[0][1] = U[1];
 t->e[0][2] = U[2];
 t->e[0][3] = U[3];

 t->e[1][0] = u[1]*U[0] / N1;
 t->e[1][1] = -(u[0]*U[0]+u[3]*U[3]) / N1;
 t->e[1][2] = 0.0;
 t->e[1][3] = u[1]*U[3] / N1;

 t->e[2][0] = u[2]*U[0] / N2;
 t->e[2][1] = u[2]*U[0] / N1;
 t->e[2][2] = (1.+u[2]*U[2]) / N2;
 t->e[2][3] = u[2]*U[3] / N2;

 t->e[3][0] = -u[0] / N3;
 t->e[3][1] = 0.0;
 t->e[3][2] = 0.0;
 t->e[3][3] = +u[3] / N3;

 t->metric = *m;
}



DEVICEFUNC
void tetrad_zamo(sim5metric *m, sim5tetrad *t)
{
 t->e[0][0] = sqrt(m->g33/(sqr(m->g03) - m->g33*m->g00));
 t->e[0][1] = 0.0;
 t->e[0][2] = 0.0;
 t->e[0][3] = -t->e[0][0] * m->g03/m->g33;

 t->e[1][0] = 0.0;
 t->e[1][1] = 1./sqrt(m->g11);
 t->e[1][2] = 0.0;
 t->e[1][3] = 0.0;

 t->e[2][0] = 0.0;
 t->e[2][1] = 0.0;
 t->e[2][2] = -1./sqrt(m->g22);
 t->e[2][3] = 0.0;

 t->e[3][0] = 0.0;
 t->e[3][1] = 0.0;
 t->e[3][2] = 0.0;
 t->e[3][3] = 1./sqrt(m->g33);

 t->metric = *m;
}



DEVICEFUNC
void tetrad_radial(sim5metric *m, double v_r, sim5tetrad *t)
{
 if (v_r==0.0) return tetrad_zamo(m, t);

 double g00 = m->g00;
 double g11 = m->g11;
 double U0 = sqrt((-1.-sqr(v_r)*g11)/g00);
 double U1 = v_r;

 t->e[0][0] = U0;
 t->e[0][1] = U1;
 t->e[0][2] = 0.0;
 t->e[0][3] = 0.0;

 double UG = U0*U0*g00 + U1*U1*g11;
 t->e[1][0] = -U1*sqrt(UG*g11*g00)*U0/(g11*UG)*g11/(U0*g00);
 t->e[1][1] = sqrt(UG*g11*g00)*U0/(g11*UG);
 t->e[1][2] = 0.0;
 t->e[1][3] = 0.0;

 t->e[2][0] = 0.0;
 t->e[2][1] = 0.0;
 t->e[2][2] = -1./sqrt(m->g22);  
 t->e[2][3] = 0.0;

 t->e[3][0] = 0.0;
 t->e[3][1] = 0.0;
 t->e[3][2] = 0.0;
 t->e[3][3] = 1./sqrt(m->g33);

 t->metric = *m;
}



DEVICEFUNC
void tetrad_azimuthal(sim5metric *m, double Omega, sim5tetrad *t)
{
 if (Omega==0.0) return tetrad_zamo(m, t);

 double g00 = m->g00;
 double g33 = m->g33;
 double g03 = m->g03;
 double U0 = sqrt(-1.0/(g00 + 2.*Omega*g03 + sqr(Omega)*g33));
 double U3 = U0*Omega;

 t->e[0][0] = U0;
 t->e[0][1] = 0.0;
 t->e[0][2] = 0.0;
 t->e[0][3] = U3;

 t->e[1][0] = 0.0;
 t->e[1][1] = sqrt(1./m->g11);
 t->e[1][2] = 0.0;
 t->e[1][3] = 0.0;

 t->e[2][0] = 0.0;
 t->e[2][1] = 0.0;
 t->e[2][2] = -sqrt(1./m->g22);  
 t->e[2][3] = 0.0;

 double k1 = (g03*U3+g00*U0);
 double k2 = (g33*U3+g03*U0);
 t->e[3][0] =  - sign(k1)*k2 / sqrt((g33*g00-g03*g03)*(g00*U0*U0+g33*U3*U3+2.0*g03*U0*U3));
 t->e[3][1] = 0.0;
 t->e[3][2] = 0.0;
 t->e[3][3] = t->e[3][0] * (-k1/k2);

 t->metric = *m;
}



DEVICEFUNC
void tetrad_surface(sim5metric *m, double Omega, double V, double dhdr, sim5tetrad *t)

{
 double g00 = m->g00;
 double g11 = m->g11;
 double g22 = m->g22;
 double g33 = m->g33;
 double g03 = m->g03;

 
 sim5metric M;
 kerr_metric_contravariant(m->a, m->r, m->m, &M);


 
 
 
 
 
 double S0r = 1.0/sqrt(g11+g22*sqr(dhdr));
 double S0h = S0r*dhdr;

 

 
 double ur = V/sqrt(1.-V*V)/sqrt(g11);
 double v = sign(V) * sqrt((sqr(ur/S0r)*(-g00-2.*Omega*g03-sqr(Omega)*g33))/(1.+sqr(ur/S0r)));

 
 
 t->e[0][0] = 1.0;
 t->e[0][1] = v*S0r;
 t->e[0][2] = v*S0h;
 t->e[0][3] = Omega;
 vector_norm_to(t->e[0], -1.0, m);
 

 
 
 
 
 
 t->e[1][0] = (v*t->e[0][0]);
 t->e[1][1] = (v*t->e[0][1] + S0r/t->e[0][0]);
 t->e[1][2] = (v*t->e[0][2] + S0h/t->e[0][0]);
 t->e[1][3] = (v*t->e[0][3]);
 vector_norm_to(t->e[1], 1.0, m);
 

 
 
 
 
 
 
 
 
 t->e[2][0] = 0.0;
 t->e[2][1] = dhdr;
 t->e[2][2] = -1.0;
 t->e[2][3] = 0.0;
 vector_norm_to(t->e[2], 1.0, m);

 
 
 t->e[3][0] = -(g03+g33*Omega)/(g00+g03*Omega);
 t->e[3][1] = 0.0;
 t->e[3][2] = 0.0;
 t->e[3][3] = 1.0;
 vector_norm_to(t->e[3], 1.0, m);
 

 
 
 
 
 
 

 t->metric = *m;
}



DEVICEFUNC
void bl2on(double Vin[4], double Vout[4], sim5tetrad* t)
{
 Vout[0] = -dotprod(t->e[0], Vin, &t->metric);
 Vout[1] = +dotprod(t->e[1], Vin, &t->metric);
 Vout[2] = +dotprod(t->e[2], Vin, &t->metric);
 Vout[3] = +dotprod(t->e[3], Vin, &t->metric);
}


DEVICEFUNC
void on2bl(double Vin[4], double Vout[4], sim5tetrad* t)
{
 int i,j;
 for (i=0;i<4;i++) {
     Vout[i] = 0.0;
     for (j=0;j<4;j++) Vout[i] += Vin[j] * t->e[j][i];
 }
 
 
 
 
}






DEVICEFUNC INLINE
double r_bh(double a)
{
 return 1. + sqrt(1.-sqr(a));
}



DEVICEFUNC INLINE
double r_ms(double a)
{
 double z1 = 1. + sqrt3(1.-sqr(a))*(sqrt3(1.+a) + sqrt3(1.-a));
 double z2 = sqrt(3.*sqr(a) + sqr(z1));
 return 3.+z2-/*+*/sqrt((3.-z1)*(3.+z1+2.*z2));
}


DEVICEFUNC INLINE
double r_mb(double a)
{
 return (2.-a) + 2.*sqrt(1.-a);
}


DEVICEFUNC INLINE
double r_ph(double a)
{
 return 2.0*(1.0+cos(2./3.*acos(-a)));
}



DEVICEFUNC INLINE
double OmegaK(double r, double a)
{
 return 1./(a + pow(r,1.5));
}


DEVICEFUNC INLINE
double ellK(double r, double a)
{
 /*
 double r2 = sqr(r);
 double a2 = sqr(a);
 double D = r2 + a2 - 2.*r;
 double gtt = -1. + 2.0/r;
 double gff = (sqr(r2 + a2) - a2*D)/r2;
 double gtf = -2.*a/r;

 double Omega = OmegaK(r,a);
 return -(gtf + Omega*gff) / (gtt + Omega*gtf);
 */

 
 return (sqr(r)-2.*a*sqrt(r)+sqr(a)) / (sqrt(r)*r-2.*sqrt(r)+a);
}


DEVICEFUNC INLINE
double omega_r(double r, double a)
{
 return OmegaK(r,a) * sqrt(1.-6./r+8.*a/sqrt(r*r*r)-3.*a*a/sqr(r));
}


DEVICEFUNC INLINE
double omega_z(double r, double a)
{
 return OmegaK(r,a) * sqrt(1.-4.*a/sqrt(r*r*r)+3.*a*a/sqr(r));
}


DEVICEFUNC INLINE
double Omega_from_ell(double ell, sim5metric *m)
{
 return  -(m->g03 + ell*m->g00) / (m->g33 + ell*m->g03);
}


DEVICEFUNC INLINE
double ell_from_Omega(double Omega, sim5metric *m)
{
 return -(m->g03 + m->g33*Omega)/(m->g00 + m->g03*Omega);
}


DEVICEFUNC INLINE
double gfactorK(double r, double a, double l)
{
 double OmegaK = 1./(a + pow(r,1.5));
 return sqrt(1. - 2./r * pow(1.-a*OmegaK,2.) - (r*r+a*a)*sqr(OmegaK)) / (1. - OmegaK*l);
}





DEVICEFUNC
void photon_momentum(double a, double r, double m, double l, double q, double r_sign, double m_sign, double k[4])
{
 double a2 = sqr(a);
 double l2 = sqr(l);
 double r2 = sqr(r);
 double m2 = sqr(m);
 double S = r2 + a2*m2;
 double D = r2 - 2.*r + a2;

 
 double R = sqr(r2+a2-a*l) - D*( sqr(l-a) + q );        
 double M = q - l2*m2/(1.-m2) + a2*m2;                  

 if ((M<0.0) && (-M<1e-8)) M = 0.0;
 if ((R<0.0) && (-R<1e-8)) R = 0.0;

 #ifndef CUDA
 if (R<0.0) error("(photon_momentum): R<0 (R=%.10e)", R);
 if (M<0.0) {
     error("(photon_momentum): M<0 (M=%.10e, l=%.4e  q=%.4e)", M, l, q);
     k[0]=k[1]=k[2]=k[3]=NAN;
     return;
 }
 #endif

 k[0] = +1/S * ( -a*(a*(1.-m2)-l) + (r2+a2)/D*(r2+a2-a*l) );
 k[1] = +1/S * sqrt(R);
 k[2] = +1/S * sqrt(M);
 k[3] = +1/S * ( -a + l/(1.-m2) + a/D*(r2+a2-a*l) );

 if (r_sign<0.0) k[1] = -k[1];
 if (m_sign<0.0) k[2] = -k[2];

 /*
 
 if (dk) {
     int i,b,c;
     double G[4][4][4];
     kerr_connection(a, r, m, G);
     for (i=0;i<4;i++) {
         
         
         dk[i] = 0.0;
         for (b=0;b<4;b++) for (c=b;c<4;c++) dk[i] += -G[i][b][c]*k[b]*k[c];
     }
 }
 */
}


DEVICEFUNC
void photon_motion_constants(double a, double r, double m, double k[4], double* L, double* Q)
{
 double a2 = sqr(a);
 double r2 = sqr(r);
 double s2 = 1.-m*m;
 double D  = r2 - 2.*r + a2;

 double l;
 double nf = k[3]/k[0];
 double nh = sqr(k[2])/sqr(k[0]);

 
 *L = l =(-a*a2 + sqr(a2)*nf + nf*sqr(r2) + a*(D-r2) + a2*nf*(2.*r2-D*s2))*s2 /
     (D - a*s2*(a-a2*nf + nf*(D-r2)));

 
 *Q = pow(a*(l-a*s2) + ((a2+r2)*(a2-a*l+r2))/D, 2.0) *
     (nh - (sqr(D*m)*(sqr(l)-a2*s2))/(-s2*pow(sqr(a2)-a*a2*l+sqr(r2)+a*l*(D-r2)+a2*(2.*r2-D*s2),2.0)));

 #ifndef CUDA
 if (isnan(*L)) warning("ERR (photon_motion_constants): L is NaN (%e, k=%e/%e/%e/%e)\n", *L, k[0], k[1], k[2], k[3]);
 if (isnan(*Q)) warning("ERR (photon_motion_constants): Q is NaN (%e, k=%e/%e/%e/%e)\n", *Q, k[0], k[1], k[2], k[3]);
 #endif
}



DEVICEFUNC
double photon_carter_const(double k[4], sim5metric *metric)
{
 double m2 = sqr(metric->m);
 double kt = k[0]*metric->g00 + k[3]*metric->g03;
 double kh = k[2]*metric->g22;
 double kf = k[3]*metric->g33 + k[0]*metric->g03;
 return sqr(kh) + sqr(kf)*m2/(1.-m2) - sqr(metric->a)*sqr(kt)*m2;
}






DEVICEFUNC INLINE
void fourvelocity_zamo(sim5metric *m, double U[4])
{
 U[0] = sqrt(m->g33/(sqr(m->g03) - m->g33*m->g00));
 U[1] = 0.0;
 U[2] = 0.0;
 U[3] = -U[0] * m->g03/m->g33;
}



DEVICEFUNC INLINE
void fourvelocity_azimuthal(double Omega, sim5metric *m, double U[4])
{
 U[0] = sqrt(-1.0/(m->g00 + 2.*Omega*m->g03 + sqr(Omega)*m->g33));
 U[1] = 0.0;
 U[2] = 0.0;
 U[3] = U[0]*Omega;
}



DEVICEFUNC INLINE
void fourvelocity_radial(double vr, sim5metric *m, double U[4])
{
 U[0] = sqrt((-1.0 - sqr(vr)*m->g11)/m->g00);
 U[1] = vr;
 U[2] = 0.0;
 U[3] = 0.0;
}



DEVICEFUNC INLINE
double fourvelocity_norm(
 double U1, double U2, double U3, sim5metric *m)
{
 double D = sqr(m->g03*U3) - m->g00*m->g11*sqr(U1) - m->g00*m->g22*sqr(U2) - m->g00*m->g33*sqr(U3) - m->g00;
 return (-m->g03*U3-sqrt(D))/m->g00;
}



DEVICEFUNC
void fourvelocity(
 double U1, double U2, double U3, sim5metric *m, double U[])
{
 double D = sqr(m->g03*U3) - m->g00*m->g11*sqr(U1) - m->g00*m->g22*sqr(U2) - m->g00*m->g33*sqr(U3) - m->g00;
 double N = (-m->g03*U3 - sqrt(D))/m->g00;
 U[0] = 1./N;
 U[1] = U1/N;
 U[2] = U2/N;
 U[3] = U3/N;
}








/*
replaced by tetrad_general

DEVICEFUNC
void ortho_tetrad_U(
 double U[4],
 double g00, double g11, double g22, double g33, double g03,
 double e0[4], double e1[4], double e2[4], double e3[4])
{

 double k1 = sqrt(fabs(g33*sqr(U[3]) + g11*sqr(U[1]) + U[0]*(2.*g03*U[3]+g00*U[0])));
 double k2 = sqrt(fabs(g33*sqr(U[3]) + g11*sqr(U[1]) + g22*sqr(U[2]) + U[0]*(2.*g03*U[3]+g00*U[0])));
 double k3 = sqrt(fabs(g33*sqr(U[3]) + U[0]*(2.*g03*U[3]+g00*U[0])));
 double f;

 e0[0] = U[0];
 e0[1] = U[1];
 e0[2] = U[2];
 e0[3] = U[3];

 f = +1./(k1*k3);
 e1[0] = f * sqrt(g11)*U[1]*U[0];
 e1[1] = f * sqr(k3)/sqrt(g11);
 e1[2] = 0.0;
 e1[3] = f * sqrt(g11)*U[1]*U[3];

 f = -1./(k1*k2);
 e2[0] = f * sqrt(g22)*U[2]*U[0];
 e2[1] = f * sqrt(g22)*U[2]*U[1];
 e2[2] = f * sqr(k1)/sqrt(g22);
 e2[3] = f * sqrt(g22)*U[2]*U[3];

 f = +1./(k3*sqrt(fabs(g33*g00-g03*g03)));
 e3[0] = f * (+g33*U[3] + g03*U[0]);
 e3[1] = 0.0;
 e3[2] = 0.0;
 e3[3] = f * (-g03*U[3] - g00*U[0]);
}
*/


/*
commented out because orientation of vecotrs is unclear

DEVICEFUNC
void ortho_tetrad_U_phi_r_motion(
 double U[4],
 double g00, double g11, double g22, double g33, double g03,
 double e0[4], double e1[4], double e2[4], double e3[4])

{

 double k1 = sqrt(fabs(g33*sqr(U[3]) + g11*sqr(U[1]) + U[0]*(2.*g03*U[3]+g00*U[0])));
 double k3 = sqrt(fabs(g33*sqr(U[3])                 + U[0]*(2.*g03*U[3]+g00*U[0])));
 double f;

 e0[0] = U[0];
 e0[1] = U[1];
 e0[2] = 0.0;
 e0[3] = U[3];

 f = +1./(k1*k3);
 e1[0] = f * sqrt(g11)*U[1]*U[0];
 e1[1] = f * sqr(k3)/sqrt(g11);
 e1[2] = 0.0;
 e1[3] = f * sqrt(g11)*U[1]*U[3];

 e2[0] = 0.0;
 e2[1] = 0.0;
 e2[2] = -1./sqrt(g22);
 e2[3] = 0.0;

 f = +1./(k3*sqrt(fabs(g33*g00-g03*g03)));
 e3[0] = f * (+g33*U[3] + g03*U[0]);
 e3[1] = 0.0;
 e3[2] = 0.0;
 e3[3] = f * (-g03*U[3] - g00*U[0]);
}
*/







/*
#ifdef CUDA
__host__ __device__ void error(char *s) {
 
}
#endif
*/

#define theta_int(x) (g->mK*jacobi_icn((x)/sqrt(g->m2p),g->mm))
#define theta_inv(x) (sqrt(g->m2p)*jacobi_cn((x)/g->mK,g->mm))

DEVICEFUNC double geodesic_priv_RR(geodesic *g, double r);
DEVICEFUNC double geodesic_priv_TT(geodesic *g, double m);
DEVICEFUNC int geodesic_priv_R_roots(geodesic *g, double r0, int *error);
DEVICEFUNC int geodesic_priv_T_roots(geodesic *g, double m0, int *error);



DEVICEFUNC
int geodesic_init_inf(double i, double a, double alpha, double beta, geodesic *g, int *error)
{

 if ((a < 0.0) || (a > 1.-1e-6)) {
     if (error) *error = GD_ERROR_SPIN_RANGE;
     return FALSE;
 }

 if ((i <= 0.0) || (i>=PI_half)) {
     if (error) *error = GD_ERROR_INCL_RANGE;
     return FALSE;
 }

 if (beta == 0.0) beta = +1e-6;

 g->a = fmax(1e-4, a);
 g->incl  = i;
 g->cos_i = cos(i);
 g->alpha = alpha;
 g->beta  = beta;

 
 g->l = -alpha*sin(i);
 g->q = sqr(beta) + sqr(cos(i))*(sqr(alpha)-sqr(a));

 if (g->q == 0.0) {
     
     
     if (error) *error = GD_ERROR_Q_RANGE;
     return FALSE;
 }

 
 if (!geodesic_priv_R_roots(g, DBL_MAX, error)) return FALSE;
 if (!geodesic_priv_T_roots(g, g->cos_i, error)) return FALSE;

 
	g->Tpp = 2.*theta_int(0.0);

 
	g->Tip = theta_int(g->cos_i);

 if (error) *error = GD_OK;
 return TRUE;
}




DEVICEFUNC
int geodesic_init_src(double a, double r, double m, double k[4], int ppc, geodesic *g, int *error)
{
 
 double l,q;
 photon_motion_constants(a, r, m, k, &l, &q);

 g->a = fmax(1e-8, a);
 g->l = l;
 g->q = q;

 
 g->cos_i = g->alpha = g->beta = NAN;

 if (!geodesic_priv_R_roots(g, r, error)) return FALSE;
 if (!geodesic_priv_T_roots(g, m, error)) return FALSE;

 
 
 
 
 
 

 
 
 if (isnan(g->cos_i) && (r > g->rp)) {
     double T, Tmp, Tpp, sign_dm;
     Tmp = theta_int(m);                         
     Tpp = 2.*theta_int(0.0);                    
     T = geodesic_P_int(g, r, ppc);              
     sign_dm = (k[2]<0.0) ? +1.0 : -1.0;           
     
     T += (sign_dm>0.0) ? Tpp-Tmp : Tmp;
     
     while (T > Tpp) {
         T -= Tpp;
         sign_dm = -sign_dm;
     }
     
     g->cos_i = -sign_dm*theta_inv(T);
     g->incl  = acos(g->cos_i);
     g->alpha = -g->l/sqrt(1.0-sqr(g->cos_i));
     g->beta  = -sign_dm * sqrt(g->q - sqr(g->cos_i)*(sqr(g->alpha)-sqr(g->a)));
 }

 
	g->Tpp = 2.*theta_int(0.0);

 
	g->Tip = theta_int(g->cos_i);

 if (error) *error = GD_OK;
 return TRUE;
}




DEVICEFUNC
double geodesic_P_int(geodesic *g, double r, int ppc)
{
 double r1,r2,r3,r4,u,v, mm, R, A, B;

 #ifdef CUDA
 if (r  < g->rp) asm("exit;");
 #else
 if (r  < g->rp) error("(geodesic_P_int): r < periastron (%.3e < %.3e; nrr=%d)", r, g->rp, g->nrr);
 #endif
 if (r == g->rp) return g->Rpc;


 switch (g->type) {
     case GEOD_TYPE_RR:
         r1 = creal(g->r1);
         r2 = creal(g->r2);
         r3 = creal(g->r3);
         r4 = creal(g->r4);
         mm = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
         R  = 2./sqrt((r1-r3)*(r2-r4)) * jacobi_isn(sqrt(((r2-r4)*(r-r1))/((r1-r4)*(r-r2))), mm);
         return (ppc) ? g->Rpc + R : g->Rpc - R;

     case GEOD_TYPE_RR_DBL:
         #ifndef CUDA
         error("(geodesic_P_int): not implemented for GEOD_TYPE_RR_DBL");
         #endif
         return NAN;

     case GEOD_TYPE_RR_BH:
         r1 = creal(g->r1);
         r2 = creal(g->r2);
         r3 = creal(g->r3);
         r4 = creal(g->r4);
         mm = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
         R  = 2./sqrt((r1-r3)*(r2-r4)) * jacobi_isn(sqrt((r1-r3)/(r2-r3)*(r2-r)/(r1-r)), mm);
         return (ppc) ? g->Rpc + R : g->Rpc - R;

     case GEOD_TYPE_RC:
         r1 = creal(g->r1);
         r2 = creal(g->r2);
         u  = creal(g->r3);
         v  = cimag(g->r3);
         A = sqrt(sqr(r1-u)+sqr(v));
         B = sqrt(sqr(r2-u)+sqr(v));
         mm = (sqr(A+B) - sqr(r1-r2)) / (4.*A*B);
         R  = 1./sqrt(A*B) * jacobi_icn(((A-B)*r+r1*B-r2*A)/((A+B)*r-r1*B-r2*A), mm);
         
         
         return g->Rpc-R;

     case GEOD_TYPE_CC:
         r1 = creal(g->r1);  
         r2 = creal(g->r3);  
         r3 = cimag(g->r1);  
         r4 = cimag(g->r3);  
         A = sqrt(sqr(r1-r2) + sqr(r3+r4));
         B = sqrt(sqr(r1-r2) + sqr(r3-r4));
         double g1 = sqrt((4.*sqr(r3)-sqr(A-B))/(sqr(A+B)-4.*sqr(r3)));
         mm = 4.*A*B/sqr(A+B);
         R = 2./(A+B) * jacobi_itn((r-r1+r3*g1)/(r3+r1*g1-g1*r), mm);
         if ((R)<0) {
             fprintf(stderr, "CC: res=%e ppc=%d Rpc=%e R=%e\n", R, ppc, g->Rpc, R);
         }
         return g->Rpc-R;
 }

 
 return NAN;
}



DEVICEFUNC
void geodesic_position(geodesic *g, double P, double x[4])
{
 return;
}



DEVICEFUNC
double geodesic_position_rad(geodesic *g, double P)
{
 double r1,r2,r3,r4,u,v;

 if ((P<=0.0)||(P>=2.*g->Rpc)) {
     #ifndef CUDA
     error("(geodesic_position_rad) P out of range (P=%e, 2Rpa=%e)", P, 2*g->Rpc);
     #endif
     return NAN;
 }
 if (P == g->Rpc) return g->rp;

 switch (g->type) {
     case GEOD_TYPE_RR:
         r1 = creal(g->r1);
         r2 = creal(g->r2);
         r3 = creal(g->r3);
         r4 = creal(g->r4);
         double m4 = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
         double x4 = 0.5*fabs    (P - g->Rpc)*sqrt((r2-r4)*(r1-r3));
         double sn2 = pow( jacobi_sn(x4,m4), 2.0);
         return ( r1*(r2-r4)-r2*(r1-r4)*sn2 ) / ( r2-r4-(r1-r4)*sn2 );

     case GEOD_TYPE_RR_DBL:
         #ifndef CUDA
         error("(geodesic_position_rad): not implemented for GEOD_TYPE_RR_DBL");
         #endif
         return NAN;

     case GEOD_TYPE_RR_BH:
         #ifndef CUDA
         error("(geodesic_position_rad): not implemented for GEOD_TYPE_RR_BH");
         #endif
         return NAN;

     case GEOD_TYPE_RC:
         
         if (P > g->Rpc) return NAN;
         r1 = creal(g->r1);
         r2 = creal(g->r2);
         u  = creal(g->r3);
         v  = cimag(g->r3);
         double A = sqrt(sqr(r1-u)+sqr(v));
         double B = sqrt(sqr(r2-u)+sqr(v));
         double m2 = (sqr(A+B) - sqr(r1-r2)) / (4.*A*B);
         double cn = jacobi_cn(sqrt(A*B)*(g->Rpc-P), m2);
         
         return (r2*A - r1*B - (r2*A+r1*B)*cn ) / ( (A-B) - (A+B)*cn );

     case GEOD_TYPE_CC:
         #ifndef CUDA
         error("(geodesic_position_rad): not implemented for GEOD_TYPE_CC");
         #endif
         return NAN;
 }

 
 return NAN;
}




DEVICEFUNC
double geodesic_position_pol(geodesic *g, double P)
{
 double sign_dm, T;

 switch (g->type) {
     
     case GEOD_TYPE_RR:
     case GEOD_TYPE_RC:
     case GEOD_TYPE_CC:
         
         
         sign_dm = (g->beta>=0.0) ? +1.0 : -1.0;
         T = (sign_dm>0.0) ? -(g->Tpp-g->Tip) : -(g->Tip);
         
         while (P > T+g->Tpp) {
             T += g->Tpp;
             sign_dm = -sign_dm;
             
         }
         return -sign_dm*theta_inv(P-T);

     case GEOD_TYPE_RR_DBL:
         #ifndef CUDA
         error("(geodesic_dm_sign): not implemented for GEOD_TYPE_RR_DBL");
         #endif
         return NAN;

     case GEOD_TYPE_RR_BH:
         #ifndef CUDA
         error("(geodesic_dm_sign): not implemented for GEOD_TYPE_BH");
         #endif
         return NAN;
 }

 
 return NAN;
}




DEVICEFUNC
double geodesic_position_pol_sign_k_theta(geodesic *g, double P)
{
 double sign_dm, T;

 switch (g->type) {
     
     case GEOD_TYPE_RR:
     case GEOD_TYPE_RC:
     case GEOD_TYPE_CC:
         
         
         sign_dm = (g->beta>=0.0) ? +1.0 : -1.0;
         T = (sign_dm>0.0) ? -(g->Tpp-g->Tip) : -(g->Tip);
         
         while (P > T+g->Tpp) {
             T += g->Tpp;
             sign_dm = -sign_dm;
         }
         
         return (sign_dm<0) ? +1 : -1;

     case GEOD_TYPE_RR_DBL:
         #ifndef CUDA
         error("(geodesic_dm_sign): not implemented for GEOD_TYPE_RR_DBL");
         #endif
         return NAN;

     case GEOD_TYPE_RR_BH:
         #ifndef CUDA
         error("(geodesic_dm_sign): not implemented for GEOD_TYPE_BH");
         #endif
         return NAN;
 }

 
 return NAN;
}




DEVICEFUNC
double geodesic_position_azm(geodesic *g, double r, double m, double P)
{
 double phi = 0.0;

 int    ppc  = (g->nrr>0) && (P > g->Rpc);              
 double a2 = sqr(g->a);
 double rp   = 1. + sqrt(1.-a2);
 double rm   = 1. - sqrt(1.-a2);
 double r1, r2, r3, r4, A, B;

 switch (g->type) {
     
     case GEOD_TYPE_RR:
         r1 = creal(g->r1);
         r2 = creal(g->r2);
         r3 = creal(g->r3);
         r4 = creal(g->r4);
         A = integral_R_rp_re_inf(r1, r2, r3, r4, rp) + (ppc?+1:-1)*integral_R_rp_re(r1, r2, r3, r4, rp, r);
         B = integral_R_rp_re_inf(r1, r2, r3, r4, rm) + (ppc?+1:-1)*integral_R_rp_re(r1, r2, r3, r4, rm, r);
         phi += 1./sqrt(1.-a2) * ( A*(g->a*rp-g->l*a2/2.) - B*(g->a*rm-g->l*a2/2.) );
         break;

     case GEOD_TYPE_RR_DBL:
         #ifndef CUDA
         error("(geodesic_position_azm): not implemented for GEOD_TYPE_RR_DBL");
         #endif
         return NAN;

     case GEOD_TYPE_RR_BH:
         #ifndef CUDA
         error("(geodesic_position_azm): not implemented for GEOD_TYPE_RR_BH");
         #endif
         return NAN;

     case GEOD_TYPE_RC:
         r1 = creal(g->r1);
         r2 = creal(g->r2);
         A = integral_R_rp_cc2_inf(r1, r2, g->r3, rp, r);
         B = integral_R_rp_cc2_inf(r1, r2, g->r3, rm, r);
         phi += 1./sqrt(1.-a2) * ( A*(g->a*rp-g->l*a2/2.) - B*(g->a*rm-g->l*a2/2.) );
         
         break;

     case GEOD_TYPE_CC:
         #ifndef CUDA
         error("(geodesic_position_azm): not implemented for GEOD_TYPE_CC");
         #endif
         return NAN;
 }

 
 
 double phi_pp = 2.0*g->l/g->a*integral_T_mp(g->m2m, g->m2p, 1.0, 0.0);
 double phi_ip =     g->l/g->a*integral_T_mp(g->m2m, g->m2p, 1.0, g->cos_i);
 double phi_mp =     g->l/g->a*integral_T_mp(g->m2m, g->m2p, 1.0, m);

 double T;
 double sign_dm = (g->beta>=0.0) ? +1.0 : -1.0;
 if (sign_dm > 0.0) {
     T = -(g->Tpp-g->Tip);
     phi -= phi_pp-phi_ip;
 } else {
     T = -g->Tip;
     phi -= phi_ip;
 }

 while (P >= T+g->Tpp) {
     T += g->Tpp;
     phi += phi_pp;
     sign_dm = -sign_dm;
     break;
 }

 phi += (sign_dm<0) ? phi_mp : phi_pp-phi_mp;
 

 
 
 return phi;
}



DEVICEFUNC
double geodesic_timedelay(geodesic *g, double P1, double r1, double m1, double P2, double r2, double m2)
{
 double time = 0.0;

 
 if (P1 > P2) {
     double tmp;
     tmp=P2; P2=P1; P1=tmp;
     tmp=r2; r2=r1; r1=tmp;
     tmp=m2; m2=m1; m1=tmp;
 }

 if (r1 == 0) {
     r1 = geodesic_position_rad(g, P1);
     m1 = geodesic_position_pol(g, P1);
 }

 if (r2 == 0) {
     r2 = geodesic_position_rad(g, P2);
     m2 = geodesic_position_pol(g, P2);
 }


 double a2 = sqr(g->a);
 double rp   = 1. + sqrt(1.-a2);
 double rm   = 1. - sqrt(1.-a2);
 double ra = creal(g->r1);
 double rb = creal(g->r2);
 double rc = creal(g->r3);
 double rd = creal(g->r4);
 double R0, R1, R2, RA, RB, A, B, s;

 
 
 

 if (r1 < g->rp) error("geodesic_timedelay: r1 < r_p (%e/%e)", r1, g->rp);
 if (r2 < g->rp) error("geodesic_timedelay: r2 < r_p (%e/%e)", r2, g->rp);

 switch (g->type) {
     
     case GEOD_TYPE_RR:
         s = (((P1 > g->Rpc)&&(P2 < g->Rpc)) || ((P1 < g->Rpc)&&(P2 > g->Rpc))) ? +1 : -1;
         R0 = integral_R_r0_re(ra, rb, rc, rd, r1)     + s*integral_R_r0_re(ra, rb, rc, rd, r2);
         R1 = integral_R_r1_re(ra, rb, rc, rd, r1)     + s*integral_R_r1_re(ra, rb, rc, rd, r2);
         R2 = integral_R_r2_re(ra, rb, rc, rd, r1)     + s*integral_R_r2_re(ra, rb, rc, rd, r2);
         RA = integral_R_rp_re(ra, rb, rc, rd, rp, r1) + s*integral_R_rp_re(ra, rb, rc, rd, rp, r2);
         RB = integral_R_rp_re(ra, rb, rc, rd, rm, r1) + s*integral_R_rp_re(ra, rb, rc, rd, rm, r2);
         A = (-g->a*g->l+4.)*rp - 2.*a2;
         B = (+g->a*g->l-4.)*rm + 2.*a2;
         time += 4.*fabs(R0) + 2.*fabs(R1) + fabs(R2) + (A*fabs(RA) + B*fabs(RB))/sqrt(1.-a2);
         
         
         
         break;

     case GEOD_TYPE_RR_DBL:
         #ifndef CUDA
         error("(geodesic_timedelay): not implemented for GEOD_TYPE_RR_DBL");
         #endif
         return NAN;

     case GEOD_TYPE_RR_BH:
         #ifndef CUDA
         error("(geodesic_timedelay): not implemented for GEOD_TYPE_RR_BH");
         #endif
         return NAN;

     case GEOD_TYPE_RC:
         R0 = integral_R_r0_cc(ra, rb, g->r3, r1) -integral_R_r0_cc(ra, rb, g->r3, r2);
         R1 = (r1<r2) ? integral_R_r1_cc(ra, rb, g->r3, r1, r2) : integral_R_r1_cc(ra, rb, g->r3, r2, r1);
         R2 = (r1<r2) ? integral_R_r2_cc(ra, rb, g->r3, r1, r2) : integral_R_r2_cc(ra, rb, g->r3, r2, r1);
         RA = (r1<r2) ? integral_R_rp_cc2(ra, rb, g->r3, rp, r1, r2) : integral_R_rp_cc2(ra, rb, g->r3, rp, r2, r1);
         RB = (r1<r2) ? integral_R_rp_cc2(ra, rb, g->r3, rm, r1, r2) : integral_R_rp_cc2(ra, rb, g->r3, rm, r2, r1);
         A = (-g->a*g->l+4.)*rp - 2.*a2;
         B = (+g->a*g->l-4.)*rm + 2.*a2;
         time += 4.*fabs(R0) + 2.*fabs(R1) + fabs(R2) + (A*fabs(RA) + B*fabs(RB))/sqrt(1.-a2);
         
         
         
         break;

     case GEOD_TYPE_CC:
         #ifndef CUDA
         error("(geodesic_timedelay): not implemented for GEOD_TYPE_CC");
         #endif
         return NAN;
 }

 /*

 
 
 double time_pp = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, 0.0)*2;        
 double time_m1 = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, fabs(m1));     
 double time_m2 = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, fabs(m2));     

 double dP = P2-P1;

 while (T >= g->Tpp) {
     T    -= g->Tpp;
     time += time_pp;
 }

 if (T > g->Tpp/2) {
     time += fabs(time_m2 - time_m1);
 } else {
     time += fabs(time_m2 - time_m1);
 }


 double sign_dm = (g->beta>=0.0) ? +1.0 : -1.0;
 if (sign_dm > 0.0) {
     T = -(g->Tpp-g->Tip);
     time -= time_pp-time_ip;
 } else {
     T = -g->Tip;
     time -= time_ip;
 }

 while (fabs(P2-P1) >= T+g->Tpp) {
     T += g->Tpp;
     time += time_pp;
     sign_dm = -sign_dm;
     break;
 }

 time += (sign_dm<0) ? time_mp : time_pp-time_mp;
*/
 return time;

/*
old:
 double Tmp = g->mK*elliptic_k(g->mm);
 double Tmo = g->mK*jacobi_icn(g->cos_i/sqrt(g->m2p),g->mm);
 double Pmp = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, 0.0);
 double Pmo = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, g->cos_i);
 double Pme = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, fabs(mu_e));
 


 x += (g->beta>=0.0) ? Tmp-Tmo : Tmo;                             
 T -= (g->beta>=0.0) ? Pmp-Pmo : Pmo;
 int k = (g->beta>=0.0) ? 3 : 0;
 while (x >= Tmp) {
     x -= Tmp;
     T += Pmp;
     k++;
 }
 switch(k%4) {
     case 0: T += Pme; break;
     case 1: T += Pmp-Pme; break;
     case 2: T += Pme; break;
     case 3: T += Pmp-Pme; break;
 }
 */
 return time;
}




DEVICEFUNC
double geodesic_dm_sign(geodesic *g, double P)
{
 double sign_dm, T;

 switch (g->type) {
     
     case GEOD_TYPE_RR:
     case GEOD_TYPE_RC:
     case GEOD_TYPE_CC:
         
         
         sign_dm = (g->beta>=0.0) ? +1.0 : -1.0;
         T = (sign_dm>0.0) ? -(g->Tpp-g->Tip) : -(g->Tip);

         
         while (P > T+g->Tpp) {
             T += g->Tpp;
             sign_dm = -sign_dm;
         }
         return sign_dm;

     case GEOD_TYPE_RR_DBL:
         #ifndef CUDA
         error("(geodesic_dm_sign): not implemented for GEOD_TYPE_RR_DBL");
         #endif
         return NAN;

     case GEOD_TYPE_RR_BH:
         #ifndef CUDA
         error("(geodesic_dm_sign): not implemented for GEOD_TYPE_BH");
         #endif
         return NAN;
 }

 
 return NAN;
}




DEVICEFUNC
void geodesic_momentum(geodesic *g, double P, double r, double m, double k[])
{
 double dm;

 
 if ((r==0.0) && (m==0.0)) {
     r = geodesic_position_rad(g, P);
     m = geodesic_position_pol(g, P);
 }

 switch (g->type) {
     
     case GEOD_TYPE_RR:
     case GEOD_TYPE_RC:
     case GEOD_TYPE_CC:
         dm = geodesic_dm_sign(g, P);
         photon_momentum(g->a, r, m, g->l, g->q, (P<g->Rpc?-1:+1), dm, k);
         return;

     case GEOD_TYPE_RR_DBL:
         #ifndef CUDA
         error("(geodesic_dm_sign): not implemented for GEOD_TYPE_RR_DBL");
         #endif
         k[0]=k[1]=k[2]=k[3]=NAN;
         return;

     case GEOD_TYPE_RR_BH:
         #ifndef CUDA
         error("(geodesic_dm_sign): not implemented for GEOD_TYPE_BH");
         #endif
         k[0]=k[1]=k[2]=k[3]=NAN;
         return;
 }

 return;
}




DEVICEFUNC
double geodesic_find_midplane_crossing(geodesic *g, int order)
{
 if (g->q<=0.0) {
     
     return NAN;
 }

 double u = g->cos_i/sqrt(g->m2p);
 if (!ensure_range(&u, -1.0, +1.0, 1e-4)) {
     #ifndef CUDA
     error("(geodesic_find_midplane_crossing): u out of range (%e)", u);
     #endif
     return NAN;
 }

 double pos;
 if (g->beta > 0.0)
     pos = g->mK*( (2.*(double)order+1.)*elliptic_k(g->mm) + jacobi_icn(u,g->mm) );
 else if (g->beta < 0.0)
     pos = g->mK*( (2.*(double)order+1.)*elliptic_k(g->mm) - jacobi_icn(u,g->mm) );
 else
     pos = g->mK*( (2.*(double)order+1.)*elliptic_k(g->mm) );

 if (pos > 2.*g->Rpc) pos = NAN;

 return pos;
}




DEVICEFUNC
void geodesic_follow(geodesic *g, double step, double *P, double *r, double *m, int *status)
{
 const double MAXSTEP_FACTOR = 5e-2;

 do {
     double truestep = step/fabs(step) * fmin(fabs(step), MAXSTEP_FACTOR*sqrt(*r));
     (*P) = (*P) + truestep/(sqr(*r)+sqr((g->a)*(*m)));   
     (*r) = geodesic_position_rad(g, *P);
     (*m) = geodesic_position_pol(g, *P);
     if ((*r) < 1.01*r_bh(g->a)) {
         if (status) *status = 0;
         return;
     }
     if ((*P < 0.0) || (*P>2.*g->Rpc)) {
         if (status) *status = 0;
         return;
     }
     step -= truestep;
 } while (fabs(step) > 1e-5);

 if (status) *status = 1;
}





































DEVICEFUNC
double geodesic_priv_RR(geodesic *g, double r)
{
 return sqr4(r) + (sqr(g->a) - sqr(g->l) - g->q)*sqr(r) + 2.*(g->q + sqr(g->l-g->a))*r - sqr(g->a)*g->q;
}


DEVICEFUNC
double geodesic_priv_TT(geodesic *g, double m)
{
 return sqr(g->a)*(g->m2m + m*m)*(g->m2p - m*m);
}


DEVICEFUNC
int geodesic_priv_R_roots(geodesic *g, double r0, int *error)
{
 double a  = g->a;
 double l  = g->l;
 double q  = g->q;
 double a2 = sqr(a);
 double l2 = sqr(l);
 double A,B,C,D,E,F,X,Z,z;

 
 C = sqr(a-l)+q;
 D = 2./3.*(q+l2-a2);
 E = 9./4.*sqr(D)-12.*a2*q;
 F = -27./4.*sqr3(D)-108.*a2*q*D+108.*sqr(C);
 X = sqr(F)-4.*sqr3(E);
 if (X >= 0) {
     A = (F>sqrt(X)?+1:-1)*1./3.*pow(fabs(F-sqrt(X))/2.,1./3.)+ (F>-sqrt(X)?+1:-1)*1./3.*pow(fabs(F+sqrt(X))/2.,1./3.);
 } else {
      Z = sqrt(pow(F/54.,2) + pow(sqrt(-X)/54.,2));
      z = atan2(sqrt(-X)/54.,F/54.);
      A = pow(Z,1./3.)*2.*cos(z/3.);
 }
 B = sqrt(A+D);
 g->r1 = +B/2. + .5*csqrt(makeComplex(-A+2.*D-4.*C/B,0.0));
 g->r2 = +B/2. - .5*csqrt(makeComplex(-A+2.*D-4.*C/B,0.0));
 g->r3 = -B/2. + .5*csqrt(makeComplex(-A+2.*D+4.*C/B,0.0));
 g->r4 = -B/2. - .5*csqrt(makeComplex(-A+2.*D+4.*C/B,0.0));
 sort_roots(&g->nrr, &g->r1, &g->r2, &g->r3, &g->r4);

 
 switch (g->nrr) {
     case 4:
         g->type = GEOD_TYPE_RR;
         
         
         if ((r0<creal(g->r3)) || ((r0>creal(g->r2)) && (r0<creal(g->r1)))) {
             if (error) *error = GD_ERROR_UNKNOWN_SOLUTION;
             return FALSE;
         }
         
         if (fabs(creal(g->r1)-creal(g->r2)) < 1e-8) {
             g->type = GEOD_TYPE_RR_DBL;
             if (error) *error = GD_ERROR_TYPE_RR_DOUBLE;
             return FALSE;
         }
         
         if ((r0>=creal(g->r3)) && (r0<=creal(g->r2))) {
             g->type = GEOD_TYPE_RR_BH;
         }
         break;
     case 2:
         g->type = GEOD_TYPE_RC;
         break;
     case 0:
         g->type = GEOD_TYPE_CC;
         break;
     default:
         if (error) *error = GD_ERROR_UNKNOWN_SOLUTION;
         return FALSE;
 }

 
 double r1,r2,r3,r4,u,v,mm;
 switch (g->type) {
     case GEOD_TYPE_RR:
         r1 = creal(g->r1);
         r2 = creal(g->r2);
         r3 = creal(g->r3);
         r4 = creal(g->r4);
         mm = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
         g->rp  = r1;
         g->Rpc = 2./sqrt((r1-r3)*(r2-r4)) * jacobi_isn(sqrt((r2-r4)/(r1-r4)), mm);
         break;

     case GEOD_TYPE_RR_BH:
         r1 = creal(g->r1);
         r2 = creal(g->r2);
         r3 = creal(g->r3);
         r4 = creal(g->r4);
         mm = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
         g->rp  = r2;
         g->Rpc = 2./sqrt((r1-r3)*(r2-r4)) * elliptic_k(mm); 
         break;

     case GEOD_TYPE_RC:
         r1 = creal(g->r1);
         r2 = creal(g->r2);
         u  = creal(g->r3);
         v  = cimag(g->r3);
         A = sqrt(sqr(r1-u)+sqr(v));
         B = sqrt(sqr(r2-u)+sqr(v));
         mm = (sqr(A+B) - sqr(r1-r2)) / (4.*A*B);
         g->rp  = r1;
         g->Rpc = 1./sqrt(A*B) * jacobi_icn((A-B)/(A+B), mm);
         break;

     case GEOD_TYPE_CC:
         r1 = creal(g->r1);   
         r2 = creal(g->r3);   
         r3 = cimag(g->r1);   
         r4 = cimag(g->r3);   
         A = sqrt(sqr(r1-r2) + sqr(r3+r4));
         B = sqrt(sqr(r1-r2) + sqr(r3-r4));
         double g1 = sqrt((4.*sqr(r3)-sqr(A-B))/(sqr(A+B)-4.*sqr(r3)));
         mm = 4.*A*B/sqr(A+B);
         g->rp  = r1 - r3*g1;
         if (g->rp > r_bh(g->a)) warning("geodesic_priv_R_roots(): g->rp > r_bh(g->a)");
         g->Rpc = 2./(A+B) * jacobi_itn(-1./g1, mm);
         break;

     default:
         return FALSE;
 }


 return TRUE;
}




DEVICEFUNC
int geodesic_priv_T_roots(geodesic *g, double m, int *error)
{
 double a  = g->a;
 double l  = g->l;
 double q  = g->q;
 double a2 = sqr(a);
 double l2 = sqr(l);

 
 
 
 
 
 #ifndef CUDA
 long double qla = q + l2 - a2;
 long double X = sqrt(sqr(qla)+4.*q*a2) + qla;
 long double dbla = a2+a2;
 long double dblq = q+q;
 g->m2m = X/dbla;
 g->m2p = dblq/X;
 #else
 double qla = q + l2 - a2;
 double X = sqrt(sqr(qla)+4.*q*a2) + qla;
 double dbla = a2+a2;
 double dblq = q+q;
 g->m2m = X/dbla;
 g->m2p = dblq/X;
 #endif

 if ((g->m2p<=0.0) || (g->m2p >= 1.0)) {
     if (error) *error = GD_ERROR_MUPLUS_RANGE;
     return FALSE;
 }

 if (q > 0.0) {
     g->mm = g->m2p/(g->m2p+g->m2m);

     if ((g->mm<0.0) || (g->mm>=1.0)) {
         if (error) *error = GD_ERROR_MM_RANGE;
         return FALSE;
     }

     if (fabs(m) > sqrt(g->m2p)) {
         if (error) *error = GD_ERROR_MU0_RANGE;
         return FALSE;
     }

     g->mK = 1./sqrt(a2*(g->m2p+g->m2m));
 } else
 if (q < 0.0) {
     g->mm = (g->m2p+g->m2m)/g->m2p;

     if ((g->mm<0.0) || (g->mm>=1.0)) {
         if (error) *error = GD_ERROR_MM_RANGE;
         return FALSE;
     }

     if ((fabs(m) > sqrt(g->m2p)) || (fabs(m) < sqrt(-g->m2m))) {
         if (error) *error = GD_ERROR_MU0_RANGE;
         return FALSE;
     }

     g->mK = 1./sqrt(a2*g->m2p);
 } else {
     
     
     
     if (error) *error = GD_ERROR_Q_RANGE;
     return FALSE;
 }

 return TRUE;
}



#undef theta_int
#undef theta_inv





#ifndef CUDA


static float disk_nt_bh_mass = 10.0;
static float disk_nt_bh_spin = 0.0;
static float disk_nt_disk_mdot  = 0.1;
static float disk_nt_disk_rms   = 6.0;
static float disk_nt_disk_alpha = 0.1;
static int   disk_nt_options = 0;



DEVICEFUNC
int disk_nt_setup(double M, double a, double mdot_or_L, double alpha, int options)
{
 disk_nt_bh_mass    = M;
 disk_nt_bh_spin    = a;
 disk_nt_disk_rms   = disk_nt_r_min();
 disk_nt_disk_alpha = alpha;
 disk_nt_options    = options;
 if (disk_nt_options & DISK_NT_OPTION_LUMINOSITY) {
     DEVICEFUNC double disk_nt_find_mdot_for_luminosity(double L0);
     disk_nt_disk_mdot = disk_nt_find_mdot_for_luminosity(mdot_or_L);
     
     
     
     
 }
 else {
     disk_nt_disk_mdot = mdot_or_L;
     
     
     
     
 }
 return 0;
}



DEVICEFUNC
void disk_nt_done()
{
}



DEVICEFUNC
double disk_nt_r_min()
{
 double a = disk_nt_bh_spin;
 double z1,z2,r0;
 double sga = (a>=0.0) ? +1. : -1.;
 z1 = 1.+pow(1.-a*a, 1./3.)*(pow(1.+a, 1./3.)+pow(1.-a, 1./3.));
 z2 = sqrt(3.*a*a+z1*z1);
 r0 = 3.+z2-sga*sqrt((3.-z1)*(3.+z1+2.*z2));
 return r0+1e-3;
}



DEVICEFUNC
double disk_nt_flux(double r)
{
 if (r <= disk_nt_disk_rms) return 0.0;
 double a = disk_nt_bh_spin;
 double x=sqrt(r);
 double x0,x1,x2,x3;
 x0=sqrt(disk_nt_disk_rms);
 x1=+2.*cos(1./3.*acos(a)-M_PI/3.);
 x2=+2.*cos(1./3.*acos(a)+M_PI/3.);
 x3=-2.*cos(1./3.*acos(a));
 double f0,f1,f2,f3,F;
 
 f0=x-x0-1.5*a*log(x/x0);
 f1=3.*sqr(x1-a)/(x1*(x1-x2)*(x1-x3))*log((x-x1)/(x0-x1));
 f2=3.*sqr(x2-a)/(x2*(x2-x1)*(x2-x3))*log((x-x2)/(x0-x2));
 f3=3.*sqr(x3-a)/(x3*(x3-x1)*(x3-x2))*log((x-x3)/(x0-x3));
 F = 1./(4.*M_PI*r) * 1.5/(x*x*(x*x*x-3.*x+2.*a)) * (f0-f1-f2-f3);

 
 
 
 

 
 

 return 9.1721376255e+28 * F * disk_nt_disk_mdot/disk_nt_bh_mass; 
}



DEVICEFUNC
double disk_nt_lumi()
{
 const float disk_rmax = 1e5;

 
 

 double func_luminosity(double log_r)
 {
     double r = exp(log_r);
     
     double gtt = -1. + 2./r;
     double gtf = -2.*disk_nt_bh_spin/r;
     double gff = sqr(r) + sqr(disk_nt_bh_spin) + 2.*sqr(disk_nt_bh_spin)/r;
     double Omega = 1./(disk_nt_bh_spin + pow(r,1.5));
     double U_t = sqrt(-1.0/(gtt + 2.*Omega*gtf + sqr(Omega)*gff)) * (gtt + Omega*gtf);
     double F = disk_nt_flux(r);
     
     return 2.*M_PI*r*2.0*(-U_t)*F * r;
 }

 double L = integrate_simpson(func_luminosity, log(disk_nt_disk_rms), log(disk_rmax), 1e-5);

 
 L *= sqr(disk_nt_bh_mass*grav_radius);

 return L/(L_Edd*disk_nt_bh_mass);
}



DEVICEFUNC
double disk_nt_mdot()
{
 return disk_nt_disk_mdot;
}



DEVICEFUNC
double disk_nt_sigma(double r)
{
 if (r < disk_nt_disk_rms) return 0.0;
 double a = disk_nt_bh_spin;

 double x=sqrt(r);
 double x0,x1,x2,x3;
 x0=sqrt(disk_nt_disk_rms);
 x1=+2.*cos(1./3.*acos(a)-M_PI/3.);
 x2=+2.*cos(1./3.*acos(a)+M_PI/3.);
 x3=-2.*cos(1./3.*acos(a));

 double xA, xB, xC, xD, xE, xL;
 xA = 1. + sqr(a)/sqr(r) + 2.*sqr(a)/sqr3(r);
 xB = 1. + a/(x*x*x);
 xC = 1. - 3./(x*x) + 2.*a/(x*x*x);
 xD = 1. - 2./r + sqr(a)/sqr(r);
 xE = 1. + 4.*sqr(a)/sqr(r) - 4.*sqr(a)/sqr3(r) + 3.*sqr4(a)/sqr4(r);

 double f0, f1, f2, f3;
 f0=x-x0-1.5*a*log(x/x0);
 f1=3.*(x1-a)*(x1-a)/(x1*(x1-x2)*(x1-x3))*log((x-x1)/(x0-x1));
 f2=3.*(x2-a)*(x2-a)/(x2*(x2-x1)*(x2-x3))*log((x-x2)/(x0-x2));
 f3=3.*(x3-a)*(x3-a)/(x3*(x3-x2)*(x3-x1))*log((x-x3)/(x0-x3));
 xL = (1.+a/(x*x*x))/sqrt(1.-3./(x*x)+2.*a/(x*x*x))/x * (f0-f1-f2-f3);

 double xMdot = disk_nt_disk_mdot*disk_nt_bh_mass*Mdot_Edd/1e17;
 double r_im = 40.*(pow(disk_nt_disk_alpha,2./21.)/pow(disk_nt_bh_mass/3.,2./3.)*pow(xMdot,16./20.)) * pow(xA,20./21.) *
 pow(xB,-36./21.) * pow(xD,-8./21.) * pow(xE,-10./21.) * pow(xL,16./21.);

 double Sigma;
 if (r < r_im)
     Sigma = 20. * (disk_nt_bh_mass/3.)/xMdot/disk_nt_disk_alpha * sqrt(r*r*r) * 1./(xA*xA) * pow(xB,3.) * sqrt(xC) * xE * 1./xL;
 else {
     Sigma = 5e4 * pow(disk_nt_bh_mass/3.,-2./5.)*pow(xMdot,3./5.)*pow(disk_nt_disk_alpha,-4./5.) * pow(r,-3./5.) * pow(xB,-4./5.) * sqrt(xC) * pow(xD,-4./5.) * pow(xL,3./5.);
 }

 return Sigma;
}



DEVICEFUNC
double disk_nt_ell(double r)
{
 double a = disk_nt_bh_spin;
 r = fmax(disk_nt_disk_rms, r);
 return (r*r-2.*a*sqrt(r)+a*a) / (sqrt(r)*r-2.*sqrt(r)+a);
}



DEVICEFUNC
double disk_nt_vr(double r)
{
 return 0.0;
}



DEVICEFUNC
double disk_nt_h(double r)
{
 return 0.0;
}



DEVICEFUNC
double disk_nt_dhdr(double r)
{
 return 0.0;
}



DEVICEFUNC
void disk_nt_dump(char* filename)
{
 FILE* stream = stdout;
 if (filename) stream = fopen(filename, "w");
 if (!stream) {
     fprintf(stderr, "disk_nt_dump: cannot open output (%s)\n", filename);
     return;
 }


 const float disk_rmax = 2000.;
 fprintf(stream, "# (sim5disk-nt) dump\n");
 fprintf(stream, "#-------------------------------------------\n");
 fprintf(stream, "# M        = %.4f\n", disk_nt_bh_mass);
 fprintf(stream, "# a        = %.4f\n", disk_nt_bh_spin);
 fprintf(stream, "# rmin     = %.4f\n", disk_nt_disk_rms);
 fprintf(stream, "# rmax     = %.4f\n", disk_rmax);
 fprintf(stream, "# alpha    = %.4f\n", disk_nt_disk_alpha);
 fprintf(stream, "# options  = %d\n",   disk_nt_options);
 fprintf(stream, "# L        = %e\n", disk_nt_lumi());
 fprintf(stream, "# mdot     = %e\n", disk_nt_mdot());
 fprintf(stream, "#-------------------------------------------\n");
 fprintf(stream, "# r   flux   sigma   ell   vr   H   dH/dr\n");
 fprintf(stream, "#-------------------------------------------\n");

 double r;
 for (r=disk_nt_disk_rms; r<disk_rmax; r*=1.05) {
     fprintf(stream,
         "%e  %e  %e  %e  %e  %e  %e\n",
         r,
         disk_nt_flux(r),
         disk_nt_sigma(r),
         disk_nt_ell(r),
         disk_nt_vr(r),
         disk_nt_h(r),
         disk_nt_dhdr(r)
     );
 }

 fflush(stream);
 if (filename) fclose(stream);
}


DEVICEFUNC
double disk_nt_find_mdot_for_luminosity(double L0)
{
 double L;

 double fce(double xmdot) {
     disk_nt_disk_mdot = xmdot;
     return L0-disk_nt_lumi();
 }

 int res = rtbis(0.0, 100.0, 1e-6, fce, &L);
 return (res) ? L : 0.0;
}

#endif


DEVICEFUNC
void polarization_vector(double k[4], sim5complex wp, sim5metric *metric, double f[4])
{
 /*
 
 double a = metric->a;
 double m = metric->m;
 double r = metric->r;
 double r2 = r*r;
 double a2 = a*a;
 double s2 = 1.0-m*m;
 double s = sqrt(s2);

 
 double d1 = r*k[0] - r*a*s2*k[3];
 double d2 = a2*s*m*k[0] - a*s*m*(r2+a2)*k[3];
 double d3 = a*r*s2*k[1] + a*s*m*(r2+a2)*k[2];
 double g1 = a*m*k[0] - a2*m*s2*k[3];
 double g2 = -a*r*s*k[0] + r*(r2+a2)*s*k[3];
 double g3 = a2*m*s2*k[1] - r*(r2+a2)*s*k[2];
 double K1 = creal(wp);
 double K2 = cimag(wp);

 
 f[0] = 0.0;
 f[1] =  (g2*K1-d2*K2)*(metric->g33*k[3]+metric->g03*k[0]) - (g3*K1-d3*K2)*metric->g22*k[2];
 f[2] = -(g1*K1-d1*K2)*(metric->g33*k[3]+metric->g03*k[0]) - (g3*K1-d3*K2)*metric->g11*k[1];
 f[3] =  (g1*K1-d1*K2)*(metric->g22*k[2])                  - (g2*K1-d2*K2)*metric->g11*k[1];

 vector_norm_to(f, 1.0, metric);
 */

 
 
 double a = metric->a;
 double m = metric->m;
 double r = metric->r;
 double s = sqrt(1.0-m*m);
 double ra2 = r*r + a*a;
 double r2 = r*r;
 double a2 = a*a;
 double s2 = 1.0-m*m;

 
 if (s < 1e-12) {
     s  = 1e-12;
     s2 = 1e-24;
     m  = 1.0-0.5*s;
 }

 
 
 double A1 = (+r*creal(wp) - a*m*cimag(wp))/(r*r + a*a*m*m);
 double A2 = (-r*cimag(wp) - a*m*creal(wp))/(r*r + a*a*m*m);

 
 
 
 
 
 
 
 
 
 
 f[0] = 0.0;
 f[3] = (
           + metric->g11*A1*k[1]*(s*r2*k[3] + s*a2*k[3] - s*a*k[0])
           + metric->g22*A2*k[2]*(k[0] - a*s2*k[3])
        ) / (
          + sqr(k[0])*metric->g33*(s*k[3]*a)
          + sqr(k[0])*metric->g03*(s*k[0]*a - s*r2*k[3] - s*a2*k[3] - a2*s*s2*k[3])
          + sqr(k[1])*metric->g11*a*s*s2*(+ r2*k[3] + a2*k[3] - a*k[0])
          + sqr(k[2])*metric->g22*(a2*a*s*s2*k[3] + r2*a*s*s2*k[3] - s*r2*k[0] - s*a2*k[0])
          + sqr(k[3])*metric->g33*s*(k[3]*a*s2*r2 + k[3]*a2*a*s2 - k[0]*r2 - k[0]*a2 - a2*s2*k[0])
          + sqr(k[3])*metric->g03*a*s*s2*(r2*k[0] + a2*k[0])
        );
 f[1] = (A1-a*s*s*k[1]*f[3])/(k[0]-a*s*s*k[3]);
 f[2] = (A2 + s*k[2]*f[3]*ra2)/(s*k[3]*ra2 - s*a*k[0]);

 vector_norm_to(f, 1.0, metric);
}



/*
DEVICEFUNC
sim5complex photon_wp_const(double k[4], double f[4], sim5metric *metric)
{
 double a = metric->a;
 double r = metric->r;
 double m = metric->m;
 double s = sqrt(1.0-m*m);

 
 
 double A1 = (k[0]*f[1]-k[1]*f[0]) + a*s*s*(k[1]*f[3]-k[3]*f[1]);
 double A2 = s*( (r*r+a*a)*(k[3]*f[2]-k[2]*f[3]) - a*(k[0]*f[2]-k[2]*f[0]) );
 double wp1 = +r*A1 - a*m*A2;
 double wp2 = -r*A2 - a*m*A1;

 return wp1 + ComplexI*wp2;
}
*/


DEVICEFUNC
sim5complex polarization_constant(double k[4], double f[4], sim5metric *metric)
{
 double a = metric->a;
 double m = metric->m;
 double r = metric->r;

 
 double A1 = (k[0]*f[1]-k[1]*f[0]) + a*(1.-m*m)*(k[1]*f[3]-k[3]*f[1]);
 double A2 = sqrt(1.-m*m)*( (r*r+a*a)*(k[3]*f[2]-k[2]*f[3]) - a*(k[0]*f[2]-k[2]*f[0]) );
 double wp1 = +r*A1 - a*m*A2;
 double wp2 = -r*A2 - a*m*A1;

 return wp1 + ComplexI*wp2;

/*
 double a = metric->a;
 double m = metric->m;
 double r = metric->r;
 double s = sqrt(1.0-m*m);
 double r2 = r*r;
 double a2 = a*a;
 double s2 = 1.0-m*m;

 
 double d1 = r*k[0] - r*a*s2*k[3];
 double d2 = a2*s*m*k[0] - a*s*m*(r2+a2)*k[3];
 double d3 = a*r*s2*k[1] + a*s*m*(r2+a2)*k[2];
 double g1 = a*m*k[0] - a2*m*s2*k[3];
 double g2 = -a*r*s*k[0] + r*(r2+a2)*s*k[3] ;
 double g3 = a2*m*s2*k[1] - r*(r2+a2)*s*k[2];

 
 double K1 = d1*f[1] + d2*f[2] + d3*f[3];
 double K2 = g1*f[1] + g2*f[2] + g3*f[3];

 return K1 + ComplexI*K2;
*/


 /*
 
 
 
 
 
 
 double a = metric->a;
 double r = metric->r;
 double m = metric->m;
 double s = sqrt(1.0-m*m);

 
 double A1 = (k[0]*f[1]-k[1]*f[0]) + a*s*s*(k[1]*f[3]-k[3]*f[1]);
 double A2 = s*( (r*r+a*a)*(k[3]*f[2]-k[2]*f[3]) - a*(k[0]*f[2]-k[2]*f[0]) );
 double wp1 = +r*A1 - a*m*A2;
 double wp2 = -r*A2 - a*m*A1;

 return wp1 + ComplexI*wp2;
 */
}



/*
DEVICEFUNC
sim5complex photon_wp_const_infinity(double a, double inc, double alpha, double beta)
{
 
 
 

 
 double wp1 = -(-alpha - a*sin(inc));
 double wp2 = -beta;
 return wp1 + ComplexI*wp2;
}
*/

DEVICEFUNC
sim5complex polarization_constant_infinity(double a, double alpha, double beta, double incl)
{
 double gamma = - alpha - a*sin(incl);

 double K1 = -gamma;
 double K2 = -beta;

 return K1 + ComplexI*K2;
}


DEVICEFUNC
double polarization_angle_rotation(double a, double inc, double alpha, double beta, sim5complex kappa)
{
 double kappa1 = creal(kappa);
 double kappa2 = cimag(kappa);
 double S = -alpha - a*sin(inc);
 double T = +beta;
 double X = (-S*kappa2 - T*kappa1)/(S*S+T*T);
 double Y = (-S*kappa1 + T*kappa2)/(S*S+T*T);
 return atan2(Y,X);
}




/*
DEVICEFUNC
void kappa_pw(double a, double r, double m, double k[4], double f[4], double *kappa1, double *kappa2)
{
 
 double A1 = (k[0]*f[1]-k[1]*f[0]) + a*(1.-m*m)*(k[1]*f[3]-k[3]*f[1]);
 double A2 = sqrt(1.-m*m)*( (r*r+a*a)*(k[3]*f[2]-k[2]*f[3]) - a*(k[0]*f[2]-k[2]*f[0]) );
 *kappa1 = +r*A1 - a*m*A2;
 *kappa2 = -r*A2 - a*m*A1;
}


DEVICEFUNC
double polarization_angle_infty(double a, double inc, double alpha, double beta, sim5complex kappa)
{
 double kappa1 = creal(kappa);
 double kappa2 = cimag(kappa);
 double S = -alpha - a*sin(inc);
 double T = +beta;
 double X = (+S*kappa2 - T*kappa1)/(S*S+T*T);
 double Y = (-S*kappa1 + T*kappa2)/(S*S+T*T);
 return atan2(Y,X);
}
*/


/*

DEVICEFUNC INLINE
void stokes_add(stokesp* dest, stokesp value)
{
 dest->i += value.i;
 dest->q += value.q;
 dest->u += value.u;
 dest->v += value.v;
}


DEVICEFUNC INLINE
double stokes_poldeg(stokesp sp)
{
 return (sp.i>0.0) ? sqrt(sqr(sp.q) + sqr(sp.u))/sp.i : 0.0;
}


DEVICEFUNC INLINE
double stokes_polang(stokesp sp)
{
 double pxang = (sp.i>0.0) ? 0.5*atan2(sp.u/sp.i,sp.q/sp.i) : 0.0;
 while (pxang < 0.0)  pxang += 2.*M_PI;
 while (pxang > M_PI) pxang -= M_PI;
 return pxang;
}


DEVICEFUNC INLINE
void lorentz_boost_x(double V, double X[])
{
 double t = X[0];
 double x = X[1];
 double gamma = 1./sqrt(1.-V*V);
 X[0] = gamma * (t + V*x);
 X[1] = gamma * (x + V*t);
}

*/











DEVICEFUNC
double blackbody_Iv(double T, double hardf, double cos_mu, double E)
{
 if (T<=0.0) return 0.0;
 
 double limbf = (cos_mu>=0.0) ? 0.5+0.75*cos_mu : 1.0;
 double freq  = kev2freq*E;
 return limbf * 2.0*planck_h*sqr3(freq)/sqr(speed_of_light)/sqr4(hardf) / expm1((planck_h*freq)/(boltzmann_k*hardf*T)) * (1./freq2kev);
 
}



DEVICEFUNC
void blackbody(double T, double hardf, double cos_mu, double E[], double Iv[], int N)
{
 if (T<=0.0) return;

 
 int i;
 double limbf = (cos_mu>=0.0) ? 0.5+0.75*cos_mu : 1.0;
 double BB1 = limbf * 2.0*planck_h/sqr(speed_of_light)/sqr4(hardf)*sqr4(kev2freq);
 double BB2 = (planck_h*kev2freq)/(boltzmann_k*hardf*T);
 for (i=0; i<N; i++) Iv[i] = BB1*sqr3(E[i])/expm1(BB2*E[i]);
}



DEVICEFUNC INLINE
double blackbody_photons(double T, double hardf, double cos_mu, double E)
{
 return blackbody_Iv(T, hardf, cos_mu, E) / (E*kev2erg);
}



DEVICEFUNC
double blackbody_photons_total(double T, double hardf)
{
 
 return M_PI * 4.808227612 * sqr3(T) * sqr3(boltzmann_k) / sqr3(planck_h) / speed_of_light2 / hardf;
}



DEVICEFUNC
double blackbody_photon_energy_random(double T)
{
 double u1 = sim5urand();
 double u2 = sim5urand();
 double u3 = sim5urand();
 double u4 = sim5urand();
 int m;
 for (m=1; ; m++) {
     double sum_j = 0.0;
     for (int i=1; i<=m; i++) sum_j += 1.0/(double)(i*i*i);
     if (1.202*u1 < sum_j) break;
 }
 return boltzmann_k*T * (-log(u2*u3*u4)) / (double)(m) * erg2kev;
}



