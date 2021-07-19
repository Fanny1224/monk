#### Executable binaries

After compilation one can obtain the following binaries for photon propagation:

- `3dcorona`: pre-calculates the geodesic from the disc to the corona
- `3dcorona_mpi`: ray-traces the photons based on the geodesics
  produced by `3dcorona`

and one post-processing program with which one can calculate the
energy/polarisation spectrum with files produced by `3dcorona_mpi`:

- `calspec`: calculate energy and polarisation spectrum

For detailed descriptions, please run `doxygen doc.config` to generate the documentation file and look
for the links to the corresponding pages in `./html/files.html`. See below for a short introduction.

#### Building
##### System Requirements
- GNU Make
- C++ compiler that supports c++14 standard; for GCC starting from GCC 5.0
- MPI library
- libstdc++fs; included in GCC starting from GCC 5.3

##### Procedure
1. build `sim5` library under ./sim5
```
cd sim5
make
```
2. modify line 17 of `electron_population.h` to set the correct value of
   `hotdir`, the path to the directory that contains ``logthetat.dat`, `logx.dat`, and `hotx.dat`, the
   data files that are used to save hot Klein-Nishina cross section. I put them
   in `./data`, but you have to use the absoluate path.
	
3. edit `Makefile` to assign the appropriate values to the following variables:
	- `SIM5INC`: path to directory that contains `sim5lib.h`
	- `SIM5OBJ`: path to `sim5lib.o`
	- `OBJDIR`: directory for object file
	- `BINDIR`: directory for binary files
	- `TRASHDIR`: directory for trash
	- `INSTALLDIR`: directory where you put executable binaries
5. `make objs` to build objects
6. `make 3dcorona 3dcorona_mpi calspec` 

#### Usage
##### AGN disc-corona system
1. Create a parameter file for `3dcorona` with `type = 0`. An exmple for
   spherical corona is as follows:
```
[physical]
# black hole spin
a = 0.998
# the geometry of the corona; 0 for stationary spherical coronae; 
gtype = 0
# the height of the corona
h = 10.
# radius of the corona
radius = 1.
# inner edge of the disc, in GM/c^2; if rin < rms, then rin = rms
rin = -1.
# outer edge of the disc, in GM/c^2
rout = 1000.

[gridsize]
# number of radial bins on the disc
nr = 100
# number of bins in the polar emission angle of the seed photon
nphi = 100
# number of bins in the azimuthal emission angle of the seed photon
ntheta = 100
# the maximum polar emission angle over PI; 0.5 means the maximum angle is half
# PI
thetamax = 0.5

[option]
# tells `3dcorona` to perform geodesic calculations
type = 0
```
Notes on the parameter file:
- Different kinds of paramters are needed for different geometries. In the
  following I list the parameter files for a slab corona and a truncated disc
  geometry:
- At this stage we are calculating the geodeisc from the disc to the surface of
  the corona, therefore only geometrical info the corona is required.

The parameter file for a slab corona:
```
[physical]
# black hole spin
a = 0.998
# the geometry of the corona; 1 for a stationary slab; and 5 for a slab
# co-rotating with the disc
gtype = 1
# the minimum height of the corona, in GM/c2
hmin = 9.
# the maximum height of the corona, in GM/c2
hmax = 11.
# radial extent of the corona, in GM/c2
s = 10.
# inner edge of the disc, in GM/c^2; if rin < rms, then rin = rms
rin = -1.
# outer edge of the disc, in GM/c^2
rout = 1000.

[gridsize]
# number of radial bins on the disc
nr = 100
# number of bins in the polar emission angle of the seed photon
nphi = 100
# number of bins in the azimuthal emission angle of the seed photon
ntheta = 100
# the maximum polar emission angle over PI; 0.5 means the maximum angle is half
# PI
thetamax = 0.5

[option]
# tells `3dcorona` to perform geodesic calculations
type = 0
```
And for a truncated disc:
```
[physical]
# black hole spin
a = 0.998
# the geometry of the corona; 2 for a stationary truncated disc
gtype = 2
# the radius of the corona, in GM/c2
radius = 30.
# the inner radius of the corona
corona_rin = 2.
# inner edge of the disc, in GM/c^2; if rin < rms, then rin = rms
rin = -1.
# outer edge of the disc, in GM/c^2
rout = 1000.

[gridsize]
# number of radial bins on the disc
nr = 100
# number of bins in the polar emission angle of the seed photon
nphi = 100
# number of bins in the azimuthal emission angle of the seed photon
ntheta = 200
# the maximum polar emission angle over PI; 1. means the maximum angle is PI
thetamax = 1.

[option]
# tells `3dcorona` to perform geodesic calculations
type = 0
```

2. Run `3dcorona` to produce `sca_params.dat` and `disk_params.dat`
3. Make a directory for disc photons; e.g., `disc`; and another directory for
   corona photons, e.g., `corona`
4. In `disc`, create a parameter file with `type = 1`, and run `3dcorona`. One
   example of the parameter is as follows:
```
[physical]
a = 0.998
# BH mass in solar mass
m = 1e7
# mass accretion rate in Eddington unit; defined in such a way that the physical mass accretion
# rate \dot{M} = 2.225475942e+18 * m * mdot [g/s]
mdot = 0.0194
# the color correction factor
fcol = 2.4
# the inner edge of the disc; if rtr < rms, rtr = rms. Note that we assume that
# the torque is zero at rtr, not rms if rtr lies outside rms
rtr = -1.

[gridsize]
# the number of photons sampled for one entry
nphoton = 1

[option]
# tell `3dcorona` to sample disc photons
type = 1
# if the polarized radiative transfer if switched on.
# if pol = 0, we assume the disc emission is isotropic and is unploarized.
# Otherwise we use the Chandrasekhar's formula for a semi-infinite scattering
# atmosphere to determine the disc photon's directionality and polarization angle/degree
pol = 0
# whether to show progress in the terminal
progress = 1
# the path to disc_params.dat
discparafile = ../disc_params.dat
```

5. In `corona`, make two directories `inf` and `disc`, create a parameter file with `type = 2` and run `3dcorona_mpi`.
   As `3dcorona_mpi` is a MPI program, it has to be executed differently from
   usual executables. One one's own laptop/desktop, it can be done like this:
   ```
   mpirun -n 4 3dcorona_mpi
   ```
   here I specify the number of threads to be 4. On clusters one needs to refer
   to the manual to know how to run a MPI program. One example of the parameter
   file:
```
[physical]
a = 0.998
m = 1e7
mdot = 0.0194
# corona temperature in kev
te = 100.
# Thomson optical depth; for spherical corona tau = ne * sigma_t * R, where ne
# is the electron density, sigma_t is the Thomson cross section, and R is the
# radius of the corona. For slab, tau = ne * sigma_t * h / 2, where h is the
# thickness of the corona
tau = 0.2
fcol = 2.4
rin = -1.

[gridsize]
# number of photons sampled for one geodesic
nphoton = 1000

[option]
# tell `3dcorona_mpi` to perform radiative transfer inside the corona
type = 2
# if the polarized radiative transfer is switched on
pol = 0
# the step size
dr = 0.001
# whether to print the progress
progress = 1
# whether to assume Klein-Nishina cross section
KN = 1
# if chandra = 0, we assume the disc emission is isotropic and is unploarized.
# Otherwise we use the Chandrasekhar's formula (Section X, Chandrasekhar 1960) for a semi-infinite scattering
# atmosphere to determine the disc photon's directionality and polarization angle/degree
chandra = 0
```
6. In `disc` and `corona`, run `calspec` to calculate spectra of disc and
   corona, respectively. `calspec` creates three/seven files depending on if
   polarized radiative transfer is performed. The syntax of `calspec` is:
   ```
   	# For photons going to all directions:
	# dir: the path to the data directory
	# ne: the number of energy bins. Linear scale if ne is positive; logscale
	# 	otherwise.
	# emin: lower boundary of the energy bin, in keV
	# emax: upper boundary of the energy bin, in keV
	calspec dir ne emin emax
   ```
   ```
	# For photons reaching a inclination between imin and imax (in degrees)
	calspec dir ne emin emax nsca
   ```
   For unpolarized cases, the files are:
   - `en.dat`: the center of the energy bins in keV, 
   - `de.dat`: the width of the energy bins in keV
   - `flux.dat`: the name is a bit misleading; it's actually the normalized
	 luminosity density. The normalization factor is 
	 ```
	 SPECNORM_SP = 4 * APERY_CONST * RG * RG / (C * C * H * H * H);
	 ```
	 where APERY_CONST is the Apery's constant; RG=G*M_sun/c^2 is the
	 gravitational radius for a solar-mass black hole; C is the speed of light;
	 and H is the Planck constant. SPECNORM_SP * `flux.dat` is the luminosity
	 density is [photons/s/keV].
	
	For polarized radiative transfer there are four additional files:
	- `qflux.dat`: the Stokes Q parameter; with the same unit and normalization
	  as `flux.dat`
	- `uflux.dat`: the Stokes U parameter; with the same unit and normalization
	  as `flux.dat`
	- `poldeg.dat`: the polarization degree
	- `polang.dat`: the polarization angle, in radians 

	One can also calculate the spectrum for a certain scattering order. For the
	spectrum of photons having scattered once:
	```
	calspec dir ne emin emax imin imax 1
	``` 
	and all the products have the `_1sca` suffix.

	For the spectrum of photons having scattered at least once:
	```
	calspec dir ne emin emax imin imax -1
	``` 



	`
