//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

/*! \file 3dcorona.cpp
The source file for binary `3dcorona`. This program first reads a parameter file. Given different values of `type`
provided in the parameter file, different kinds of job are performed.

### Syntax 
    - `3dcorona`: reads default parameter file `params.txt`
    - `3dcorona parafile`: now the parameter file is `parafile`

In this case `3dcorona` calls \ref grid3dcorona_sp_writeparam. Basically this function writes geodesic info from equatorial plane disc
into the following two binary files:
    - `disc_params.dat`: info of photons escaping to infinity or striking at the disc without entering the corona
    - `sca_params.dat`: info of photons entering the corona
For details of these two files see \ref grid3dcorona_sp_writeparam.

An example of the parameter file:
\code
[physical]
# BH spin
a = 0.998
# disc outer edge radius in GM/c^2
rout = 1000.
# corona type
ctype = 0
h = 10.
radius = 1.

[gridsize]
nr = 100
ntheta = 100
nphi = 100

[option]
# whether to calculate polarisation; the directionarity of disc emission also depends on this option
pol = 0
\endcode
For `ctype`, see \ref grid3dcorona_sp_writeparam; for gridsize parameters: see \ref calmuobsarr; 
for `pol`: see \ref register_sp_nt and \ref register_sp_tridsca.
*/
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include <random>
#include <experimental/filesystem>
#include "const.h"
#include "scatter.h"
#include "gridgeod.h"
#include "tridgeo.h"
#include "utils.h"
#include "./external/INIReader.h"
#include "electron_population.h"

namespace fs = std::experimental::filesystem;

int main(int argc, char * argv[]) {
    std::vector<std::string> args(argv, argv + argc);
    std::string parafile;
    
    long gtype; // geometry, 0: sphere; 
    double a, rin, rout, radius, h, s, hmin, hmax, thetamax,
        corona_rin, velocity, smin, smax, maxmu, thickness, opening_angle;
    long nr, ntheta, nphi;
    
    if (argc == 1) {
        parafile = "params.txt";
    } else {
        parafile = args[1];
    }
    
    if (!fs::exists(parafile)) {
        std::cerr << parafile + " does not exist!" << std::endl;
        return 1;
    }
    
    INIReader parfile(parafile);
        
    std::ofstream logfile("3dcorona.log");
    std::random_device rd;
    std::mt19937_64 gen(rd());

    if (!parfile.GetReal("physical", "a", a)) {
        std::cout << "missing parameter a!" << std::endl;
        return 1;
    }
    if (!parfile.GetReal("physical", "rin", rin)) {
        std::cout << "missing parameter rin!" << std::endl;
        return 1;
    }
    nt disc(a, 1., 1.);
    if (rin < disc.rms)
        rin = disc.rms;
    if (!parfile.GetReal("physical", "rout", rout)) {
        std::cout << "missing parameter rout!" << std::endl;
        return 1;
    }
    //if (!parfile.GetReal("gridsize", "thetamax", thetamax)) {
    //    std::cout << "missing parameter thetamax!" << std::endl;
    //    return 1;
    //}
	thetamax = 0.5;
    if (!parfile.GetInteger("physical", "gtype", gtype)) {
        std::cout << "missing parameter gtype!" << std::endl;
        return 1;
    }
    if (gtype == 0) {
        if (!parfile.GetReal("physical", "h", h)) {
            std::cout << "missing parameter h!" << std::endl;
            return 1;
        }
        if (!parfile.GetReal("physical", "radius", radius)) {
            std::cout << "missing parameter radius!" << std::endl;
            return 1;
        }
    } else if (gtype == 1 || gtype == 5) {
        if (!parfile.GetReal("physical", "s", s)) {
            std::cout << "missing parameter s!" << std::endl;
            return 1;
        }
        if (!parfile.GetReal("physical", "hmin", hmin)) {
            std::cout << "missing parameter hmin!" << std::endl;
            return 1;
        }
        if (!parfile.GetReal("physical", "hmax", hmax)) {
            std::cout << "missing parameter hmax!" << std::endl;
            return 1;
        }
    } else if (gtype == 2) {
        if (!parfile.GetReal("physical", "radius", radius)) {
            std::cout << "missing parameter h!" << std::endl;
            return 1;
        }
        if (!parfile.GetReal("physical", "corona_rin", corona_rin)) {
            std::cout << "missing parameter radius!" << std::endl;
            return 1;
        }
        if (corona_rin <= disc.rms)
            corona_rin = disc.rms;
        if (corona_rin > radius) {
            std::cerr << "corona rin > r!" << std::endl;
            return 1;
        }
    } else if (gtype == 3) {
        if (!parfile.GetReal("physical", "h", h)) {
            std::cout << "missing parameter h!" << std::endl;
            return 1;
        }
        if (!parfile.GetReal("physical", "radius", radius)) {
            std::cout << "missing parameter radius!" << std::endl;
            return 1;
        }
        if (!parfile.GetReal("physical", "velocity", velocity)) {
            std::cout << "missing parameter velocity!" << std::endl;
            return 1;
        }
    } else if (gtype == 6 || gtype == 7 || gtype == 8) {
        if (!parfile.GetReal("physical", "smin", smin)) {
            std::cout << "missing parameter smin!" << std::endl;
            return 1;
        }
        if (smin < disc.rms) {
            smin = disc.rms;
        }
        if (!parfile.GetReal("physical", "smax", smax)) {
            std::cout << "missing parameter smax!" << std::endl;
            return 1;
        }
        if (!parfile.GetReal("physical", "maxmu", maxmu)) {
            std::cout << "missing parameter maxmu!" << std::endl;
            return 1;
        }
    } else if (gtype == 9) {
        if (!parfile.GetReal("physical", "hmin", hmin)) {
            std::cout << "missing parameter hmin!" << std::endl;
            return 1;
        }
        if (!parfile.GetReal("physical", "thickness", thickness)) {
            std::cout << "missing parameter thickness!" << std::endl;
            return 1;
        }
        if (!parfile.GetReal("physical", "opening_angle", opening_angle)) {
            std::cout << "missing parameter opening_angle!" << std::endl;
            return 1;
        }
        if (!parfile.GetReal("physical", "velocity", velocity)) {
            std::cout << "missing parameter velocity!" << std::endl;
            return 1;
        }
	} else if (gtype == 11 || gtype == 12) {
        if (!parfile.GetReal("physical", "s", s)) {
            std::cout << "missing parameter s!" << std::endl;
            return 1;
        }
        if (!parfile.GetReal("physical", "smin", smin)) {
            std::cout << "missing parameter smin!" << std::endl;
            return 1;
        }
        if (!parfile.GetReal("physical", "hmin", hmin)) {
            std::cout << "missing parameter hmin!" << std::endl;
            return 1;
        }
        if (!parfile.GetReal("physical", "hmax", hmax)) {
            std::cout << "missing parameter hmax!" << std::endl;
            return 1;
        }
    }
        
    if (!parfile.GetInteger("gridsize", "nr", nr)) {
        std::cout << "missing parameter nr!" << std::endl;
        return 1;
    }
    if (!parfile.GetInteger("gridsize", "ntheta", ntheta)) {
        std::cout << "missing parameter ntheta!" << std::endl;
        return 1;
    }
    if (!parfile.GetInteger("gridsize", "nphi", nphi)) {
        std::cout << "missing parameter nphi!" << std::endl;
        return 1;
    }
    
    std::unique_ptr<tridgeo> tridptr(nullptr);
    std::vector<double> gparams;
    switch (gtype) {
        case 0 : {
            tridptr.reset(new zamosphere3d(a, h, radius, 0.1));
            gparams = std::vector<double> ({h, radius});
            break;
        }
        case 1 : {
            tridptr.reset(new zamoslab3d(a, s, hmin, hmax, 0.1));
            gparams = std::vector<double> ({s, hmin, hmax});
            break;
        }
        case 2 : {
            tridptr.reset(new zamosphere3d_truncated(a, radius, corona_rin, 0.1));
            gparams = std::vector<double> ({radius, corona_rin});
            break;
        }
        case 3 : {
            tridptr.reset(new rotsphere3d(a, h, radius, velocity, 0.1));
            gparams = std::vector<double> ({h, radius, velocity});
            break;
        }
        
        case 5 : {
            tridptr.reset(new kepslab3d(a, s, hmin, hmax, 0.1));
            gparams = std::vector<double> ({s, hmin, hmax});
            break;
        }
        
        case 6 : {
            tridptr.reset(new wedge3d(a, smin, smax, maxmu, 0.1));
            gparams = std::vector<double> ({smin, smax, maxmu});
            break;
        }
        case 7 : {
            tridptr.reset(new wedge3d_zamo_uniform(a, smin, smax, maxmu, 0.1));
            gparams = std::vector<double> ({smin, smax, maxmu});
            break;
        }
        
        case 8 : {
            tridptr.reset(new wedge3d_corotation_uniform(a, smin, smax, maxmu, 0.1));
            gparams = std::vector<double> ({smin, smax, maxmu});
            break;
        }
        
        case 9 : {
            tridptr.reset(new conical(a, hmin, thickness, opening_angle, velocity, 0.1));
            gparams = std::vector<double> ({hmin, thickness, opening_angle, velocity});
            break;
        }
		case 11 : {
            tridptr.reset(new zamosandwich(a, s, smin, hmin, hmax, 1.));
            gparams = std::vector<double> ({s, smin, hmin, hmax});
            break;
        }

        case 12 : {
            tridptr.reset(new kepsandwich(a, s, smin, hmin, hmax, 1.));
            gparams = std::vector<double> ({s, smin, hmin, hmax});
            break;
        }

        default : {
            std::cerr << "3dcorona: Invalid value of gtype!" << std::endl;
            return 1;
        }
    }
    wdoublevec(gparams, "gparams.dat");
    
    //truncated_nt disc(a, m, mdot, rin);
    grid3dcorona_sp_writeparam(rin, rout, thetamax, disc, *tridptr, nr, ntheta, nphi, "disc_params.dat", "selfirr_params.dat", "sca_params.dat");
    
    // log
    logfile << "a = " << a << std::endl;
    logfile << "# rin: if rin < rms, then rin = rms." << std::endl;
    logfile << "rin = " << rin << std::endl;
    logfile << "rout = " << rout << std::endl;
    logfile << "thetamax = " << thetamax << " PI" << std::endl;
    logfile << "geometry = " << tridptr->name << std::endl;
    
    switch (gtype) {
        case 0 : {
            logfile << "h = " << h << std::endl;
            logfile << "radius = " << radius << std::endl;
            break;
        }
        case 1 : {
            logfile << "hmin = " << hmin << std::endl;
            logfile << "hmax = " << hmax << std::endl;
            logfile << "s = " << s << std::endl;
            break;
        }
        case 2 : {
            logfile << "radius = " << radius << std::endl;
            logfile << "# corona_rin: if corona_rin < rms, then corona_rin = rms." << std::endl;
            logfile << "corona_rin = " << corona_rin << std::endl;
            break;
        }
        case 3 : {
            logfile << "h = " << h << std::endl;
            logfile << "radius = " << radius << std::endl;
            logfile << "velocity = " << velocity << std::endl;
            break;
        }
        case 5 : {
            logfile << "hmin = " << hmin << std::endl;
            logfile << "hmax = " << hmax << std::endl;
            logfile << "s = " << s << std::endl;
            break;
        }
        
        case 6 : {
            logfile << "smin = " << smin << std::endl;
            logfile << "smax = " << smax << std::endl;
            logfile << "maxmu = " << maxmu << std::endl;
            break;
        }
        case 7 : {
            logfile << "smin = " << smin << std::endl;
            logfile << "smax = " << smax << std::endl;
            logfile << "maxmu = " << maxmu << std::endl;
            break;
        }
        
        case 8 : {
            logfile << "smin = " << smin << std::endl;
            logfile << "smax = " << smax << std::endl;
            logfile << "maxmu = " << maxmu << std::endl;
            break;
        }
        
        case 9 : {
            logfile << "hmin = " << hmin << std::endl;
            logfile << "thickness = " << thickness << std::endl;
            logfile << "opening_angle = " << opening_angle << std::endl;
            logfile << "velocity = " << velocity << std::endl;
            break;
        }
		case 11 : {
            logfile << "hmin = " << hmin << std::endl;
            logfile << "hmax = " << hmax << std::endl;
            logfile << "smin = " << smin << std::endl;
            logfile << "s = " << s << std::endl;
            break;
        }

        case 12 : {
            logfile << "hmin = " << hmin << std::endl;
            logfile << "hmax = " << hmax << std::endl;
            logfile << "smin = " << smin << std::endl;
            logfile << "s = " << s << std::endl;
            break;
        }
        
        default : {
            std::cerr << "gtype input not valid!" << std::endl;
            return 1;
        }
    }
    logfile << "nr = " << nr << std::endl; 
    logfile << "ntheta = " << ntheta << std::endl; 
    logfile << "nphi = " << nphi << std::endl; 
	logfile.close();
}
        
