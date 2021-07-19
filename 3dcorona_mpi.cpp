//-----------------------------
// Author: Wenda Zhang (zhang@asu.cas.cz)
// Astronomical Institute, Czech Academy of Sciences
// Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//-----------------------------

/*! \file 3dcorona_mpi.cpp
The source file for binary `3dcorona_mpi`. This program first reads a parameter file. Given different values of `type`
provided in the parameter file, different kinds of job are performed.

### `type == 2`
In this case `3dcorona_mpi` progates the photons that enters that corona recorded in `sca_params.dat`. 
Similarly as `type == 1`, `3dcorona` always assumes that `sca_params.dat` is located in the parent directory.
For detail see also \ref register_sp_tridsca.

#### Output files:
The results are written into several binary files, all containing data in double precision.
One can use *calspec* provided by \ref calspec.cpp to calculate the spectrum.
    - `en0.dat`: dimensionless energies of the superphotons
    - `weight.dat`: statistical weights of the superphotons
    - `muinf.dat`: 
        - photon inclination angle at infinity, if the photon arrives at infinity
        - disc radius, if the photon strikes at the disc
    - `l.dat`: constant of motion
    - `q.dat`: constant of motion
    - `krobssign`: sign of \f$k^r\f$ at infinity or disc
    - `kmuobssign`: sign of \f$k^\theta\f$ at infinity or disc
    
Two additional files will be outputted if `pol` is on:
    - `qweight.dat`: Stokes parameter Q
    - `uweight.dat`: Stokes parameter U


Example of paramfile:
\code
[physical]
# BH mass in solar mass
m = 1e7
# mass accretion rate in Eddington rate
mdot = 0.0194
# electron temperature in keV
te = 100.

# scattering mean free path
tau = 20.
# color correction factor
fcol = 2.4

[gridsize]
# number of photon per geodesic
nphoton = 100

[option]
type = 2
# whether to do polarisation calculation; this option also affects scattering and disc photon angular distribution
pol = 0
# step size
dr = 0.01
# whether to print progress
progress = 0
\endcode
For `pol` option, see \ref register_sp_tridsca.

*/
#include <mpi.h>
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include <random>
#include <exception>
#include <experimental/filesystem>
#include "const.h"
#include "scatter.h"
//#include "gridgeod.h"
#include "tridgeo.h"
#include "utils.h"
#include "kerr.h"
#include "./external/INIReader.h"
#include "electron_population.h"

namespace fs = std::experimental::filesystem;

int main(int argc, char * argv[]) {
    // define MPI related variables 
    int size, rank;
    // Rank Tag
    const int RANK_ID_TAG = 2;
    // Tag for sending superphoton info
    const int RUN_PARAMS_TAG = 3;
    const int DIE_TAG = 4;
    int id_from;
    const int BPARAMS_LENGTH = 15;
    const int RUNPARAMS_LENGTH = 14;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Info info = MPI_INFO_NULL;
    MPI_File fh;
    std::array<double, BPARAMS_LENGTH> bparams;
    std::array<double, RUNPARAMS_LENGTH> run_params;
    
    // vector that save
    std::vector<std::vector<double>> vecarr_disc(10), vecarr_inf(8);
    std::vector<std::string> vecname_inf = {"en0", "weight", "muinf", "l", "q", "ktheta", "qweight", "uweight"};
    std::vector<std::string> vecname_disc = {"en0", "weight", "rhit", "l", "q", "kr", "ktheta", "P", "K1", "K2"};
    std::vector<int> nsca_inf, nsca_disc;
    
    // attributes
    std::vector<std::string> attnames = {"pol", "corona_type", "a", "m", "mdot", "fcol", "corona_par1", "corona_par2", "corona_par3", 
        "te", "tau", "dr"};
    std::string filename;
        
    std::string pol_str;
    std::unique_ptr<tridgeo> tridptr(nullptr);
    
    bool pol, KN, chandra;
    long gtype;
    std::vector<unsigned long> ndiscarr(size), ninfarr(size);
    unsigned long ndisc, ndisc_all, ninf, ninf_all, offset_inf = 0, offsetnow_inf = 0, offset_disc = 0, offsetnow_disc = 0, nfile;

	std::string parafile;
	std::vector<std::string> args(argv, argv + argc);
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
        
	double a, rin;
    if (!parfile.GetReal("physical", "rin", rin)) {
        std::cout << "missing parameter rin!" << std::endl;
        return 1;
    }

    if (rank == 0) {
        bool progress, percent_20_printed = false, percent_80_printed = false, percent_50_printed = false;
        std::vector<double> params, theta, muarr, darkarr, poldegarr;
        double teff, tcol_emass, flimb, dtheta, weightnorm, poldegnow, K1, K2, ftheta, fphi, qfac, ufac, cos2psi, sin2psi,
            theta1, dit, ennow;

        long ntheta, n = 0, nparam = 16, it0, it1, nseg, offset = 6;
        double m, te, te_emass, tau, fcol, mdot, dr;
        long nphoton;
        long type; // geometry, 0: sphere; 
        std::vector<double> gparams;
        
        // check if output directories exist
        /*
        std::vector<std::string> dirs = {"./inf", "./disc"};
        for (int i = 0; i < 2; ++i) {
            if (!fs::exists(dirs[i])) {
                fs::create_directory(dirs[i]);
            } else {
                if (!fs::is_directory(dirs[i])) {
                    std::cerr << dirs[i] << " is not a directory!" << std::endl;
                    return -1;
                }
            }
        }
        */
        
        // check if the two directories `inf` and `disc` exist
        if (!fs::exists("./inf") || !fs::is_directory("./inf") || !fs::exists("./disc") || !fs::is_directory("./disc")) {
            std::cerr << "./inf or ./disc does not exist or is not directory" << std::endl;
            return 1;
        }
            
        std::ofstream logfile("3dcorona.log");
        
        if (!parfile.GetInteger("option", "type", type)) {
            std::cout << "missing parameter type!" << std::endl;
            return 1;
        }
        
        if (type != 2) {
            std::cerr << "Type should be equal to 2!" << std::endl;
            return 1;
        }
        
        if (!parfile.GetReal("physical", "fcol", fcol)) {
            std::cout << "missing parameter fcol!" << std::endl;
            return 1;
        }
        if (!parfile.GetReal("physical", "m", m)) {
            std::cout << "missing parameter m!" << std::endl;
            return 1;
        }
        if (!parfile.GetReal("physical", "mdot", mdot)) {
            std::cout << "missing parameter mdot!" << std::endl;
            return 1;
        }
        if (!parfile.GetInteger("gridsize", "nphoton", nphoton)) {
            std::cout << "missing parameter nphoton!" << std::endl;
            return 1;
        }
        if (!parfile.GetBoolean("option", "pol", pol)) {
            std::cout << "missing parameter pol!" << std::endl;
            return 1;
        }
        if (!parfile.GetBoolean("option", "KN", KN)) {
            std::cout << "missing parameter KN!" << std::endl;
            return 1;
        }
        if (!parfile.GetBoolean("option", "chandra", chandra)) {
            std::cout << "missing parameter chandra!" << std::endl;
            return 1;
        }
                
        logfile << "pol = " << pol << std::endl;
        logfile << "KN = " << KN << std::endl;
        logfile << "m = " << m << " Msun" << std::endl;
        logfile << "mdot = " << mdot << " MEdd" << std::endl;
        logfile << "fcol = " << fcol << std::endl;
            
        if (!parfile.GetReal("physical", "te", te)) {
            std::cout << "missing parameter te!" << std::endl;
            return 1;
        }
        if (!parfile.GetReal("physical", "tau", tau)) {
            std::cout << "missing parameter tau!" << std::endl;
            return 1;
        }
        if (!parfile.GetReal("option", "dr", dr)) {
            std::cout << "missing parameter dr!" << std::endl;
            return 1;
        }
        if (!parfile.GetBoolean("option", "progress", progress)) {
            std::cout << "missing parameter progress!" << std::endl;
            return 1;
        }
		std::string scafile, gparamfile;
        if (!parfile.GetString("option", "scafile", scafile)) {
            std::cout << "missing parameter option/scafile!" << std::endl;
            return 1;
        }
		if (!fs::exists(scafile)) {
			std::cout << scafile + " does not exist!" << std::endl;
			return 1;
		}
        if (!parfile.GetString("option", "gparamfile", gparamfile)) {
            std::cout << "missing parameter option/gparamfile!" << std::endl;
            return 1;
        }
		if (!fs::exists(gparamfile)) {
			std::cout << gparamfile + " does not exist!" << std::endl;
			return 1;
		}

        te_emass = te / ME_KEV;
        logfile << "te = " << te << " keV" << std::endl;
        logfile << "tau = " << tau << std::endl;
        logfile << "dr = " << dr << std::endl;
        
        rdoublevec(gparams, gparamfile);
        while (gparams.size() < 4)
            gparams.push_back(0.);
        
        logfile << "nphoton = " << nphoton << std::endl;
        logfile.close();
        
        rdoublevec(params, scafile);
        a = params[0];
        gtype = (long)(std::round(params[2]));
        nseg = (params.size() - offset) / nparam;
        
        ntheta = 501;
        dtheta = 0.5 * M_PI / (double)(ntheta - 1);
        
        theta.resize(ntheta);
        std::generate(theta.begin(), theta.end(), [&n, dtheta]{return (double)(n++) * dtheta;});
        muarr.resize(ntheta);
        darkarr.resize(ntheta);
        std::transform(theta.cbegin(), theta.cend(), muarr.begin(), [](double x){return std::abs(std::cos(x));});
        poldegarr.resize(ntheta);
        calpoldeg(muarr, poldegarr);
        calchandra_darkening(muarr, darkarr);
        
        bparams = std::array<double, BPARAMS_LENGTH>({(double)(pol), (double)(gtype), a, m, mdot, fcol, gparams[0], gparams[1], 
            gparams[2], gparams[3], te_emass, tau, dr, (double)(KN), (double)(chandra)});
        
        // read hot cross section table and broadcast
        //MPI_Bcast(HOTX::LOGTHETATARR.data(), HOTX::NTHETAT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //MPI_Bcast(HOTX::LOGXARR.data(), HOTX::NX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //MPI_Bcast(HOTX::HOTXARR.data(), HOTX::NX * HOTX::NTHETAT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(bparams.data(), BPARAMS_LENGTH, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        std::array<double, 4> kmunow;
        std::complex<double> wpnow;
        std::random_device rd;
        std::mt19937_64 gen(rd());
        nt disk(a, m, mdot);

        //for (int i = 1; i < size; ++i) {
        //    MPI_Send(bparams.data(), BPARAMS_LENGTH, MPI_DOUBLE, i, BCAST_TAG, MPI_COMM_WORLD);
        //}
        std::cout << "nseg = " << nseg << std::endl;
            
        for (long i = 0; i < nseg; ++i) {
            auto pos = i * nparam + offset;
            std::copy(params.cbegin() + pos + 6, params.cbegin() + pos + 10, kmunow.begin());
        
            teff = disk.teff_kev(params[pos]);
        //std::cout << "(T / Tin)^3 = " << teff / tin << std::endl;
            tcol_emass = fcol * teff / ME_KEV * params[pos + 10];;
            
            theta1 = (params[pos + 1] <= 0.5 * M_PI) ? params[pos + 1] : M_PI - params[pos + 1];
            dit = (theta1 - 0.5 * dtheta) / dtheta;
            it0 = (long)(std::floor(dit));
            it1 = it0 + 1;
            
            if (chandra) {
                flimb = (darkarr[it0] * (theta1 - theta[it0]) + darkarr[it1] * (theta[it1] - theta1)) / dtheta;
                poldegnow = (poldegarr[it0] * (theta1 - theta[it0]) + poldegarr[it1] * (theta[it1] - theta1)) / dtheta;
            } else {
                flimb = 1.;
                poldegnow = 0.;
            }
            //std::cout << "param[i * nparam + 1] = " << params[i * nparam + 1] << std::endl;
            //std::cout << "theta[it] = "  << theta[it] << std::endl;
            
            K1 = -1. * params[pos + 13] * params[pos + 12];
            K2 = -1. * params[pos + 13] * params[pos + 11];
            ftheta = -1. * params[pos + 15];
            fphi = params[pos + 14];
            sin2psi = 2. * ftheta * fphi;
            cos2psi = 2. * ftheta * ftheta - 1.;
            qfac = poldegnow * cos2psi;
            ufac = poldegnow * sin2psi;
           
            weightnorm = params[pos + 3] * disk.m * disk.m  * flimb * teff * teff * teff / fcol / (double)(nphoton);
        
            for (long ip = 0; ip < nphoton; ++ip) {
                ennow = sample_bb(tcol_emass, gen);
                    
                run_params = std::array<double, RUNPARAMS_LENGTH>({params[pos + 4], params[pos + 5], 
                    kmunow[0], kmunow[1], kmunow[2], kmunow[3], ennow, poldegnow,
                    K1, K2, qfac, ufac, weightnorm, params[pos + 2]});
                        
                MPI_Recv(&id_from, 1, MPI_INT, MPI_ANY_SOURCE, RANK_ID_TAG, MPI_COMM_WORLD, &status);
                MPI_Send(run_params.data(), RUNPARAMS_LENGTH, MPI_DOUBLE, id_from, RUN_PARAMS_TAG, MPI_COMM_WORLD);
                    
            }
            if (progress) {
                showprogress((double)(i) / (double)(nseg));
            }
            if (!percent_20_printed && i > nseg * 2 / 10) {
                std::cout << "20% completed" << std::endl;
                percent_20_printed = true;
            }
            
            if (!percent_50_printed && i > nseg / 2) {
                std::cout << "50% completed" << std::endl;
                percent_50_printed = true;
            }
            if (!percent_80_printed && i > nseg * 8 / 10) {
                std::cout << "80% completed" << std::endl;
                percent_80_printed = true;
            }
        }
        std::cout << "Scattering calculation finished" << std::endl;
        for (int i = 1; i < size; ++i) {
            MPI_Send(0, 0, MPI_DOUBLE, i, DIE_TAG, MPI_COMM_WORLD);
        }
        ndisc = 0;
        ninf = 0;
        //std::cout << "unsigned int size = " << sizeof(unsigned int) << std::endl;
        //std::cout << "double size = " << sizeof(double) << std::endl;

    } else {
        //HOTX::LOGTHETATARR.resize(HOTX::NTHETAT);
        //HOTX::LOGXARR.resize(HOTX::NX);
        //HOTX::HOTXARR.resize(HOTX::NX * HOTX::NTHETAT);
        
        //MPI_Bcast(HOTX::LOGTHETATARR.data(), HOTX::NTHETAT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //MPI_Bcast(HOTX::LOGXARR.data(), HOTX::NX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //MPI_Bcast(HOTX::HOTXARR.data(), HOTX::NX * HOTX::NTHETAT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(bparams.data(), BPARAMS_LENGTH, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        //std::array<double, BPARAMS_LENGTH> bparams;
        //std::array<double, RUNPARAMS_LENGTH> run_params;
        nt disk(bparams[2], 1., 1.);
        pol = (bool)(bparams[0]);
        KN = (bool)(bparams[13]);
        gtype = (long)(bparams[1]);

		if (rin < disk.rms) {
			rin = disk.rms;
		}
		std::cout << "rin = " << rin << std::endl;
        
        std::random_device rd;
        std::mt19937_64 gen(rd());
        
        switch(gtype) {
            case 0 : {
                tridptr.reset(new zamosphere3d(disk.a, bparams[6], bparams[7], bparams[11]));
                break;
            }
            case 1 : {
                tridptr.reset(new zamoslab3d(disk.a, bparams[6], bparams[7], bparams[8], bparams[11]));
                break;
            }
            case 2 : {
                tridptr.reset(new zamosphere3d_truncated(disk.a, bparams[6], bparams[7], bparams[11]));
                break;
            }
            case 3 : {
                tridptr.reset(new rotsphere3d(disk.a, bparams[6], bparams[7], bparams[8], bparams[11]));
                break;
            }
            case 5 : {
                tridptr.reset(new kepslab3d(disk.a, bparams[6], bparams[7], bparams[8], bparams[11]));
                break;
            }
            case 6 : {
                tridptr.reset(new wedge3d(disk.a, bparams[6], bparams[7], bparams[8], bparams[11]));
                break;
            }
            case 7 : {
                tridptr.reset(new wedge3d_zamo_uniform(disk.a, bparams[6], bparams[7], bparams[8], bparams[11]));
                break;
            }
            case 8 : {
                tridptr.reset(new wedge3d_corotation_uniform(disk.a, bparams[6], bparams[7], bparams[8], bparams[11]));
                break;
            }
            case 9 : {
                tridptr.reset(new conical(disk.a, bparams[6], bparams[7], bparams[8], bparams[9], bparams[11]));
                break;
            }
			case 11 : {
                tridptr.reset(new zamosandwich(disk.a, bparams[6], bparams[7], bparams[8], bparams[9], bparams[11]));
                break;
            }
            case 12 : {
                tridptr.reset(new kepsandwich(disk.a, bparams[6], bparams[7], bparams[8], bparams[9], bparams[11]));
                break;
            }
            default : {
                std::cerr << "Invalid gtype value of " << gtype << std::endl;
            }
        }
        
        if (rank == 1) {
            std::ofstream logfile1("input_rank1.log", std::ios::out);
            logfile1 << "Corona type = " << tridptr->name << std::endl << 
                "rmin = " << tridptr->rmin << std::endl << "rmax = " << tridptr->rmax << std::endl;
            logfile1 << "te = " << bparams[10] << std::endl << "tau = " << bparams[11] << std::endl;
            logfile1 << "pol = " << pol << ", KN = " << KN << std::endl << "a = " << tridptr->a << std::endl << "dr = " << bparams[12] << std::endl;
            logfile1.close();
        }
        //isotropic_thermal epop(bparams[10] * ME_KEV);
        uniform_thermal edist(bparams[10] * ME_KEV);
        
        while(true) {
            MPI_Send(&rank, 1, MPI_INT, 0, RANK_ID_TAG, MPI_COMM_WORLD);
            MPI_Recv(&run_params[0], RUNPARAMS_LENGTH, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (status.MPI_TAG == DIE_TAG) {
                std::cout << "rank " << rank << " receives die tag" << std::endl;
                break;
            }
            //run_params = std::array<double, 13>({params[pos + 4], params[pos + 5], 
            //    kmunow[0], kmunow[1], kmunow[2], kmunow[3], ennow, poldegnow,
            //    K1, K2, qfac, ufac, weightnorm});
                        
            std::array<double, 4> kmu1 = {run_params[2], run_params[3], run_params[4], run_params[5]};
            superphoton sp;
            sp.kmu = kmu1;
            sp.en = run_params[6];
            sp.rnow = run_params[0];
            sp.munow = run_params[1];
            sp.poldeg = run_params[7];
            sp.wp.real(run_params[8]);
            sp.wp.imag(run_params[9]);
            sp.nscatter = 0;
            sp.weight = 1.;
            //std::cout << "rank " << rank << ", a = " << bparams[0] << ", sp.rnow = " << sp.rnow << std::endl;
            
            tridsca_sp(pol, KN, disk.a, rin, edist, run_params[12], run_params[10], run_params[11],
                run_params[13], bparams[12], sp, *tridptr, vecarr_inf, nsca_inf, vecarr_disc, nsca_disc, gen);
            //std::cout << "run_params = " << run_params << std::endl;
            //std::cout << "rnow = " << sp.rnow << ", munow = " << sp.munow << ", pdeg = " << sp.poldeg << ", muinf = " <<
            //    run_params[13] << std::endl;
        
        }
        std::cout << rank << " scattering calculation finished" << std::endl;
        /*
        for (int i = 0; i <= rank; ++i) {
            for (int j = 0; j < vecarr_inf.size(); ++j)
                vecarr_inf[j].push_back((double)(i));
            for (int j = 0; j < vecarr_disc.size(); ++j)
                vecarr_disc[j].push_back((double)(i));
        }
        */
            
        ndisc = vecarr_disc[0].size();
        ninf = vecarr_inf[0].size();
        
        //for (size_t i = 0; i < nsca_inf.size(); ++i) {
        //    if (nsca_inf[i] > 10) {
        //        std::cout << "nsca_inf[i] = " << nsca_inf[i] << std::endl;
        //    }
        //}
        //ninf = 1;
        //ndisc = rank;
    }
    // calculate total number of superphotons to infinity and disc, and broadcast the numbers; we are going to write data collectively
    MPI_Barrier(MPI_COMM_WORLD);
        
    if (rank >= 1) {
        if (rank > 1) {
            MPI_Recv(&offset_inf, 1, MPI_UNSIGNED_LONG, rank - 1, 10, MPI_COMM_WORLD, &status);
        }
        offsetnow_inf = offset_inf;
        offset_inf += ninf;
    }
    if (rank >= 1 && rank < size - 1) {
        MPI_Send(&offset_inf, 1, MPI_UNSIGNED_LONG, rank + 1, 10, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank >= 1) {
        if (rank > 1) {
            MPI_Recv(&offset_disc, 1, MPI_UNSIGNED_LONG, rank - 1, 11, MPI_COMM_WORLD, &status);
        }
        offsetnow_disc = offset_disc;
        offset_disc += ndisc;
    }
    if (rank >= 1 && rank < size - 1) {
        MPI_Send(&offset_disc, 1, MPI_UNSIGNED_LONG, rank + 1, 11, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == size - 1) {
        ninf_all = offset_inf;
        ndisc_all = offset_disc;
    }
    
    MPI_Bcast(&ndisc_all, 1, MPI_UNSIGNED_LONG, size - 1, MPI_COMM_WORLD);
    MPI_Bcast(&ninf_all, 1, MPI_UNSIGNED_LONG, size - 1, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    //std::cout << "offsetnow_inf = " << offsetnow_inf << std::endl;
    std::cout << "ninf = " << ninf << std::endl;
    
    nfile = pol ? 8 : 6;
    for (unsigned int i = 0; i < nfile; ++i) {
        filename = "./inf/" + vecname_inf[i] + ".dat";
        MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_RDWR, info, &fh);
        MPI_File_set_view(fh, offsetnow_inf * sizeof(double), MPI_DOUBLE, MPI_DOUBLE, "native", info);
        MPI_File_write_at_all(fh, 0, vecarr_inf[i].data(), vecarr_inf[i].size(), MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    filename = "./inf/nsca.dat";
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_RDWR, info, &fh);
    MPI_File_set_view(fh, offsetnow_inf * sizeof(int), MPI_INT, MPI_INT, "native", info);
    MPI_File_write_at_all(fh, 0, nsca_inf.data(), nsca_inf.size(), MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
    MPI_Barrier(MPI_COMM_WORLD);

    nfile = pol ? 10 : 7;
    for (unsigned int i = 0; i < nfile; ++i) {
        filename = "./disc/" + vecname_disc[i] + ".dat";
        MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_RDWR, info, &fh);
        MPI_File_set_view(fh, offsetnow_disc * sizeof(double), MPI_DOUBLE, MPI_DOUBLE, "native", info);
        MPI_File_write_at_all(fh, 0, vecarr_disc[i].data(), vecarr_disc[i].size(), MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    filename = "./disc/nsca.dat";
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_RDWR, info, &fh);
    MPI_File_set_view(fh, offsetnow_disc * sizeof(int), MPI_INT, MPI_INT, "native", info);
    MPI_File_write_at_all(fh, 0, nsca_disc.data(), nsca_disc.size(), MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
