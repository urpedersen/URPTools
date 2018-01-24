//============================================================================
// Author      : Ulf R. Pedersen
// Build       : g++ -O3 lattice.h lattice.cc make_lattice.cc -lboost_iostreams -o make_lattice
//============================================================================

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cfloat>
#include <unistd.h>
#include <getopt.h>
#include <vector>
#include <string>
#include <sstream>

#include <boost/iostreams/filtering_stream.hpp>   // Linker -lboost_iostreams
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include "lattice.h"
#include "split.h"    // Split a string into a vector<string>

using namespace std;

int main(int argc, char **argv) {

	// Set default values
	string lattice_type = "sc";
	unsigned nx=5;
	unsigned ny=nx;
	unsigned nz=nx;
	vector<unsigned> num_par;
	vector<double> mass;
	bool mix_positions=false;
	bool quiet = false;
	string filename = "start.xyz.gz";
	double density=1.0;
	bool   density_is_given_as_input=false;
	double min_distance=-1.0;
	double Lx=-1.0;
	double Ly=-1.0;
	double Lz=-1.0;
	double temperature=1.0;
	unsigned seed=0;

	// Handle command line options
	vector<string> vecstr;
	int c;
	while(1) {
		static struct option long_options[] = {
				{"help",	no_argument, 0, 'h'},
				{"quiet",	no_argument, 0, 'q'},
				{"lattice",	optional_argument, 0, 'l'},
				{"cells",	optional_argument, 0, 'c'},
				{"rho",	optional_argument, 0, 'r'},
				{"minimum_distance", optional_argument, 0, 'd'},
				{"length", optional_argument, 0, 'L'},
				{"num_par", optional_argument, 0, 'N'},
				{"mass", optional_argument, 0, 'u'},
				{"mix_positions", no_argument, 0, 'm'},
				{"temperature", optional_argument, 0, 'T'},
				{"seed", optional_argument, 0, 's'},
				{"output",	optional_argument, 0, 'o'},
				{0, 0, 0, 0}
		};
		int option_index = 0;
		c = getopt_long (argc, argv, "hql:c:N:u:r:d:L:mT:s:o:",long_options,&option_index);
		if(c==-1) break;
		switch (c) {
		/*case 0:
			if (long_options[option_index].flag != 0) break;
			printf ("option %s", long_options[option_index].name);
			if (optarg) printf (" with arg %s", optarg);
			printf ("\n");
			break;*/
		case 'h':
			cout << endl;
			cout << "    This program generate a crystal lattice configuration" << endl;
			cout << "      by Ulf R. Pedersen (http://urp.dk/tools)."<< endl << endl;
			cout << "Usage examples:" << endl;
			cout << argv[0] << " --lattice=fcc --cells=6,6,6" << endl;
			cout << argv[0] << " -lfcc -c6,6,6" << endl;
			cout << argv[0] << " -c10 -N800,200 -r1.2 -m -o KobAndersen.xyz.gz" << endl << endl;
			cout << "Optional flags:" << endl;
			cout << " -h, --help                    Print information on using this program." << endl;
			cout << " -q, --quiet                   Hide program output." << endl;
			cout << " -l, --lattice=STR             Set lattice type. Default: sc." << endl;
			cout << "                               Allowed values (n=\"particles in unit cell\"):" << endl;
			cout << "                                 sc     Simple cubic (n=1)." << endl;
			cout << "                                 rp     Randomly packed (n=1). Warning: Slow for large systems." << endl;
			cout << "                                 bcc    Body centered cubic (n=2)." << endl;
			cout << "                                 fcc    Face centered cubic (n=4)." << endl;
			cout << "                                 hcp    Hexagonal close packed (n=4)." << endl;
			cout << "                                 hex    Hexagonal layers in xy-planes (n=2)." << endl;
			cout << "                                 dc     Diamond cubic lattice (n=8)." << endl;
			cout << "                                 NaCl   Rock salt lattice. (n=2x4)" << endl;
			cout << "                                 CsCl   Cesium Chloride lattice (n=2x1)." << endl;
			cout << " -l, --lattice=FILE            Read unit-cell from a *.xyz or *.xyz.gz file." << endl;
			cout << " -c, --cells=INT               Set number of unit cells. Default: 5." << endl;
			cout << "     --cells=INT,INT,INT       " << endl;
			cout << " -N, --num_par=INT,INT,...     Reset particle types in lattice." << endl;
			cout << "                                 Note: The sum of particles must be < or = lattice sites." << endl;
			cout << " -u, --mass=NUM,NUM,...        Set masses of types. Default: 1." << endl;
			cout << " -m, --mix_positions           Make positions of particles on lattice sites random." << endl;
			cout << " -r, --rho=NUM                 Set number density. Default: 1." << endl;
			cout << " -L, --length=NUM              Change length of box vectors. Skip resetting when <0."<<endl;
			cout << "     --length=NUM,NUM,NUM        Set one value =0 for a 2D configuration." << endl;
			cout << " -d, --minimum_distance=NUM    Do not allow distances shorter than min_dist. Warning: Slow for large systems." << endl;
			cout << " -T, --temperature=NUM         Temperature of random velocity vectors." << endl;
			cout << " -s, --seed=INT                Seed for pseudo random numbers." << endl;
			cout << " -o, --output=FILE             *.xyz or *.xyz.gz or data.*.gz " << endl;
			cout << "                                 The latter is a data file for LAMMPS (http://lammps.sandia.gov)."<< endl;
			exit(0);
			break;
		case 'q':
			quiet=true;
			break;
		case 'l':
			lattice_type = optarg;
			break;
		case 'c':
			vecstr = split(optarg,',');
			if(vecstr.empty()){
				cerr << "error: unknown input for -c, --cells.\nTry -h or --help for more information." << endl;
				abort();
			}
			if( vecstr.size()==1 ) {
				nx = atoi(optarg);
				ny = nx;
				nz = nx;
			} else if ( vecstr.size()==3 ) {
				nx = atoi(vecstr.at(0).c_str());
				ny = atoi(vecstr.at(1).c_str());
				nz = atoi(vecstr.at(2).c_str());
			} else {
				cerr << "error: unknown input for -c, --cells.\nTry -h or --help for more information." << endl;
				abort();
			}
			break;
		case 'r':
			density=atof(optarg);
			density_is_given_as_input=true;
			break;
		case 'd':
			min_distance=atof(optarg);
			break;
		case 'L':
			vecstr = split(optarg,',');
			if( vecstr.size()==1 ) {
				Lx = atof(optarg);
				Ly = Lx;
				Lz = Lx;
			} else if (vecstr.size()==3 ) {
				Lx = atof(vecstr.at(0).c_str());
				Ly = atof(vecstr.at(1).c_str());
				Lz = atof(vecstr.at(2).c_str());
			} else {
				cerr << "error: unknown input for -s, --size.\nTry -h or --help for more information." << endl;
				abort();
			}
			break;
		case 'N':
			vecstr = split(optarg,',');
			for(unsigned i=0;i<vecstr.size();i++)
				num_par.push_back(atoi(vecstr.at(i).c_str()));
			break;
		case 'u':
			vecstr = split(optarg,',');
			for(unsigned i=0;i<vecstr.size();i++)
							mass.push_back(atof(vecstr.at(i).c_str()));
			break;
		case 'm':
			mix_positions=true;
			break;
		case 'T':
			temperature = atof(optarg);
			break;
		case 's':
			seed=atoi(optarg);
			break;
		case 'o':
			filename = optarg;
			break;
		/*case '?':
			break; */
		default:
			cerr << "Try -h or --help for more information." << endl;
			abort();
		}
	}

	if(!quiet)
			cout << argv[0] << ": Use  -h or --help for more information." << endl;

	// Apply user input to lattice object
	Lattice lattice(lattice_type,nx,ny,nz,seed);
	if(mix_positions)
		lattice.mix_positions();
	if(!num_par.empty())
		lattice.reset_particle_types(num_par);

	// Set masses of types
	if(mass.empty())
		for(unsigned i=0;i<lattice.number_of_types();i++)
			mass.push_back(1.0);
	if(lattice.number_of_types()!=mass.size()) {
		cerr << "error: number of masses given is not the same as number of types. Try -h or --help for more information." << endl;
		abort();
	}
	lattice.reset_mass_of_types(mass);

	// Reset box volume
	if(density_is_given_as_input)
		lattice.set_density(density);
	if(!(Lx<0.0))
		lattice.scale_x_coordinates(Lx);
	if(!(Ly<0.0))
		lattice.scale_y_coordinates(Ly);
	if(!(Lz<0.0))
		lattice.scale_z_coordinates(Lz);
	if(min_distance>0.0)	// TODO lattice.get_min_distance() scales badly, and could make the program slow
		if(min_distance>lattice.get_min_distance())
			lattice.set_min_distance(min_distance);
	if(!quiet)
		cout << lattice.info();


	{	// Write configuration to file
		using namespace boost::iostreams;
		vector<string> fnames = split(filename,'.');
		filtering_ostream out;
		if(fnames.size()<2){
			cerr << "error: Incompatible name of output file. Should be *.xyz, *.xyz.gz or data.*.gz (lammps data file)" << endl;
			abort();
		}
		if(fnames.back()=="xyz"){
			out.push(file_sink(filename));
			lattice.write_xyz(out,temperature);
		} else if ( fnames.back()=="gz" && fnames.at(fnames.size()-2)=="xyz" )  {
			out.push(gzip_compressor());
			out.push(file_sink(filename));
			lattice.write_xyz(out,temperature);
		} else if (fnames.back()=="gz" && fnames.front()=="data"){
			out.push(gzip_compressor());
			out.push(file_sink(filename));
			lattice.write_lammps(out,temperature);
		}
		else {
			cerr << "error: Incompatible name of output file. Should be *.xyz or *.xyz.gz" << endl;
			abort();
		}

	}
	if(!quiet)
		cout << "Write configuration to " << filename  << " with temperature " << temperature << endl;


	return 0;
}
