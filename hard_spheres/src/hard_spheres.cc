//============================================================================
// Name        : hard_spheres.cc
// Author      : Ulf R. Pedersen (May 2018)
// Build       : Note, require the Boost libary, boost.org
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

#include "HSSim.h"
#include "split.h"

using namespace std;

int main(int argc, char **argv) {

	// Set default values
	bool quiet=false;
	string ifilename = "none";
	unsigned time_steps = 20;
	unsigned frames = 250;
	double step_length = 0.05;
	double Lx=16,Ly=16,Lz=16;
	double neighbour_cutoff=1.7;
	string ofilename = "traj.xyz";
	double pressure = 1.0;
	double volume_step = 0.0;

	// Handle command line options
	vector<string> vecstr;
	int c;
	while(1) {
		static struct option long_options[] = {
				{"help",	no_argument      , 0, 'h'},
				{"quiet",	no_argument      , 0, 'q'},
				{"input",	optional_argument, 0, 'i'},
				{"time_steps", optional_argument, 0, 't'},
				{"frames",  optional_argument, 0, 'f'},
				{"step_length",optional_argument, 0, 's'},
				{"Lengths",	optional_argument, 0, 'L'},
				{"pressure",optional_argument, 0, 'p'},
				{"volume_step",optional_argument, 0, 'v'},
				{"rcut",	optional_argument, 0, 'r'},
				{"output",	optional_argument, 0, 'o'},
				{0, 0, 0, 0}
		};
		int option_index = 0;
		c = getopt_long (argc, argv, "hq:i:t:f:s:L:p:v:r:o:",long_options,&option_index);
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
			cout << "      Simulation of hard-spheres, by Ulf R. Pedersen (http://urp.dk)." << endl << endl;
			cout << " -h, --help                       Prints this help." << endl;
			cout << " -q, --quiet                      Hide program output." << endl;
			cout << " -r, --rcut=NUM        [1.7]      Neighbour cut-off distance." << endl;
			cout << " -i, --input=FILE      [none]     Input file (*.xyz, *.xyz.gz or *.atom)." << endl;
			cout << "                                    Default [none]: an ideal gas configuration." << endl;
			cout << " -t, --time_steps=INT  [20]       Time steps per frame." << endl;
			cout << " -f, --frames=INT      [250]      Number of frames." << endl;
			cout << " .s, --step_length=NUM [0.05]      The maximum particle step length." << endl;
			cout << " -L, --Lenghts=NUM     [10]          Size of periodic box," << endl;
			cout << "     --Lenghts=NUM,NUM,NUM           used when not provided in the input file." << endl;
			cout << " -p, --pressure=NUM    [1.0]       Pressure for barostat (if applied)." << endl;
			cout << " -v, --volume_step=NUM [0.0]       Volume step for barostat (if applied)." << endl;
			cout << "                                     The default is 0.0 resulting in a NVT simulation" << endl;
			cout << " -o, --output=FILE     [traj.xyz]  Output file for trajectory (*.xyz or *.xyz.gz)." << endl << endl;
			cout << "  Documentation and source code on github: https://github.com/urpedersen/URPTools " << endl << endl;
			exit(0);
			break;
		case 'q':
			quiet=true;
			break;
		case 'r':
			neighbour_cutoff = atof(optarg);
			break;
		case 'i':
			ifilename = optarg;
			break;
		case 't':
			time_steps = atoi(optarg);
			break;
		case 'f':
			frames = atoi(optarg);
			break;
		case 's':
			step_length = atof(optarg);
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
				cerr << "error: unknown input for -L, --Lengths.\nTry -h or --help for more information." << endl;
				abort();
			}
			break;
		case 'p':
			pressure = atof(optarg);
			break;
		case 'v':
			volume_step = atof(optarg);
			break;
		case 'o':
			ofilename = optarg;
			break;
		/*case '?':
			break; */
		default:
			cerr << "Error in command line parameters. Try -h or --help for more information." << endl;
			abort();
		}
	}
	if(!quiet)
			cout << "  Simulation of Hard Spheres. Type " << endl 
                             << "    " << argv[0] << " -h" << endl 
			     << "for more information." << endl;

	// Create object to compute rotational order
	HSSim sim;
	sim.set_neighbour_cutoff(neighbour_cutoff);
	int load_frame = 0;
	if(ifilename=="none")
	{
    	sim.generate_ideal_gas_positions(1000,Lx,Ly,Lz);
	  	cout << sim.info() << endl;
	} else {
		bool sucessfull_load = sim.load_xyz(ifilename,load_frame,Lx,Ly,Lz);
		if(sucessfull_load){
	  		cout << sim.info() << endl;
		} else {
	  		cout << "error: could not load " << ifilename << endl;
		}
	}
	
	// Make MC steps
	//sim.monte_carlo_NVT(50,0.1,250);
	sim.monte_carlo_NpT(time_steps,frames,step_length,pressure,volume_step);

	// TODO impliment Event-driven simulation  https://algs4.cs.princeton.edu/61event/

	// Write info to user
	sim.write_xyz("final.xyz");
	cout << sim.info() << endl;

	return 0;
}
