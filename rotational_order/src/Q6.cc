//============================================================================
// Name        : Q6.cc
// Author      : Ulf R. Pedersen (Oct. 2016)
// Build       : g++ -O3 Rotational_order.h Rotational_order.cc Q6.cc -lboost_iostreams -o Q6
//               Require the Boost libary, boost.org
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

#include "Rotational_order.h"

using namespace std;


/**
 * Split a string into a vector string
 */

vector<string> split(string str, char delimiter) {
  vector<string> output;stringstream ss(str);string substr;
  while(getline(ss, substr, delimiter)) output.push_back(substr);
  return output;
}

int main(int argc, char **argv) {

	// Set default values
	bool quiet=false;
	unsigned degree = 6;
	string ifilename = "traj.xyz";
	unsigned frame=0;
	double Lx=10,Ly=10,Lz=10;
	double neighbour_cutoff=1.4;
	string ofilename = "none";
	double Qmin=0.0;
	double Qmax=1.0;
	double Sij_min=-1.0;

	// Handle command line options
	vector<string> vecstr;
	int c;
	while(1) {
		static struct option long_options[] = {
				{"help",	no_argument, 0, 'h'},
				{"quiet",	no_argument, 0, 'q'},
				{"degree",	optional_argument, 0, 'l'},
				{"input",	optional_argument, 0, 'i'},
				{"frame",	optional_argument, 0, 'f'},
				{"Lengths",	optional_argument, 0, 'L'},
				{"rcut",	optional_argument, 0, 'r'},
				{"QminQmax",	optional_argument, 0, 'Q'},
				{"Sij",		optional_argument, 0, 'S'},
				{"output",	optional_argument, 0, 'o'},
				{0, 0, 0, 0}
		};
		int option_index = 0;
		c = getopt_long (argc, argv, "hql:i:f:L:r:Q:S:o:",long_options,&option_index);
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
			cout << "      Compute the Q6 rotational bond order-parameters by Steinhard (ql)," << endl;
			cout << "              and the Lechner-Dellago averaged versions (qlAvl)." << endl;
			cout << " Local bond-orders parameters are printed to the last columns of an xyz-file." << endl << endl;
			cout << "              Written by Ulf R. Pedersen (2016), www.urp.dk." << endl << endl;
			cout << "Usage examples:" << endl;
			cout << argv[0] << " --input=input.xyz --Lenghts=15.0,15.0.30.0" << endl;
			cout << argv[0] << " -i bcc.xyz -L 15.0 -r 1.5 -l 4 -o Q4.xyz" << endl << endl;
			cout << "Optional flags:     [default]"      << endl;
			cout << " -h, --help                     Print information on using this program." << endl;
			cout << " -q, --quiet                    Hide program output." << endl;
			cout << " -l, --l=INT        [6]         The degree of the bond rotational-order." << endl;
			cout << " -r, --rcut=NUM     [1.4]       Neighbour cutoff distance." << endl;
			cout << " -i, --input=FILE   [input.xyz] Input file (*.xyz, *.xyz.gz or *.atom)." << endl;
			cout << " -f, --frame=INT    [0]         Frame of input file (0=first frame)." << endl;
			cout << " -L, --Lenghts=NUM  [10]        Size of periodic box," << endl;
			cout << "     --Lenghts=NUM,NUM,NUM        unless it is provided in the input file." << endl; 
			cout << " -Q  --QminQmax=NUM,NUM         Minimum and maximum Q6 limits" << endl;
			cout << "                    [0.0,1.0]     Default values includes all particles." << endl;
			cout << " -S, --Sij=NUM      [-1.0]      Min threshold values for the Sij connection matrix." << endl;
			cout << "                                  Default value -1.0 skip computation." << endl;
			cout << "                                  A sparse matrix is written to node_connections.dat and" << endl;
			cout << "                                  the largest cluster is written to largest_cluster.xyz." << endl;
			cout << "                                  The Lechner-Dellago vectors are for Sij." << endl; 
			cout << " -o, --output=FILE  [none]      Output file (*.xyz or *.xyz.gz)." << endl;
			cout << "                                  Default value none does not produce a file." << endl;
			cout << "                                  The 5th column gives number of neighbours. " << endl;
			cout << "                                  The 6th column gives Steinhard bond-order. " << endl;
			cout << "                                  The 7th column gives Lechner-Dellago bond-order. " << endl;
			cout << "                                  The 8th column gives the number of Sij connections, " << endl;
			cout << "                                  if the -S flag is used. " << endl;
			exit(0);
			break;
		case 'q':
			quiet=true;
			break;
		case 'l':
			degree = atoi(optarg);
		case 'r':
			neighbour_cutoff = atof(optarg);
			break;
		case 'i':
			ifilename = optarg;
			break;
		case 'f':
			frame = atoi(optarg);
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
				cerr << "error: unknown input for -L, --Lengths.\nTry -h or --help for more information." << endl;
				abort();
			}
			break;
		case 'Q':
			vecstr = split(optarg,',');
			if( vecstr.size()==2 ) {
				Qmin = atof(vecstr.at(0).c_str());
				Qmax = atof(vecstr.at(1).c_str());
			} else {
				cerr << "error: unknown input for -Q, --QminQmax.\nTry -h or --help for more information." << endl;
			}
		case 'S':
			Sij_min=atof(optarg);
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
			cout << "Computing bond orientational order. Type " << endl 
                             << "    " << argv[0] << " -h" << endl 
			     << "for more information." << endl;

	// Create object to compute rotational order
	Rotational_order rot;
	rot.load_xyz(ifilename,frame,Lx,Ly,Lz,neighbour_cutoff);
	rot.compute_ql(degree);
	if(Sij_min>-1.0) {
		rot.compute_Sij(Sij_min,"node_connections.dat","largest_cluster.xyz");
		if(!quiet)
			cout << "Wrote Sij matrix to node_connections.dat and largest cluster of connected particles to largest_cluster.xyz." << endl;
	}
	if(ofilename!="none"){
		rot.write_xyz(ofilename,Qmin,Qmax);
		cout << "Wrote " << ofilename <<" with local bond-orderparamters of particles within Qmin and Qmax values." << endl;
	}

	// Say goodby to the nice user (unless you are asked to be quiet). 
	if(!quiet){
		cout << rot.info(Qmin,Qmax) << endl << endl;
	}
	
	return 0;
}
