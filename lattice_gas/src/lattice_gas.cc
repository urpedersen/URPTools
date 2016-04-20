 /**
 *      Author: Ulf R. Pedersen
 */

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
#include <string>
#include <sstream>
#include <fstream>
#include <cassert>
using namespace std;

#include "lattice.h"

int main(int ac, char* av[]) {

	//system("touch afile");

	// Handle program setting using Boost library
    string config_file ;
    string save_conf ;
    string save_conf_restart ;
    string load_conf ;
    unsigned int verbose ;
    double beta ;
    unsigned int types ;
    vector<double> chem_pot ;
    vector<double> umbrella_center ;
    vector<unsigned int> umbrella_center_halfperiod;
    vector<double> umbrella_kappa;
	unsigned int Lx ;
	unsigned int Ly ;
	unsigned int steps ;
	unsigned int freq_ener ;
	unsigned int freq_conf ;
	string ini_lattice_filling ;
	double weight_particle_insertion ;
	double weight_neighbour_swop     ;
	double weight_switch_identity    ;
	unsigned int seed				 ;
	try {
		po::options_description generic("Program options");
		generic.add_options()
				("help,h", "Print this help message and exit.")
	            ("input,i", po::value<string>(&config_file)->default_value("input.cfg"),"Path to input  file with simulation parameters.")
	            ("verbose,v", po::value<unsigned int>(&verbose)->default_value(5),  "Integer [0,9], 0=silent and 9=loud.")
	    ;

		po::options_description config("Simulation parameters (from input file or command line)");
		config.add_options()
						("inverse_temperature,b",
								po::value<double>(&beta)->default_value(1.0),
								"Inverse temperature, beta = 1/kT.")
						("number_of_types,t",
								po::value<unsigned int>(&types)->default_value(2),
								"Number of chemical species.")
						("chemical_potential,m",
								po::value< vector<double> >(),
								"Set chemical potentials.\n  (Must be used once for each species.)")
						("umbrella_kappa",
								po::value< vector<double> >(),
								"Set spring constant, kappa, of umbrella potentials.\n  (Must be used once for each species.)")
								("umbrella_center",
								po::value< vector<double> >(),
								"Set anchor point, a, of umbrella potentials.\n  (Must be used once for each species.)")
						("umbrella_center_period",
								po::value< vector<unsigned int> >(),
								"Set period of center of umbrella potentials to sweep over all number of particles.\n  Set to 0 to disable.\n  umbrella_center value is note used\n  (Must be used once for each species.)")
						("lattice_size_in_X,X",
								po::value<unsigned int>(&Lx)->default_value(8),
								"Lattice size in x direction.")
						("lattice_size_in_Y,Y",
								po::value<unsigned int>(&Ly)->default_value(16),
								"Lattice size in y direction.")
						("steps,t",
								po::value<unsigned int>(&steps)->default_value(10000),
								"Number of MC steps.")
						("write_frequency_ener,n",
								po::value<unsigned int>(&freq_ener)->default_value(50),
								"Write frequency of energy stats.")
						("write_frequency_conf,c",
								po::value<unsigned int>(&freq_conf)->default_value(1000),
								"Write frequency of configuration stats.")
						("weight_particle_insertion",
								po::value<double>(&weight_particle_insertion)->default_value(1.0),
								"Weight of particle insertion/removal MC move")
						("weight_neighbour_swop",
								po::value<double>(&weight_neighbour_swop)->default_value(0.0),
								"Weight of particle MC swop move")
						("weight_switch_identity",
								po::value<double>(&weight_switch_identity)->default_value(0.0),
								"Weight of particle change identity MC move")
						("seed,s",
										po::value<unsigned int>(&seed)->default_value(2013),
										"Seed for the pseudo random number generator.")
						("ini_lattice_filling,f",
								po::value<string>(&ini_lattice_filling)->default_value("empty"),
								"Set the initial filling of the lattice:\n   empty = empty lattice,\n   load  = load from disk,\n   half  = half filled lattice.")
						("load_conf,l",
								po::value<string>(&load_conf)->default_value("ini.grid"),
								"Where to load initial configuration.")
			            ("save_conf,e",
			            		po::value<string>(&save_conf)->default_value("end.grid"),
			            		"Where to save final configuration.")
			            ("save_conf_restart",
			            		po::value<string>(&save_conf_restart)->default_value("restart.grid"),
			            		"Where to save restart configurations.")
		;

		po::options_description cmdline_options;
		cmdline_options.add(generic).add(config);

		po::options_description config_file_options;
		config_file_options.add(config);

		po::options_description visible("Input options");
		visible.add(generic).add(config);

		po::positional_options_description p;
		p.add("config", -1);

		po::variables_map vm;
		store(po::command_line_parser(ac, av).
			  options(cmdline_options).positional(p).run(), vm);
		notify(vm);

        if (vm.count("help")) {
        	cout << "Program for simulating a 2D square lattice gas: " << endl;
        	cout <<	"  - attractions of 8 neighbour particles" << endl;
        	cout <<	"  - harmonic potential on number of particles" << endl;
        	cout << "Usage example:\n  mc_lattice -X 5 -Y 5 -i input.cfg\n\n";
            cout << visible << "\n";
            return 0;
        }

		ifstream ifs(config_file.c_str());
		if (!ifs)
		{
			cout << "error: can not open input file ( " << config_file << " ). Use --help for program details." <<  endl;
			return 0;
		} else {
            store(parse_config_file(ifs, config_file_options), vm);
            notify(vm);
        }

		// Set chemical potential parameters
        chem_pot.clear();
        if (vm.count("chemical_potential"))
        	chem_pot=vm["chemical_potential"].as< vector<double> >();
        if ( chem_pot.size()==0 )
        	for(unsigned int type = 0 ; type < types ; type++)
        		chem_pot.push_back(0.0);
        if(chem_pot.size()!=types){
        	cout << "error: There are " << types << " chemical species, but " << chem_pot.size() << " chemical_potential(s). Use --help for program details.\n";
        	return 1;
        }

        // Set umbrella parameters
        umbrella_center.clear();
        if (vm.count("umbrella_center"))
        	umbrella_center=vm["umbrella_center"].as< vector<double> >();
        if ( umbrella_center.size()==0 )
        	for(unsigned int type = 0 ; type < types ; type++)
        		umbrella_center.push_back(0.0);
        if(umbrella_center.size()!=types){
        	cout << "error: There are " << types << " chemical species, but " << umbrella_center.size() << " umbrella_center(s). Use --help for program details.\n";
        	return 1;
        }

        umbrella_center_halfperiod.clear();
         if (vm.count("umbrella_center_period"))
         	umbrella_center_halfperiod=vm["umbrella_center_period"].as< vector<unsigned int> >();
         if ( umbrella_center_halfperiod.size()==0 )
         	for(unsigned int type = 0 ; type < types ; type++)
         		umbrella_center_halfperiod.push_back(0);
         if(umbrella_center_halfperiod.size()!=types){
         	cout << "error: There are " << types << " chemical species, but " << umbrella_center_halfperiod.size() << " umbrella_center_period(s). Use --help for program details.\n";
         	return 1;
         }


        umbrella_kappa.clear();
        if (vm.count("umbrella_kappa"))
        	umbrella_kappa=vm["umbrella_kappa"].as< vector<double> >();
        if ( umbrella_kappa.size()==0 )
        	for(unsigned int type = 0 ; type < types ; type++)
        		umbrella_kappa.push_back(0.0);
        if(umbrella_kappa.size()!=types){
        	cout << "error: There are " << types << " chemical species, but " << umbrella_kappa.size() << " umbrella_kappa(s). Use --help for program details.\n";
        	return 1;
        }

	} catch(std::exception& e) {

        cout << "error: " << e.what() << ". Use --help for program details." << endl;
        return 1;

    }

	if(verbose > 4) cout << "Use --help for program details." << endl;
	assert(verbose<10);




    // Setup lattice and etc. and perform computations
	Lattice lattice;
	lattice.setup(Lx,Ly,types,seed);
	lattice.set_beta(beta);
	lattice.set_weights(
			weight_particle_insertion ,
			weight_neighbour_swop     ,
			weight_switch_identity    );
	for (unsigned int type=0;type<types;type++)
		lattice.set_chemical_potential(type,chem_pot.at(type));
	for (unsigned int type=0;type<types;type++)
		lattice.set_umbrella_parameters(type,umbrella_kappa.at(type),umbrella_center.at(type));


	//Fill lattice sites with particles for the initial configurations
	if( ini_lattice_filling == "empty" ) {

		// Do nothing

	} else if (ini_lattice_filling == "load" ) {
		// Read configuration from input file

		ifstream ifile;
		ifile.open(load_conf.c_str(),std::ifstream::in);

		if (!ifile) {
			cerr << "error: " << load_conf << " could not be opened for reading. Exit program." << endl;
			return 1;
		} else {
			// Helper variables
			string line;
			char * pEnd;

			getline(ifile,line); // Check lattice size
			Lx = strtol(line.c_str(),&pEnd,10);
			Ly = strtol(pEnd,&pEnd,10);
			if( Lx!=lattice.getLx() || Ly!=lattice.getLy()) cout << "Error: The input file has a lattice size that is inconsistent with the simulations parameters.";
			assert(Lx==lattice.getLx());assert(Ly==lattice.getLy());

			getline(ifile,line); // Comment line
			if(verbose > 4) cout << "Comment line in input file: " << line << endl;

			unsigned int x = 0;
			while( ifile.good() ) {
				unsigned int y = 0;
				getline(ifile,line);
				int type = strtol(line.c_str(),&pEnd,10);
				for ( unsigned int i = 0 ; i<Ly ; i++ ) {
					assert(x<Lx);assert(y<Ly);
					if(type>0) lattice.create_particle(x,y,type);
					type = strtol(pEnd,&pEnd,10);
					y++;
				}
				if(type>0) lattice.create_particle(x,y,type);
				x++;
			}
		}

		ifile.close();

	} else if (ini_lattice_filling == "half") {
		for(unsigned int y=0 ; y < Ly/2 ; y++ )
			for(unsigned int x=0 ; x < Lx ; x++ )
				lattice.create_particle(x,y,1);
	} else {
		cerr << "error: Unknown value of ini_lattice_filling = " << ini_lattice_filling << endl;
		return 1;
	}


	if(verbose > 4) cout << lattice.str() << endl;


	// Initialize file for writing restart configurations
	if(verbose > 4) cout << "Writing restart configurations to "<< save_conf_restart << ".\n";
	ofstream ofile_restart;
	ofile_restart.open (save_conf_restart.c_str());


	if(verbose > 4) cout << "Make " << steps << " Monte Carlo steps.\n";
	for(unsigned int t = 0 ; t < steps ; t++ ) {

		// Reset umbrella_center if the user want to scan it parameters
		for (unsigned int type=0;type<types;type++){
			if(umbrella_center_halfperiod.at(type)>0){
				unsigned int hp = umbrella_center_halfperiod.at(type) ;
				double new_center;
				unsigned int thp = t%hp ;
				if((t/hp)%2==0) {
					new_center = (double)lattice.volume() *         (double)(thp)/(double)hp   ;
				} else {
					new_center = (double)lattice.volume() * ( 1.0 - (double)(thp)/(double)hp );
				}
				lattice.set_umbrella_parameters(
						type,
						umbrella_kappa.at(type),
						new_center) ;
				//cout << new_center << endl ;
			}
		}


		if(t%freq_ener == 0 ) {
			if(verbose>4) cout << "ener: " << t << " " << lattice.energy_str() << endl;
		}

		if(t%freq_conf == 0 ) {
			if(verbose>4) cout << "conf: " << t << endl << lattice.configuration_str();
			ofile_restart << lattice.grid_str() << endl;

			stringstream filename;
			char ifmt[4];
			sprintf (ifmt,"%04d",t/freq_conf) ;
			filename << "conf" << ifmt << ".ppm";
			if(verbose>4) cout << "Write image to "<< filename.str() << endl;
			ofstream ofile;
			ofile.open (filename.str().c_str());
			ofile << lattice.ppm_str();
			ofile.close();

		}

		for(unsigned int i = 0 ; i < Lx*Ly ; i++ ) {
			lattice.mc_step();
		}

	}

	if(verbose > 4) cout << lattice.str() << endl;

	{
		// Write final configuration to file
		if(verbose > 4) cout << "Writing final configuration to "<< save_conf << ".\n";
		ofstream ofile;
		ofile.open (save_conf.c_str());
		ofile << lattice.grid_str();
		ofile.close();
	}
	{
		// Write ppm image
		string save_img = "end.ppm" ;
 		if(verbose > 4) cout << "Writing final configuration to image = "<< save_img << ".\n";
		ofstream ofile;
		ofile.open (save_img.c_str());
		ofile << lattice.ppm_str();
		ofile.close();


		//lattice.write_image ("final.png") ;
	}

	// Clean up
	ofile_restart.close();

	return 0;

}
