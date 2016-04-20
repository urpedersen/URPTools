/**
 * Lattice.h
 *
 *  Created on: Aug 15, 2012
 *      Author: Ulf R. Pedersen
 */

#ifndef LATTICE_H_
#define LATTICE_H_

#include <vector>
#include <string>
#include <sstream>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

using namespace std;

class Lattice {
private:
	// Thermodynamic data
	double beta;
	vector<double> chem_pot;

	// Grid data
	unsigned int Lx,Ly;
	vector<unsigned int> grid;

	// Particle data
	unsigned int types;
	vector<vector<unsigned int> > position_x; // type --> particle index
	vector<vector<unsigned int> > position_y;

	// Umbrella parameters
	vector<double> umbrella_center ;
	vector<double> umbrella_kappa  ;

	// MC
	double weight_particle_insertion ;
	double weight_neighbour_swop     ;
	double weight_switch_identity    ;
	unsigned int tries_particle_insertion ;
	unsigned int tries_neighbour_swop ;
	unsigned int tries_switch_identity ;
	unsigned int success_particle_insertion  ;
	unsigned int success_neighbour_swop ;
	unsigned int success_switch_identity  ;

	boost::random::mt19937 gen;

public:
	Lattice();
	virtual ~Lattice();

	// Manipulations
	void setup ( unsigned int , unsigned int , unsigned int , unsigned int ) ;
	void set_beta( double ) ;
	void set_chemical_potential( unsigned int , double ) ;
	void create_particle(unsigned int,unsigned int,unsigned int);
	void abolish_particle(unsigned int,unsigned int);
	void swop_particles(unsigned int,unsigned int,unsigned int,unsigned int);

	// MC stuff
	void set_weights(double,double,double);

	unsigned int mc_step();
	unsigned int mc_step_insert();
	unsigned int mc_step_swop ( );
	unsigned int mc_step_type ( );
	// Return stuff

	unsigned int volume();
	unsigned int area();
	unsigned int getLx();
	unsigned int getLy();


	double site_energy(unsigned int,unsigned int);
	double site_energy_8neighbours(unsigned int,unsigned int);
	double site_energy_4neighbours(unsigned int,unsigned int);
	double site_energy_type(unsigned int,unsigned int,unsigned int);
	double total_energy();

	void set_umbrella_parameters ( unsigned int , double , double ) ;
	double umbrella_energy(unsigned int,unsigned int) ;
	double umbrella_energy(unsigned int) ;

	unsigned int type_at(unsigned int,unsigned int);

	string str();
	string grid_str();
	string configuration_str();
	string energy_str();
	string ppm_str( );

	unsigned int number_of_particles(unsigned int) ;
	unsigned int number_of_particles() ;
	unsigned int site ( int , int ) ;
	unsigned int x_img (int) ;
	unsigned int y_img (int) ;

};

#endif /* LATTICE_H_ */
