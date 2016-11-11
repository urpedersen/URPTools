#ifndef LATTICE_H
#define LATTICE_H
/*
 * Lattice.h
 *
 *  Created on: Nov 9, 2015
 *      Author: urp
 */

#include <vector>
#include <string>
using namespace std;

class Lattice {
private:
	string lattice_type;
	int number_of_sites_on_lattice;
	vector<unsigned> type;
	vector<double> x;
	vector<double> y;
	vector<double> z;
	double Lx;
	double Ly;
	double Lz;
	vector<double> mass_of_types;

	double random_velocity(double temperature,double mass);

public:
	Lattice(string lattice_type,unsigned nx,unsigned ny,unsigned nz,unsigned seed);
	virtual ~Lattice();

	void addParticle(unsigned in_type,double in_x,double in_y,double in_z);

	void mix_positions();
	void set_density(double new_density);
	double get_min_distance();
	void set_min_distance(double new_min_distance);
	void scale_coordinates(double sf);
	void scale_x_coordinates(double new_Lx);
	void scale_y_coordinates(double new_Ly);
	void scale_z_coordinates(double new_Lz);
	void translate_all_particles(double dx,double dy,double dz);
	void reset_particle_types(vector<unsigned>);
	void reset_mass_of_types(vector<double>);
	bool reset_particle_type(unsigned particle_index,unsigned new_type);
	void reset_number_of_particles(unsigned num_remove);
	void write_xyz(ostream& out,double temperature);

	unsigned number_of_particles();
	unsigned number_of_types();
	unsigned number_of_particles_of_type(unsigned test_type);
	double volume();
	string info();
};

#endif // LATTICE_H
