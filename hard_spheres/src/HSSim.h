#ifndef HSSIM_H
#define HSSIM_H
/*
 * HSSim.h
 *
 *  Created on: May 2, 2018
 *      Author: Ulf R. Pedersen, www.urp.dk
 */

#include <vector>
#include <string>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>


#include "Cell_list.h"

using namespace std;

class HSSim {
private:
	vector<unsigned> type;
	vector<double> x;
	vector<double> y;
	vector<double> z;
	Cell_list cell_list;
	double Lx;
	double Ly;
	double Lz;
	double neighbour_cutoff;
	
	boost::random::mt19937 gen;

public:
	HSSim();
	virtual ~HSSim();
	
	void clear();  /// Clear memory for a new computation
	
	void addParticle(unsigned in_type,double in_x,double in_y,double in_z);

  void generate_ideal_gas_positions(unsigned num_particles,double in_Lx,double in_Ly,double in_Lz);
  bool load_xyz(string ifilename);
	bool load_xyz(string ifilename,unsigned frame,double Lx,double Ly,double Lz);
	void build_neighbour_list();
	void wrap_into_box(double xO,double yO,double zO);
	void write_xyz();
	void write_xyz(string ofilename);

	unsigned number_of_particles();
	unsigned number_of_types();
	unsigned number_of_particles_of_type(unsigned test_type);
	double volume();

	void monte_carlo(unsigned steps,double stepSize,unsigned frames);

	string info();
};

#endif // HSSIM_H
