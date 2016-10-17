#ifndef ROTATIONAL_ORDER_H
#define ROTATIONAL_ORDER_H
/*
 * Rotational_order.h
 *
 *  Created on: Oct 13, 2016
 *      Author: Ulf R. Pedersen, www.urp.dk
 */

#include <vector>
#include <string>
using namespace std;

class Rotational_order {
private:
	vector<unsigned> type;
	vector<double> x;
	vector<double> y;
	vector<double> z;
	vector<unsigned> number_of_neighbors;
	vector<unsigned> number_of_S_connections;
	vector<double> qlm;
	vector<double> qi;
	vector<double> qlmAvg;
	vector<double> qiAvg;
	//vector< vector<double> > qlm_imag;
	double Lx;
	double Ly;
	double Lz;
	unsigned degree;	/// Typically refered to index l
	double neighbour_cutoff;

public:
	Rotational_order();
	virtual ~Rotational_order();

	void addParticle(unsigned in_type,double in_x,double in_y,double in_z);

	void compute_ql(unsigned in_degree);
	void compute_Sij(double S_min,string filename);

	//void translate_all_particles(double dx,double dy,double dz);
	//void load_xyz(ifstream& in);
	void  load_xyz(string ifilename,double Lx,double Ly,double Lz,double neighbour_cutoff);
	//void write_xyz(ostream& out);
	void write_xyz(string ofilename,double Qmin, double Qmax);

	unsigned number_of_particles();
	unsigned number_of_types();
	unsigned number_of_particles_of_type(unsigned test_type);
	double volume();
	string info(double Qmin,double Qmax);
};

#endif // ROTATIONAL_ORDER_H
