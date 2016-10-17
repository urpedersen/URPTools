#ifndef CELL_LIST_H
#define CELL_LIST_H
/*
 * Rotational_order.h
 *
 *  Created on: Oct 16, 2016
 *      Author: Ulf R. Pedersen, www.urp.dk
 */


#include <vector>
#include <string>
using namespace std;


class Cell_list {
private:
	vector<vector<unsigned> > cells;
	vector<unsigned> index2cell;
	unsigned nx;
	unsigned ny;
	unsigned nz;
	
public: 
	Cell_list();
	virtual ~Cell_list();
	
	void build(
		vector<double> x,
		vector<double> y,
		vector<double> z,
		double X,
		double Y,
		double Z,
		double r
	);
	void neighbors(unsigned particle_index, vector<unsigned>& out);
	string info();
};

#endif // CELL_LIST_H

