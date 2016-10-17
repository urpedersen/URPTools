/*
 * Cell_list.cc
 *
 *  Created on: Nov 16, 2016
 *      Author: Ulf R. Pedersen, www.urp.dk
 */

#include "Cell_list.h"

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

/**
 * Constructor
 */
Cell_list::Cell_list() :
	cells(),
	index2cell(),
	nx(0),
	ny(0),
	nz(0)
	{
	// Do nothing in constructor
}

/**
* Destructor
*/
Cell_list::~Cell_list() {
	cells.clear();
}


/**
*  Build the cell list of a periodic box of size (X,Y,Z)
*  r gives a minimum lengths of the cell
*/
void Cell_list::build(
	vector<double> x,
	vector<double> y,
	vector<double> z,
	double X,
	double Y,
	double Z,
	double r
){
	assert(x.size()==y.size());
	assert(x.size()==z.size());
	
	nx=floor(X/r);
	ny=floor(Y/r);
	nz=floor(Z/r);
	double lx = X/(double)nx;
	double ly = Y/(double)ny;
	double lz = Z/(double)nz;
	
	//cout << "r: " << r << " Box: " << X << Y << Z  <<  " Cell size: " << lx << " " << ly << " " << lz << endl;
	
	cells.resize(nx*ny*nz);
	
	for ( unsigned i = 0 ; i<x.size(); i++) {
		double xx=x.at(i)-X*floor(x.at(i)/X);
		double yy=y.at(i)-Y*floor(y.at(i)/Y);
		double zz=z.at(i)-Z*floor(z.at(i)/Z);
		unsigned ix = floor(xx/lx);
		unsigned iy = floor(yy/ly);
		unsigned iz = floor(zz/lz);
		unsigned index = iz*ny*nx+iy*nx+ix;
		cells.at(index).push_back(i);
		index2cell.push_back(index);
	}
}

/**
*  Return particle indexes of the 27 nearby images.
*/
vector<unsigned> Cell_list::neighbors(unsigned index){
	vector<unsigned> out;
	// TODO look in 27 neighbors cells
	return out;
}


/**
*  Print information about the cell list
*/
string Cell_list::info(){
	stringstream out;
	out << "    ..:: Cell list information ::.. ";
	out << "index2cell.size() = " << index2cell.size() << endl;
	out << "nx = " << nx << endl;
	out << "ny = " << ny << endl;
	out << "nz = " << nz << endl;
	out << "nx*ny*nz = " << nx*ny*nz << endl;
	out << "cells.size() = " << cells.size() << endl;
	out << "cells.at(0).size() = " << cells.at(0).size() << endl;

	for(unsigned iz=0;iz<nz;iz++) 
		for(unsigned iy=0;iy<ny;iy++)
			for(unsigned ix=0;ix<nx;ix++) {
				unsigned index = iz*ny*nx+iy*nx+ix;
				out << "Size at cell " << index << " : " << cells.at(index).size() << endl;
	}
				
			
	return out.str();
} 

