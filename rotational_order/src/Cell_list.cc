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
	const vector<double>& x,
	const vector<double>& y,
	const vector<double>& z,
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
void Cell_list::neighbors(unsigned particle_index,vector<unsigned>& out){
	out.clear();
	unsigned index = index2cell.at(particle_index); // Convert particle index to a lattice index
	
	int ix=index%nx;
	int iy=(index/nx)%ny;
	int iz=index/(nx*ny);
	//cout << " index = " << index << " ix = " <<  ix << " iy =  " << iy << " iz " << iz << endl; 
	//cout << "   cells.size() = " << cells.size() << endl; 
	//cout << "     nx = " << nx << " ny =  " << ny << " nz " << nz << endl; 
	//cout << "       cells.at(index).size() = " << cells.at(index).size() << endl;
		
	for(int dix=-1;dix<2;dix++){
		for(int diy=-1;diy<2;diy++){
			for(int diz=-1;diz<2;diz++){
				unsigned iix=(ix+dix+nx)%nx;
				unsigned iiy=(iy+diy+ny)%ny;
				unsigned iiz=(iz+diz+nz)%nz;
				//cout << " (" << iix  << "," << iiy << "," << iiz << ") ";
				unsigned nindex = iiz*ny*nx+iiy*nx+iix;
				//cout << "nindex: " << nindex << endl;
				out.insert(out.end(),cells.at(nindex).begin(),cells.at(nindex).end() );
				/*for(unsigned i = 0; i < cells.at(nindex).size(); i++){
					out.push_back( cells.at(nindex).at(i) );
				}*/
			}
		}
	}	
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
	//out << "cells.at(0).size() = " << cells.at(0).size() << endl;
	
	for(unsigned i = 0;i<cells.size();i++)
		out << "cells.at( " << i << " ).size() = " << cells.at(i).size() << endl;
	

	return out.str();
} 

