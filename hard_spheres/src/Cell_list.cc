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
	index2cell.clear();
	neighbor_list.clear();
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
	assert(r<0.5*X);	// TODO handle special cases where this requerement is not meet
	assert(r<0.5*Y);
	assert(r<0.5*Z);
	
	nx=floor(X/r);
	ny=floor(Y/r);
	nz=floor(Z/r);
	double lx = X/(double)nx;
	double ly = Y/(double)ny;
	double lz = Z/(double)nz;
	
	//cout << "r: " << r << " Box: " << X << Y << Z  <<  " Cell size: " << lx << " " << ly << " " << lz << endl;
	cells.clear();
	cells.resize(nx*ny*nz);
	index2cell.clear();
	
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
	
	// Construct neighbour list
	neighbor_list.clear();
	neighbor_list.resize(x.size());
	for( unsigned i = 0 ; i<x.size(); i++ ) {
		vector<unsigned> n27;
		neighbors27(i,n27);
		for ( unsigned in=0; in<n27.size() ; in++) {
			unsigned j=n27.at(in);
			double dx=x.at(j)-x.at(i);
			double dy=y.at(j)-y.at(i);
			double dz=z.at(j)-z.at(i);
			dx-=X*nearbyint(dx/X);
			dy-=Y*nearbyint(dy/Y);
			dz-=Z*nearbyint(dz/Z);
			double d2=dx*dx+dy*dy+dz*dz;
			if(d2<r*r) {
				neighbor_list.at(i).push_back(j);
			}
		}
	}
}

/**
*  Return neighbor list
*/
void Cell_list::neighbors(unsigned particle_index,vector<unsigned>& out){
	out = neighbor_list.at(particle_index);
}

/**
*  Return particle indexes of the 27 nearby images.
*/
void Cell_list::neighbors27(unsigned particle_index,vector<unsigned>& out){
	out.clear();
	unsigned index = index2cell.at(particle_index); // Convert particle index to a lattice index
	
	int ix=index%nx;
	int iy=(index/nx)%ny;
	int iz=index/(nx*ny);
	//cout << " index = " << index << " ix = " <<  ix << " iy =  " << iy << " iz " << iz << endl; 
	//cout << "   cells.size() = " << cells.size() << endl; 
	//cout << "     nx = " << nx << " ny =  " << ny << " nz " << nz << endl; 
	//cout << "       cells.at(index).size() = " << cells.at(index).size() << endl;
	
	// Loop particles
	// (Ensure that particles are not included twize for small systems)
	int min_ix=-1;
	int max_ix=1;
	int min_iy=-1;
	int max_iy=1;
	int min_iz=-1;
	int max_iz=1;
	if(nx==1){min_ix=0;max_ix=0;}
	if(ny==1){min_iy=0;max_iy=0;}
	if(nz==1){min_iz=0;max_iz=0;}
	if(nx==2){min_ix=0;max_ix=1;}
	if(ny==2){min_iy=0;max_iy=1;}
	if(nz==2){min_iz=0;max_iz=1;}

	for(int dix=min_ix;dix<max_ix+1;dix++){
		for(int diy=min_iy;diy<max_iy+1;diy++){
			for(int diz=min_iz;diz<max_iz+1;diz++){
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
	out << "    ..:: Cell list information ::.. " << endl;
	out << "index2cell.size() = " << index2cell.size() << endl;
	out << "nx = " << nx << endl;
	out << "ny = " << ny << endl;
	out << "nz = " << nz << endl;
	out << "nx*ny*nz = " << nx*ny*nz << endl;
	out << "cells.size() = " << cells.size() << endl;
	out << "cells.at(0).size() = " << cells.at(0).size() << endl;	
	
	for(unsigned i = 0;i<cells.size();i++)
		out << "cells.at( " << i << " ).size() = " << cells.at(i).size() << endl;
	
	out << "neighbor_list.at(0).size() = " << neighbor_list.at(0).size() << endl;
	
	vector<unsigned> n;
	neighbors27(0,n);
	out << "this.neighbors27(0,n): n.size() = " << n.size() << ", elements = ";
	for(unsigned i = 0;i<n.size();i++) out << " " << n.at(i);
	out << endl;
	
	n.clear();
	neighbors(0,n);
	out << "this.neighbors(0,n): n.size() = " << n.size() << ", elements = ";
	for(unsigned i = 0;i<n.size();i++) out << " " << n.at(i);
	out << endl;

	return out.str();
} 

