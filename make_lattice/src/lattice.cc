/*
 * Lattice.cc
 *
 *  Created on: Nov 9, 2015
 *      Author: urp
 */

#include "lattice.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>

using namespace std;

/**
 * Return random component to velocity vector TODO See if this is the right distribution
 */
double random_velocity(double temperature,double mass){
	double tmp = 0.0;
	for(unsigned i=0;i<12;i++) tmp+=((double)rand()/(double)RAND_MAX-0.5); // Normal distribution with sigma=1
	return tmp*sqrt(temperature/mass);
}

/**
 *
 * Useful link:
 * http://www.geocities.jp/ohba_lab_ob_page/Structure.html
 */
Lattice::Lattice(string in_lattice_type,unsigned nx, unsigned ny,unsigned nz,unsigned seed) :
		lattice_type(in_lattice_type),
		number_of_sites_on_lattice(0),
		type(),
		x(),
		y(),
		z(),
		Lx(nx),
		Ly(ny),
		Lz(nz),
		mass_of_types()
	{
	srand(seed);
	if(lattice_type=="sc"){			// Simple cubic
		for(unsigned ix=0;ix<nx;ix++){
			for(unsigned iy=0;iy<ny;iy++){
				for(unsigned iz=0;iz<nz;iz++){
					addParticle(0,ix+0.5,iy+0.5,iz+0.5);
				}
			}
		}
	}
	else if(lattice_type=="rp"){	// Randomly packed
		for(unsigned ix=0;ix<nx;ix++){
			for(unsigned iy=0;iy<ny;iy++){
				for(unsigned iz=0;iz<nz;iz++){
					addParticle(0,Lx*(double)rand()/(double)RAND_MAX,Lx*(double)rand()/(double)RAND_MAX,Lx*(double)rand()/(double)RAND_MAX);
					//addParticle(0,ix,iy,iz);
				}
			}
		}
		// TODO Below MC simulation of hard spheres scales badly (N**2). Can be fixed with a cell list.
		for(double size=0.0;size<1.1;size+=0.001){
			for(unsigned n=0;n<x.size();n++){
				double x_old=x.at(n);
				double y_old=y.at(n);
				double z_old=z.at(n);
				double x_move=0.1*((double)rand()/(double)RAND_MAX-0.5);
				double y_move=0.1*((double)rand()/(double)RAND_MAX-0.5);
				double z_move=0.1*((double)rand()/(double)RAND_MAX-0.5);
				x.at(n)+=x_move;
				y.at(n)+=y_move;
				z.at(n)+=z_move;
				x.at(n)-=floor(x.at(n)/Lx)*Lx;
				y.at(n)-=floor(y.at(n)/Ly)*Ly;
				z.at(n)-=floor(z.at(n)/Lz)*Lz;
				bool overlap = false;
				for(unsigned m=0;m<x.size();m++){
					if(n!=m){
						double xx=x.at(n)-x.at(m);
						xx-=round(xx/Lx)*Lx;xx*=xx;
						double yy=y.at(n)-y.at(m);
						yy-=round(yy/Ly)*Ly;yy*=yy;
						double zz=z.at(n)-z.at(m);
						zz-=round(zz/Lz)*Lz;zz*=zz;
						if(xx+yy+zz<size*size)
							overlap = true;
					}
				}
				if(overlap){
					x.at(n)=x_old;
					y.at(n)=y_old;
					z.at(n)=z_old;
				}
			}
		}
	}
	else if(lattice_type=="bcc"){	// Body centre cubic
		for(unsigned ix=0;ix<nx;ix++){
			for(unsigned iy=0;iy<ny;iy++){
				for(unsigned iz=0;iz<nz;iz++){
					addParticle(0,ix+0.0,iy+0.0,iz+0.0);
					addParticle(0,ix+0.5,iy+0.5,iz+0.5);
				}
			}
		}
	} else if (lattice_type=="fcc") {	// Face centre cubic
		for(unsigned ix=0;ix<nx;ix++){
			for(unsigned iy=0;iy<ny;iy++){
				for(unsigned iz=0;iz<nz;iz++){
					//addParticle(0,ix+0.0,iy+0.0,iz+0.0);
					//addParticle(0,ix+0.5,iy+0.5,iz+0.0);
					//addParticle(0,ix+0.5,iy+0.0,iz+0.5);
					//addParticle(0,ix+0.0,iy+0.5,iz+0.5);
					addParticle(0,ix+0.25,iy+0.25,iz+0.25);
					addParticle(0,ix+0.75,iy+0.75,iz+0.25);
					addParticle(0,ix+0.75,iy+0.25,iz+0.75);
					addParticle(0,ix+0.25,iy+0.75,iz+0.75);
				}
			}
		}
	} else if (lattice_type=="hcp") {	// Hexagonal close packed
		for(unsigned ix=0;ix<nx;ix++){
			for(unsigned iy=0;iy<ny;iy++){
				for(unsigned iz=0;iz<nz;iz++){
					addParticle(0,ix+0.0,iy+0.0,iz+0.0);
					addParticle(0,ix+0.5,iy+0.5,iz+0.0);
					addParticle(0,ix+0.5,iy+5.0/6.0,iz+0.5);
					addParticle(0,ix+0.0,iy+1.0/3.0,iz+0.5);
				}
			}
		}
		scale_y_coordinates(ny*sqrt(3.0));
		scale_z_coordinates(nz*sqrt(8.0/3.0));
	} else if (lattice_type=="hex") {	// Hexagonal layers
			for(unsigned ix=0;ix<nx;ix++) {
				for(unsigned iy=0;iy<ny;iy++) {
					for(unsigned iz=0;iz<nz;iz++) {
						addParticle(0,ix+0.0,iy+0.0,iz+0.0);
						addParticle(0,ix+0.5,iy+0.5,iz+0.0);
					}
				}
			}
			scale_y_coordinates(ny*sqrt(3.0));
	} else if (lattice_type=="dc") {	// Diamond cubic lattice
			for(unsigned ix=0;ix<nx;ix++){
				for(unsigned iy=0;iy<ny;iy++){
					for(unsigned iz=0;iz<nz;iz++){
						addParticle(0,ix+0.00,iy+0.00,iz+0.00);
						addParticle(0,ix+0.00,iy+0.50,iz+0.50);
						addParticle(0,ix+0.50,iy+0.00,iz+0.50);
						addParticle(0,ix+0.50,iy+0.50,iz+0.00);
						addParticle(0,ix+0.75,iy+0.75,iz+0.75);
						addParticle(0,ix+0.75,iy+0.25,iz+0.25);
						addParticle(0,ix+0.25,iy+0.75,iz+0.25);
						addParticle(0,ix+0.25,iy+0.25,iz+0.75);
					}
				}
			}
	} else if (lattice_type=="NaCl") {	// Rock salt lattice
		for(unsigned ix=0;ix<nx;ix++){
			for(unsigned iy=0;iy<ny;iy++){
				for(unsigned iz=0;iz<nz;iz++){
					addParticle(0,ix+0.0,iy+0.0,iz+0.0);
					addParticle(0,ix+0.5,iy+0.5,iz+0.0);
					addParticle(0,ix+0.5,iy+0.0,iz+0.5);
					addParticle(0,ix+0.0,iy+0.5,iz+0.5);
					addParticle(1,ix+0.5,iy+0.5,iz+0.5);
					addParticle(1,ix+0.0,iy+0.0,iz+0.5);
					addParticle(1,ix+0.0,iy+0.5,iz+0.0);
					addParticle(1,ix+0.5,iy+0.0,iz+0.0);
				}
			}
		}
	} else if (lattice_type=="CsCl") {	// Cesium Chloride lattice (AuCd or TiNi)
		for(unsigned ix=0;ix<nx;ix++){
			for(unsigned iy=0;iy<ny;iy++){
				for(unsigned iz=0;iz<nz;iz++){
					addParticle(0,ix+0.00,iy+0.00,iz+0.00);
					addParticle(1,ix+0.50,iy+0.50,iz+0.50);
				}
			}
		}
	} else if (lattice_type=="CuZr2") {
		// Mail from Toby Hudson to Ulf R. Pedersen:
		// CuZr2 structure type optimized at rA/rB = 6.0/5.0 (NB - rB<rA) Lx=1.134695 Ly=1.134695 Lz=5.04776
		// A 0.000000 0.000000 0.000000
		// A 0.567348 0.567348 2.523880
		// B 0.000000 0.000000 1.738161
		// B 0.000000 0.000000 3.309599
		// B 0.567348 0.567348 0.785719
		// B 0.567348 0.567348 4.262041
		//    MoSi2, C11b type structure. Pearson symbol is tI6. Space group is I4/mmm.
		double sx=1.134695;double sy=1.134695;double sz=5.04776;
		for(unsigned ix=0;ix<nx;ix++){
			for(unsigned iy=0;iy<ny;iy++){
				for(unsigned iz=0;iz<nz;iz++){
					addParticle(0,ix+0.000000/sx,iy+0.000000/sy,iz+0.000000/sz);
					addParticle(0,ix+0.567348/sx,iy+0.567348/sy,iz+2.523880/sz);
					addParticle(1,ix+0.000000/sx,iy+0.000000/sy,iz+1.738161/sz);
					addParticle(1,ix+0.000000/sx,iy+0.000000/sy,iz+3.309599/sz);
					addParticle(1,ix+0.567348/sx,iy+0.567348/sy,iz+0.785719/sz);
					addParticle(1,ix+0.567348/sx,iy+0.567348/sy,iz+4.262041/sz);
				}
			}
		}
		scale_x_coordinates(ny*sx);scale_z_coordinates(ny*sy);scale_z_coordinates(nz*sz);
	}
	else {
		cerr << "error: unknown lattice type " << lattice_type << "." << endl;
		abort();
	}
	for(unsigned i=0;i<number_of_types();i++)
		mass_of_types.push_back(1.0);
	// translate lattice to origin at center of box
	translate_all_particles(-0.5*Lx,-0.5*Ly,-0.5*Lz);
	number_of_sites_on_lattice=number_of_particles();
}

Lattice::~Lattice() {
	type.clear();
	mass_of_types.clear();
	x.clear();
	y.clear();
	z.clear();
}

/**
 * Add a particle of in_type to the lattice
 */
void Lattice::addParticle(unsigned in_type,double in_x,double in_y,double in_z){
	type.push_back(in_type);
	x.push_back(in_x);
	y.push_back(in_y);
	z.push_back(in_z);
}

/**
 * Randomly swop position of particles
 */
void Lattice::mix_positions(){
	for(unsigned i=0;i<x.size();i++){
		double this_x=x.at(i);
		double this_y=y.at(i);
		double this_z=z.at(i);
		unsigned j=rand()*x.size()/RAND_MAX;
		x.at(i)=x.at(j);
		y.at(i)=y.at(j);
		z.at(i)=z.at(j);
		x.at(j)=this_x;
		y.at(j)=this_y;
		z.at(j)=this_z;
	}
}

/**
 * Reset the density of the lattice
 */
void Lattice::set_density(double new_density){
	double sf=pow( (number_of_particles()/volume())/new_density , 1.0/3.0 );
	scale_coordinates(sf);
}

void Lattice::scale_coordinates(double sf){
	Lx=sf*Lx;
	Ly=sf*Ly;
	Lz=sf*Lz;
	for (unsigned i=0;i<number_of_particles();i++) {
		x.at(i)=sf*x.at(i);
		y.at(i)=sf*y.at(i);
		z.at(i)=sf*z.at(i);
	}
}

void Lattice::scale_x_coordinates(double new_Lx){
	double sf=new_Lx/Lx;
	Lx=sf*Lx;
	for (unsigned i=0;i<number_of_particles();i++) {
		x.at(i)=sf*x.at(i);
	}
}

void Lattice::scale_y_coordinates(double new_Ly){
	double sf=new_Ly/Ly;
	Ly=sf*Ly;
	for (unsigned i=0;i<number_of_particles();i++) {
		y.at(i)=sf*y.at(i);
	}
}

void Lattice::scale_z_coordinates(double new_Lz){
	double sf=new_Lz/Lz;
	Lz=sf*Lz;
	for (unsigned i=0;i<number_of_particles();i++) {
		z.at(i)=sf*z.at(i);
	}
}

void Lattice::translate_all_particles(double dx,double dy,double dz){
	for(unsigned n=0;n<x.size();n++){
		x.at(n)+=dx;
		y.at(n)+=dy;
		z.at(n)+=dz;
	}
}

void Lattice::reset_particle_types(vector<unsigned> num_par){
	mass_of_types.clear();
	// Rename particle types
	unsigned c=0;
	for(unsigned i=0;i<num_par.size();i++){
		mass_of_types.push_back(1.0);
		for(unsigned n=0;n<num_par.at(i);n++)
			if(reset_particle_type(c,i))
				c++;
			else {
				cerr<<"error: Could not reset particle types."<<endl;
				abort();
			}
	}
	// Remove extra particles
	reset_number_of_particles(c);
}

void Lattice::reset_mass_of_types(vector<double> mass){
	mass_of_types.clear();
	for(unsigned i=0;i<mass.size();i++)
		mass_of_types.push_back(mass.at(i));
	if(number_of_types()!=mass_of_types.size()) {
		cerr << "error: number of masses given is not the same as number of types." << endl;
		abort();
	}
}

/**
 * Reset the type of a particle. Returns false if unsuccessful.
 */
bool Lattice::reset_particle_type(unsigned i,unsigned new_type){
	bool output=false;
	if(i<type.size()){
		type.at(i)=new_type;
		output=true;
	}
	return output;
}


void Lattice::reset_number_of_particles(unsigned num_par) {
	while(x.size()>num_par){
		type.pop_back();
		x.pop_back();
		y.pop_back();
		z.pop_back();
	}
}

/**
 * Return the minimum distance found in lattice
 * TODO The get_min_distance() function scales badly, and makes the program slow for large systems.
 */
double Lattice::get_min_distance(){
	cout << "get_min_distance()" << endl;
	double min_distance=Lx+Lz+Lz;
	for (unsigned i=0;i<number_of_particles()-1;i++) {
		for (unsigned j=i+1;j<number_of_particles();j++) {
			double xx = x.at(i)-x.at(j);
			xx-=round(xx/Lx)*Lx;
			xx*=xx;
			double yy = y.at(i)-y.at(j);
			yy-=round(yy/Ly)*Ly;
			yy*=yy;
			double zz = z.at(i)-z.at(j);
			zz-=round(zz/Lz)*Lz;
			zz*=zz;
			if(xx+yy+zz<min_distance*min_distance){
				min_distance = sqrt(xx+yy+zz);
			}
		}
	}
	return min_distance;
}

/**
 * Scale coordinates so that the minimum distance is given by new_min_distance
 */
void Lattice::set_min_distance(double new_min_distance){
	scale_coordinates(new_min_distance/get_min_distance());
}

/**
 * Return the number of particles
 */
unsigned Lattice::number_of_particles(){
	return x.size();
}


/**
 * Return the number of particles types
 */
unsigned Lattice::number_of_types(){
	unsigned out=0;
	if( type.size()>0 )
		for(unsigned i=0;i<type.size();i++)
			if(type.at(i)>=out) out=type.at(i)+1;
	return out;
}

unsigned Lattice::number_of_particles_of_type(unsigned test_type){
	unsigned out=0;
	for(unsigned n=0;n<type.size();n++)
		if(type.at(n)==test_type)
			out++;
	return out;
}

/**
 * Return volume of lattice
 */
double Lattice::volume(){
	return Lx*Ly*Lz;
}


/**
 * Write coordinates of particles to xyz-file.
 */
void Lattice::write_xyz(ostream& out,double temperature){

	//stringstream out;
	out << this->number_of_particles()<<endl;
	out << "ioformat=2 numTypes=" << number_of_types();
	out << setprecision(16);
	out << " sim_box=RectangularSimulationBox," << Lx << "," << Ly << "," << Lz;
	out << " mass=";
	for(unsigned i=0;i<mass_of_types.size();i++){
		if(i<mass_of_types.size()-1)
			out << mass_of_types.at(i) << ",";
		else
			out << mass_of_types.at(i);
	}
	out << " columns=type,x,y,z,imx,imy,imz,vx,vy,vz " << endl;
	for (unsigned i=0;i<number_of_particles();i++) {
		out << type.at(i) << " " << x.at(i) << " " << y.at(i) << " " << z.at(i);
		out << " 0 0 0";
		// for(unsigned i=0;i<3;i++)
		if(Lx>0) {
			out << " " << random_velocity(temperature,mass_of_types.at(type.at(i)));
		}else{
			out << " 0.0";
		}
		if(Ly>0) {
			out << " " << random_velocity(temperature,mass_of_types.at(type.at(i)));
		}else{
			out << " 0.0";
		}
		if(Lz>0) {
			out << " " << random_velocity(temperature,mass_of_types.at(type.at(i)));
		}else{
			out << " 0.0";
		}
		out << endl;
	}
}

/**
 * Return a string with various information about the lattice
 */
string Lattice::info(){
	stringstream out;
	out << "Lattice type:                 " << lattice_type << endl;
	out << "Number of lattice sites:      " << number_of_sites_on_lattice << endl;
	out << "Total number of particles:    " << number_of_particles() << endl;
	out << "Number of types:              " << number_of_types() << endl;
	out << "Number of particles of types:";
	for(unsigned i=0;i<number_of_types();i++)
		out << " " << number_of_particles_of_type(i);
	out << endl;
	out << "Mass of types:               ";
	for(unsigned i=0;i<mass_of_types.size();i++)
		out << " " << mass_of_types.at(i);
	out << endl;
	out << "Lengths of box vectors:       " << Lx << " " << Ly << " " << Lz << endl;
	out << "Box volume:                   " << volume() << endl;
	out << "Number density:               " << number_of_particles()/volume() << endl;
	// out << "Shortest distance in lattice: " << get_min_distance() << endl;  // TODO A call to get_min_distance() is removed since it scales badly
	return out.str();
}
