/*
 * Rotational_order.cc
 *
 *  Created on: Nov 13, 2016
 *      Author: Ulf R. Pedersen, www.urp.dk
 */

#include "Rotational_order.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>


#include <boost/iostreams/filtering_stream.hpp>   // Linker -lboost_iostreams
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include "Cell_list.h"


using namespace std;

/**
 * Split a string into a vector string
 */
vector<string> rotational_order_split(string str, char delimiter) {
  vector<string> output;stringstream ss(str);string substr;
  while(getline(ss, substr, delimiter)) output.push_back(substr);
  return output;
}


/**
 * Constructor
 */
Rotational_order::Rotational_order() :
		type(),
		x(),
		y(),
		z(),
		number_of_neighbors(),
		qlm(),
		qi(),
		qlmAvg(),
		qiAvg(),
		Lx(1.0),
		Ly(1.0),
		Lz(1.0),
		degree(6),
		neighbour_cutoff(1.5)
	{
	// Do nothing in constructor
}

Rotational_order::~Rotational_order() {
	type.clear();
	x.clear();
	y.clear();
	z.clear();
}


/**
 * Add a particle of in_type to the lattice
 */
void Rotational_order::addParticle(unsigned in_type,double in_x,double in_y,double in_z){
	type.push_back(in_type);
	x.push_back(in_x);
	y.push_back(in_y);
	z.push_back(in_z);
}

/**
* Compute the qlm vector for each particle
*/
void Rotational_order::compute_ql(unsigned in_degree){
	degree=in_degree;
	qlm.clear();
	qi.clear();
	qlmAvg.clear();
	qiAvg.clear();
	double pi=acos(-1.);
	
	// Build a cell list
	Cell_list cell_list;
	cell_list.build(x,y,z,Lx,Ly,Lz,neighbour_cutoff);	
	//cout << cell_list.info();
	
	// Compute complex bond vector qlm(i)=sum_j=neighbour { Ylm(r_ij) } where Ylm is a speherical harmonic
	for(unsigned i=0;i<number_of_particles();i++){
		for(int m=(-1)*(int)degree;m<(int)degree+1;m++){
			unsigned num_bonds=0;
			double this_qlm_real=0;
			double this_qlm_imag=0;
			for(unsigned j=0;j<number_of_particles();j++){
				if(i!=j){
					double dx=x.at(j)-x.at(i);
					double dy=y.at(j)-y.at(i);
					double dz=z.at(j)-z.at(i);
					dx-=Lx*nearbyint(dx/Lx);
					dy-=Ly*nearbyint(dy/Ly);
					dz-=Lz*nearbyint(dz/Lz);
					double r=sqrt(dx*dx+dy*dy+dz*dz);
					double phi=atan2(dy,dx);       // goes from 0 to 2pi
					if(phi<0) phi+=2.*pi;   // acos(-1.)=pi=3.1415...
					double theta=acos(dz/r);       // goes from 0 to pi
					if(r<neighbour_cutoff){
						num_bonds++;
	//if(m==0)
	//	cout << "i=" << i << " j=" << j << " r= " << r << "  " << "phi= " << phi << " theta= " << theta << endl; 
						using namespace boost::math;
						this_qlm_real+=spherical_harmonic_r(degree,m,theta,phi);
						this_qlm_imag+=spherical_harmonic_i(degree,m,theta,phi);
					}
				}
			}
			this_qlm_real/=(double)num_bonds;
			this_qlm_imag/=(double)num_bonds;
			qlm.push_back(this_qlm_real);
			qlm.push_back(this_qlm_imag);
			if(m==0)
				number_of_neighbors.push_back(num_bonds);
			//cout << "i= " << i << " m= " << m << " qlm_real=  " << " this_qlm_real= " << this_qlm_real << " this_qlm_imag=" << this_qlm_imag << endl;
		}
		{ // Compute the Steinhard Order parameters
			double this_qi=0;
			for(unsigned c=qlm.size()-(2*degree+1)*2;c<qlm.size();c++)
				this_qi+=qlm.at(c)*qlm.at(c);
			this_qi*=4.*pi/(2*degree+1);
			this_qi=sqrt(this_qi);
			qi.push_back(this_qi);
			//cout << " qi= " << this_qi << endl;
		}
	}// Done looping particles
	// Compute the Lechner-Dellago averaged versions of the Steinhard order-parameter*
	for(unsigned i=0;i<number_of_particles();i++){
		for(int m=(-1)*(int)degree;m<(int)degree+1;m++){
			unsigned num_bonds=0;
			double this_qlm_real=0;
			double this_qlm_imag=0;
			for(unsigned j=0;j<number_of_particles();j++){
				double dx=x.at(j)-x.at(i);
				double dy=y.at(j)-y.at(i);
				double dz=z.at(j)-z.at(i);
				dx-=Lx*nearbyint(dx/Lx);
				dy-=Ly*nearbyint(dy/Ly);
				dz-=Lz*nearbyint(dz/Lz);
				double r=sqrt(dx*dx+dy*dy+dz*dz);
				int index = (int)j*(2*(int)degree+1)*2;
				index+=2*(m+(int)degree);
				if(r<neighbour_cutoff && i!=j){
					num_bonds++;
					this_qlm_real+=qlm.at(index);
					this_qlm_imag+=qlm.at(index+1);
				}
			}
			this_qlm_real/=(double)num_bonds;
			this_qlm_imag/=(double)num_bonds;	
			qlmAvg.push_back(this_qlm_real);
			qlmAvg.push_back(this_qlm_imag);
		}
		{ // The Lechner-Dellago order-parameter
			double this_qi=0.;
			for(unsigned c=qlmAvg.size()-(2*degree+1)*2;c<qlmAvg.size();c++)
				this_qi+=qlmAvg.at(c)*qlmAvg.at(c);
			this_qi*=4.*pi/(2*degree+1);
			this_qi=sqrt(this_qi);
			qiAvg.push_back(this_qi);
		}
	}// Done looping particles
}


/**
 * Return the number of particles
 */
unsigned Rotational_order::number_of_particles(){
	return x.size();
}


/**
 * Return the number of particles types
 */
unsigned Rotational_order::number_of_types(){
	unsigned out=0;
	if( type.size()>0 )
		for(unsigned i=0;i<type.size();i++)
			if(type.at(i)>=out) out=type.at(i)+1;
	return out;
}

unsigned Rotational_order::number_of_particles_of_type(unsigned test_type){
	unsigned out=0;
	for(unsigned n=0;n<type.size();n++)
		if(type.at(n)==test_type)
			out++;
	return out;
}

/**
 * Return volume of lattice
 */
double Rotational_order::volume(){
	return Lx*Ly*Lz;
}

/*
void Rotational_order::load_xyz(ifstream& ifile){
		string line;
		getline(ifile,line);
		unsigned num_atoms=atoi(line.c_str());
		getline(ifile,line); // Comment line
		unsigned line_counter=0;
		for(unsigned i=0;i<num_atoms;i++){
				getline(ifile,line);
				char * pEnd;
				int    i0;
				double d0, d1, d2;
				i0 = strtol (line.c_str(),&pEnd,10);
				d0 = strtod (pEnd,&pEnd);
				d1 = strtod (pEnd,&pEnd);
				d2 = strtod (pEnd,&pEnd);
				type.push_back(i0);
				x.push_back(d0);
				y.push_back(d1);
				z.push_back(d2);
		}
	
	
}*/


/**
* Load coordinates of particles from xyz-file or an zipped xyz.gz file.
*/
void Rotational_order::load_xyz(string ifilename,double in_Lx,double in_Ly,double in_Lz,double in_neighbour_cutoff){
	using namespace boost::iostreams;
	
	Lx=in_Lx;
	Ly=in_Ly;
	Lz=in_Lz;
	neighbour_cutoff=in_neighbour_cutoff;

	filtering_istream in;
	vector<string> fnames = rotational_order_split(ifilename,'.');
	//in.push(file_source(ifilename));
	if(fnames.size()<2){
		cerr << "error: Incompatible name of input file. Should be an *.xyz or *.xyz.gz file." << endl;
		abort();
	}
	if(fnames.back()=="xyz"){
		in.push(file_source(ifilename));
	} else if ( fnames.back()=="gz" && fnames.at(fnames.size()-2)=="xyz" ) {
		in.push(gzip_decompressor());
		in.push(file_source(ifilename));
	} else {
		cerr << "error: Incompatible name of output file. Should be an *.xyz or *.xyz.gz file." << endl;
		abort();
	}
	
	string line;
	getline(in,line);
	unsigned num_atoms=atoi(line.c_str());
	getline(in,line); // Comment line
	unsigned line_counter=0;
	for(unsigned i=0;i<num_atoms;i++){
		getline(in,line);
		char * pEnd;
		int    i0;
		double d0, d1, d2;
		i0 = strtol (line.c_str(),&pEnd,10);
		d0 = strtod (pEnd,&pEnd);
		d1 = strtod (pEnd,&pEnd);
		d2 = strtod (pEnd,&pEnd);
		type.push_back(i0);
		x.push_back(d0);
		y.push_back(d1);
		z.push_back(d2);
	}
	
	
	//cout << line;
	
	//abort();
	/*
	//filtering_istream in;
	//in.push(file_sink(ifilename));
	
	ifstream ifile;
	ifile.open(ifilename.c_str());
	filtering_streambuf<input> in;
	
	
	if(ifile.is_open()){
		load_xyz(ifile);
	}else{
		cerr << "error: Could not load input file name " << ifilename << endl;
		abort();
	}*/
}


/**
 * Write coordinates of particles to xyz-file.
 */
//void Rotational_order::write_xyz(ostream& out){}
void Rotational_order::write_xyz(string ofilename){
	using namespace boost::iostreams;
	vector<string> fnames = rotational_order_split(ofilename,'.');
	filtering_ostream out;
	if(fnames.size()<2){
		cerr << "error: Incompatible name of output file. Should be *.xyz or *.xyz.gz" << endl;
		abort();
	}
	if(fnames.back()=="xyz"){
		out.push(file_sink(ofilename));
	} else if ( fnames.back()=="gz" && fnames.at(fnames.size()-2)=="xyz" )  {
		out.push(gzip_compressor());
		out.push(file_sink(ofilename));
	} else {
		cerr << "error: Incompatible name of output file. Should be *.xyz or *.xyz.gz" << endl;
		abort();
	}
	out << this->number_of_particles()<<endl;
	out << " numTypes=" << number_of_types();
	out << " sim_box=RectangularSimulationBox," << Lx << "," << Ly << "," << Lz;
	out << " columns=type,x,y,z,Ni,ql,qlAvg " << endl;
	for (unsigned i=0;i<number_of_particles();i++) {
		out << type.at(i); 
		out << " " << x.at(i)  << " " << y.at(i) << " " << z.at(i);
		out << " " << number_of_neighbors.at(i);
		out << " " << qi.at(i) << " " << qiAvg.at(i);
		out << endl;
	}
}



/**
 * Return a string with various information about the lattice
 */
string Rotational_order::info(){
	double avg;
	
	stringstream out;
	out << "Total number of particles:    " << number_of_particles() << endl;
	out << "Number of types:              " << number_of_types() << endl;
	out << "Number of particles of types:";
	for(unsigned i=0;i<number_of_types();i++)
		out << " " << number_of_particles_of_type(i);
	out << endl;
	out << "Lengths of box vectors:       " << Lx << " " << Ly << " " << Lz << endl;
	out << "Box volume:                   " << volume() << endl;
	out << "Number density:               " << number_of_particles()/volume() << endl;
	out << "degree:                    l= " << degree << endl;
	out << "rcut:                     rc= " << neighbour_cutoff << endl;
	avg = 0;for(unsigned i=0;i<number_of_particles();i++) avg += (double)number_of_neighbors.at(i);avg/=(double)number_of_particles();
	out << "Average number of neighbors:  " << avg << endl;
	avg = 0;for(unsigned i=0;i<number_of_particles();i++) avg += qi.at(i);avg/=(double)number_of_particles();
	out << "Average qi:               Ql= " << avg << endl;
	avg = 0;for(unsigned i=0;i<number_of_particles();i++) avg += qiAvg.at(i);avg/=(double)number_of_particles();
	out << "Average qiAvg:         QlAvg= " << avg << endl;
	out << "           DEBUG INFO" << endl;
	out << " qi.size()= " << qi.size() << " qiAvg.size()= " << qiAvg.size() <<endl;
	out << " qlm.size()= " << qlm.size() << " qlmAvg.size()= " << qlmAvg.size() <<endl;
	out << " number_of_neighbors.size()= " << number_of_neighbors.size() <<endl;

	return out.str();
}
