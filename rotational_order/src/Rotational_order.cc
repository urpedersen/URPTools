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

#include "split.h"
#include "Cell_list.h"
#include "cluster_analysis.h"

using namespace std;

/**
 * Split a string into a vector string
vector<string> rotational_order_split(string str, char delimiter) {
  vector<string> output;stringstream ss(str);string substr;
  while(getline(ss, substr, delimiter)) output.push_back(substr);
  return output;
}
 */

/**
 * Constructor
 */
Rotational_order::Rotational_order() :
		type(),
		x(),
		y(),
		z(),
		number_of_neighbors(),
		number_of_S_connections(),
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
			vector<unsigned> n;
			cell_list.neighbors(i,n);
			//for(unsigned j=0;j<number_of_particles();j++){ // N^2 loop
			for(unsigned in=0;in<n.size();in++){
				unsigned j=n.at(in);
				if(i!=j){
					double dx=x.at(j)-x.at(i);
					double dy=y.at(j)-y.at(i);
					double dz=z.at(j)-z.at(i);
					dx-=Lx*nearbyint(dx/Lx);
					dy-=Ly*nearbyint(dy/Ly);
					dz-=Lz*nearbyint(dz/Lz);
					double r=sqrt(dx*dx+dy*dy+dz*dz);
					if(r<neighbour_cutoff){
						double phi=atan2(dy,dx);        // goes from 0 to 2pi
							if(phi<0) phi+=2.*pi;
						double theta=acos(dz/r);        // goes from 0 to pi
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
			vector<unsigned> n;
			cell_list.neighbors(i,n);
			//for(unsigned j=0;j<number_of_particles();j++){ // N^2 loop
			for(unsigned in=0;in<n.size();in++){
				unsigned j=n.at(in);
				double dx=x.at(j)-x.at(i);
				double dy=y.at(j)-y.at(i);
				double dz=z.at(j)-z.at(i);
				dx-=Lx*nearbyint(dx/Lx);
				dy-=Ly*nearbyint(dy/Ly);
				dz-=Lz*nearbyint(dz/Lz);
				double r=sqrt(dx*dx+dy*dy+dz*dz);
				unsigned index = (int)j*(2*(int)degree+1)*2;
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

void Rotational_order::compute_Sij(double S_min,string filename,string filename_cluster){
	
	ofstream out;
	out.open(filename.c_str());
	if(!out.is_open()){
		cout << "error: could not open " << filename << " for writing Sij matrix";
		abort();
	}
	out << number_of_particles() << endl;
	out << "List of connected particles using Sij > " << S_min << endl;
	
	// Build a cell list
	Cell_list cell_list;
	cell_list.build(x,y,z,Lx,Ly,Lz,neighbour_cutoff);	
	
	for(unsigned i=0;i<number_of_particles();i++){
		unsigned num_connections=0;
		vector<unsigned> n;
		cell_list.neighbors(i,n);
		//for(unsigned j=0;j<number_of_particles();j++){ // N^2 loop
		for(unsigned in=0;in<n.size();in++){
			unsigned j=n.at(in);
			double dx=x.at(j)-x.at(i);
			double dy=y.at(j)-y.at(i);
			double dz=z.at(j)-z.at(i);
			dx-=Lx*nearbyint(dx/Lx);
			dy-=Ly*nearbyint(dy/Ly);
			dz-=Lz*nearbyint(dz/Lz);
			double r=sqrt(dx*dx+dy*dy+dz*dz);
			if(i!=j && r<neighbour_cutoff){
				double Sij_real=0;
				double Sij_imag=0;
				for(int m=(-1)*(int)degree;m<(int)degree+1;m++){
					unsigned index = (int)i*(2*(int)degree+1)*2;
					index+=2*(m+(int)degree);
					unsigned jndex = (int)j*(2*(int)degree+1)*2;
					jndex+=2*(m+(int)degree);
					/*double qlm_i_real=qlm.at(index);
					double qlm_i_imag=qlm.at(index+1);
					double qlm_j_real=qlm.at(jndex);
					double qlm_j_imag=qlm.at(jndex+1);*/
					double qlm_i_real=qlmAvg.at(index);
					double qlm_i_imag=qlmAvg.at(index+1);
					double qlm_j_real=qlmAvg.at(jndex);
					double qlm_j_imag=qlmAvg.at(jndex+1);
					Sij_real+=qlm_i_real*qlm_j_real+qlm_i_imag*qlm_j_imag;
					Sij_imag+=-qlm_i_real*qlm_j_imag+qlm_i_imag*qlm_j_real;
				}
				if(Sij_real>S_min){
					//cout << "S( " << i << " , " << j << " ) = " << Sij_real << " + i " << Sij_imag << endl;
					out << i << " " << j << " # S = " << Sij_real << "  r = " << r << endl;
					num_connections++;
				}
			}
		}
		number_of_S_connections.push_back(num_connections);
	}
	out.close();

	// Perform cluster analysis
	Cluster_analysis clu;
	clu.verbose=6;
	clu.load_data_from_file(filename.c_str());
	clu.assign_nodes_to_clusters();
	//int clu_index = clu.get_index_of_largest_cluster();
	//cout << " clu_index = " << clu_index << " mass = " <<  clu.get_mass_of_largest_cluster();
	vector<unsigned> largest;
	clu.get_nodes_in_cluster(clu.get_index_of_largest_cluster(),largest);
	
	ofstream out_cluster;
	out_cluster.open(filename_cluster.c_str());
	if(!out_cluster.is_open()){
		cout << "error: could not open " << filename_cluster << " for writing coordinates of the largest cluster";
		abort();
	}
	out_cluster << largest.size() << endl;
	out_cluster << " Largest cluster of Sij connections ";
	out_cluster << " numTypes=" << number_of_types();
	out_cluster << " sim_box=RectangularSimulationBox," << Lx << "," << Ly << "," << Lz;
	out_cluster << " columns=type,x,y,z" << endl;
	for(unsigned i=0;i<largest.size();i++){
		unsigned p = largest.at(i);
		out_cluster << type.at(p) << " " << x.at(p) << " " << y.at(p) << " " << z.at(p) << endl;
	}
	out_cluster.close();
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
	
	
}
*/


/**
* Load coordinates of particles from xyz-file or an zipped xyz.gz file or a LAMMPS *.atom file
*/
void Rotational_order::load_xyz(string ifilename,unsigned frame,double in_Lx,double in_Ly,double in_Lz,double in_neighbour_cutoff){
	using namespace boost::iostreams;
	
	string fileformat="";

	Lx=in_Lx;
	Ly=in_Ly;
	Lz=in_Lz;
	neighbour_cutoff=in_neighbour_cutoff;

	filtering_istream in;
	vector<string> fnames = split(ifilename,'.');
	//in.push(file_source(ifilename));
	if(fnames.size()<2){
		cerr << "error: Incompatible name of input file. Should be an *.xyz or *.xyz.gz file." << endl;
		abort();
	}
	if(fnames.back()=="atom"){
		in.push(file_source(ifilename));
		fileformat="lammps";
	}
	else if(fnames.back()=="xyz"){
		in.push(file_source(ifilename));
		fileformat="xyz";
	} else if ( fnames.back()=="gz" && (fnames.at(fnames.size()-2)=="xyz" ) ) {
		in.push(gzip_decompressor());
		in.push(file_source(ifilename));
		fileformat="xyz";
	} else {
		cerr << "error: Incompatible name of input file. Should be an *.xyz, *.xyz.gz or *.atom file." << endl;
		abort();
	}

	
	string line;
	if(fileformat=="xyz"){ // READ xyz file 
		if(!in.good()){cerr << "error: while reading " << ifilename << endl;abort();}
		// Begin to read header of selected frame
		getline(in,line); // Number of atoms
		unsigned num_atoms=atoi(line.c_str());
		getline(in,line); // Comment line
		
		// Skip frames
		for(unsigned f=0;f<frame;f++){
			for(unsigned i = 0;i<num_atoms;i++){
				if(!in.good()){cerr << "error: while reading " << ifilename << endl;abort();}
				getline(in,line);
			}
			if(!in.good()){cerr << "error: while reading " << ifilename << endl;abort();}
			getline(in,line); // Number of atoms
			num_atoms=atoi(line.c_str());
			getline(in,line); // Comment line
		}
		
		// Attempt to find box vectors in header
		vector<string> sections = split(line,' ');
		for(unsigned i = 0 ; i < sections.size() ; i++){
			vector<string> elements = split(sections.at(i),'=');
			if(elements.size()>0 && elements.at(0)=="sim_box" && elements.size()==2){
				vector<string> vars = split(elements.at(1),',');
				if(elements.size()>0 && vars.at(0)=="RectangularSimulationBox" && vars.size()==4){
					Lx = atof(vars.at(1).c_str()); 
					Ly = atof(vars.at(2).c_str()); 
					Lz = atof(vars.at(3).c_str()); 
				}
			}
		}
	
		// Read atom positions in xyz file
		for(unsigned i=0;i<num_atoms;i++){
			if(!in.good()){cerr << "error: while reading " << ifilename << endl;abort();}
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
	} else if( fileformat=="lammps" ) {
	
		unsigned current_frame = -1;
		unsigned num_atoms=0;
		bool not_done=true;
		while(not_done){
			if(!in.good()){cerr << "error: while reading " << ifilename << endl;abort();}
			getline(in,line);
			vector<string> sections = split(line,' ');
			if(sections.size()>3 && sections.at(0)=="ITEM:" && sections.at(1)=="NUMBER" && sections.at(2)=="OF" && sections.at(3)=="ATOMS"){
				current_frame++;
				if(!in.good()){cerr << "error: while reading " << ifilename << endl;abort();}
				getline(in,line);
				num_atoms=atoi(line.c_str());
				//cout << " num_atoms = " << num_atoms << " current_frame = " << current_frame << endl;
				if(frame==current_frame){
					type.resize(num_atoms,0);
					x.resize(num_atoms,0.0);
					y.resize(num_atoms,0.0);
					z.resize(num_atoms,0.0);
				}
			}
			else if(sections.size()>2 && sections.at(0)=="ITEM:" && sections.at(1)=="BOX" && sections.at(2)=="BOUNDS" ){
				getline(in,line);
				vector<string> vars=split(line,' ');
				if(vars.size()>1)
					Lx=atof(vars.at(1).c_str())-atof(vars.at(0).c_str());
				getline(in,line);
				vars=split(line,' ');
				if(vars.size()>1)
					Ly=atof(vars.at(1).c_str())-atof(vars.at(0).c_str());
				getline(in,line);
				vars=split(line,' ');
				if(vars.size()>1)
					Lz=atof(vars.at(1).c_str())-atof(vars.at(0).c_str());
				//cout << " Lx = " << Lx << " Ly = " << Ly << " Lz = " << Lz << endl;
			}else if(sections.size()>1 && sections.at(0)=="ITEM:" && sections.at(1)=="ATOMS" && current_frame==frame){
				if(!in.good()){cerr << "error: while reading " << ifilename << endl;abort();}
				for(unsigned n=0;n<num_atoms;n++){
					if(!in.good()){cerr << "error: while reading " << ifilename << endl;abort();}
					getline(in,line);
					vector<string> vars=split(line,' ');
					if(vars.size()>4){
						unsigned i = atoi(vars.at(0).c_str())-1;
						//cout << i << endl;
						type.at(i) = atoi(vars.at(1).c_str());
						x.at(i)    = atof(vars.at(2).c_str())*Lx;
						y.at(i)    = atof(vars.at(3).c_str())*Ly;
						z.at(i)    = atof(vars.at(4).c_str())*Lz;
						//cout << i << " " <<  type.at(i) << " " << x.at(i) << " " << y.at(i) << " " << z.at(i) << endl;
					}
				}
				not_done=false;
			}
		}
	} else {
		cerr << "error: Unknown input file format. filename = " <<  ifilename;
	}
	//cout << "Done reading " << ifilename << " Lx = " << Lx << " Ly = " << Ly << " Lz = " << Lz <<  " num_atoms = " << x.size() << endl;
}


/**
 * Write coordinates of particles to xyz-file.
 */
//void Rotational_order::write_xyz(ostream& out){}
void Rotational_order::write_xyz(string ofilename,double Qmin,double Qmax){
	using namespace boost::iostreams;
	// Count number of particles in the allowed interval
	unsigned numSelected = 0;
	for (unsigned i=0;i<number_of_particles();i++) 
		if(qiAvg.at(i)>Qmin && qiAvg.at(i)<Qmax)
			numSelected++;
	double Navg   = 0;for(unsigned i=0;i<number_of_particles();i++) Navg += (double)number_of_neighbors.at(i);Navg/=(double)number_of_particles();
	double Ql    = 0;for(unsigned i=0;i<number_of_particles();i++) Ql += qi.at(i);Ql/=(double)number_of_particles();
	double QlAvg = 0;for(unsigned i=0;i<number_of_particles();i++) QlAvg += qiAvg.at(i);QlAvg/=(double)number_of_particles();
	
	vector<string> fnames = split(ofilename,'.');
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
	out << numSelected <<endl;
	out << " numTypes=" << number_of_types();
	out << " sim_box=RectangularSimulationBox," << Lx << "," << Ly << "," << Lz;
	out << " columns=type,x,y,z,Ni,ql,qlAvg,Ns";
	out << " Navg=" << Navg << " Ql=" << Ql << " QlAvg=" << QlAvg << endl; 
	for (unsigned i=0;i<number_of_particles();i++) {
		if(qiAvg.at(i)>Qmin && qiAvg.at(i)<Qmax){
			out << type.at(i); 
			out << " " << x.at(i)  << " " << y.at(i) << " " << z.at(i);
			out << " " << number_of_neighbors.at(i);
			out << " " << qi.at(i) << " " << qiAvg.at(i);
			if(number_of_S_connections.size()==number_of_particles())
			out << " " << number_of_S_connections.at(i);
			out << endl;
		}
	}
}



/**
 * Return a string with various information about the lattice
 */
string Rotational_order::info(double Qmin,double Qmax){
	// Count number of particles in the allowed interval
	unsigned numSelected = 0;
	for (unsigned i=0;i<number_of_particles();i++) 
		if(qiAvg.at(i)>Qmin && qiAvg.at(i)<Qmax)
			numSelected++;
	double Navg   = 0;for(unsigned i=0;i<number_of_particles();i++) Navg += (double)number_of_neighbors.at(i);Navg/=(double)number_of_particles();
	double Ql    = 0;for(unsigned i=0;i<number_of_particles();i++) Ql += qi.at(i);Ql/=(double)number_of_particles();
	double QlAvg = 0;for(unsigned i=0;i<number_of_particles();i++) QlAvg += qiAvg.at(i);QlAvg/=(double)number_of_particles();
	
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
	out << "Average number of neighbors:  " << Navg << endl;
	out << "Average qi:               Ql= " << Ql << endl;
	out << "Average qiAvg:         QlAvg= " << QlAvg << endl;
	out << "Q minimum value:        Qmin= " << Qmin << endl;
	out << "Q maximum value:        Qmax= " << Qmax << endl;
	out << "Number of selected par.:      " << numSelected << endl;
	/*
	out << "           DEBUG INFO" << endl;
	out << " qi.size()= " << qi.size() << " qiAvg.size()= " << qiAvg.size() <<endl;
	out << " qlm.size()= " << qlm.size() << " qlmAvg.size()= " << qlmAvg.size() <<endl;
	out << " number_of_neighbors.size()= " << number_of_neighbors.size() <<endl;
	out << " number_of_S_connections.size()= " << number_of_S_connections.size() <<endl;
	*/
	return out.str();
}
