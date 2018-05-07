/*
 * HSSim.cc
 *
 *  Created on: May 2, 2018
 *      Author: Ulf R. Pedersen, www.urp.dk
 */

#include "HSSim.h"

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

#include <boost/random/mersenne_twister.hpp>      // For random numbers
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "split.h"
#include "Cell_list.h"

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
HSSim::HSSim() :
		type(),
		x(),
		y(),
		z(),
		cell_list(),
		Lx(1.0),
		Ly(1.0),
		Lz(1.0),
		neighbour_cutoff(1.5)
	{
	// Do nothing in constructor
}

HSSim::~HSSim() {
	this->clear();
}

/**
* Clear memory for a new computation
*/ 
void HSSim::clear() {
	type.clear();
	x.clear();
	y.clear();
	z.clear();
}


/**
 * Add a particle of in_type to the lattice
 */
void HSSim::add_particle(unsigned in_type,double in_x,double in_y,double in_z){
	type.push_back(in_type);
	x.push_back(in_x);
	y.push_back(in_y);
	z.push_back(in_z);
}
	

void HSSim::set_neighbour_cutoff(double in_neighbour_cutoff){
  neighbour_cutoff = in_neighbour_cutoff;
}


/**
 * Return the number of particles
 */
unsigned HSSim::number_of_particles(){
	return x.size();
}


/**
 * Return the number of particles types
 */
unsigned HSSim::number_of_types(){
	unsigned out=0;
	if( type.size()>0 )
		for(unsigned i=0;i<type.size();i++)
			if(type.at(i)>=out) out=type.at(i)+1;
	return out;
}

unsigned HSSim::number_of_particles_of_type(unsigned test_type){
	unsigned out=0;
	for(unsigned n=0;n<type.size();n++)
		if(type.at(n)==test_type)
			out++;
	return out;
}

/**
 * Return volume of lattice
 */
double HSSim::volume(){
	return Lx*Ly*Lz;
}


/**
  Generate an ideal gas configuration
*/
void HSSim::generate_ideal_gas_positions(unsigned num_particles,double in_Lx,double in_Ly,double in_Lz){
	type.clear();
	x.clear();
	y.clear();
	z.clear();
	
	Lx=in_Lx;
	Ly=in_Ly;
	Lz=in_Lz;
	
	boost::random::uniform_real_distribution<> dice( 0.0 , 1.0 );

	for(unsigned i = 0;i<num_particles;i++){
	  float random_x = dice(gen)*Lx;
	  float random_y = dice(gen)*Ly;
	  float random_z = dice(gen)*Lz;
	  add_particle(0,random_x,random_y,random_z);
	}	
	
	cell_list.build(x,y,z,Lx,Ly,Lz,neighbour_cutoff);	
}

/**
* Load coordinates of particles from xyz-file or an zipped xyz.gz file or a LAMMPS *.atom file
*/
bool HSSim::load_xyz(string ifilename){
  return load_xyz(ifilename,0,10.0,10.0,10.0);
}

/**
* Load coordinates of particles from xyz-file or an zipped xyz.gz file or a LAMMPS *.atom file
*/
bool HSSim::load_xyz(string ifilename,unsigned frame,double in_Lx,double in_Ly,double in_Lz){
	using namespace boost::iostreams;

	type.clear();
	x.clear();
	y.clear();
	z.clear();
	
	bool sucessfull_load_of_frame = true;
	string fileformat="";

	Lx=in_Lx;
	Ly=in_Ly;
	Lz=in_Lz;

	filtering_istream in;
	vector<string> fnames = split(ifilename,'.');
	//in.push(file_source(ifilename));
	if(fnames.size()<2){
		cerr << "error: Incompatible name of input file. Should be an *.xyz, *.xyz.gz or atom file." << endl;
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
		// Begin to read header of selected frame
		getline(in,line); // Number of atoms
		if(!in.good()) sucessfull_load_of_frame = false;
		unsigned num_atoms=atoi(line.c_str());

		getline(in,line); // Comment line
		if(!in.good()) sucessfull_load_of_frame = false;
		
		// Skip frames
		for(unsigned f=0;f<frame;f++){

			for(unsigned i = 0;i<num_atoms;i++){
				getline(in,line);
				if(!in.good()) sucessfull_load_of_frame = false;
			}
			getline(in,line); // Number of atoms
			if(!in.good())	sucessfull_load_of_frame = false;
			num_atoms=atoi(line.c_str());
			getline(in,line); // Comment line
			if(!in.good())	sucessfull_load_of_frame = false;
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
			if(!in.good())	sucessfull_load_of_frame = false;
			
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
			getline(in,line);
			if(!in.good()){
				sucessfull_load_of_frame = false;
				not_done=false;
			}
			if(!in.good()){
				sucessfull_load_of_frame = false;
				not_done=false;
			}
			vector<string> sections = split(line,' ');
			if(sections.size()>3 && sections.at(0)=="ITEM:" && sections.at(1)=="NUMBER" && sections.at(2)=="OF" && sections.at(3)=="ATOMS"){
				current_frame++;
				getline(in,line);
				if(!in.good()){
					sucessfull_load_of_frame = false;
					not_done=false;
				}
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
					getline(in,line);
					if(!in.good()){cerr << "error: while reading " << ifilename << endl;abort();}
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

	cell_list.build(x,y,z,Lx,Ly,Lz,neighbour_cutoff);	

	//cout << "Done reading " << ifilename << " Lx = " << Lx << " Ly = " << Ly << " Lz = " << Lz <<  " num_atoms = " << x.size() << endl;
	return sucessfull_load_of_frame;
}

void HSSim::wrap_into_box(double xO,double yO,double zO){
	for (unsigned p=0;p<number_of_particles();p++){
		x.at(p)-=xO;
		y.at(p)-=yO;
		z.at(p)-=zO;
		x.at(p)-=Lx*round(x.at(p)/Lx);
		y.at(p)-=Ly*round(y.at(p)/Ly);
		z.at(p)-=Lz*round(z.at(p)/Lz);
	}
}

/**
  Write coordinates to traj.xyz
*/
void HSSim::write_xyz(){
  write_xyz("traj.xyz");
}

/**
 * Write coordinates of particles to xyz-file.
 */
void HSSim::write_xyz(string ofilename){
	using namespace boost::iostreams;
	
	vector<string> fnames = split(ofilename,'.');
	filtering_ostream out;

	if(fnames.size()<2){
		cerr << "error: Incompatible name of output file. Should be *.xyz or *.xyz.gz" << endl;
		abort();
	}
	if(fnames.back()=="xyz"){
		out.push(file_sink(ofilename,ios_base::app));
	} else if ( fnames.back()=="gz" && fnames.at(fnames.size()-2)=="xyz" )  {
		out.push(gzip_compressor());
		out.push(file_sink(ofilename,ios_base::app));
	} else {
		cerr << "error: Incompatible name of output file. Should be *.xyz or *.xyz.gz" << endl;
		abort();
	}

	out << number_of_particles() << endl;
	out << "Hard-spheres numTypes=" << number_of_types();
	out << " sim_box=RectangularSimulationBox," << Lx << "," << Ly << "," << Lz;
	out << " columns=type,x,y,z" << endl;
	for (unsigned i=0;i<number_of_particles();i++) {
		out << type.at(i); 
		out << " " << x.at(i)  << " " << y.at(i) << " " << z.at(i);
		out << endl;
	}
}

/**
  Run a monte_carlo (MC) simulation
*/
void HSSim::monte_carlo(unsigned steps,double step_size,unsigned frames){
	boost::random::uniform_real_distribution<> dice_m1to1( -1.0 , 1.0 );
	boost::random::uniform_real_distribution<> dice_0to1( 0.0 , 1.0 );
	boost::random::uniform_int_distribution<> dice_particle( 0 , number_of_particles()-1 );
	cout << "  Perform " << frames << " frames of " << steps << " MC steps" 
	  << ", stepsize = " << step_size << endl;

	unsigned rejected = 0;
	unsigned attempts = 0;

	for(unsigned frame=0;frame<frames;frame++){
	  if      (frame%50==0 ) cout << endl << frame << "\t|";
	  else if (frame%10==0) cout << "|";
	  else if (frame%5==0) cout << ":";
	  else if (frame%50==49) cout << ".|";
	  else cout << ".";
	  cout << flush;


	  wrap_into_box(0,0,0);
	  write_xyz("traj.xyz");

	  for(unsigned s=0;s<steps*number_of_particles();s++){
		// TODO only build neighbour list when nessesary, and safely!
		if(s%number_of_particles()*5==0)
	    	cell_list.build(x,y,z,Lx,Ly,Lz,neighbour_cutoff);	
		
		// Pick a random particle
      	unsigned p = dice_particle(gen);

      	// Random step (in sphere)
      	double r2=2.0;
      	double xold=x.at(p);
	  	double yold=y.at(p);
	  	double zold=z.at(p);
		bool old_is_overlapping = is_overlapping(p);
	  	double dx,dy,dz;
	    while(r2>1) {
		  dx = dice_m1to1(gen);
		  dy = dice_m1to1(gen);
		  dz = dice_m1to1(gen);
		  r2=dx*dx+dy*dy+dz*dz;
	  	}
	  	// Move particles, move back if the move result in overlap
	  	x.at(p)+=dx*step_size;
	  	y.at(p)+=dy*step_size;
	  	z.at(p)+=dz*step_size;
	    attempts++;
		if( is_overlapping(p) && !old_is_overlapping ){
	      rejected++;
		  // Reject move
		  x.at(p)=xold;
		  y.at(p)=yold;
		  z.at(p)=zold;
		}
	  }
	}
	cout << endl << "Rejected " << rejected << " of " << attempts 
	  << " MC move attempts (" << 100.0*(double)rejected/(double)attempts << "%)" << endl;
}

/**
 Return true if particle p is overlapping wit
 one of the particle in the neighbour list
*/
bool HSSim::is_overlapping(unsigned i){
  bool found_overlap=false;
  vector<unsigned> n;
  cell_list.neighbors(i,n);
  
  for(unsigned in=0;in<n.size();in++){
	unsigned j=n.at(in);
	if(!found_overlap && i!=j) {
	  double dx=x.at(j)-x.at(i);
	  double dy=y.at(j)-y.at(i);
	  double dz=z.at(j)-z.at(i);
	  dx-=Lx*nearbyint(dx/Lx);
	  dy-=Ly*nearbyint(dy/Ly);
	  dz-=Lz*nearbyint(dz/Lz);
	  double r2=dx*dx+dy*dy+dz*dz;
	  static double r2cut=1.0*1.0;
	  if(r2<r2cut){
		found_overlap = true;
	  }
	}
  }
  return found_overlap;
}

/**
  Return true if two or more particles are overlapping
  */
bool HSSim::is_overlapping(){
  bool found_overlap=false;
  for(unsigned i=0;i<x.size();i++)
	if(is_overlapping(i)) found_overlap=true;
  return found_overlap;
}

/**
 * Return a string with various information about the lattice
 */
string HSSim::info(){

  cout << number_of_particles() << endl;
	// Compute mean values

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
	out << "Is overlapping:               " << is_overlapping() << endl;
	out << "rcut:                     rc= " << neighbour_cutoff << endl;
	
	return out.str();
}
