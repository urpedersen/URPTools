/**
 * Lattice.cpp
 *
 *  Created on: Aug 15, 2012
 *      Author: Ulf R. Pedersen
 */

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cassert>
#include <cmath>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

/**
#define png_infopp_NULL (png_infopp)NULL
#define int_p_NULL (int*)NULL
#include <boost/gil/gil_all.hpp>
#include <boost/gil/extension/io/png_dynamic_io.hpp>
using namespace boost::gil;
*/

#include "lattice.h"

using namespace std;

Lattice::Lattice() {

	beta = 0.0 ;
	chem_pot.clear() ;

	Lx=0 ;
	Ly=0 ;
	grid.clear() ;

	types=0 ;
	position_x.clear() ;
	position_y.clear() ;

	tries_particle_insertion = 0 ;
	tries_neighbour_swop = 0 ;
	tries_switch_identity = 0 ;
	success_particle_insertion = 0 ;
	success_neighbour_swop = 0 ;
	success_switch_identity = 0 ;

}

Lattice::~Lattice() {
	// TODO Auto-generated destructor stub
}

/*
 * Set up a Lx times Ly lattice of vacancies. Define two types of particles + vacancies (types=3)
 */
void Lattice::setup(
		unsigned int in_Lx,
		unsigned int in_Ly,
		unsigned int in_types,
		unsigned int seed
		) {

	// Set-up seed for pseudo random number generator
	//unsigned int seed = 5;
	gen.seed(seed);
	//cout << "Mersenne twister 19937: gen = " << gen << endl;

	types = in_types ; // (add vacancy and input types of particles)
	types++;

	chem_pot.clear();
	chem_pot.push_back(0.0);
	for (unsigned int type = 0 ; type < in_types ; type++ )
		chem_pot.push_back(0.0);

	// Set up lattice
	Lx = in_Lx ;
	Ly = in_Ly ;

	grid.clear();
	grid.assign ( Lx*Ly , 0 ) ;

	position_x.clear();
	position_y.clear();


	// Add vacuum and allocate memory for position vector.
	{
		vector<unsigned int> tmpX;
		position_x.push_back ( tmpX ) ;
		for(unsigned int i = 0 ; i < Lx*Ly ; i++ )
			position_x.at(0).push_back(i%Lx);

		vector<unsigned int> tmpY;
		position_y.push_back ( tmpY ) ;
		for(unsigned int i = 0 ; i < Lx*Ly ; i++ )
			position_y.at(0).push_back(i/Lx);
	}
	for ( unsigned int i = 1 ; i < types ; i++ ) {

		vector<unsigned int> tmpX   ;
		position_x.push_back ( tmpX );

		vector<unsigned int> tmpY    ;
		position_y.push_back ( tmpY );
	}

	// Set MC parameters
	weight_particle_insertion = 1/3.0 ;
	weight_neighbour_swop     = weight_particle_insertion ;
	weight_switch_identity    = weight_particle_insertion ;

	// Initialize umbrella parameters
	umbrella_center.clear();
	umbrella_center.assign( types , 0.0 );
	umbrella_kappa.clear();
	umbrella_kappa.assign( types , 0.0 );

}

/**
 * (Re)Set the inverse temperature
 */
void Lattice::set_beta(double in_beta){
	beta=in_beta;
}

/**
 * (Re)Set the chemical potential of a given spices
 */
void Lattice::set_chemical_potential(unsigned int type,double in_chem_pot ){
	assert(type+1<types);
	chem_pot.at(type+1)=in_chem_pot;
}

/**
 * Set the weights that determines how often veries moves are made
 */
void Lattice::set_weights(
		double in_weight_particle_insertion ,
		double in_weight_neighbour_swop     ,
		double in_weight_switch_identity    ){

	weight_particle_insertion = in_weight_particle_insertion;
	weight_neighbour_swop = in_weight_neighbour_swop;
	weight_switch_identity = in_weight_switch_identity;

	double weight_sum  = weight_particle_insertion;
	        weight_sum += weight_neighbour_swop ;
	        weight_sum += weight_switch_identity ;

    weight_particle_insertion /= weight_sum ;
    weight_neighbour_swop /= weight_sum ;
    weight_switch_identity /= weight_sum ;

}



/**
 * Add particle.
 */
void Lattice::create_particle(unsigned int x, unsigned int y, unsigned int type) {
	assert (x<Lx) ; // Particle is on lattice.
	assert (y<Ly) ;
	assert (type<types) ; // Particle is a allowed type.
	assert (type!=0) ;

	// Confirm that lattice site is empty and add particle to site grid vector
	assert(type_at(x,y)==0) ;
	grid.at(site(x,y))=type ;

	position_x.at(type).push_back(x) ;
	position_y.at(type).push_back(y) ;

	// TOOO Move a vacancy
	//unsigned int x_back = position_x.at(type).back() ;
	//unsigned int y_back = position_y.at(type).back() ;

	//position_x.at(0).at(index) = x_back ;
	//position_y.at(0).at(index) = y_back ;

	position_x.at(0).pop_back();
	position_y.at(0).pop_back();
}



/**
 * Remove particle. The last particle is put in its place in the position_x and position_y arrays (and the last element is eliminated).
 */
void Lattice::abolish_particle(unsigned int type, unsigned int index ) {

	assert( type  < types ) ;
	assert( index < position_x.at(type).size() ) ;
	assert( index < position_y.at(type).size() ) ;

	unsigned int x = position_x.at(type).at(index) ;
	unsigned int y = position_y.at(type).at(index) ;

	assert( grid.at(site(x,y))==type ); // Consistency test of particle type.
	grid.at(site(x,y)) = 0 ;

	unsigned int x_back = position_x.at(type).back() ;
	unsigned int y_back = position_y.at(type).back() ;

	position_x.at(type).at(index) = x_back ;
	position_y.at(type).at(index) = y_back ;

	position_x.at(type).pop_back();
	position_y.at(type).pop_back();

	position_x.at(0).push_back(x);
	position_y.at(0).push_back(y);
}



/**
 * Swop the position of two particles.
 */
void Lattice::swop_particles(
	unsigned int type, unsigned int index,
	unsigned int type2, unsigned int index2){

	unsigned int x = position_x.at(type).at(index);
	unsigned int y = position_y.at(type).at(index);
	unsigned int x2 = position_x.at(type2).at(index2);
	unsigned int y2 = position_y.at(type2).at(index2);

	grid.at(site(x,y))  =type2;
	grid.at(site(x2,y2))=type;

	position_x.at(type).at(index)=x2   ;
	position_y.at(type).at(index)=y2   ;
	position_x.at(type2).at(index2)=x  ;
	position_y.at(type2).at(index2)=y  ;
}



/**
 * Do a Monte Carlo (MC) attempt
 */
unsigned int Lattice::mc_step ( ) {

/*
	unsigned int tries_particle_insertion = 0 ;
	unsigned int tries_neighbour_swop = 0 ;
	unsigned int tries_switch_identity = 0 ;
	unsigned int success_particle_insertion = 0 ;
	unsigned int success_neighbour_swop = 0 ;
	unsigned int success_switch_identity = 0 ;
*/

	boost::random::uniform_real_distribution<> dist_test( 0.0 , 1.0 );
	double test = dist_test(gen);

	unsigned int success = 0;

	if( test < weight_particle_insertion ){
		success = mc_step_insert (  ) ;
		tries_particle_insertion++ ;
		success_particle_insertion += success;
	} else if ( test < weight_particle_insertion + weight_neighbour_swop  ) {
		success = mc_step_swop (  ) ;
		tries_neighbour_swop++ ;
		success_neighbour_swop += success ;
	} else {
		success = mc_step_type (  ) ;
		tries_switch_identity++ ;
		success_switch_identity += success ;
	}

	return success ;
}


/**
 * Swop the position of two particles.
 */
unsigned int Lattice::mc_step_swop ( ) { // TODO

	if( number_of_particles() == 0 || types < 3 ){
		return 0;
	}

	// First particle
	boost::random::uniform_int_distribution<> dist_type( 1 , types - 1  );
	unsigned int type = dist_type(gen);

	unsigned int x;
	unsigned int y;
	unsigned int index;

	if( number_of_particles(type) > 0 ){
		boost::random::uniform_int_distribution<> dist_index( 0 , number_of_particles(type) - 1 );
		index = dist_index(gen);
		x = position_x.at(type).at(index);
		y = position_y.at(type).at(index);
	} else {
		return 0;
	}

	// 2nd particle
	unsigned int x2;
	unsigned int y2;
	unsigned int type2 ;
	unsigned int index2;
	type2 = type;
	while(type2 == type) {
		boost::random::uniform_int_distribution<> dist_type2( 1 , types - 1  );
		type2 = dist_type2(gen);

		if( number_of_particles(type2) > 0 ) {
			boost::random::uniform_int_distribution<> dist_index2( 0 , number_of_particles(type2) - 1 );
			index2 = dist_index2(gen);
			x2 = position_x.at(type2).at(index2);
			y2 = position_y.at(type2).at(index2);
		} else {
			return 0;
		}
	}


	double current_energy =
			 site_energy( x ,  y  )
			+site_energy( x2 , y2 );

	double new_energy =
			 site_energy_type( x ,  y  , type2 )
			+site_energy_type( x2 , y2 , type );

	double dU = new_energy - current_energy ;

	double  acceptance_rule  = exp ( -1.0 * beta * dU );

	boost::random::uniform_real_distribution<> dist_test( 0.0 , 1.0 );
	double test = dist_test(gen);

	if ( test < acceptance_rule ) {
		swop_particles(type,index,type2,index2);
		return 1;
	} else {
		return 0;
	}

	return 0; // The program should never reach this
}



unsigned int Lattice::mc_step_type ( ) { // TODO
	return 0;
}



/**
 * Perform a Monte Carlo step where a particle is inserted or removed.
 */
unsigned int Lattice::mc_step_insert ( ) {

	boost::random::uniform_int_distribution<> dist_type( 1 , types - 1  );
	int type = dist_type(gen);

	boost::random::uniform_int_distribution<> dist_insert_or_remove( 0 , 1 );
	int insert_or_remove = dist_insert_or_remove(gen);

	if( insert_or_remove==1 ) {	// Attempt to insert particle

		boost::random::uniform_int_distribution<> dist_x( 0 , Lx-1 );
		int x = dist_x(gen);

		boost::random::uniform_int_distribution<> dist_y( 0 , Ly-1 );
		int y = dist_y(gen);

		if( grid.at(site(x,y)) == 0 ) {

			double current_energy  = site_energy(x,y);
					current_energy += umbrella_energy( type , number_of_particles(type)     );

			double new_energy      = site_energy_type(x,y,type);
				    new_energy     += umbrella_energy( type , number_of_particles(type) + 1 );

			double dU = new_energy - current_energy ;

			double  acceptance_rule  = (double)volume();
					acceptance_rule /= (double)(number_of_particles(type)+1);
			        acceptance_rule *= exp ( -1.0 * beta * ( dU  - chem_pot.at(type) ) );

			boost::random::uniform_real_distribution<> dist_test( 0.0 , 1.0 );
			double test = dist_test(gen);

			//cout << endl << "Insert: " << dU << " " <<  acceptance_rule << " " << test << endl ;

			//cout << " " << number_of_particles(type) ;
			if ( test < acceptance_rule ) {
				create_particle(x,y,type);
				return 1;
			} else {
				return 0;
			}
			//cout << " " << number_of_particles(type) ;

		}

	} else { // Attempt to remove a particle
		if ( number_of_particles(type)>0 ) {

			boost::random::uniform_int_distribution<> dist_index( 0 , number_of_particles(type) - 1 );
			int index = dist_index(gen);

			unsigned int x = position_x.at(type).at(index);
			unsigned int y = position_y.at(type).at(index);

			double current_energy  = site_energy(x,y);
					current_energy += umbrella_energy( type , number_of_particles(type) );

			double new_energy  = site_energy_type(x,y,0);
				    new_energy += umbrella_energy( type , number_of_particles(type) - 1 );

			double dU = new_energy - current_energy ;

			double   acceptance_rule  = (double)number_of_particles(type);
					 acceptance_rule /= (double)volume();
			         acceptance_rule *= exp ( - 1.0 * beta * ( dU  + chem_pot.at(type) ) );

			boost::random::uniform_real_distribution<> dist_test( 0.0 , 1.0 );
			double test = dist_test(gen);

			//cout << endl << "Remove: " << dU << " " <<  acceptance_rule << endl ;

			if ( test < acceptance_rule ) {
				abolish_particle(type,index);
				return 1;
			} else {
				return 0;
			}
		}
	}

	return 0;

}


/**
 * Return the volume of the system
 */
unsigned int Lattice::volume(){
	return Lx*Ly;
}

/**
 * Return the area of the system
 */
unsigned int Lattice::area(){
	return Lx*Ly;
}


/**
 * Return x width of the system
 */
unsigned int Lattice::getLx(){
	return Lx;
}

/**
 * Return x hight of the system
 */
unsigned int Lattice::getLy(){
	return Ly;
}


/**
 * Return the lattice site energy .
 * Note, the lattice model is determined here.
 * TODO Let the user select between different kind of energy models, i.e. 4 or 8 neighbors etc.
 */
double Lattice::site_energy(unsigned int x,unsigned int y) {
	return site_energy_8neighbours(x,y);
}




/**
 * Return energy of lattice site (x,y) of model with binding energy to 8 neighbors.
 */
double Lattice::site_energy_8neighbours(unsigned int x,unsigned int y) {
	assert(x<Lx);
	assert(y<Ly);

	// Interaction parameters
	static int range = 1;
	static double same_kind = -1.0;
	static double another_kind = -2.0;

	double out=0.0;

	if ( type_at(x,y)==0 ) {
		return 0.0 ;
	} else {
		for ( int dx=-range ; dx<range+1 ; dx++ ) {
			int xx = (x-dx+Lx)%Lx;
			for ( int dy=-range ; dy<range+1 ; dy++ ) {
				int yy = (y-dy+Ly)%Ly;
				if(type_at(xx,yy)==0.0){
					out+=0.0;
				}else if(type_at(x,y)==type_at(xx,yy)){
					out+=same_kind;
				}else{
					out+=another_kind;
				}
				// cout << "("<< x << "," same_kind<< y << ")  (" << xx << "," << yy << ") " << out << endl;
			}
		}
		out-=same_kind; // We have (accidently) added a self-bond that should be subtracted
	}

	return out;

}



/**
 * Return energy of lattice site (x,y).
 */
double Lattice::site_energy_4neighbours(unsigned int x,unsigned int y) {
	assert(x<Lx);
	assert(y<Ly);

	// Interaction parameters
	static double same_kind = -1.0;
	static double another_kind = -2.0;

	double out=0.0;

	if ( type_at(x,y)==0 ) {
		return 0.0 ;
	} else {

		{	// 1st
			int dx=-1;
			int dy=0;
			int xx = (x-dx+Lx)%Lx;
			int yy = (y-dy+Ly)%Ly;

			if(type_at(xx,yy)==0.0){
				out+=0.0;
			}else if(type_at(x,y)==type_at(xx,yy)){
				out+=same_kind;
			}else{
				out+=another_kind;
			}
		}

		{   // 2nd
			int dx=1;
			int dy=0;
			int xx = (x-dx+Lx)%Lx;
			int yy = (y-dy+Ly)%Ly;

			if(type_at(xx,yy)==0.0){
				out+=0.0;
			}else if(type_at(x,y)==type_at(xx,yy)){
				out+=same_kind;
			}else{
				out+=another_kind;
			}
		}

		{   // 3rd
			int dx=0;
			int dy=-1;
			int xx = (x-dx+Lx)%Lx;
			int yy = (y-dy+Ly)%Ly;

			if(type_at(xx,yy)==0.0){
				out+=0.0;
			}else if(type_at(x,y)==type_at(xx,yy)){
				out+=same_kind;
			}else{
				out+=another_kind;
			}
		}

		{   // 4th
			int dx=0;
			int dy=1;
			int xx = (x-dx+Lx)%Lx;
			int yy = (y-dy+Ly)%Ly;

			if(type_at(xx,yy)==0.0){
				out+=0.0;
			}else if(type_at(x,y)==type_at(xx,yy)){
				out+=same_kind;
			}else{
				out+=another_kind;
			}
		}
	}

	return out;

}

/**
 * Return energy of lattice site (x,y)
 *   AS IF it was populated with a particle of input type
 */
double Lattice::site_energy_type(unsigned int x,unsigned int y,unsigned int type) {

	assert(x<Lx);
	assert(y<Ly);
	assert(type<types);

	unsigned int old_type = type_at(x,y);
	grid.at(site(x,y)) = type ;
	double out = site_energy(x,y);
	grid.at(site(x,y)) = old_type ;

	return out ;

}


/**
 * Return the total energy of the system
 */
double Lattice::total_energy(){
	double out=0.0;
	for ( unsigned int x = 0 ; x < Lx ; x++ ) {
		for ( unsigned int y = 0 ; y < Ly ; y++ ) {
			out+=0.5*site_energy(x,y);
		}
	}

	return out;

}


/**
 * Set the umbrella parameter of a given particle type.
 */
void Lattice::set_umbrella_parameters ( unsigned int type , double kappa , double center ){
	umbrella_kappa.at(type+1) = kappa   ;
	umbrella_center.at(type+1) = center ;
}


/**
 * Return the umbrella-energy for the type'th umbrella.
 */
double Lattice::umbrella_energy ( unsigned int type , unsigned int num_par ) {
	assert(type<types);
	assert(umbrella_center.size()>0);
	assert(umbrella_kappa.size()>0);

	double delta = num_par - umbrella_center.at(type) ;
	return 0.5 * umbrella_kappa.at(type) * delta * delta ;
}

double Lattice::umbrella_energy ( unsigned int type ) {
	return umbrella_energy( type , (double)number_of_particles(type) );
}

unsigned int Lattice::type_at(unsigned int x, unsigned int y ){
	assert(x<Lx);
	assert(y<Ly);

	return grid.at(y*Lx+x);

}


/**
 * Return string with info about the lattice
 */
string Lattice::str(){
	stringstream out;

	out << "Types = " << types << endl;
	out << "beta = " << beta << endl;

	out << "Chemical potentials: " ;
	for ( unsigned int type = 0 ; type < types ; type++ ) out << chem_pot.at(type) << " ";
	out << endl;

	out << "Umbrella kappas: " ;
	for ( unsigned int type = 0 ; type < types ; type++ ) out << umbrella_kappa.at(type) << " ";
	out << endl;

	out << "Umbrella centers: " ;
	for ( unsigned int type = 0 ; type < types ; type++ ) out << umbrella_center.at(type) << " ";
	out << endl;

	out << "weight_particle_insertion = " << weight_particle_insertion << endl ;
	out << "weight_neighbour_swop = " << weight_neighbour_swop << endl ;
	out << "weight_switch_identity = " << weight_switch_identity << endl ;

	out << "Particle insertion success: "
			<< success_particle_insertion << "/" << tries_particle_insertion
			<< " = " << (double)success_particle_insertion /(double)tries_particle_insertion  << endl;

	out << "Neighbor swop success: "
			<< success_neighbour_swop << "/" << tries_neighbour_swop
			<< " = " << (double)success_neighbour_swop /(double)tries_neighbour_swop  << endl;

	out << "Switch identity success: "
			<< success_switch_identity << "/" << tries_switch_identity
			<< " = " << (double)success_switch_identity /(double)tries_switch_identity  << endl;

	out << "X = " << Lx << endl;
	out << "Y = " << Ly << endl;
	out << "Volume = " << volume() << endl;
	out << "Bond_energy = " << total_energy() << endl;

	out << "Type population: " ;
	for ( unsigned int type = 0 ; type < types ; type++ ) out << number_of_particles(type) << " ";
	out << endl;

	out << configuration_str();

	if(false) { // TODO make it a function to print particle positions
		for(unsigned int type = 1 ; type < types ; type++ ) {
			out << "Positions of " <<  position_x.at(type).size() << " type " << type << " particles:\n";
			for( unsigned int index = 0 ; index<position_x.at(type).size() ; index++ ) {
				out << position_x.at(type).at(index) << " " << position_y.at(type).at(index) << endl;
			}
		}
	}

	return out.str();
}

/**
 * Return string with lattice decoration in simple text.
 */
string Lattice::grid_str(){
	stringstream out;

	out << Lx << " " << Ly << endl;
	out << "# Comment line";
	for (unsigned int x = 0 ; x < Lx ; x++ ) {
		out << endl ;
		for (unsigned int y = 0 ; y < Ly ; y++ )
			out << grid.at(site(x,y)) << " ";
	}

	return out.str();
}



/**
 * Return string with lattice decoration in the the "portable pixmap format" (ppm) format.
 */
string Lattice::ppm_str (  ) {
	stringstream out;

	out << "P3" << endl;
	out << Lx << " " << Ly << endl;
	out << "255" << endl;

	for (unsigned int y = 0 ; y < Ly ; y++ ) {
		for (unsigned int x = 0 ; x < Lx ; x++ )
			if(grid.at(site(x,y))==0){
				out << "80 80 200" << "   ";
			} else {
				out << "80 " << (int)(-1.*site_energy(x,y)/8.*200.) << " 80" << "   ";
			}
		out << endl ;
	}

	return out.str();
}


string Lattice::configuration_str(){

	stringstream out;
	out << "Lattice decoration:\n" ;
	for ( unsigned int y = Ly-1 ; y+1 > 0 ; y-- ) {
		if(y/10==0) out << " ";
		out << y << " [ ";
		for ( unsigned int x = 0 ; x < Lx ; x++ ) {
			out << type_at(x,y) << " ";
		}
		out << "]" << endl;
	};
	out	<< "    ";
	for ( unsigned int x = 0 ; x < Lx ; x++ ) out << " " << x ;
	out << endl;

	out << "Lattice energies times -1:\n";
	for ( unsigned int y = Ly-1 ; y+1 > 0 ; y-- ) {
		if(y/10==0) out << " ";
		out << y << " [ ";
		for ( unsigned int x = 0 ; x < Lx ; x++ ) {
			unsigned int tmp = (unsigned int)(-1.0*site_energy(x,y));
			out << tmp << " ";
			if(tmp/10==0) out << " ";
		}
		out << "]" << endl;
	};
	out	<< "    ";
	for ( unsigned int x = 0 ; x < Lx ; x++ ){
		out << " " << x ;
		if(x/10==0) out << " ";
	}
	out << endl;

	return out.str();
}

string Lattice::energy_str(){

	stringstream out;

	out << total_energy();

	for ( unsigned int type = 0 ; type < types ; type++ )
		out << " " << number_of_particles(type);

	for ( unsigned int type = 0 ; type < types ; type++ )
		out << " " << umbrella_energy(type);

	for ( unsigned int type = 0 ; type < types ; type++ )
		out << " " << umbrella_center.at(type);

	return out.str() ;
}

unsigned int Lattice::number_of_particles(unsigned int type){
	return position_x.at(type).size();
}

unsigned int Lattice::number_of_particles(){
	double sum = 0.0 ;
	for (unsigned int type = 1 ; type < types ; type++)
		sum += number_of_particles(type);
	return sum ;
}


/**
 * Return lattice position integer (in grid vector) in the the primary image.
 */
unsigned int Lattice::site ( int x , int y ) {

	while( x <  0 ) x += Lx ;	// TODO use modulo function for the periodic boundaries
	while( y <  0 ) y += Ly ;
	while( x >= Lx) x -= Lx ;
	while( y >= Ly) y -= Ly ;

	return Lx*y+x ;
}


unsigned int Lattice::x_img ( int x ) {

	while( x <  0 ) x += Lx ;	// TODO use modulo function for the periodic boundaries
	while( x >= Lx) x -= Lx ;

	return x ;
}

unsigned int Lattice::y_img ( int y ) {

	while( y <  0 ) y += Ly ;	// TODO use modulo function for the periodic boundaries
	while( y >= Ly) y -= Ly ;

	return y;
}
