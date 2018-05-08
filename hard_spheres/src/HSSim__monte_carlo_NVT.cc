
#include "HSSim.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>

#include <boost/random/mersenne_twister.hpp>      // For random numbers
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

using namespace std;

/**
  Run a monte_carlo (MC) simulation
*/
void HSSim::monte_carlo_NVT(unsigned steps,double step_size,unsigned frames){
	boost::random::uniform_real_distribution<> dice_m1to1( -1.0 , 1.0 );
	boost::random::uniform_real_distribution<> dice_0to1( 0.0 , 1.0 );
	boost::random::uniform_int_distribution<> dice_particle( 0 , number_of_particles()-1 );
	cout << "  Perform " << frames << " frames of " << steps << " MC steps (constant NVT):" 
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

