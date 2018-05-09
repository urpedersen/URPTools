
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
void HSSim::monte_carlo_NpT(unsigned steps,double step_size,unsigned frames,double pressure,double volume_step_size){
	boost::random::uniform_real_distribution<> dice_m1to1( -1.0 , 1.0 );
	boost::random::uniform_real_distribution<> dice_0to1( 0.0 , 1.0 );
	boost::random::uniform_int_distribution<> dice_particle( 0 , number_of_particles()-1 );
	cout << "  Perform " << frames << " frames of " << steps << " MC steps (constant NpT):" 
	  << ", stepsize = " << step_size << ", pressure = " << pressure << endl;

	unsigned rejected = 0;
	unsigned attempts = 0;

	cout << "ener: frame volume is_overlapping" << endl;
	for(unsigned frame=0;frame<frames;frame++){
	  cout << "enr: " << frame << " " << Lx*Ly*Lz << " " << is_overlapping() << endl;
      //cell_list.build(x,y,z,Lx,Ly,Lz,neighbour_cutoff);
      //	  cout << is_overlapping() <<  endl;
	  /* Write pretty progress 
	  if      (frame%50==0 ) cout << endl << frame << "\t|";
      else if (frame%10==0 ) cout << "|";
	  else if (frame%5==0  ) cout << ":";
	  else if (frame%50==49) cout << ".|";
	  else cout << ".";
	  cout << flush;
	  */

	  // wrap_into_box(0,0,0); try and see if removing this fixes it?
	  cell_list.build(x,y,z,Lx,Ly,Lz,neighbour_cutoff);
	  write_xyz("traj.xyz");

	  for(unsigned s=0;s<steps*number_of_particles();s++){
		
		// TODO only build neighbour list when nessesary, and safely!
		if(s%number_of_particles()*5==0) cell_list.build(x,y,z,Lx,Ly,Lz,neighbour_cutoff);

		if(s%number_of_particles()*1==0 and !is_overlapping()){
			// Make a volume MC step
			double Vold=Lx*Ly*Lz;
			double dV = dice_m1to1(gen)*volume_step_size;
			double Vnew=Vold+dV;
			double dL = pow(Vnew/Vold,1./3.);
			// Scale positions of particles and box (if the configuration is not overlapping)
			for(unsigned p = 0;p<x.size();p++){
				x.at(p)*=dL;
				y.at(p)*=dL;
		  	    z.at(p)*=dL;
			}
			Lx*=dL;
			Ly*=dL;
			Lz*=dL;

			double randf = dice_0to1(gen);
			// See if move should be rejected TODO ....
			// use implicit beta=1
			double arg = pressure*dV-number_of_particles()*log(Vnew/Vold);
			if( randf<exp(arg) ){ // Restore old state
			  for(unsigned p = 0;p<x.size();p++){
				  x.at(p)/=dL;
				  y.at(p)/=dL;
		    	z.at(p)/=dL;
			  }
			  Lx/=dL;
			  Ly/=dL;
			  Lz/=dL;
				}
		}
		
		// Pick a random particle, and make a single particle step
      	unsigned p = dice_particle(gen);

      	// Random step (in sphere)
      	double r2=2.0;
      	double xold=x.at(p);
	  	double yold=y.at(p);
	  	double zold=z.at(p);
		bool old_is_overlapping = is_overlapping(p);
	  	double dx,dy,dz;
	    while(r2>1.0) {
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

