
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
void HSSim::monte_carlo_NpT(unsigned steps,unsigned frames,double step_size,double pressure,double volume_step_size){
	boost::random::uniform_real_distribution<> dice_m1to1( -1.0 , 1.0 );
	boost::random::uniform_real_distribution<> dice_0to1( 0.0 , 1.0 );
	boost::random::uniform_int_distribution<> dice_particle( 0 , number_of_particles()-1 );
	cout << "  Perform " << frames << " frames of " << steps << " MC steps (constant NpT):" 
	  << ", stepsize = " << step_size << ", pressure = " << pressure << endl;

	unsigned rejected = 0;
	unsigned attempts = 0;
	unsigned rejected_volume = 0;
	unsigned attempts_volume = 0;
	unsigned overlap_tests = 0; // Count overlapping configurations
	unsigned overlap_hits_before = 0;
	unsigned overlap_hits_after = 0;

	cout << "ener: frame volume packing_fraction" << endl;
	for(unsigned frame=0;frame<frames;frame++){
      //	  cout << is_overlapping() <<  endl;
	  /* Write pretty progress 
	  if      (frame%50==0 ) cout << endl << frame << "\t|";
      else if (frame%10==0 ) cout << "|";
	  else if (frame%5==0  ) cout << ":";
	  else if (frame%50==49) cout << ".|";
	  else cout << ".";
	  */
	  double volume = Lx*Ly*Lz;
	  double packing_fraction = (double)x.size()*atan(1)*4/6/volume; 
	  cout << "ener: " << frame << " " << volume << " " << packing_fraction << " "  << endl;
	  cout << flush;

	  wrap_into_box(0,0,0);
	  //cell_list.build(x,y,z,Lx,Ly,Lz,neighbour_cutoff);
	  write_xyz(ofilename,true);

	  for(unsigned s=0;s<steps*x.size();s++){
		
		// TODO only build neighbour list when nessesary, and safely (using a skin).
		if(s%number_of_particles()*10==0){
		  overlap_tests++;
		  if(is_overlapping()) overlap_hits_before++;
		  cell_list.build(x,y,z,Lx,Ly,Lz,neighbour_cutoff);
		  if(is_overlapping()) overlap_hits_after++;
		}

		if(s%number_of_particles()*1==0 and volume_step_size>0.0 and !is_overlapping()){
			// Make a volume MC step
			double lnV_old = log(Lx*Ly*Lz);
			double lnV_new = lnV_old + dice_m1to1(gen)*volume_step_size;
			double dV = exp(lnV_new)-exp(lnV_old);
			double scale_factor = exp((lnV_new-lnV_old)/3.0);
			// Scale positions of particles and box (if the configuration is not overlapping)
			for(unsigned p = 0;p<x.size();p++){
				x.at(p)*=scale_factor;
				y.at(p)*=scale_factor;
		  	    z.at(p)*=scale_factor;
			}
			Lx*=scale_factor;
			Ly*=scale_factor;
			Lz*=scale_factor;

			// Accept move at random
			// implicit use temperature one, i.e. beta=1
			attempts_volume++;
			double arg = -pressure*dV+(number_of_particles()+1)*(lnV_new-lnV_old);
			double randf = dice_0to1(gen);
			if( randf>exp(arg) || is_overlapping() ){ // Reject, and restore old state
			  rejected_volume++;
			  for(unsigned p = 0;p<x.size();p++){
				x.at(p)/=scale_factor;
				y.at(p)/=scale_factor;
				z.at(p)/=scale_factor;
			  }
			  Lx/=scale_factor;
			  Ly/=scale_factor;
			  Lz/=scale_factor;
			}
		}
		
		// Pick a random particle, and make a single particle step
      	unsigned p = dice_particle(gen);
      	//unsigned p = s/x.size();

      	// Random step (in sphere)
      	double r2=2.0;
      	double xold=x.at(p);
	  	double yold=y.at(p);
	  	double zold=z.at(p);
		
		//bool old_is_overlapping = is_overlapping(p);
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
		if( is_overlapping(p) ){
	      rejected++;
		  // Reject move
		  x.at(p)=xold;
		  y.at(p)=yold;
		  z.at(p)=zold;
		}
	  }
	}
	cout << "Tested " << overlap_tests << " configutations and found " << overlap_hits_before << " before neighbour list updates (NLU)" << endl;
	cout <<         "  and " << overlap_hits_after << " after NLU." << endl;
	cout << "Accepted " << attempts-rejected << " of " << attempts 
	  << " MC move attempts (" << 100.0*(1.0-(double)rejected/(double)attempts) << "%)" << endl;
	cout << "Accepted " << attempts_volume-rejected_volume << " of " << attempts_volume
	  << " volume change attempts (" << 100.0*(1-(double)rejected_volume/(double)attempts_volume) << "%)" << endl;
}

