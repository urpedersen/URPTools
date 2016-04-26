/*
 * traj.cpp
 *
 *  Created on: Sep 23, 2010
 *      Author: Ulf R. Pedersen, www.urp.dk
 *
 *  Minimal example of program to analyze trajectories. This program only load trajectory.
 */


#include "traj.h"

int main(int narg, char **arg) {

	// Initialize.
	string tmp=set_variable_string(narg,arg);
	Traj traj(tmp);

	// Welcome user.
	if(traj.verbose>4){
		cout << "Minimal example of program to analyze trajectories. This program only load trajectory." << endl
			<< "Usage: " << arg[0] << " --num_frames=1" << endl
			<< "See documentation and output for variables that can and should be set with --variable=value flags." << endl;
	}

	// Load trajectory.
	traj.load_data_from_disk();

	// PUT STUFF HERE THAT DOES ANALYSIS

	// Say goodbye to the kind user.
	if(traj.verbose>4) cout << endl << "Variable string: " << traj.variables << endl;
	if(traj.verbose>4) cout << "Happy ending. " << arg[0] <<endl;
}

