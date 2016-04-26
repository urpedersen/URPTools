//============================================================================
// Name        : traj_gr.cpp
// Author      : Ulf R. Pedersen
//               Made in 2010
//============================================================================

#include <iostream>
using namespace std;

#include "traj.h"

int main(int narg, char **arg) {
	// Initialize.
	string tmp=set_variable_string(narg,arg);
	Traj traj(tmp);

	// Welcome user.
	if(traj.verbose>4){
		cout << "Calculate the radial distribution function" << endl
			<< "Usage: " << arg[0] << " --num_frames=1" << endl
			<< "See documentation and output for variables that can and should be set with --variable=value flags." << endl;
	}

	// Load trajectory
	traj.load_data_from_disk();

	// Analysis
	traj.radial_distribution("gr.dat");

	// Say goodbye to the kind user.
	if(traj.verbose>4) cout << endl << "Variable string: " << traj.variables << endl;
	if(traj.verbose>4) cout << "Happy ending. " << arg[0] <<endl;

	return 0;
}
