//============================================================================
// Name        : traj_Fs.cpp
// Author      : Ulf R. Pedersen
//============================================================================

#include <iostream>

#include "traj.h"

using namespace std;

int main(int narg, char **arg) {

	// Initialize.
	string tmp=set_variable_string(narg,arg);
	Traj traj(tmp);


	// Welcome user.
	if(traj.verbose>4){
		cout << "Calculate mean the self intermediate scattering function F_s(t) for a given q-vector" << endl
			<< "Usage: " << arg[0] << " --num_frames=1" << endl
			<< "See documentation and output for variables that can and should be set with --variable=value flags." << endl;
	}

	// Load trajectory.
	traj.load_data_from_disk();

	// Self intermediate scattering function of particles
	double q_vector = get_variable(traj.variables,"q_vector",6.28318530717959);
	traj.incoherent_self_scattering("Fs.dat",q_vector);

	// Say goodbye to the kind user.
	if(traj.verbose>4) cout << endl << "Variable string: " << traj.variables << endl;
	if(traj.verbose>4) cout << "Happy ending. " << arg[0] <<endl;
}

