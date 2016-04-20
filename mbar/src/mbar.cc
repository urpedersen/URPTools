/**
 *  mbar.cpp
 * 
 *  Source file for the program mbar 
 * 
 *  Version: 2.0
 *  Created on: Sep 24, 2010
 *  Modified on: Aug 6, 2015
 *  Author: Ulf R. Pedersen (http://urp.dk)
 *  Copyright: GNU GPL 3.0
 *
 *  Implementation of the MBAR algorithm for reweighing of probability distributions.
 *  Please read [Shirts and Chodera, J. Chem. Phys. 129, 124105 (2008); doi:10.1063/1.2978177] for details on the MBAR algorithm.
 *
 */


#include "mbar.h"
#include "time.h"

using namespace std;

int main(int narg, char **arg) {

	// Construct Mbar object and load data.
	string variables = set_variable_string(narg,arg);
	Mbar mbar(variables);

	// Welcome the kind user.
	if( mbar.verbose>4 ) {
		cout << "mbar v. 2.0 by Ulf R. Pedersen (www.urp.dk)" << endl;
		cout << "  - an implementation of the MBAR algorithm for reweighing of probability distributions." << endl;
		cout << "Please read [Shirts and Chodera, J. Chem. Phys. 129, 124105 (2008); doi:10.1063/1.2978177] for details." << endl << endl;
	}

	// Timing of program
	clock_t start_clock = clock();

	// Load data from disk
	mbar.load_data();

	// Write histograms of raw inputs
	mbar.write_histogram_of_umbrellas();

	// Get free energies by iterations using the MBAR algorithm (this is the heavy part of the program)
	mbar.iterate_free_energies();

	// Write the natural logarithm of probability along observable(s).
	mbar.write_log_probability_of_observables();

	// Write the natural logarithm of probability of the input from umbrellas.
	mbar.write_log_probability_of_umbrellas();
	mbar.write_log_probability_of_umbrella_contributions();

	// Write moments of the computed distribution
	mbar.write_moments();
	mbar.write_central_moments();

	// Say goodbye to the kind user.
	if( mbar.verbose>4 ) { 
		cout << endl << "Variable string: " << variables << endl;
		cout << "Usage time: " << (float)(clock()-start_clock)/(float)CLOCKS_PER_SEC << " secounds." << endl;
		cout << "Happy ending. " << arg[0] <<endl;
	}

}
