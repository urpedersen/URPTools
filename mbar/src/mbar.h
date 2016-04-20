/**
 * mbar.h
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
 *  Write reduced potential function of microstate x_n
 *  at macro state i as (later referred to as window i)
 *
 *     u_i(x_n) = beta_i [ U(x_n)+v_i(x_n)]
 *
 *  where beta_i=1/kT_i is the inverse temperature.
 *
 *  The essential functions of this program and the MBAR algorithm
 *  are call in the function iterate_free_energies. They are
 *	recalc_free_energies
 *	recalc_log_sample_weights
 *  that uses
 *	log_sum_exp
 *	
 */

#ifndef MBAR_H_
#define MBAR_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

using namespace std;

string set_variable_string(int, char **); 		/// Let user set variables
int get_variable(string,string,int);
double get_variable(string,string,double);

double log_sum_exp(double *,int,double);		/// Sum exponentials without loss of precision

class Mbar {
public:
	int verbose;					/// Determine how loud the program is (0 for silent and 9 loud)
	string variables;				/// String of variables from user
	string mbar_window_info_filename;		/// File containing information about where windows

	int number_of_window_dimensions;		/// Dimensions of umbrella windows
	int number_of_umbrella_parameters;		/// Number of parameters for the umbrella energy (2 for harmonic)
	double *umbrella_parameters;			/// List of parameters of the umbrella energy

	int number_of_windows;				/// Number of umbrella windows
	int *number_of_samples_in_window;		/// List of number of samples in windows
	double *log_number_of_samples_in_window; 	/// ... and log of that
	int *first_sample_in_window;			/// First index of sample in given window
	string *window_filename;			/// File names

	int number_of_samples;				/// number_of_samples
	int sample_stride;				/// Gives the option to skip point on load to save memory and computation time
	double *samples;				/// Value of all elements
	double *log_sample_weights;			/// The logarithm of the sample weights (when there is no umbrella)
 	double *potential_energy;			/// The potential energy of sample TODO this is not implemented, but important if temperatures of windows differ
 
	bool save_energies_to_mem;
	double *harmonic_window_energy_of_sample;

 	double log_sum_exp_tolerance;			/// Tolerance for needed in the log_sum_exp function
	
	double* sample_terms;				/// Arrays for inner-loop in MBAR computations
	double* window_terms;

	double *window_beta;				/// Inverse temperature of window, 1/k_BT, in inverse energy units.
	double *free_energies;				/// free energy of windows

	double *sample_observables;			/// Observable to analyse
	int number_of_observable_dimensions;		/// Dimensions of observable

	Mbar(string);					/// Constructor
	~Mbar();					/// Destructor

	void load_data();				/// Loading data from disk and allocating arrays
	int get_number_of_windows();
	void load_window_variables();
	int load_window_metadata(int,string);
	void load_window(int);

	void iterate_free_energies();			// Perform iterations to get free energies
	void recalc_free_energies();
	void recalc_free_energies_normalize();
	void recalc_log_sample_weights();
	double log_sum_of_sample_weights();
	double get_harmonic_window_energy_of_sample(int,int);
	void set_harmonic_window_energy_of_sample();

	void write_log_probability_of_observables();	// For analysis
	void write_histogram_of_umbrellas();
	void write_log_probability_of_umbrella_contributions();
	void write_log_probability_of_umbrellas();
	void write_moments();
	void write_central_moments();
};









/**
 *   Return a ''command line string'' with variables.
 *   This is send to the get_variable functions
 *
 *  @param variable_string
 *  @param variable_name
 *  @param default_value
 */
string set_variable_string(int narg, char **arg){
	string variables="";
	for(int i=0;i<narg;i++){
		variables.append(arg[i]);
		variables.append(" ");
	}
	variables.append(" ");

	return variables;
}







/**
*   Return variable value given in variable_string.
*      ... this could look some like "var1=146 var2=2e5 and a mandatory comment var3=-200"
*      If the variable is set more than ones, the first one is used.
*
* @param variable_string
* @param variable_name
* @param default_value
*/
int get_variable(string variable_string,string variable_name,int default_value){
	int verbose=5;

	size_t pos;
	variable_name.insert(0," --");   // make sure that variable name is surrounded by " " and "="
	variable_name.append("=");
	pos = variable_string.find(variable_name,0);

	if(pos==-1){ // No match. NB: comparing signed and unsigned integer
		if(verbose>4) cout << "Warning: Missing variable. Use default"<< variable_name << default_value << "\n" ;
	}
	else{
		string tmp_str=variable_string.substr(pos + variable_name.size());
		pos = tmp_str.find(" ");
		tmp_str=tmp_str.substr(0,pos);
		default_value=atol(tmp_str.c_str());
		if(verbose>4) cout << "Use "<< variable_name << default_value << "\n" ;
	}
	return default_value;
}






/**
*   Same as int get_variable(string variable_string,string variable_name,int default_value)
*     but returns a double
* @param variable_string
* @param variable_name
* @param default_value
*/
double get_variable(string variable_string,string variable_name,double default_value){
	//int verbose=DEFAULT_VERBOSE;
	int verbose=5;

	size_t pos;
	variable_name.insert(0," --");   // make sure that variable name is surrounded by white space
	variable_name.append("=");     //   and "="
	pos = variable_string.find(variable_name,0);

	if(pos==-1){  // No match. NB: comparing signed and unsigned integer
		if(verbose>4) cout << "Warning: Missing variable. Use default "<< variable_name << default_value << " \n" ;
	}
	else{
		string tmp_str=variable_string.substr(pos + variable_name.size());
		pos = tmp_str.find(" ");
		tmp_str=tmp_str.substr(0,pos);
		default_value=atof(tmp_str.c_str());
		if(verbose>4) cout << "Use "<< variable_name << default_value << " \n";
	}
	return default_value;
}








/**
 * Returns log of sum of exponentials
 *
 * ln( Sum{ exp(s_i) }  )
 *
 * and num_lnTerms is the number of terms
 *
 * @param double* lnTerms   list of s_i's ()
 * @param int default_value   length of lnTerms
 */
double log_sum_exp(double lnTerms[],int num_lnTerms,double log_sum_exp_tolerance)
{
	// Find the largest term
	double maxTerm = lnTerms[0];
	for(int i=0;i<num_lnTerms;i++)
		if(lnTerms[i]>maxTerm) maxTerm=lnTerms[i];

	// Make sum of exp(s_i)/exp(s_max)
	double threshold = maxTerm-log(num_lnTerms)+log(log_sum_exp_tolerance);
	double result = 0.0;
	for(int i=0;i<num_lnTerms;i++)
		if(lnTerms[i]>=threshold)
			result += exp( lnTerms[i] - maxTerm );
		
	result = log( result );
	result += maxTerm;

	return result;
}





/**
 * (Standard) Constructor
 * Input: variables_in is a string with variables as taken from the commandline
 */
Mbar::Mbar(string variables_in){

	// Set variables string
	variables=variables_in;
	verbose = get_variable(variables,"verbose",5);

}



/**
*  Destructor
*/
Mbar::~Mbar(){
	// Free memory
	delete[] samples;
}



/**
 * Load data from files and allocate (all) variables
 */
void Mbar::load_data(){

	mbar_window_info_filename="mbar_window_info.txt"; // TODO read from variables string

	number_of_windows=get_number_of_windows();

	number_of_samples_in_window = new int[number_of_windows];
	log_number_of_samples_in_window = new double[number_of_windows];
	first_sample_in_window = new int[number_of_windows];

	number_of_window_dimensions=get_variable(variables,"number_of_window_dimensions",1);

	number_of_umbrella_parameters=2;  // <-- This program can only do 2/kappa*(X-X_0)^2 umbrellas thus 2 parameters

	umbrella_parameters = new double[number_of_window_dimensions
	                                 *number_of_umbrella_parameters
	                                 *number_of_windows];

	window_beta = new double[number_of_windows];  // <-- TODO I'm not sure that this program will give right answer if beta's are not the same

	window_filename = new string[number_of_windows];

	free_energies = new double[number_of_windows];

	// Nothing is loaded yet, thus 0.
	number_of_samples=0;

	// Set the stride when loading samles
	sample_stride = get_variable(variables,"sample_stride",1);

	// Load and set (above) window variables
	load_window_variables();

	// Allocate arrays for inner-loop terms 
	sample_terms = new double[number_of_samples];
	window_terms = new double[number_of_windows];

	// Allocate array with all elements, and load them
	samples = new double[number_of_samples*number_of_window_dimensions];

	number_of_observable_dimensions = get_variable(variables,"number_of_observable_dimensions",1);
	sample_observables = new double[number_of_samples*number_of_observable_dimensions];

	if(verbose>4) cout << "number_of_samples=" << number_of_samples << " number_of_observable_dimensions=" << number_of_observable_dimensions << endl;

	for(int window=0;window<number_of_windows;window++){
		load_window(window);
	}

	// Allocate sample weights
	log_sample_weights = new double[number_of_samples];

	// Allow the user to save window energies to memory instead of computing them every time
	if(verbose>4) cout << "If save_energies_to_mem=1 then window_energy_of_sample are saved to memory, " << endl
			   << "   the default is that they are computed every time (saves memory). " << endl ;
	int LOAD_save_energies_to_mem = get_variable(variables,"save_energies_to_mem",0);
	save_energies_to_mem = (LOAD_save_energies_to_mem==1);
	if(save_energies_to_mem){
		harmonic_window_energy_of_sample = new double[number_of_samples*number_of_windows];
		set_harmonic_window_energy_of_sample();	// Make list with umbrella energies of samples 
	} else {
		harmonic_window_energy_of_sample = NULL;
	}

	
	


	log_sum_exp_tolerance = get_variable(variables,"log_sum_exp_tolerance",1e-18);
}






/**
 * Return and set number_of_windows by counting the number of lines in the mbar_window_info_filename input file.
 */
int Mbar::get_number_of_windows(){
	ifstream info_file;
	number_of_windows=0;
	info_file.open(mbar_window_info_filename.c_str());
	if (info_file.is_open()) {
		string line;
		getline(info_file,line);
		while(! info_file.eof() ){
			number_of_windows++;
			getline(info_file,line);
		}
	}else{
		cout << "Error: Unable to open " << mbar_window_info_filename << ". Exit." << endl;
		exit(0);
	}
	info_file.close();
	return number_of_windows;
}



/**
 * Return and set number_of_windows by counting the number of lines in the mbar_window_info_filename input file.
 */
void Mbar::load_window_variables(){

	if(verbose>4) cout << endl << " ..:: Read meta-data of windows (reading " << mbar_window_info_filename << ") ::.." << endl;
	if(verbose>6) cout << "Lines should look like this: datafile.dat 2.5 36 1.2" << endl;
	if(verbose>6) cout << "  where the inverse temperature is bata=2.5, and x0=36 kappa=1.2 are parameters of umbrella potential, 0.5*kappa(x-x0)^2." << endl;
	if(verbose>4) cout << "  number_of_windows=" << number_of_windows << " number_of_window_dimensions=" << number_of_window_dimensions << endl;

	ifstream info_file;
	info_file.open(mbar_window_info_filename.c_str());
	if (info_file.is_open()) {
		for(int window=0;window<number_of_windows;window++){
			string line;
			getline(info_file,line);
			size_t pos;
			pos = line.find(" ");
			window_filename[window]=line.substr(0,pos);
			number_of_samples+=load_window_metadata(window,window_filename[window]);

			line=line.substr(pos+1);
			pos = line.find(" ");
			window_beta[window]=atof(line.substr(0,pos).c_str());
			if(verbose>4) cout << "window " << window << ": " << "window_filename=" << window_filename[window] << " beta=" << window_beta[window];
			for(int dim=0;dim<number_of_window_dimensions;dim++){
				if(verbose>4) cout << endl << "   Umbrella dimension " << dim << ": ";
				for(int parameter=0;parameter<number_of_umbrella_parameters;parameter++){
					line=line.substr(pos+1);
					pos = line.find(" ");
					int ndx = window*number_of_window_dimensions*number_of_umbrella_parameters+dim*number_of_umbrella_parameters+parameter;
					umbrella_parameters[ndx]=atof(line.substr(0,pos).c_str());
					if(verbose>4) cout << " umbrella_parameters[" << ndx << "]=" << umbrella_parameters[ndx];
				}
			}
			if(verbose>4) cout << endl << "   number_of_samples_in_window[" << window << "]=" << number_of_samples_in_window[window] << " first_sample_in_window[" << window << "]=" << first_sample_in_window[window] << endl;
		}
	}else{
		cout << "Error: Unable to open " << mbar_window_info_filename << ". Exit." << endl;
		exit(0);
	}
	info_file.close();
}


/**
 * Set first_element_in_window number_of_elements_in_window
 */
int Mbar::load_window_metadata(int window,string filename){

	// Set list of positions in sample array
	if(window==0){
		first_sample_in_window[window]=0;
	}else{
		first_sample_in_window[window]=first_sample_in_window[window-1]+number_of_samples_in_window[window-1];
	}
	number_of_samples_in_window[window]=0;

	// Read number of lines in input file
	ifstream ifile;
	ifile.open(filename.c_str());
	if (ifile.is_open()) {
		string line;
		int line_count=0;

		getline(ifile,line);
		while(! ifile.eof() ){
			if(line_count%sample_stride==0)	number_of_samples_in_window[window]++;
			getline(ifile,line);
			line_count++;
		}
	}else{
		cout << "Error: Unable to open " << filename << ". Exit." << endl;
		exit(0);
	}
	ifile.close();

	log_number_of_samples_in_window[window]=log(number_of_samples_in_window[window]);
	return number_of_samples_in_window[window];
}









/**
 * Load window from window_filename[window]
 *
 * One line for each sample, and one field per dimension, separated by a single white space.
 * Thus, it could look something like this for three samples in two dimensions:
 *    34 53
 *    66 23
 *    23 53
 */
void Mbar::load_window(int window){
	ifstream ifile;
	ifile.open(window_filename[window].c_str());
	if (ifile.is_open()) {
		string line;

		for(int sample=first_sample_in_window[window];sample<first_sample_in_window[window]+number_of_samples_in_window[window];sample++){
			getline(ifile,line);
			if(verbose>8) cout << "line: " << line << endl;

			// Load sample umbrella variables
			for(int dim=0;dim<number_of_window_dimensions;dim++){
				size_t pos = line.find(" ");
				samples[sample*number_of_window_dimensions+dim]=atof(line.substr(0,pos).c_str());
				if(verbose>8) cout << "samples[" << sample*number_of_window_dimensions+dim << "]=" << samples[sample*number_of_window_dimensions+dim] << endl;
				line=line.substr(pos+1);
			}

			// Load sample observable
			for(int dim=0;dim<number_of_observable_dimensions;dim++){
				size_t pos = line.find(" ");
				sample_observables[sample*number_of_observable_dimensions+dim]=atof(line.substr(0,pos).c_str());
				if(verbose>8) cout << "sample_observables[" << sample*number_of_observable_dimensions+dim << "]=" << sample_observables[sample*number_of_observable_dimensions+dim] << endl;
				line=line.substr(pos+1);
			}

			// To stride of input data
			for(int step=0;step<sample_stride-1;step++) getline(ifile,line);
		}
	}else{
		cout << "Error: Unable to open " << window_filename[window] << ". Exit." << endl;
		exit(0);
	}

	// Print the first line of the loaded data
	if(verbose>4){
		cout << "First sample in window " << window << ":" ;
		for(int dim=0;dim<number_of_window_dimensions;dim++) cout << " " << "umbrella_value(dim=" << dim << ")=" << samples[first_sample_in_window[window]*number_of_window_dimensions+dim];
		for(int dim=0;dim<number_of_observable_dimensions;dim++) cout << " " << "observable(dim=" << dim << ")=" << sample_observables[first_sample_in_window[window]*number_of_observable_dimensions+dim];
		cout << endl;
	}

}








/**
 * Perform iterations of recalc_free_energies() and recalc_log_sample_weights()
 * until change of all free energies is smaller that threshold
 *
 * TODO Now this algorithm scales as N*K^2 where K is number of windows. However, if sums are only performed over overlapping window, this the scaling would be NK.
 *
 */
void Mbar::iterate_free_energies(){
	if(verbose>4) cout << endl << " ..:: Perform iterations to calculate free_energies of umbrella windows ::.." << endl;

	double threshold = get_variable(variables,"threshold",1e-4);

	double *free_energies_old = new double[number_of_windows];
	double largest_change = threshold+1.0;

	// Timing of iterations
	clock_t start_clock = clock();	

	// Try to load file free_energies.dat with free_energies as initial guesses.
	// The free energy of i'th window should be on i'th line
	// Otherwise, set initial guesses to zero
	ifstream ifile;
	ifile.open("free_energies.dat");
	if (ifile.is_open()) {
		string line;
		for(int i=0;i<number_of_windows;i++){
			getline(ifile,line);
			free_energies[i]=atof(line.c_str());
		}
	}else{
		for(int i=0;i<number_of_windows;i++) free_energies[i]=0.0;
	}
	ifile.close();
	if(verbose>4){
		cout << "Initial free energies: f =";
		for(int i=0;i<number_of_windows;i++) cout << " " << free_energies[i];
		cout << endl;
	}
	recalc_log_sample_weights();

	// Do iterations until convergence
	int counter = 0;
	int next_write_out=1;

	while( largest_change > threshold ) {
		counter++;
		for(int i=0;i<number_of_windows;i++) free_energies_old[i]=free_energies[i];

		// Do one iteration
		recalc_free_energies();
		recalc_log_sample_weights();

		// Calculate largest change
		largest_change=sqrt((free_energies_old[0]-free_energies[0])*(free_energies_old[0]-free_energies[0]));
		for(int i=0;i<number_of_windows;i++){
			double change = sqrt((free_energies_old[i]-free_energies[i])*(free_energies_old[i]-free_energies[i]));
			if(change>largest_change) largest_change=change;
		}

		// Inform user about progress
		if(verbose>4 && (counter==next_write_out || largest_change < threshold)){
			next_write_out*=2;
			cout << "Iteration " << counter << ": largest_change=" << largest_change << " clock=" << (float)(clock()-start_clock)/(float)CLOCKS_PER_SEC << endl;
			cout << "f =";
			for(int i=0;i<number_of_windows;i++) cout << " " << free_energies[i];
			cout << endl;
			cout << "f/beta =";
			for(int i=0;i<number_of_windows;i++) cout << " " << free_energies[i]/window_beta[i];
			cout << endl;

			// Make a restart file with free energies
			FILE * restartfile = fopen("restart_free_energies.dat","w");
			if(restartfile!=NULL) {
				for(int i=0;i<number_of_windows;i++) fprintf (restartfile, "%16.16f\n",free_energies[i]);
				fclose(restartfile);
			}
		}
	}

	// Print final free_energies
	//recalc_free_energies_normalize();
	if(verbose>4){
		cout << "Final free energies:";
		for(int i=0;i<number_of_windows;i++) cout << " " << free_energies[i];
		cout << endl;
	}

	// Write free energies to ASCII file
	FILE * ofile = fopen("free_energies.dat","w");
	if(ofile!=NULL) {
		for(int i=0;i<number_of_windows;i++) fprintf (ofile, "%16.16f\n",free_energies[i]);
		fclose(ofile);
	}
}






/**
 * Recalculate free energies using equation (11) in
 *    [Shirts & Chodera, J. Chem. Phys. 129, 123105 (2008)]
 * That is, iterating this with recalc_log_sample_weights will yield self consistent values of free energies.
 */
void Mbar::recalc_free_energies(){

	// double* sample_terms;
	// sample_terms = new double[number_of_samples];

	int jn;
	for(int i=0;i<number_of_windows;i++){
		for(int j=0;j<number_of_windows;j++){
				if(save_energies_to_mem){
					for(int n=0;n<number_of_samples_in_window[j];n++) {
						jn=first_sample_in_window[j]+n;
						sample_terms[jn]=-harmonic_window_energy_of_sample[jn*number_of_windows+i]+log_sample_weights[jn];
					}
				} else {
					for(int n=0;n<number_of_samples_in_window[j];n++) {
						jn=first_sample_in_window[j]+n;
						sample_terms[jn]=-get_harmonic_window_energy_of_sample(jn,i)+log_sample_weights[jn];
					}
					
				}
		}
		free_energies[i]=-log_sum_exp(sample_terms,number_of_samples,log_sum_exp_tolerance);
	}

	recalc_free_energies_normalize();
	
}


/**
 * Normalize distribution by shifting free energies
 */
void Mbar::recalc_free_energies_normalize(){
	double shift=log_sum_of_sample_weights();
	for(int i=0;i<number_of_windows;i++) free_energies[i]+=shift;
}









/**
 * Recalculate sample weights, sample_weight[sample]
 *
 */
void Mbar::recalc_log_sample_weights(){

	//double* terms;
	//terms = new double[number_of_windows];

	// Sum over all samples, though a double sum over windows and samples in windows
	for(int j=0;j<number_of_windows;j++){
		for(int n=0;n<number_of_samples_in_window[j];n++){
			int jn=first_sample_in_window[j]+n;		// index of sample

			for(int k=0;k<number_of_windows;k++){
				if(save_energies_to_mem){
					window_terms[k]=log_number_of_samples_in_window[k]+free_energies[k]-harmonic_window_energy_of_sample[jn*number_of_windows+k];
				} else {
					window_terms[k]=log_number_of_samples_in_window[k]+free_energies[k]-get_harmonic_window_energy_of_sample(jn,k);
				}
			}
			log_sample_weights[jn]=-log_sum_exp(window_terms,number_of_windows,log_sum_exp_tolerance);
		}
	}
	
	//delete[] terms;
}




/**
 * Return the natural logarithm of sum of log_sample_weights.
 * When distribution is normalizes, this function return zero.
 */
double Mbar::log_sum_of_sample_weights(){
	//double* terms = new double[number_of_samples];
	for(int jn=0;jn<number_of_samples;jn++) sample_terms[jn]=log_sample_weights[jn];
	double out = log_sum_exp(sample_terms,number_of_samples,log_sum_exp_tolerance);
	//delete[] terms;
	return out;
}






/**
 * Return the reduced harmonic energy of sample in window.
 */
double Mbar::get_harmonic_window_energy_of_sample(int jn,int window){
	double energy=0.0;

	// Print a warning to user
	if( ( ! number_of_umbrella_parameters==2 ) && verbose>4) cout << "Warning: number_of_umbrella_parameters=" << number_of_umbrella_parameters << " but expected 2 in get_harmonic_window_energy_of_sample(int,int)" << endl;

	for(int dim=0;dim<number_of_window_dimensions;dim++){
		double a=    umbrella_parameters[window*number_of_window_dimensions*number_of_umbrella_parameters+dim*number_of_umbrella_parameters+0];
		double kappa=umbrella_parameters[window*number_of_window_dimensions*number_of_umbrella_parameters+dim*number_of_umbrella_parameters+1];
		double distance=samples[jn*number_of_window_dimensions+dim]-a;
		energy+=0.5*kappa*distance*distance;
	}
	return window_beta[window]*energy;
}

/*
 * Set array with energy of sample in a given window.
 */
 void Mbar::set_harmonic_window_energy_of_sample(){
	 for(int jn=0;jn<number_of_samples;jn++)
		 for(int k=0;k<number_of_windows;k++)
			 harmonic_window_energy_of_sample[jn*number_of_windows+k]=get_harmonic_window_energy_of_sample(jn,k);
 }








/**
 * Write natural logarithm of probability distribution of observable log( P( observable ) ) to file
 */
void Mbar::write_log_probability_of_observables(){

	// Allocate arrays with bin parameters
	double *bin_min = new double[number_of_observable_dimensions];
	double *bin_width = new double[number_of_observable_dimensions];
	int *bin_num = new int[number_of_observable_dimensions];
	int bin_num_total = 1;

	// Get bin limits from user
	for(int dim=0;dim<number_of_observable_dimensions;dim++){

		char *dim_str = new char[16];
		sprintf(dim_str,"%d",dim);

		string str;

		str="bin_min";
		str.append(dim_str);
		bin_min[dim] = get_variable(variables,str,0.0);

		str="bin_width";
		str.append(dim_str);
		bin_width[dim] = get_variable(variables,str,1.0);

		str="bin_num";
		str.append(dim_str);
		bin_num[dim] = get_variable(variables,str,100);

		bin_num_total*=bin_num[dim];
	}
	if(verbose>4) cout << "bin_num_total=" << bin_num_total << endl;


	// Determine what bin samples are in (put -1 if sample is outside histogram).
	// Note, if binX is the bin in the X'th dimension, then
	// bin = bin0 + bin1*num_bin0 + bin2*num_bin1*num_bin0 ...
	int *bin_of_sample = new int[number_of_samples];
	bool a_sample_is_outside_bins = false;
	for(int jn=0;jn<number_of_samples;jn++){

		bin_of_sample[jn]=0;
		int num_bin_product=1;
		bool outside=false;
		for(int dim=0;dim<number_of_observable_dimensions;dim++){
			int bin_in_this_dimention=(int)floor((sample_observables[jn*number_of_observable_dimensions+dim]-bin_min[dim])/bin_width[dim]+0.5);
			if( bin_in_this_dimention<0 || bin_in_this_dimention>bin_num[dim]-1) outside=true;
			bin_of_sample[jn]+=bin_in_this_dimention*num_bin_product;
			num_bin_product*=bin_num[dim];
		}
		if(outside) bin_of_sample[jn]=-1;
		if(outside) a_sample_is_outside_bins=true;
	}
	if(a_sample_is_outside_bins && verbose>4) cout << "Warning: One or more samples are outside bins. It is recommended to increase limits of bins." << endl;


	// Write distribution to log_probability.dat
	FILE * ofile = fopen("log_probability.dat","w");
	// TODO write below header so that it is specific for the number of dimensions
	fprintf (ofile, "# [central value in dimension 0 ; central value in dimension 1 ; ... ; natural logarithm of probability ; global bin index ; bin index in dimension 0 ; bin index in dimension 1 ; ... ;  number of samples in bin ]\n");
	for(int bin=0;bin<bin_num_total;bin++){

		int num_terms=0;
		for(int jn=0;jn<number_of_samples;jn++) if(bin_of_sample[jn]==bin) num_terms++;

		if(num_terms>0){

			// Retrieve/print central bin value
			int tmp_int=bin;
			for(int dim=0;dim<number_of_observable_dimensions;dim++){
				int bin_in_this_dimention=tmp_int%bin_num[dim];
				tmp_int-=bin_in_this_dimention;
				tmp_int/=bin_num[dim];
				//fprintf (ofile, " %d", bin_in_this_dimention);
				fprintf (ofile, "%f ", bin_min[dim]+bin_width[dim]*(double)bin_in_this_dimention);
			}

			double *terms = new double[num_terms];
			int counter=0;

			for(int jn=0;jn<number_of_samples;jn++) if(bin_of_sample[jn]==bin) terms[counter++]=log_sample_weights[jn];

			fprintf (ofile, "%16.16f ", log_sum_exp(terms,num_terms,log_sum_exp_tolerance));

			// Retrieve/print central bin indexes
			fprintf (ofile, "%d", bin);
			tmp_int=bin;
			for(int dim=0;dim<number_of_observable_dimensions;dim++){
				int bin_in_this_dimention=tmp_int%bin_num[dim];
				tmp_int-=bin_in_this_dimention;
				tmp_int/=bin_num[dim];
				fprintf (ofile, " %d", bin_in_this_dimention);
				//fprintf (ofile, " %f", bin_min[dim]+bin_width[dim]*(double)bin_in_this_dimention);
			}


			fprintf (ofile, " %d\n", num_terms);

			delete [] terms;
		}else{
			// No sampling in this bin. Do nothing.
		}
	}
}







/**
 * Calculate and print histogram of umbrellas
 */
void Mbar::write_histogram_of_umbrellas(){

	// Allocate, e.i. arrays with bin parameters
	double *bin_min = new double[number_of_observable_dimensions];
	double *bin_width = new double[number_of_observable_dimensions];
	int *bin_num = new int[number_of_observable_dimensions];
	int bin_num_total = 1;

	// Get bin limits from user
	for(int dim=0;dim<number_of_observable_dimensions;dim++){

		char *dim_str = new char[16];
		sprintf(dim_str,"%d",dim);

		string str;

		str="bin_min";
		str.append(dim_str);
		bin_min[dim] = get_variable(variables,str,0.0);

		str="bin_width";
		str.append(dim_str);
		bin_width[dim] = get_variable(variables,str,1.0);

		str="bin_num";
		str.append(dim_str);
		bin_num[dim] = get_variable(variables,str,100);

		bin_num_total*=bin_num[dim];
	}
	if(verbose>4) cout << "bin_num_total=" << bin_num_total << endl;

	// Allocate histogram
	int *histogram = new int[bin_num_total];

	// Output file
	FILE * ofile = fopen("umbrella_histograms.dat","w");
	fprintf (ofile, "# [central value in dimension 0 ; central value in dimension 1 ; ... ; counts in bin ; global bin index ; bin index in dimension 0 ; bin index in dimension 0 ; ... ]\n");

	// Loop of windows and write histograms after each other
	int jn = 0;
	for(int j=0;j<number_of_windows;j++){
		fprintf (ofile, "# Window %d \n",j);

		// Reset arrays
		for(int bin=0;bin<bin_num_total;bin++) histogram[bin]=0;

		bool a_sample_is_outside_bins = false;

		for(int n=0;n<number_of_samples_in_window[j];n++){

			// Determine what bin samples are in.
			// Note, if binX is the bin in the X'th dimension, then
			// bin = bin0 + bin1*num_bin0 + bin2*num_bin1*num_bin0 ... and so on
			// Shift variable is used in above as product of num_binX'es calculation.
			int bin=0;
			int num_bin_product=1;
			bool outside=false;
			for(int dim=0;dim<number_of_observable_dimensions;dim++){
				int bin_in_this_dimention=(int)floor((sample_observables[jn*number_of_observable_dimensions+dim]-bin_min[dim])/bin_width[dim]+0.5);
				if( bin_in_this_dimention<0 || bin_in_this_dimention>bin_num[dim]-1) outside=true;
				bin+=bin_in_this_dimention*num_bin_product;
				num_bin_product*=bin_num[dim];
			}

			if(outside) a_sample_is_outside_bins=true;

			// Add to histogram if not ouside the bin range (then the stored bin value makes no sense)
			if(!outside) histogram[bin]++;

			// Increase jn double loop counter.
			jn++;
		}

		if(a_sample_is_outside_bins && verbose>4) cout << "Warning: In window " << j << ", one or more samples were outside bin range. It is recommended to increase limits of bins." << endl;

		// Write to file
		for(int bin=0;bin<bin_num_total;bin++){

			if(histogram[bin]>0){

				// Retrieve/print  central bin value
				int tmp_int=bin;
				for(int dim=0;dim<number_of_observable_dimensions;dim++){
					int bin_in_this_dimention=tmp_int%bin_num[dim];
					tmp_int-=bin_in_this_dimention;
					tmp_int/=bin_num[dim];
					fprintf (ofile, "%f ", bin_min[dim]+bin_width[dim]*(double)bin_in_this_dimention);
				}

				// Print bin value
				fprintf (ofile, "%d ", histogram[bin]);

				// Print bin index
				fprintf (ofile, "%d", bin);

				// Retrieve/print bin index in various dimensions
				tmp_int=bin;
				for(int dim=0;dim<number_of_observable_dimensions;dim++){
					int bin_in_this_dimention=tmp_int%bin_num[dim];
					tmp_int-=bin_in_this_dimention;
					tmp_int/=bin_num[dim];
					fprintf (ofile, " %d", bin_in_this_dimention);
					//fprintf (ofile, " %f", bin_min[dim]+bin_width[dim]*(double)bin_in_this_dimention);
				}

				fprintf (ofile, "\n", histogram[bin]);

			}
		}

		fprintf (ofile, " \n");
	}

	fclose(ofile);
}



/**
 * Write natural logarithm of probability distribution of observable log( P( observable ) ) to file
 *
 *  TODO document this function (or delete)
 *
 */
void Mbar::write_log_probability_of_umbrella_contributions(){

	// Allocate arrays with bin parameters
	double *bin_min = new double[number_of_observable_dimensions];
	double *bin_width = new double[number_of_observable_dimensions];
	int *bin_num = new int[number_of_observable_dimensions];
	int bin_num_total = 1;

	// Get bin limits from user
	for(int dim=0;dim<number_of_observable_dimensions;dim++){

		char *dim_str = new char[16];
		sprintf(dim_str,"%d",dim);

		string str;

		str="bin_min";
		str.append(dim_str);
		bin_min[dim] = get_variable(variables,str,0.0);

		str="bin_width";
		str.append(dim_str);
		bin_width[dim] = get_variable(variables,str,1.0);

		str="bin_num";
		str.append(dim_str);
		bin_num[dim] = get_variable(variables,str,100);

		bin_num_total*=bin_num[dim];
	}
	if(verbose>4) cout << "bin_num_total=" << bin_num_total << endl;


	// Determine what bin samples are in (put -1 if sample is outside histogram).
	// Note, if binX is the bin in the X'th dimension, then
	// bin = bin0 + bin1*num_bin0 + bin2*num_bin1*num_bin0 ...
	int *bin_of_sample = new int[number_of_samples];
	bool a_sample_is_outside_bins = false;

	for(int jn=0;jn<number_of_samples;jn++){

		bin_of_sample[jn]=0;
		int num_bin_product=1;
		bool outside=false;
		for(int dim=0;dim<number_of_observable_dimensions;dim++){
			int bin_in_this_dimention=(int)floor((sample_observables[jn*number_of_observable_dimensions+dim]-bin_min[dim])/bin_width[dim]+0.5);
			if( bin_in_this_dimention<0 || bin_in_this_dimention>bin_num[dim]-1) outside=true;
			bin_of_sample[jn]+=bin_in_this_dimention*num_bin_product;
			num_bin_product*=bin_num[dim];
		}
		if(outside) bin_of_sample[jn]=-1;
		if(outside) a_sample_is_outside_bins=true;
	}
	if(a_sample_is_outside_bins && verbose>4) cout << "Warning: One or more samples are outside bins. It is recommended to increase limits of bins." << endl;


	// Write result to file
	FILE * ofile = fopen("umbrella_log_probability_contributions.dat","w");
	// TODO write below header so that it is specific for the number of dimensions
	fprintf (ofile, "# [central value in dimension 0 ; central value in dimension 1 ; ... ; natural logarithm of probability ; global bin index ; bin index in dimension 0 ; bin index in dimension 1 ; ... ;  number of samples in bin ]");

	int jn = 0;
	for(int j=0;j<number_of_windows;j++){

		// Write distribution to log_probability.dat
		fprintf (ofile, "\n# Window %d \n",j);

		for(int bin=0;bin<bin_num_total;bin++){

			int num_terms=0;
			for(int jn=first_sample_in_window[j];jn<first_sample_in_window[j]+number_of_samples_in_window[j];jn++){
				if(bin_of_sample[jn]==bin) num_terms++;
			}

			if(num_terms>0){

				// Retrieve/print central bin value
				int tmp_int=bin;
				for(int dim=0;dim<number_of_observable_dimensions;dim++){
					int bin_in_this_dimention=tmp_int%bin_num[dim];
					tmp_int-=bin_in_this_dimention;
					tmp_int/=bin_num[dim];
					//fprintf (ofile, " %d", bin_in_this_dimention);
					fprintf (ofile, "%f ", bin_min[dim]+bin_width[dim]*(double)bin_in_this_dimention);
				}

				double *terms = new double[num_terms];
				int counter=0;

				for(int jn=first_sample_in_window[j];jn<first_sample_in_window[j]+number_of_samples_in_window[j];jn++){
					if(bin_of_sample[jn]==bin) terms[counter++]=log_sample_weights[jn];
				}

				fprintf (ofile, "%16.16f ", log_sum_exp(terms,num_terms,log_sum_exp_tolerance));

				// Retrieve/print central bin indexes
				fprintf (ofile, "%d", bin);
				tmp_int=bin;
				for(int dim=0;dim<number_of_observable_dimensions;dim++){
					int bin_in_this_dimention=tmp_int%bin_num[dim];
					tmp_int-=bin_in_this_dimention;
					tmp_int/=bin_num[dim];
					fprintf (ofile, " %d", bin_in_this_dimention);
					//fprintf (ofile, " %f", bin_min[dim]+bin_width[dim]*(double)bin_in_this_dimention);
				}

				fprintf (ofile, " %d\n", num_terms);

				delete [] terms;
			}else{
				// No sampling in this bin. Do nothing.
			}
		}
	}

	fclose (ofile ) ;
}






/**
 * Calculate and print histogram of umbrellas
 *
 * TODO the major part of this was copied from Mbar::write_histogram_of_umbrellas(). This part of the code should not be copied but should rather be in a separate function.
 *
 */
void Mbar::write_log_probability_of_umbrellas(){

	// Allocate, e.i. arrays with bin parameters
	double *bin_min = new double[number_of_observable_dimensions];
	double *bin_width = new double[number_of_observable_dimensions];
	int *bin_num = new int[number_of_observable_dimensions];
	int bin_num_total = 1;

	// Get bin limits from user
	for(int dim=0;dim<number_of_observable_dimensions;dim++){

		char *dim_str = new char[16];
		sprintf(dim_str,"%d",dim);

		string str;

		str="bin_min";
		str.append(dim_str);
		bin_min[dim] = get_variable(variables,str,0.0);

		str="bin_width";
		str.append(dim_str);
		bin_width[dim] = get_variable(variables,str,1.0);

		str="bin_num";
		str.append(dim_str);
		bin_num[dim] = get_variable(variables,str,100);

		bin_num_total*=bin_num[dim];
	}
	if(verbose>4) cout << "bin_num_total=" << bin_num_total << endl;

	// Allocate histogram
	int *histogram = new int[bin_num_total];

	// Output file
	FILE * ofile = fopen("umbrella_log_probability.dat","w");
	fprintf (ofile, "# [central value in dimension 0 ; central value in dimension 1 ; ... ; log probability ; global bin index ; bin index in dimension 0 ; bin index in dimension 0 ; ... ]\n");

	// Loop of windows and write histograms after each other
	int jn = 0;
	for(int j=0;j<number_of_windows;j++){
		fprintf (ofile, "# Window %d \n",j);

		// Reset arrays
		for(int bin=0;bin<bin_num_total;bin++) histogram[bin]=0;

		bool a_sample_is_outside_bins = false;

		for(int n=0;n<number_of_samples_in_window[j];n++){

			// Determine what bin samples are in.
			// Note, if binX is the bin in the X'th dimension, then
			// bin = bin0 + bin1*num_bin0 + bin2*num_bin1*num_bin0 ... and so on
			// Shift variable is used in above as product of num_binX'es calculation.
			int bin=0;
			int num_bin_product=1;
			bool outside=false;
			for(int dim=0;dim<number_of_observable_dimensions;dim++){
				int bin_in_this_dimention=(int)floor((sample_observables[jn*number_of_observable_dimensions+dim]-bin_min[dim])/bin_width[dim]+0.5);
				if( bin_in_this_dimention<0 || bin_in_this_dimention>bin_num[dim]-1) outside=true;
				bin+=bin_in_this_dimention*num_bin_product;
				num_bin_product*=bin_num[dim];
			}

			if(outside) a_sample_is_outside_bins=true;

			// Add to histogram if not ouside the bin range (then the stored bin value makes no sense)
			if(!outside) histogram[bin]++;

			// Increase jn double loop counter.
			jn++;
		}

		if(a_sample_is_outside_bins && verbose>4) cout << "Warning: In window " << j << ", one or more samples were outside bin range. It is recommended to increase limits of bins." << endl;

		// Write to file
		for(int bin=0;bin<bin_num_total;bin++){

			if(histogram[bin]>0){

				// Retrieve/print  central bin value
				int tmp_int=bin;
				double umbrella_energy = 0.0;
				for(int dim=0;dim<number_of_observable_dimensions;dim++){
					int bin_in_this_dimention=tmp_int%bin_num[dim];
					tmp_int-=bin_in_this_dimention;
					tmp_int/=bin_num[dim];
					double central_bin_value = bin_min[dim]+bin_width[dim]*(double)bin_in_this_dimention ;
					fprintf (ofile, "%f ",central_bin_value );

					// TODO NOTE here we assume that umbrella energy is spring-like
					int window = j;
					double a=    umbrella_parameters[window*number_of_window_dimensions*number_of_umbrella_parameters+dim*number_of_umbrella_parameters+0];
					double kappa=umbrella_parameters[window*number_of_window_dimensions*number_of_umbrella_parameters+dim*number_of_umbrella_parameters+1];
					double distance=central_bin_value-a;
					umbrella_energy+=0.5*kappa*distance*distance;
				}

				// Print bin value
				//cout << log(2) << endl;
				//fprintf (ofile, "%f ", log((double)histogram[bin]) - free_energies[j] + 0.8*window_beta[j]*umbrella_energy - log(number_of_samples_in_window[j]));
				fprintf (ofile, "%f ", log((double)histogram[bin]) - free_energies[j] + window_beta[j]*umbrella_energy - log(number_of_samples_in_window[j]));

				// Print bin index
				fprintf (ofile, "%d", bin);

				// Retrieve/print bin index in various dimensions
				tmp_int=bin;
				for(int dim=0;dim<number_of_observable_dimensions;dim++){
					int bin_in_this_dimention=tmp_int%bin_num[dim];
					tmp_int-=bin_in_this_dimention;
					tmp_int/=bin_num[dim];
					fprintf (ofile, " %d", bin_in_this_dimention);
					//fprintf (ofile, " %f", bin_min[dim]+bin_width[dim]*(double)bin_in_this_dimention);
				}

				fprintf (ofile, "\n", histogram[bin]);

			}
		}

		fprintf (ofile, " \n");
	}

	fclose(ofile);
}





/**
 * Calculate and print moments
 * Integer number_of_moments indicated the number of moments to be calculated
 *
 *   TODO write also non-central moments
 *
 */
void Mbar::write_moments(){
	if(verbose>4) cout << endl << "   ..:: Moments of distribution ::.." << endl;

	int number_of_moments=9;
	double* moments = new double[number_of_moments*number_of_observable_dimensions];

	//double* terms = new double[number_of_samples];

	//double* central_value = new double[number_of_observable_dimensions];

	// Find largest weight
	double largest_weight = log_sample_weights[0];
	for(int jn=0;jn<number_of_samples;jn++){
		if(log_sample_weights[jn]>largest_weight) largest_weight=log_sample_weights[jn];
	}

	for(int dim=0;dim<number_of_observable_dimensions;dim++){

			// Calculate central value
			//central_value[dim]=0.0;
			//for(int jn=0;jn<number_of_samples;jn++) central_value[dim]+=exp(log_sample_weights[jn]-largest_weight)*sample_observables[jn*number_of_observable_dimensions+dim];
			//central_value[dim]*=exp(largest_weight);
			//if(verbose>4) cout << "Observable dimension " << dim << " with central value, x0, " << central_value[dim] << endl;

			// Calculate higher order moments
			for(int i=0;i<number_of_moments;i++){
				int moment = i;
				moments[i*number_of_observable_dimensions+dim]=0.0;
				for(int jn=0;jn<number_of_samples;jn++) moments[i*number_of_observable_dimensions+dim]+=exp(log_sample_weights[jn]-largest_weight)*pow(sample_observables[jn*number_of_observable_dimensions+dim],moment);
				moments[i*number_of_observable_dimensions+dim]*=exp(largest_weight);
				if(verbose>4) cout <<  "<x^" << moment << "> = " << moments[i*number_of_observable_dimensions+dim] << endl;
			}

			// Non-Gaussian and relative screw parameter"
			if(verbose>4) cout << "Relative screw (m3/m2^(3/2)) in dimension " << dim << " : " << moments[3*number_of_observable_dimensions+dim]/pow(moments[2*number_of_observable_dimensions+dim],3.0/2.0) << endl;
			if(verbose>4) cout << "Non-Gauss parameter (m4/(3*m2^2) - 1) in dimension " << dim << " : " << moments[4*number_of_observable_dimensions+dim]/( 3.0 * moments[2*number_of_observable_dimensions+dim]*moments[2*number_of_observable_dimensions+dim] ) - 1.0 << endl;
			if(verbose>4) cout << "2D |R| Non-Gauss parameter (m4/(2*m2^2) - 1) in dimension " << dim << " : " << moments[4*number_of_observable_dimensions+dim]/( 2.0 * moments[2*number_of_observable_dimensions+dim]*moments[2*number_of_observable_dimensions+dim] ) - 1.0 << endl;
	}

	// Calculate co-variance matrix
	if(verbose>4) cout << endl << "Covariance matrix:" << endl;
	for(int dimA=0;dimA<number_of_observable_dimensions;dimA++){
		for(int dimB=0;dimB<number_of_observable_dimensions;dimB++){
			double cov = 0.0;
			for(int jn=0;jn<number_of_samples;jn++) cov+=exp(log_sample_weights[jn]-largest_weight)*(sample_observables[jn*number_of_observable_dimensions+dimA])*(sample_observables[jn*number_of_observable_dimensions+dimB]);
			if(verbose>4) cout << " " << cov*exp(largest_weight);
		}
		if(verbose>4) cout << endl;
	}
}



/**
 * Calculate and print moments
 * Integer number_of_moments indicated the number of moments to be calculated
 *
 */
void Mbar::write_central_moments(){
	if(verbose>4) cout << endl << "   ..:: Central moments of distribution ::.." << endl;

	int number_of_moments=7;
	double* moments = new double[number_of_moments*number_of_observable_dimensions];

	//double* terms = new double[number_of_samples];

	double* central_value = new double[number_of_observable_dimensions];

	// Find largest weight
	double largest_weight = log_sample_weights[0];
	for(int jn=0;jn<number_of_samples;jn++){
		if(log_sample_weights[jn]>largest_weight) largest_weight=log_sample_weights[jn];
	}

	for(int dim=0;dim<number_of_observable_dimensions;dim++){

			// Calculate central value
			central_value[dim]=0.0;
			for(int jn=0;jn<number_of_samples;jn++) central_value[dim]+=exp(log_sample_weights[jn]-largest_weight)*sample_observables[jn*number_of_observable_dimensions+dim];
			central_value[dim]*=exp(largest_weight);
			if(verbose>4) cout << "Observable dimension " << dim << " with central value, x1 = " << central_value[dim] << endl;

			// Calculate higher order moments
			for(int i=0;i<number_of_moments;i++){
				int moment = i+2;
				moments[i*number_of_observable_dimensions+dim]=0.0;
				for(int jn=0;jn<number_of_samples;jn++) moments[i*number_of_observable_dimensions+dim]+=exp(log_sample_weights[jn]-largest_weight)*pow(sample_observables[jn*number_of_observable_dimensions+dim]-central_value[dim],moment);
				moments[i*number_of_observable_dimensions+dim]*=exp(largest_weight);
				if(verbose>4) cout <<  "<(x-x1)^" << moment << "> = " << moments[i*number_of_observable_dimensions+dim] << " = ( " << moments[i*number_of_observable_dimensions+dim]/pow(moments[dim],0.5*(double)moment) << " )<(x-x1)^2>^(" << moment << "/2)"<< endl;
			}

			// Non-Gaussian and relative screw parameter"
			if(verbose>4) cout << "Relative screw (m3/m2^(3/2)) in dimension " << dim << " : " << moments[(3-2)*number_of_observable_dimensions+dim]/pow(moments[(2-2)*number_of_observable_dimensions+dim],3.0/2.0) << endl;
			if(verbose>4) cout << "Non-Gauss parameter (m4/(3*m2^2) - 1) in dimension " << dim << " : " << moments[(4-2)*number_of_observable_dimensions+dim]/( 3.0 * moments[(2-2)*number_of_observable_dimensions+dim]*moments[(2-2)*number_of_observable_dimensions+dim] ) - 1.0 << endl;
	}

	// Calculate co-variance matrix
	if(verbose>4) cout << endl << "Covariance matrix:" << endl;
	for(int dimA=0;dimA<number_of_observable_dimensions;dimA++){
		for(int dimB=0;dimB<number_of_observable_dimensions;dimB++){
			double cov = 0.0;
			for(int jn=0;jn<number_of_samples;jn++) cov+=exp(log_sample_weights[jn]-largest_weight)*(sample_observables[jn*number_of_observable_dimensions+dimA]-central_value[dimA])*(sample_observables[jn*number_of_observable_dimensions+dimB]-central_value[dimB]);
			if(verbose>4) cout << " " << cov*exp(largest_weight);
		}
		if(verbose>4) cout << endl;
	}
}




#endif /* MBAR_H_ */
