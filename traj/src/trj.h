/*
 * trj.h
 *  by Ulf R. Pedersen
 *     http://www.urp.dk
 * This is a bunch of helper functions for the Traj Class. 
 * They are, however, not specific to this class and can be used 
*    in other programs/classes for convenience.
 */

#ifndef TRJ_H_
#define TRJ_H_

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
using namespace std;

#ifndef PI
#define PI 3.14159265358979323846	// ratio of a circles circumference to its diameter
#endif

#ifndef NUM_DIM
#define NUM_DIM 3					// number of spatial dimensions
#endif

#ifndef DEFAULT_VERBOSE
#define DEFAULT_VERBOSE 5					// default_value for how verbose functions are (9 for loud, 5 for comfortable and 0 for silent)
#endif

// IO help functions
string set_variable_string(int, char **);
int get_variable(string,string,int);
double get_variable(string,string,double);
// TODO string get_commandline_variable(string,string,string);
double first_x_where_y_equals_y0(double *y,int N,double y0);

// Reading trajectory files
void get_metadata_traj(string ,int *,int *);
int load_traj(string , int ,int , double * , int *, int * , double *);
void get_num_frames_from_user (string , int *);
int make_trajectory_dummy(int,int,double *,int *,int *,double *);
void get_metadata_lammps_atomfile (string ,int *, int *);
int load_lammps_atomfile ( string , int ,int , double * , int *, int * , double * );
void get_metadata_gromacs_grofile (string ,int *, int *);
int load_gromacs_grofile ( string , int ,int , double * , int *, int * , double * );
void get_metadata_xyzfile (string ,int *, int *);
int load_xyzfile( string , int ,int , double * , int *, int * , double * );

// Getting atom related thing from trajectories
double dx_img(int, int, int, int, int, int ,       double*, int*, double*);
double dr_img(int, int, int, int, int,       double* ,int*, double*);
double dx(    int, int, int, int, int, int , int, double*, int*, double*);
double dr(    int, int, int, int, int, int , double*, int*, double*);

// Get molecule related
double get_geometric_center_img(int ,int ,int , int , int , int  , double* , double*);


// Make new data sets
void set_neighbour_list(double *,int,int,int,double *,int *,double *,int ,int *,int *);
void get_vectors(int,double*,double*,double*);

void get_cross_products (	int ,double *,double *,double *);
void get_dot_products (	int ,double *,double *,double *);
void get_cos_thetas (	int ,double *,double *,double *);

double get_dot_product (	double *,double *);
double get_cos_theta (	double *,double *);

void get_mean_square_displacement(int,int,double *, int *, double *,double *,int);
void get_mean_square_displacement_of_selected(int,int,double *, int *, double *,bool *,double *,int);

void get_self_intermediate_scattering_function(int,int,double *, int *, double *,double *,double,int);
void get_self_intermediate_scattering_function_of_selected (int,int ,double * ,int *,double *,bool *,double *,double,int);

void get_radial_distribution_function(int,int,double*,double*,double*,double,int,int);
void get_radial_distribution_function_of_selected(int,int,double*,double*,bool *,double*,double,int,int);

void get_rotational_auto_correlation(int,int,double*,double*,double*,int);
void get_collective_rotational_auto_correlation(int,int,double*,double*,double*,int);

void get_qubatic_orderparameter(int ,int ,double *,double *);






/* return a ''command line string'' with variables. This is send to the get_variable functions
* ... let the user change values of variables
*/
string set_variable_string(int narg, char **arg){
	string variables="";

	//
	for(int i=0;i<narg;i++){
		variables.append(arg[i]);
		variables.append(" ");
	}
	variables.append(" ");

	// Read variables from variables.txt
	string line;
	ifstream vars_file ( "variables.txt" );
	if ( vars_file.is_open() ){
		while ( ! vars_file.eof() ){
			getline ( vars_file , line );

			// Skip stuff after a hash
			line = line.substr(0,line.find("#"));

			// Skip if line is blank
			if(line.size()>0){
				variables.append (" ");
				variables.append ( line );
			}
		}
	}
	vars_file.close();


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
	//int verbose=DEFAULT_VERBOSE;
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
 *	Return values where y(x)=y0
 *	using a linear interpolation between point.
 *	Assuming that x={0, 1, 2, ..., N-1}
 *
 *	Return -1.0 if y(x) does not travel though y0
 */
double first_x_where_y_equals_y0(double *y,int N,double y0){

	for (int i = 0 ; i < N - 1 ; i++ )
		if( ( y[i] >= y0 ) ^ ( y[i+1] > y0 ) )
			return (double)i + ( y0 - y[i] ) / ( y[i+1] - y[i] );

	return -1.0;

}



/**
 *  Return number of particles (ptr_num_frames) and number of frames (ptr_num_atoms) of (any known) trajectory type
 *
 */
void get_metadata_traj(string commandline, int *ptr_num_frames, int *ptr_num_atoms ){
	int verbose=DEFAULT_VERBOSE;

	string filename="unknown_filename";	// name of input file with coordinates
	int load_traj_type=get_variable(commandline,"load_traj_type",1);

	// Print help to user about trajectories that can be used
	if(verbose>4){
		cout << "The load_traj_type=" << load_traj_type << " can take the following values:" << endl;
		cout << "  load_traj_type=0 : Make dummy trajectory (for testing and debugging)." << endl;
		cout << "  load_traj_type=1 : LAMMPS traj.atom file." << endl;
		cout << "     Have something like this in your LAMMPS input file:" << endl;
		cout << "     dump dumpATOM all atom 100 traj.atom" << endl;
		cout << "     dump_modify dumpATOM image yes" << endl;
		cout << "  load_traj_type=2 : GROMACS traj.gro file." << endl;
		cout << "      Use the GROMACS trjconv program like this to get image positions correct:" << endl;
		cout << "      echo 0 | trjconv -o traj.gro -pbc nojump" << endl;
		cout << "  load_traj_type=3 : XYZ (simple) traj.xyz file." << endl;
		cout << "      Note, use the flags bboxX, bboxY and bboxZ to set periodic boundaries." << endl;
	}

	switch (load_traj_type) {
		case 0:
			if(verbose>4) cout << "Do not load a input file (debugging only)" << endl;
			*ptr_num_frames=3;
			*ptr_num_atoms=3;
		break;
		case 1:
			filename="traj.atom";
			if(verbose>4) cout << "Load LAMMPS atom trajectory file (with images) " << filename << endl;
			get_metadata_lammps_atomfile( filename , ptr_num_frames , ptr_num_atoms );
		break;
		case 2:
			filename="traj.gro";
			if(verbose>4) cout << "Load GROMACS *.gro trajectory file " << filename << endl;
			get_metadata_gromacs_grofile( filename , ptr_num_frames , ptr_num_atoms );
		break;
		case 3:
			filename="traj.xyz";
			if(verbose>4) cout << "Load XYZ trajectory file " << filename << endl;
			get_metadata_xyzfile( filename , ptr_num_frames , ptr_num_atoms );
		break;
		default:
			if(verbose>0) cout << "Unknown trajectory type (get_metadata_traj). load_traj_type=" << load_traj_type << ". Exit." << endl;
			exit(0);
		break;
	}

	// Reset number of frames, if requested by user

	int num_frames_from_user=get_variable(commandline,"num_frames",-1);

	if(num_frames_from_user>(*ptr_num_frames) || num_frames_from_user<0){
		// *ptr_num_frames=num_frames_from_user;
		if(verbose>4) cout << "Warning: num_frames=" << *ptr_num_frames << " is not set by user." << endl;
	}else{
		*ptr_num_frames=num_frames_from_user;
		if(verbose>4) cout << "num_frames=" << *ptr_num_frames << " set by user" << endl;
	}
}






/**
 *
 *  Load a configuration(s) of a LAMMPS "atom" ascii output file with images. In LAMMPS, use something like:
 *
 *      dump dumpATOM all atom 5000 traj.atom
 *      dump_modify dumpATOM image yes
 *
 *  Assume same number of atoms in each frames
 *
 *  Return the number of frames (actually) loaded and -1 if some horrifying error happen.
 */
int load_traj (string commandline,
			int num_frames,
			int num_atoms,
			double bbox[],
			int type[],
			int img_coords[],
			double coords[]) {

	int verbose=DEFAULT_VERBOSE;

	// Set load_traj_type and write help to user about the trajectories that can be loaded
	int load_traj_type=get_variable(commandline,"load_traj_type",1);
	string filename="unknown_filename";

	// Uses these default values if not specified in file
	double bboxX = get_variable(commandline,"bboxX",50.0);
	double bboxY = get_variable(commandline,"bboxY",bboxX);
	double bboxZ = get_variable(commandline,"bboxZ",bboxY);

	// Erase old values
	for(int i=0;i<num_frames*num_atoms*NUM_DIM;i++) img_coords[i]=0;
	for(int i=0;i<num_frames*num_atoms*NUM_DIM;i++) coords[i]=0.0;

	// Fill arrays with data from trajectory files files
	int load_traj_out=-1;
	switch (load_traj_type) {
		case 0:
			if(verbose>4) cout << "Construct dummy trajectory. Use bboxX, bboxY and bboxZ to set system size." << endl;
			for(int frame=0;frame<num_frames;frame++){	// set boundary box
				bbox[frame+0]=bboxX;
				bbox[frame+1]=bboxY;
				bbox[frame+2]=bboxZ;
			}
			load_traj_out=make_trajectory_dummy(num_frames,num_atoms,bbox,type,img_coords,coords);
		break;
		case 1:
			filename="traj.atom";
			load_traj_out = load_lammps_atomfile(filename,num_frames,num_atoms,bbox,type,img_coords,coords);
		break;
		case 2:
			filename="traj.gro";
			load_traj_out = load_gromacs_grofile(filename,num_frames,num_atoms,bbox,type,img_coords,coords);
		break;
		case 3:
			filename="traj.xyz";
			// Get info about boundary box from user
			if(verbose>4) cout << "Note, XYZ files have no information about boundary box. We use bboxX=" << bboxX << " bboxY=" << bboxY << " bboxZ=" << bboxZ << endl;
			for(int frame=0;frame<num_frames;frame++){
				bbox[frame*NUM_DIM+0] = bboxX;
				bbox[frame*NUM_DIM+1] = bboxY;
				bbox[frame*NUM_DIM+2] = bboxZ;
			}
			load_traj_out = load_xyzfile(filename,num_frames,num_atoms,bbox,type,img_coords,coords);
		break;
		default:
			if(verbose>0) cout << "Unknown trajectory type. load_traj_type=" << load_traj_type << ". Exit." << endl;
			//exit(0);
		break;
	}

	if(load_traj_out==-1){cout<<"Error: failed to load trajectory. Exit.\n" << endl ;exit (0);}
	if( ! (load_traj_out==num_frames) ){
		if(verbose>0) cout << "Warning: reset num_frames from " << num_frames << " to " << load_traj_out<<endl;
		num_frames=load_traj_out;
	}

	return num_frames;
}


/**
 *  Allow user to set number of frames
 */
void get_num_frames_from_user (string commandline, int *ptr_num_frames){
	int verbose=DEFAULT_VERBOSE;

	// Allow user to set number of frames
	int num_frames_from_user=get_variable(commandline,"num_frames",-1);
	if(num_frames_from_user>*ptr_num_frames || num_frames_from_user<0){
		if(verbose>4) cout << "Warning: num_frames=" << *ptr_num_frames << " is not set by user." << endl;
	}else{
		*ptr_num_frames=num_frames_from_user;
	}
	if(verbose>4) cout << " num_frames=" << *ptr_num_frames << endl;
}



/**
 *    TODO this function should be tested
 *  Make a random trajectory, e.i. and ideal gas.
 *  Note expects bbox to be set.
 *  Put particles in one of the 27 central images
 *
 *  Return num_frames integer
 */
int make_trajectory_dummy(
		int num_frames,
		int num_atoms,
		double bbox[],
		int type[],
		int img_coords[],
		double coords[]){

	int verbose=DEFAULT_VERBOSE;
	int num_types=3;
	//int rand_lenght=10000;

	srand(2010);	// Set pseudo-random number seed

	//double random_number;
	for(int frame=0;frame<num_frames;frame++){
		for(int atom=0;atom<num_atoms;atom++){
			for(int dim=0;dim<NUM_DIM;dim++){
				int index_coords=frame*num_atoms*NUM_DIM+atom*NUM_DIM+dim;
				//random_number=(double)(rand()%rand_lenght)/(double)rand_lenght;
				//coords[index_coords]=random_number*bbox[frame*NUM_DIM+dim];
				coords[index_coords]=(double)(index_coords/3);
				//random_number=(double)(rand()%rand_lenght)/(double)rand_lenght;
				img_coords[index_coords]=index_coords/3;
				type[num_atoms+atom]=atom%num_types;
				if(verbose>8) cout << "i=" << index_coords << " coords[i]=" << coords[index_coords] << " img_coords[i]=" << img_coords[index_coords] << " (type "<< type[num_atoms+atom] << ")" << endl;
			}
		}
	}

	return num_frames;
}




/**
 *  Return number of particles (ptr_num_frames) and number of frames (ptr_num_atoms)
 */
void get_metadata_lammps_atomfile (string filename,int *ptr_num_frames, int *ptr_num_atoms) {
	int verbose=DEFAULT_VERBOSE;

	ifstream ifile;
	ifile.open(filename.c_str());

	int frame_counter=0;

	if (ifile.is_open()) {
		string line;
		getline(ifile,line);
		while(! ifile.eof() ){
			if( line.compare("ITEM: NUMBER OF ATOMS")==0 ){
				getline (ifile,line);
				//cout << line << endl;
				frame_counter++;
				*ptr_num_atoms=atoi(line.c_str());
			}
			getline(ifile,line);
		}
	}else{
		cout << "Error: Unable to open " << filename << ". Exit." << endl;
		exit(0);
	}
	ifile.close();

	if(verbose>4) cout << "Counted " << frame_counter << " frames in " << filename;
	if(verbose>4) cout << " and " << *ptr_num_atoms << " atoms" << endl;
	*ptr_num_frames=frame_counter;
}







/**
 *
 *  Load a configuration(s) of a LAMMPS "atom" ascii output file with images. In LAMMPS, use something like:
 *
 *      dump dumpATOM all atom 5000 traj.atom
 *      dump_modify dumpATOM image yes
 *
 *  Assume same number of atoms in each frames
 *
 *
 */
int load_lammps_atomfile (  string filename,
							int num_frames,
							int num_atoms,
							double bbox[],
							int type[],
							int img_coords[],
							double coords[]) {

	int verbose=DEFAULT_VERBOSE;
	int num_dim=NUM_DIM;		// number of dimensions

	//double *_bbox = new double[num_frames,num_atoms,num_dim];	// size of boundary box
	//int *_type = (int) type;
	//int *_img_coords = new int[num_frames,num_atoms,num_dim];    // image coordinates
	//double *_coords = new double[num_frames,num_atoms,num_dim];  // coordinates

	int frame=-1;   // number of frames, will be read from file
	ifstream ifile (filename.c_str());

	bool found_all_atomes=true;

	if (ifile.is_open()) {
		string line;
		int atom_counter=0;

		bool stopReadingInFile=false;
		while ( ! (ifile.eof())  && !(stopReadingInFile) ) {
			getline (ifile,line);
			//cout << line << endl;
			if( line.compare("ITEM: TIMESTEP")==0 ){ //This is the first line of a new configuration
				if( (! found_all_atomes) && verbose>5 ) cout<<"Warning: did not find all atoms in frame " << frame <<". atom_counter=" << atom_counter << " num_atoms=" << num_atoms <<  endl;
				found_all_atomes=false;
				atom_counter=0;

				getline (ifile,line);
				frame++;	// count number of frames
				if(verbose>6) cout << "Frame: "<< frame <<". Time Step: " << line << endl;
				if(frame==num_frames){
					// Got all frames, so stop reading in file
					frame--;
					stopReadingInFile=true;
					if(verbose>6) cout << "Note, got all frames and stop reading; stopReadingInFile=true" << endl;
				}
				//if(verbose>4 && ( frame>=num_frames )) cout << "note: Got all frames (frame>=num_frames)." << endl;
			}
			else if( line.compare("ITEM: NUMBER OF ATOMS")==0 ){
				getline (ifile,line);
				int num_atoms_temp=atoi(line.c_str());
				if(verbose>6) cout << "Number of atoms: " << line << endl;
				if( ! num_atoms_temp==num_atoms ){
					if(verbose>0) cout << "Warning: unexpected number of atoms. Read " << num_atoms_temp << ", but expected " << num_atoms << line << endl;
				}
			}
			else if( line.compare("ITEM: BOX BOUNDS")==0 || line.compare("ITEM: BOX BOUNDS pp pp pp")==0  ){
				for(int d=0;d<num_dim;d++){
					getline (ifile,line);
					if(verbose>6) cout << "Box boundaries: " << line << endl;
					char * pEnd;
					double d1,d2;
					d1 = strtod (line.c_str(),&pEnd);
					d2 = strtod (pEnd,NULL);
					bbox[frame*num_dim+d]=d2-d1;
				}
			}
			else if(line.compare("ITEM: ATOMS id type xs ys zs ix iy iz")==0 ){ //Last line before atom coordinates. Make new configuration
				if(verbose>6) cout << "Found expected header: " << line << endl;
				// TODO read in particles here and make an error in the else
				for(int iline=0;iline<num_atoms;iline++){
					getline (ifile,line);
					char * pEnd;
					int    i0, i1;
					double d0, d1, d2;
					int   ii0,ii1,ii2;
					i0 = strtol (line.c_str(),&pEnd,10);
					int atom=i0-1;
					i1 = strtol (pEnd,&pEnd,10);
					type[atom]=i1 - 1; // NB: note that indexes are not LAMMPS, but c-convention (starting with 0)
					d0 = strtod (pEnd,&pEnd);
					coords[frame*num_atoms*num_dim+atom*num_dim+0]=d0*bbox[frame*num_dim+0];
					d1 = strtod (pEnd,&pEnd);
					coords[frame*num_atoms*num_dim+atom*num_dim+1]=d1*bbox[frame*num_dim+1];
					d2 = strtod (pEnd,&pEnd);
					coords[frame*num_atoms*num_dim+atom*num_dim+2]=d2*bbox[frame*num_dim+2];
					ii0 = strtol (pEnd,&pEnd,10);
					img_coords[frame*num_atoms*num_dim+atom*num_dim+0]=ii0;
					ii1 = strtol (pEnd,&pEnd,10);
					img_coords[frame*num_atoms*num_dim+atom*num_dim+1]=ii1;
					ii2 = strtol (pEnd,&pEnd,10);
					img_coords[frame*num_atoms*num_dim+atom*num_dim+2]=ii2;

					// Unwrap coordinates
					// TODO Some old version of this function did not do the out the images, and, thus, some (old) programs may then give an error
					for (int dim = 0 ; dim < num_dim ; dim++){
						coords[frame*num_atoms*num_dim+atom*num_dim+dim] += bbox[frame*num_dim+dim]*(double)img_coords[frame*num_atoms*num_dim+atom*num_dim+dim];
					}
					//coords[frame*num_atoms*num_dim+atom*num_dim+0]+=

					atom_counter++;
					if(verbose>8){
						cout << "Particle: " << line << endl;
						cout << " i0= " << i0 ;
						cout << " i1= " << i1 ;
						cout << " d0= " << d0 ;
						cout << " d1= " << d1 ;
						cout << " d2= " << d2 ;
						cout << " ii0= " << ii0 ;
						cout << " ii1= " << ii1 ;
						cout << " ii2= " << ii2 ;
						cout << endl;
					}
				}
			if(atom_counter>=num_atoms){
				found_all_atomes=true;
				if(verbose>6)  cout<<"Got all atoms of frame " << frame << " ( num_frames=" << num_frames << " )" <<endl;
				if(atom_counter>num_atoms){
					if(verbose>0) cout<<"Warning: atoms="<<atom_counter<<" > num_atoms="<<num_atoms << endl << "... line: " << line << endl;
				}
			}
			//getline (ifile,line);
	        //cout << line << endl;
		}
	  }
	}else{
		  if(verbose>0) cout << "Warning: Unable to open file " << filename << endl;
	 }

	ifile.close();

	frame++;

	//if ( ! ( frame>num_frames || frame<0 ) ) {
	//	if (verbose>0) cout << "Warning: loaded " << frame << " frames, however, got num_frames=" << num_frames << " from input" << endl;
	//}

	return frame;
}






/** TODO
 *  Return number of particles (ptr_num_frames) and number of frames (ptr_num_atoms)
 */
void get_metadata_gromacs_grofile (string filename,int *ptr_num_frames, int *ptr_num_atoms) {
	int verbose=DEFAULT_VERBOSE;
	ifstream ifile;
	ifile.open(filename.c_str());

	int line_counter=0;

	if (ifile.is_open()) {
		string line;
		getline(ifile,line);  // read (first) header
		line_counter++;
		getline(ifile,line);  // read line with number of atoms
		line_counter++;
		*ptr_num_atoms=atoi(line.c_str());

		while(! ifile.eof() ){
			getline(ifile,line);
			line_counter++;
		}
		line_counter--; // The last "line" was end of file
	}else{
		cout << "Error: Unable to open " << filename << ". Exit." << endl;
		exit(0);
	}
	ifile.close();
	*ptr_num_frames=line_counter/(*ptr_num_atoms+3);

	if(verbose>6) cout << "Counted " << line_counter << " lines in " << filename << ", num_frames=lines/(num_atoms+3)=" << *ptr_num_frames << endl;
	int remainder=line_counter%(*ptr_num_atoms+3);
	if(verbose>6 && ! ( remainder==0 )) cout << "warning: unexpected number of lines. " << remainder << " to many lines found." << endl;
}









/**
 *  Load gromacs *.gro file
 */
int load_gromacs_grofile (string filename,int num_frames, int num_atoms, double bbox[], int type[], int img_coords[], double coords[]) {
	int verbose=DEFAULT_VERBOSE;
	ifstream ifile (filename.c_str());

	if (ifile.is_open()) {
		for (int frame=0;frame<num_frames;frame++){
			string line;
			getline (ifile,line); // Header
			getline (ifile,line); // Number of atoms TODO make consistency check
			for (int atom=0;atom<num_atoms;atom++){
				getline (ifile,line);	// Line with atom coordinates  atof(tmp_str.c_str())
				if(verbose>8) cout << "atom line: " << line << endl;  // TODO also load atom types
				double x=atof(line.substr(21,8).c_str());
				double y=atof(line.substr(29,8).c_str());
				double z=atof(line.substr(37,8).c_str());
				coords[frame*num_atoms*NUM_DIM+atom*NUM_DIM+0]=x;
				coords[frame*num_atoms*NUM_DIM+atom*NUM_DIM+1]=y;
				coords[frame*num_atoms*NUM_DIM+atom*NUM_DIM+2]=z;
				string type_str=line.substr(5,1);
				if( type_str.find('B')==0 ){     // TODO Change this, this is not really generic (I assume that there are A, and B particles).
					type[atom]=1;
				}else{
					type[atom]=0;
				}

				if(verbose>8) cout << "atom=" << atom <<  " type=" << type[atom] << " x=" << x << " y=" << y << " z=" << z << endl;
			}

			// Get periodic boundary box
			getline (ifile,line);
			char * pEnd;
			bbox[frame*NUM_DIM+0] = strtod (line.c_str(),&pEnd);
			bbox[frame*NUM_DIM+1] = strtod (pEnd,&pEnd);
			bbox[frame*NUM_DIM+2] = strtod (pEnd,&pEnd);
			if(verbose>7) cout << "bbox: " << bbox[frame*NUM_DIM+0] << " " << bbox[frame*NUM_DIM+1] << " " << bbox[frame*NUM_DIM+0] << endl;
		}
	}else{
		if (verbose>4) cout << "Unable to open " << filename << ". Exit." << endl;
		exit (0);
	}

	// TODO put coords in primary image, and set the img_coords with the image coordinates.

	return num_frames;
}






/** TODO
 *  Return number of particles (ptr_num_frames) and number of frames (ptr_num_atoms)
 */
void get_metadata_xyzfile (string filename,int *ptr_num_frames, int *ptr_num_atoms) {
	int verbose=DEFAULT_VERBOSE;
	ifstream ifile;
	ifile.open(filename.c_str());

	int line_counter=0;

	if (ifile.is_open()) {
		string line;
		getline(ifile,line);  // read (first) header
		line_counter++;
		*ptr_num_atoms=atoi(line.c_str());
		getline(ifile,line);  // read line with number of atoms
		line_counter++;

		while(! ifile.eof() ){
			getline(ifile,line);
			line_counter++;
		}
		line_counter--; // The last "line" was end of file
	}else{
		cout << "Error: Unable to open " << filename << ". Exit." << endl;
		exit(0);
	}
	ifile.close();
	*ptr_num_frames=line_counter/(*ptr_num_atoms+2);

	if(verbose>6) cout << "Counted " << line_counter << " lines in " << filename << ", num_frames=lines/(num_atoms+3)=" << *ptr_num_frames << endl;
	int remainder=line_counter%(*ptr_num_atoms+2);
	if(verbose>6 && ! ( remainder==0 )) cout << "warning: unexpected number of lines. " << remainder << " to many lines found." << endl;
}






/**
 *  Load xyz trajectory file
 */
int load_xyzfile (
		string filename,
		int num_frames,
		int num_atoms,
		double bbox[],
		int type[],
		int img_coords[],
		double coords[])
	{
	int verbose=DEFAULT_VERBOSE;

	ifstream ifile (filename.c_str());
	if (ifile.is_open()) {
		for (int frame=0;frame<num_frames;frame++){
			string line;
			getline (ifile,line); // Number of atoms TODO make consistency check
			getline (ifile,line); // Header

			if(verbose>6) cout << line << endl;
			for (int atom=0;atom<num_atoms;atom++){
				getline (ifile,line);	// Line with atom coordinates  atof(tmp_str.c_str())
				if(verbose>8) cout << "atom line: " << line << endl;  // TODO also load atom types
				char * pEnd;
				int    i0;
				double d0, d1, d2;
				i0 = strtol (line.c_str(),&pEnd,10);
				type[atom]=i0;
				d0 = strtod (pEnd,&pEnd);
				coords[frame*num_atoms*NUM_DIM+atom*NUM_DIM+0]=d0;
				d1 = strtod (pEnd,&pEnd);
				coords[frame*num_atoms*NUM_DIM+atom*NUM_DIM+1]=d1;
				d2 = strtod (pEnd,&pEnd);
				coords[frame*num_atoms*NUM_DIM+atom*NUM_DIM+2]=d2;
			}
			//cout << "test " << frame << endl;
			// Set periodic boundary box (nb XYZ files dont have such information)
//			bbox[frame*NUM_DIM+0] = bboxX;
//			bbox[frame*NUM_DIM+1] = bboxY;
//			bbox[frame*NUM_DIM+2] = bboxZ;
//			if(verbose>7) cout << "bbox: " << bbox[frame*NUM_DIM+0] << " " << bbox[frame*NUM_DIM+1] << " " << bbox[frame*NUM_DIM+0] << endl;
		}
	}else{
		if (verbose>4) cout << "Unable to open " << filename << ". Exit." << endl;
		exit (0);
	}
	cout << "test";
	return num_frames;
}







/**
 * Return minimum image distance in x (y or z) direction
 */
double dx_img(int atom,int atom2,int dim,int frame,int num_frames, int num_atoms, double bbox[], double coords[]) {
	double dx=0.0;
	dx=coords[frame*num_atoms*NUM_DIM+atom2*NUM_DIM+dim]-coords[frame*num_atoms*NUM_DIM+atom*NUM_DIM+dim];
	dx-=bbox[frame*NUM_DIM+dim]*floor(dx/bbox[frame*NUM_DIM+dim]+0.5);
	return dx;
}






/**
 * Return minimum image distance in frame between atom and atom2
 */
double dr_img(int atom,int atom2,int frame,int num_frames, int num_atoms, double bbox[], double coords[]) {
	double dr=0.0;
	for (int dim=0;dim<NUM_DIM;dim++){
		//double dx=0.0;
		//dx=coords[frame*num_atoms*num_dim+atom2*num_dim+dim]-coords[frame*num_atoms*num_dim+atom*num_dim+dim];
		//dx-=bbox[frame*num_dim+dim]*floor(dx/bbox[frame*num_dim+dim]+0.5);
		double dxtmp=dx_img(atom,atom2,dim,frame,num_frames,num_atoms,bbox,coords);
		dr+=dxtmp*dxtmp;
	}
	return pow(dr,0.5);
}






/**
 * Return distance in x (y or z) direction and between frame and frame2
 */
double dx(int atom,int atom2,int dim,int frame,int frame2,int num_frames, int num_atoms, double bbox[], double coords[]) {
	// TODO Move origen

	double dxtmp=0.0;
	dxtmp=coords[frame2*num_atoms*NUM_DIM+atom2*NUM_DIM+dim]-coords[frame*num_atoms*NUM_DIM+atom*NUM_DIM+dim];
	return dxtmp;
}
/**   OLD version, August 9th 2011
double dx(int atom,int atom2,int dim,int frame,int frame2,int num_frames, int num_atoms, double bbox[], double coords[]) {
	double dxtmp=0.0;
	dxtmp=coords[frame2*num_atoms*NUM_DIM+atom2*NUM_DIM+dim]-coords[frame*num_atoms*NUM_DIM+atom*NUM_DIM+dim];
	return dxtmp;
}*/



/**
 * Return (total) displacement between atom and atom2 in frame and frame2
 */
double dr(int atom,int atom2,int frame, int frame2, int num_frames, int num_atoms, double bbox[], double coords[]) {
	double dr=0.0;
	for (int dim=0;dim<NUM_DIM;dim++){
		//double dx=0.0;
		//dx=coords[frame*num_atoms*num_dim+atom2*num_dim+dim]-coords[frame*num_atoms*num_dim+atom*num_dim+dim];
		//dx-=bbox[frame*num_dim+dim]*floor(dx/bbox[frame*num_dim+dim]+0.5);

		double dxtmp=dx(atom,atom2,dim,frame,frame2,num_frames,num_atoms,bbox,coords);
		dr+=dxtmp*dxtmp;
	}
	return pow(dr,0.5);
}



/**
 *   Return the geometric center of molecule (calculated from minimum image distance to first atom in molecule)
 */
double get_geometric_center_img(int molecule,int dim,int frame, int num_frames, int num_molecules, int atoms_per_molecule , double bbox[] , double coords[]){
	double dx=0;
	if(atoms_per_molecule>1){
		for(int atom=1;atom<atoms_per_molecule;atom++){
			  //dx_img(int atom,                     int atom2,                int dim,int frame,int num_frames, int num_atoms, double bbox[], double coords[])
			dx+=dx_img(molecule*atoms_per_molecule+0,molecule*atoms_per_molecule+atom,dim,frame,num_frames, num_molecules*atoms_per_molecule,bbox, coords);
			//dx+=coords[frame*num_molecules*atoms_per_molecule*NUM_DIM+molecule*atoms_per_molecule*NUM_DIM+atom*NUM_DIM+dim];
		}
	}
	dx/=(double)atoms_per_molecule;

	return coords[frame*num_molecules*atoms_per_molecule*NUM_DIM+molecule*atoms_per_molecule*NUM_DIM+0*NUM_DIM+dim]+dx;
}








/**
 *   Construct a neighbor list using cut-offs given in rcutMatrix matrix
 */
void set_neighbour_list(double rcutMatrix[],
						int num_frames,
						int num_atoms,
						int num_types,
						double bbox[],
						int type[],
						double coords[],
						int maxnn,
						int num_neighbours[],
						int neighbour[]){

	int verbose=DEFAULT_VERBOSE;

	if(verbose>4) cout << "Construct neighbor list  :" << endl << "frame";

	for(int frame=0;frame<num_frames;frame++) {
		if(verbose>4) cout << " " <<  frame;
		for(int atom=0;atom<num_atoms;atom++) {
			num_neighbours[frame*num_atoms+atom]=0;
			for(int atom2=0;atom2<num_atoms;atom2++){
				//cout << dr_img(atom,atom2,frame,num_frames,num_atoms,bbox,coords) << endl;
				if(type[atom2]>=num_types || type[atom]>=num_types) cout<<"Warning:  While constructing neighbor list, atom type larger than num_types=" << num_types << endl;
				if( !(atom==atom2) && dr_img(atom,atom2,frame,num_frames,num_atoms,bbox,coords)<rcutMatrix[type[atom]*num_types+type[atom]] ){
					neighbour[frame*num_atoms*maxnn+atom*maxnn+num_neighbours[frame*num_atoms+atom]]=atom2;
					num_neighbours[frame*num_atoms+atom]++;
					if(num_neighbours[frame*num_atoms+atom]>maxnn){cout<< "Warning: While constructing neighbor list, num_neighbours=" << num_neighbours[frame*num_atoms+atom] << " in frame "<<frame<<" of atom "<<atom<<" exceeds maxnn="<< maxnn <<". Please increase maxnn, if possible."<<endl ;}
				}
			}
			if(verbose>8){
					cout << "frame " << frame << " atom " << atom << " have " << num_neighbours[frame*num_atoms+atom] << " neighbors: ";
					for(int n=0;n<num_neighbours[frame*num_atoms+atom];n++) cout << " " << neighbour[frame*num_atoms*maxnn+atom*maxnn+n];
					cout << endl;
			}
		}
	}

	if(verbose>4) cout << endl;

}






/**
 *   Calculate displacement from x0 to x1 and return to vec.
 *   The length of x0, x1 and vec should be num_vectors*NUM_DIM.
 */
void get_vectors(	int num_vectors,
					double x0[],
					double x1[],
					double vec[]){

	for ( int i = 0 ; i<num_vectors ; i++ ){
		for ( int d = 0 ; d<NUM_DIM ; d++ ){
			vec[i*NUM_DIM+d]=x1[i*NUM_DIM+d]-x0[i*NUM_DIM+d];
		}
	}
}



/**
 *   Return the cross products of num_vectors vectors in vecB and vecC to vec.
 *       A = B X C
 *   The length of the vector arrays should be num_vectors*NUM_DIM.
 */
void get_cross_products(	int num_vectors,
					double vecB[],
					double vecC[],
					double vecA[]){

	if(!NUM_DIM==3){
		cout << "WARNING: !NUM_DIM==3 and get_cross_product(int num_vectors,double vec0[],double vec1[],double vec[]) may return what was intended.";
	}

	for ( int i = 0 ; i<num_vectors ; i++ ){ // TODO with the use of some modumuls, this can be calculated by looping dimentions

			// x element of A is a_x=b_y*c_z-b_z*c_y
			// Use x -> y -> z -> x permutation for y- and a-coordinate
			for(int d = 0 ; d < NUM_DIM ; d++){
				vecA[i*NUM_DIM+d]
				     =vecB[i*NUM_DIM+(d+1)%NUM_DIM]*vecC[i*NUM_DIM+(d+2)%NUM_DIM]
				     -vecB[i*NUM_DIM+(d+2)%NUM_DIM]*vecC[i*NUM_DIM+(d+1)%NUM_DIM];
			}
	}

}





/**
 *   Return the dot products of vectors in vecA and vecB to vec.
 *   The length of the vector arrays should be num_vectors*NUM_DIM,
 *   and the length of dot_products should be num_vectors.
 */
void get_dot_products(	int num_vectors,
						double vecA[],
						double vecB[],
						double dot_products[]){

	for ( int i = 0 ; i<num_vectors ; i++ ){
			dot_products[i]=0.0;

			// dot_product = a_x*b_x + a_y*b_y + a_z*b_z + ...
			for(int d = 0 ; d < NUM_DIM ; d++){
				dot_products[i]+=vecA[i*NUM_DIM+d]*vecB[i*NUM_DIM+d];
			}
	}

}

/**
 *  Return dot product of vector A and B
 */
double get_dot_product (double vecA[],
						double vecB[] ) {
	double dot = 0.0;
	for(int d = 0 ; d < NUM_DIM ; d++){
		dot += vecA[d]*vecB[d];
	}
	return dot;
}


/**
 *  Return cos(theta) of vector A and B
 */
double get_cos_theta   (double vecA[],
						double vecB[] ) {

	double AA = 0.0;
	double BB = 0.0;
	double AB = 0.0;

	for(int d = 0 ; d < NUM_DIM ; d++){
		AA += vecA[d]*vecA[d];
		BB += vecB[d]*vecB[d];
		AB += vecA[d]*vecB[d];
	}

	return AB/sqrt(AA*BB);
}


/**
 *   Return the cos(theta) of vectors in vecA and vecB to cos_theta.
 *   The length of the vector arrays should be num_vectors*NUM_DIM,
 *   and the length of cos_theta should be num_vectors.
 */
void get_cos_thetas(	int num_vectors,
						double vecA[],
						double vecB[],
						double cos_theta[]){

	for ( int i = 0 ; i<num_vectors ; i++ ){

			double dot_product = 0.0; // dot_product = a_x*b_x + a_y*b_y + a_z*b_z + ...
			double AA = 0.0;
			double BB = 0.0;

			for(int d = 0 ; d < NUM_DIM ; d++){
				dot_product+=vecA[i*NUM_DIM+d]*vecB[i*NUM_DIM+d];
				AA+=vecA[i*NUM_DIM+d]*vecA[i*NUM_DIM+d];
				BB+=vecB[i*NUM_DIM+d]*vecB[i*NUM_DIM+d];
			}

			//   cos(theta) = A.B/(|A|*|B|)
			cos_theta[i]=dot_product/sqrt(AA*BB);
	}

}


/**
 *   Calculate the mean square displacement
 */
void get_mean_square_displacement(
		int num_frames,
		int num_atoms,
		double bbox[],
		int img_coords[],
		double coords[],
		double msd[],
		int trestart){

	int *counter = new int[num_frames];
	for(int f=0;f<num_frames-1;f++){
		msd[f]=0.0;
		counter[f]=0;
	}

	for(int a=1;a<num_atoms;a++){
		for(int f=0;f<num_frames-1;f+=trestart){
			for(int f2=f;f2<num_frames;f2++){
				int df=f2-f;
				double tmp=dr(a,a,f,f2,num_frames,num_atoms,bbox,coords);
				msd[df]+=tmp*tmp;
				counter[df]++;
			}
		}
	}

	for(int f=0;f<num_frames-1;f++){
		if(!(counter[f]==0)){
			msd[f]/=(double)counter[f];
		}
	}
}







/**
 *   Calculate the mean square displacement
 */
void get_mean_square_displacement_of_selected(
		int num_frames,
		int num_atoms,
		double bbox[],
		int img_coords[],
		double coords[],
		bool selected_atom[],
		double msd[],
		int trestart){

	int *counter = new int[num_frames];
	for(int f=0;f<num_frames-1;f++){
		msd[f]=0.0;
		counter[f]=0;
	}

	for(int a=1;a<num_atoms;a++){
		for(int f=0;f<num_frames-1;f+=trestart){
			if(selected_atom[f*num_atoms+a]){
				for(int f2=f;f2<num_frames;f2++){
					int df=f2-f;
					double tmp=dr(a,a,f,f2,num_frames,num_atoms,bbox,coords);
					msd[df]+=tmp*tmp;
					counter[df]++;
				}
			}
		}
	}

	for(int f=0;f<num_frames-1;f++){
		if(!(counter[f]==0)){
			msd[f]/=(double)counter[f];
		}
	}
}








/**
 *   Calculate self intermediate scattering function, F_s(q,t).
 */
void get_self_intermediate_scattering_function(
		int num_frames,
		int num_atoms,
		double bbox[],
		int img_coords[],
		double coords[],
		double Fs[],
		double q_vector,
		int trestart){

	int *counter = new int[num_frames];
	for(int f=0;f<num_frames;f++){
		Fs[f]=0.0;
		counter[f]=0;
	}

	for(int a=1;a<num_atoms;a++){
		for(int f=0;f<num_frames-1;f+=trestart){
			for(int f2=f;f2<num_frames;f2++){
				int df=f2-f;
				//double tmp=0.0;
				for(int d=0;d<NUM_DIM;d++){
					Fs[df]+=cos(q_vector*dx(a,a,d,f,f2,num_frames,num_atoms,bbox,coords));
					counter[df]++;
				}
				//double tmp=dr(a,a,f,f2,num_frames,num_atoms,bbox,coords);
				//Fs[df]+=tmp*tmp;
			}
		}
	}

	for(int f=0;f<num_frames-1;f++){
		if(!(counter[f]==0)){
			Fs[f]/=(double)counter[f];
		}
	}

}







/**
 *    Calculate self intermediate scattering function of selected atoms (at t=0), F_s(q,t).
 */
void get_self_intermediate_scattering_function_of_selected (
		int num_frames,
		int num_atoms,
		double bbox[],
		int img_coords[],
		double coords[],
		bool selected_atom[],
		double Fs[],
		double q_vector,
		int trestart ) {

	int *counter = new int[num_frames];
	for(int f=0;f<num_frames;f++){
		Fs[f]=0.0;
		counter[f]=0;
	}

	for(int a=1;a<num_atoms;a++){
		for(int f=0;f<num_frames-1;f+=trestart){
			if(selected_atom[f*num_atoms+a]){
				for(int f2=f;f2<num_frames;f2++){
					int df=f2-f;
					//double tmp=0.0;
					for(int d=0;d<NUM_DIM;d++){
						Fs[df]+=cos(q_vector*dx(a,a,d,f,f2,num_frames,num_atoms,bbox,coords));
						counter[df]++;
					}
					//double tmp=dr(a,a,f,f2,num_frames,num_atoms,bbox,coords);
					//Fs[df]+=tmp*tmp;
				}
			}
		}
	}

	for(int f=0;f<num_frames-1;f++){
		if(!(counter[f]==0)){
			Fs[f]/=(double)counter[f];
		}
	}

}





/**
 *   Calculate the self radial distribution function rdf
 */
void get_radial_distribution_function(
		int num_frames,
		int num_atoms,
		double bbox[],
		double coords[],
		double rdf[],
		double r_max,
		int num_bins,
		int frame_stride){

	int *hr = new int[num_bins];			// histogram of distances
	//double *gr = new double[num_bins];	// radial distribution function
	//double L=r_max; 						// Largest length considered
	double V=0.0;  		  					// (average) Box volume
	int frames_used=0;
	for (int b=0;b<num_bins;b++) hr[b]=0;
	for (int frame=0;frame<num_frames;frame+=frame_stride){
		frames_used++;
		V+=bbox[frame*NUM_DIM+0]*bbox[frame*NUM_DIM+1]*bbox[frame*NUM_DIM+2];
		for (int atom=0;atom<num_atoms-1;atom++){	// j>i doublet sum over atoms
			for(int atom2=atom+1 ; atom2<num_atoms ; atom2++){
				double dr=dr_img(atom,atom2,frame,num_frames,num_atoms,bbox,coords)  ;
				if(dr<r_max){
					hr[(int) floor( (double)num_bins * dr/r_max )]+=2;
				}
			}
		}
	}
	V/=(double)frames_used;

	// Make normalized radial distribution
	double dr=r_max/num_bins;	// bin size
	for (int b=0;b<num_bins;b++){
		double r=((double)b*r_max)/(double)num_bins;
		rdf[b] = (double)hr[b] / (double)num_atoms / (double)frames_used;
		rdf[b]/=4.0/3.0*PI*(pow(r+dr,3)-pow(r,3));
		rdf[b]/= (double)num_atoms/V;
	}
}







/**
 *   Calculate the self radial distribution function (rdf) of selected atoms
 */
void get_radial_distribution_function_of_selected(
		int num_frames,
		int num_atoms,
		double bbox[],
		double coords[],
		bool selected_atoms[],
		double rdf[],
		double r_max,
		int num_bins,
		int frame_stride){

	int *hr = new int[num_bins];			// histogram of distances
	//double *gr = new double[num_bins];	// radial distribution function
	//double L=r_max; 						// Largest length considered
	double V=0.0;  		  					// (average) Box volume
	int frames_used=0;
	int num_selected_atoms=0;
	for (int b=0;b<num_bins;b++) hr[b]=0;
	for (int frame=0;frame<num_frames;frame+=frame_stride){
		frames_used++;
		V+=bbox[frame*NUM_DIM+0]*bbox[frame*NUM_DIM+1]*bbox[frame*NUM_DIM+2];
		for (int atom=0;atom<num_atoms;atom++){	// j,i doublet sum over atoms
			if(selected_atoms[frame*num_atoms+atom]){
				num_selected_atoms++;
				for(int atom2=0 ; atom2<num_atoms ; atom2++){
					double dr=dr_img(atom,atom2,frame,num_frames,num_atoms,bbox,coords)  ;
					if(dr<r_max and dr>0.0){
						hr[(int) floor( (double)num_bins * dr/r_max )]+=1;
					}
				}
			}
		}
	}
	V/=(double)frames_used;

	// Make normalized radial distribution.
	double dr=r_max/num_bins;	// bin size
	for (int b=0;b<num_bins;b++){
		double r=((double)b*r_max)/(double)num_bins;
		rdf[b] = (double)hr[b] / (double)num_selected_atoms / (double)frames_used;
		rdf[b]/=4.0/3.0*PI*(pow(r+dr,3)-pow(r,3));
		rdf[b]/= (double)num_atoms/V;
	}
}







/**
 *   Self part of van Hove correlation function.
 *
 *   Result is returned to Gs.
 *
 *   Only frame differences on a log2 scale is returned (1,2,4,8,...), thus
 *   the length of Gs should be num_bins*num_log2_df where num_log2_df may be compute from num_frames as
 *
 *      for (int df = 1  ; df < num_frames ; df*=2) num_log2_df++;
 *
 *  and result is prited using something like:
 *
 *  for(int bin=0;bin<num_bins-1;bin++){
 *		fprintf (ofile_rdf, "%f" , (double)bin*r_max/(double)num_bins );
 *		for (int df_index = 0  ; df < num_log2_df ; df_index++){
 *			fprintf ( ofile_Gs, " %f" , Gs [df_index*num_bins+bin] ) ;
 *		}
 *	    fprintf (ofile_rdf, "\n");
 *	}
 *
 */
void get_van_Hove_self_log2(
		int num_frames,
		int num_atoms,
		double bbox[],
		double coords[],
		double Gs[],
		double r_max,
		int num_bins,
		int frame_stride){

	int num_log2_df = 0;
	for (int df = 1  ; df < num_frames ; df*=2) num_log2_df++;

	double V=0.0;  		  					// (Average) Box volume.
	int frames_used=0;

	int *hr = new int[num_bins];// histogram of distances (temporary storage)


	int df_index = 0;
	for (int df = 1  ; df < num_frames ; df*=2){

		// Reset arrays
		for (int bin=0;bin<num_bins;bin++) hr[bin]=0;
		V=0.0;
		frames_used=0;

		// Make histogram for this df frame distance
		for (int frame=0;frame<num_frames-df;frame+=frame_stride){

			frames_used++;
			V+=bbox[frame*NUM_DIM+0]*bbox[frame*NUM_DIM+1]*bbox[frame*NUM_DIM+2];
			for (int atom=0;atom<num_atoms;atom++){

				// TODO is this the jumpdistance, e.i is images done right.
				double jump_distance=dr(atom,atom,frame,frame+df,num_frames,num_atoms,bbox,coords);

				// Distance particle have travel between df frames.
				/*double jump_distance=0.0;
				for (int dim=0;dim<NUM_DIM;dim++){
					double dx=coords[frame*num_atoms*NUM_DIM+atom*NUM_DIM+dim];
					dx-=coords[(frame+df)*num_atoms*NUM_DIM+atom*NUM_DIM+dim];
					jump_distance+=dx*dx;
				}
				jump_distance=sqrt(jump_distance);
				cout << "f=" << frame << " df=" << df << " a=" << atom << " jump_distance=" << jump_distance << endl;
				*/

				if(jump_distance<r_max){
					hr[ (int) floor( (double)num_bins * jump_distance/r_max )]++;
				}
			}
		}

		// Make normalization and move result to Gs
		V/=(double)frames_used;

		double dr=r_max/num_bins;	// bin size
		for (int bin=0;bin<num_bins;bin++){
			//double r=((double)bin*r_max)/(double)num_bins;
			Gs[df_index*num_bins+bin] = (double)hr[bin] / (double)num_atoms / (double)frames_used / dr;
			//Gs[df_index*num_bins+bin]/=4.0/3.0*PI*(pow(r+dr,3)-pow(r,3));
			//Gs[df_index*num_bins+bin]/= (double)num_atoms/V;
		}

		df_index++;
	} // end of df loop

}



/**
 *   Return collective rotational autocorrelation using the first,
 *
 *     C1(t) = < cos[theta(0)] * cos[theta(t)] >,
 *
 *   and second Legendre polynomial,
 *
 *	   C2(t) = < (3/2*cos^2[theta(0)]-1/2) * (3/2*cos^2[theta(t)]-1/2) >,
 *
 * 	 where theta is the angle between summed vector and either x, y or z vectors.
 *   C1(t) returned to rotacf, and
 *   C2(t) to rotacf2.
 *
 */
void get_collective_rotational_auto_correlation(
		int num_frames,
		int num_vectors,
		double vector[],
		double rotacf[],
		double rotacf2[],
		int trestart){

}



/**
 *   Return rotational autocorrelation using the first,
 *
 *     C1(t) = < cos[theta(0)] * cos[theta(t)] >,
 *
 *   and second Legendre polynomial,
 *
 *	   C2(t) = < (3/2*cos^2[theta(0)]-1/2) * (3/2*cos^2[theta(t)]-1/2) >,
 *
 * 	 where theta is the angle between vector and either x, y or z vectors.
 *   C1(t) is returned to rotacf, and
 *   C2(t) to rotacf2.
 *
 */
void get_rotational_auto_correlation(
		int num_frames,
		int num_vectors,
		double vector[],
		double rotacf[],
		double rotacf2[],
		int trestart){

	int *frame_counter = new int[num_frames];
	for(int i=0;i<num_frames;i++) frame_counter[i]=0;

	for(int i=0;i<num_frames;i++) rotacf[i] = 0.0;
	for(int i=0;i<num_frames;i++) rotacf2[i] = 0.0;

	for (int f0 = 0; f0<num_frames ; f0+=trestart){
		for (int f1 = f0; f1<num_frames ; f1++){
			for (int v = 0;v<num_vectors;v++){

				// TODO use get_costheta_of_between_vectors function

				// Length of vector 0
				double length0=0;
				for (int d = 0; d < NUM_DIM;d++){
					length0 += vector[f0*num_vectors*NUM_DIM+v*NUM_DIM+d]
					          *vector[f0*num_vectors*NUM_DIM+v*NUM_DIM+d];
				}
				length0=sqrt(length0);

				// Length of vector 1
				double length1=0;
				for (int d = 0; d < NUM_DIM;d++){
					length1 += vector[f1*num_vectors*NUM_DIM+v*NUM_DIM+d]
					          *vector[f1*num_vectors*NUM_DIM+v*NUM_DIM+d];
				}
				length1=sqrt(length1);

				// Loop over x,y and z
				frame_counter[f1-f0]++;
				for (int d = 0 ; d < NUM_DIM ; d++ ){

					// cos(theta) = z/sqrt(x^2+y^2+z^2)
					double cos_theta0 = vector[f0*num_vectors*NUM_DIM+v*NUM_DIM+d]/length0;
					double cos_theta1 = vector[f1*num_vectors*NUM_DIM+v*NUM_DIM+d]/length1;

					rotacf[f1-f0]+=cos_theta0*cos_theta1;
					rotacf2[f1-f0]+=0.25*(3.0*cos_theta0*cos_theta0-1.0)*(3.0*cos_theta1*cos_theta1-1.0);
				}
			}
		}
	}

	// Divide to get average
	for ( int i = 0 ; i < num_frames ; i++ ){
		rotacf[i]/=(double)frame_counter[i];
		rotacf2[i]/=(double)frame_counter[i];
	}
}

/**
 *   Return the qubatic orderparameter
 *      Q_i = \frac{1}{N-1} \sum_{j\neq i}^N \frac{\cos^2(2\theta_{ij})-\frac{7}{15}}{1-\frac{7}{15}}
 */
void get_qubatic_orderparameter(int num_frames,int num_vectors,double vectors[],double Qi[]){

	// Reset Q_i's
	for (int i = 0 ; i<num_frames*num_vectors ; i++ ) Qi[i]=0.0;

	// Sum of cos(2*theta)^2;
	for (int f = 0; f<num_frames ; f++){
		for (int i = 0;i<num_vectors;i++){
			for (int j = 0;j<num_vectors;j++){

				if(!(i==j)){

					// TODO put calculation of vector angle to function, since this may be used elsewhere
					// Calculate angle = acos ( A.B/sqrt(A.A * B.B) ) where . is the dot product.
					double AA = get_dot_product (&vectors[f*num_vectors*NUM_DIM+i*NUM_DIM],&vectors[f*num_vectors*NUM_DIM+i*NUM_DIM]);
					double BB = get_dot_product (&vectors[f*num_vectors*NUM_DIM+j*NUM_DIM],&vectors[f*num_vectors*NUM_DIM+j*NUM_DIM]);
					double AB = get_dot_product (&vectors[f*num_vectors*NUM_DIM+i*NUM_DIM],&vectors[f*num_vectors*NUM_DIM+j*NUM_DIM]);
					double theta = acos( AB/sqrt(AA*BB) );

					//  acos( AB/sqrt(AA*BB) ) can return nan when theta is 0 or pi (due to round off error). Fixed below;
					if( AB/sqrt(AA*BB)>1.0  ) theta = 0;
					if( AB/sqrt(AA*BB)<-1.0  ) theta = -PI;
					// This should never be written:
					if(theta!=theta) cout << "Warning, 'not a number' angle produced: f= "<< f <<" j= " << j << " i= " << i << " AB= " << AB << " AA= " << AA << " BB= " << BB  << " AB/sqrt(AA*BB)= " << AB/sqrt(AA*BB) <<  endl;

					// Add to sum
					Qi[f*num_vectors+i] += (cos(2.0*theta)*cos(2.0*theta)-7.0/15.0)/(1.0-7.0/15.0)/((double)num_vectors-1.0);


					//if(i == 221 && f == 2) cout << "f=  " << f << " j= " <<  j  << " theta = " << theta << " sum Qi= " << Qi[f*num_vectors+i] << " AB = " << AB << endl;
				}
			}
		}
	}


	// Normalize Qi's; [15/8 cos^2(2 theta)-7/8]/(N-1)
	//for (int i=0;i<num_frames*num_vectors;i++){
		// Qi[i] = (Qi[i]-7.0/15.0)/(1.0-7.0/15.0)/((double)num_vectors-1.0);
		//Qi[i] /= ;
	//}

}



/**
 * Return
 *   cos(theta) of the
 */
/*double get_cos_theta_between_vectors(int num_frames,int num_vectors,double base_vectors[],int vector0,frame0,vector1,frame1){

	// Length of vector 0
	double length0=0.0;
	for (int d = 0; d < NUM_DIM;d++){
		length0 += vector[frame0*num_vectors*NUM_DIM+vector0*NUM_DIM+d]*vector[frame0*num_vectors*NUM_DIM+vector0*NUM_DIM+d];
	}
	length0=sqrt(length0);

	// Length of vector 1
	double length1=0.0;
	for (int d = 0; d < NUM_DIM;d++){
		length0 += vector[frame1*num_vectors*NUM_DIM+vector1*NUM_DIM+d]*vector[frame1*num_vectors*NUM_DIM+vector1*NUM_DIM+d];
	}
	length1=sqrt(length1);

	// cos(theta) = z/sqrt(x^2+y^2+z^2)
	return vector[f0*num_vectors*NUM_DIM+v*NUM_DIM+d]/length0;
}*/


/**
 *   TODO Calculate structure factor by integrating the radial distribution function
 */
/*void get_structure_factor_from_rdf(
		int num_frames,
		int num_atoms,
		double bbox[],
		double coords[],
		double rdf[],
		double Sq[],
		double k_max,
		int num_bins){

	//FILE * ofile_Sk = fopen("Sk.dat","w");
	//fprintf (ofile_Sk, "# k S(k) \n");
	//double kmin=get_variable(commandline,"kmin",0.01);
	//double kmax=get_variable(commandline,"kmax",50.01);

	double V=0.0;
	for(int frame=0;frame<num_frames;frame++){
		V+=bbox[frame*NUM_DIM+0]*bbox[frame*NUM_DIM+1]*bbox[frame*NUM_DIM+2]/(double)num_frames;
	}
	double L=0.5*pow(V,1/3);

	//for (double k=kmin;k<kmax;k+=(kmax-kmin)/(double)num_bins){
	double dk=k_max/double(num_bins);
	for (double i=0;i<num_bins;i++){
		k=dk*(double)i;
		double Sq[i]=0.0;
		for (int b=0;b<num_bins;b++){
			double r=(b*L)/num_bins+0.5*dr;
			Sq+=4*PI*num_atoms/V/k*(gr[b]-1.0)*r*sin(k*r)*dr;
		}
		Sq+=1;
		fprintf (ofile_Sk, "%f %f\n", k , Sq);
	}
	//fclose(ofile_Sk);
}
*/






#endif /* TRJ_H_ */
