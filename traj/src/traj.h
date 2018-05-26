/*
 * traj.h
 *  by Ulf R. Pedersen
 *     http://www.urp.dk
 */

#ifndef TRAJ_H_
#define TRAJ_H_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include "trj.h"
using namespace std;

#ifndef PI
#define PI 3.14159265358979323846264338328 // 3.14159265358979323846 ratio of a circles circumference to its diameter
#endif

#ifndef NUM_DIM
#define NUM_DIM 3					// number of spatial dimensions
#endif

#ifndef DEFAULT_VERBOSE
#define DEFAULT_VERBOSE 5					// default_value for how verbose functions are (9 for loud, 5 for comfortable and 0 for silent)
#endif



class Traj {
public:
		// Simple data types
		string variables;
		int verbose;	// Between 0-9, 0 makes the program silent, while 9 prints A LOT of stuff
		int num_atoms;	// Number of atoms
		int num_frames; // Number of frames
		int num_types;  // Number of types

		double *bbox;	// size of boundary box		// TODO use STL vector instead of these primitive types
		int *type;		// Atom type
		int *img_coords;// image coordinates
		double *coords;	// Coordinates
		bool *selected; // Selected atom

		double time_of_first_frame;// Time of first frame
		double frame_dt;// Time interval between frames

		double kT;		// Thermal energy

		// Neighbor list information
		vector<int> nl;			// Neighbor list
		vector<int> nl_first;	// i'th index is first element in nl of i'th particles (-1 if no neighbors)
		vector<int> nl_last;	// i'th index is last element in nl of i'th particles (-1 if no neighbors)
		vector<int> nl_num;    	// ... and, likewise, nl_num in number of neighbors .
		vector<int> n_type;		// What kind of neighbor is this (e.g. 0 for h-bond donor and 1 for accepter)
		vector<bool> nl_selected;

		// Constructors
		Traj(string&);

		// Loading data and allocate arrays
		void load_data_from_disk();

		// Return stuff about the simulation
		double get_average_box_volume();

		// Trajectory manipulations
		void smooth_by_time_integral();

		// Calculate/write stuff
		void van_Hove_self ( string ) ;
		void mean_squared_displacement ( string ) ;
		void radial_distribution ( string ) ;
		void rhok_statistics ( string , int ) ;
		void incoherent_self_scattering ( string , double ) ;
		void incoherent_self_scattering_of_selected ( string , double ) ;

		// Atoms displacement
		double get_dr(int,int,int,int);
		double get_dr_img(int,int,int);
		double get_dr_img(int,int,int,int);
		double get_x_img(int,int,int);
		double get_dx_img(int,int,int,int);
		double get_dx_img(int,int,int,int,int);

		// More stuff
		double get_time(int);
		void deselect_all();

		// Analysis that have to do with WATER
		void construct_hbond_network (  ) ;
		void enumerate_nl_hbonds ( ) ;
		void enumerate_nl_hbonds_by_occurrence ( ) ;
		void nl_bond_time_autocorrelation ( string ) ;
		void hbond_dynamic_exitations_ts (  ) ;
		void hbond_dynamic_exitations_sojourn ( ) ;
		bool do_n_type_exist_in_frame ( int  , int , int ) ;
		int n_type_in_frame ( int , int , int ) ;
		void water_dipole_rotational_correlation ( ) ;

		// Kink analysis
		void select_displacement_kinks ( double , int , int , int ) ;

};


/**
* Constructor that set
*    string variables
* and
*    int verbose
*/
Traj::Traj(string& variables_in){
	variables = variables_in;
	verbose=get_variable(variables,"verbose",5);
}


/**
*      Return the average volume of the the average volume of the boundary box
*/
double Traj::get_average_box_volume () {

	double volume = 0.0;
	for ( int f = 0 ; f < num_frames ; f++ ) {

		double volume_frame = 1.0;
		for ( int d = 0 ; d < NUM_DIM ; d++ ){
			volume_frame *= bbox[f*NUM_DIM+d];
		}
		volume += volume_frame;
	}

	return volume/(double)num_frames;
}


/**
 * Load data from files and allocate (all) variables
 */
void Traj::load_data_from_disk(){
	if(verbose>4) cout << endl << " ..:: Load trajectory ::.." << endl;

	// Load meta data
	get_metadata_traj(variables,&num_frames,&num_atoms);

	// Allocate arrays
	bbox = new double[num_frames*NUM_DIM];	// size of boundary box
	type = new int[num_atoms];					// Atom type
	img_coords = new int[num_frames*num_atoms*NUM_DIM];    // image coordinates
	coords = new double[num_frames*num_atoms*NUM_DIM];  // coordinates
	selected = new bool[num_frames*num_atoms];

	// Load data from disk
	num_frames=load_traj(variables,num_frames,num_atoms,bbox,type,img_coords,coords);

	// TODO Read frame_dt from input file (if possible)
	time_of_first_frame =  get_variable(variables,"time_of_first_frame",0.0);
	frame_dt = get_variable(variables,"frame_dt",1.0) ;

	kT = get_variable(variables,"kT",1.0);
}


/**
 * TODO  Smooth trajectory by performing a time integral of previous times.
 */
void Traj::smooth_by_time_integral(){
	// TODO
}


/**
 * Write van Hoove correlation function
 */
void Traj::van_Hove_self(string ofilename){

	if(verbose>4) cout << endl << " ..:: Self part of van Hove correlation function (" << ofilename << ") ::.." << endl;


	// Initialize variables
	int num_bins=get_variable(variables,"num_bins",5000);
	double r_max=get_variable(variables,"r_max",5);

	int frame_stride_Gs=get_variable(variables,"frame_stride_Gs",1);

	int num_log2_df = 0; // Number of possible frame differences on log2 scale (1,2,4,8 ...)
	for (int df = 1  ; df < num_frames ; df*=2) num_log2_df++;

	double *Gs = new double[num_bins*num_log2_df];

	if(verbose>4) cout << "num_log2_df=" << num_log2_df << endl;

	// Calculate self-part of van Hooves correlation function
	get_van_Hove_self_log2(num_frames,num_atoms,bbox,coords,Gs,r_max,num_bins,frame_stride_Gs);

	// Write to file
	FILE * ofile_Gs = fopen(ofilename.c_str(),"w");

	fprintf (ofile_Gs, "# Self part of van Hove correlation function, Gs, of particles centers (4*pi*r^2*Gs): [Jump distance]");
	for ( int df = 1 ; df < num_frames ; df*=2 ) fprintf(ofile_Gs," [t=%f]",(double)df*frame_dt);
	fprintf (ofile_Gs, "\n");

	for(int bin=0;bin<num_bins-1;bin++){
		fprintf (ofile_Gs, "%f" , (double)bin*r_max/(double)num_bins );
		for (int df_index = 0  ; df_index < num_log2_df ; df_index++){
			fprintf ( ofile_Gs, " %f" , Gs [df_index*num_bins+bin] );
		}
		fprintf (ofile_Gs, "\n");
	}

	fclose(ofile_Gs);

}


/**
 * Calculate and write the mean squared displacement
 */
void Traj::mean_squared_displacement(string ofilename){

	if(verbose>4) cout << endl << " ..:: Mean squared displacement (" << ofilename << ") ::.. " << endl;

	double * msd = new double[num_frames];
	int msd_trestart=get_variable(variables,"msd_trestart",10);
	get_mean_square_displacement(
			num_frames,
			num_atoms,
			bbox,
			img_coords,
			coords,
			msd,
			msd_trestart);

	// Calculate diffusion constant from as
	double diffusion_constant = 0.0;
	int diff_estimate_first=get_variable(variables,"diff_estimate_first",num_frames/20);
	int diff_estimate_last=get_variable(variables,"diff_estimate_last",num_frames/4);
	if(verbose>4) cout << "Use slope of mean squared displacement from t_first = " << frame_dt*(double)diff_estimate_first << " to t_last = " << frame_dt*(double)diff_estimate_last << " to estimate a diffusion constant of " << endl ;
	if(diff_estimate_first>=0 && diff_estimate_last<num_frames){
		diffusion_constant  = (msd[diff_estimate_last] - msd[diff_estimate_first]);
		diffusion_constant /= frame_dt*(double)(diff_estimate_last-diff_estimate_first);
		diffusion_constant /= 2.0 * (double)NUM_DIM;
		if(verbose>4) cout << "D = " << diffusion_constant << endl;
	}else{
		if(verbose>4) cout << "Warning: diff_estimate_first and/or diff_estimate_last is out of range. Diffusion constant not calculated." << endl;
	}

	FILE * ofile = fopen(ofilename.c_str(),"w");
	fprintf (ofile, "# Mean squared displacement [ time ; msd ; 6Dt ] ( D = %f )\n",diffusion_constant);
	for ( int frame = 0 ; frame < num_frames - 1  ; frame++ ) fprintf (ofile, "%f %f %f\n", frame_dt*(double)frame,msd[frame],6.0*diffusion_constant*frame_dt*(double)frame);
	fclose(ofile);



}


/**
 * Write radial distribution function
 */
void Traj::radial_distribution(string ofilename){

	if(verbose>4) cout << endl << " ..:: Radial distribution function (" << ofilename << ") ::.." << endl;

	int num_bins=get_variable(variables,"num_bins",5000);
	double r_max=get_variable(variables,"r_max",5);
	double bin_width = r_max/(double)num_bins;

	int *bins = new int[num_bins];
	for (int i = 0 ; i < num_bins ; i++) bins[i]=0;

	int frame_stride=get_variable(variables,"frame_stride_gr",1);

	double dr;
	double sum_volume = 0.0;
	int sum_atoms = 0;
	for ( int frame = 0 ; frame < num_frames ; frame+=frame_stride ) {
		sum_volume += bbox[frame*3+0]*bbox[frame*3+1]*bbox[frame*3+1];
		sum_atoms +=  num_atoms;
		for ( int atom0 = 0 ; atom0 < num_atoms-1 ; atom0++ ) {
			for ( int atom1 = atom0+1 ; atom1 < num_atoms ; atom1++ ) {
				dr = get_dr_img(atom0,atom1,frame);
				if( dr < r_max ) {
					bins[ (int) floor ( dr / r_max * (double)num_bins ) ]+=2;
				}
			}
		}
	}


	FILE * ofile = fopen(ofilename.c_str(),"w");
	fprintf (ofile, "# Radial distribution function [center of bin; radial distribution in units of number density]\n");

	double V_shell;
	double avg_density = (double)sum_atoms/sum_volume;
	for(int bin=0;bin<num_bins;bin++){
		dr = bin_width*(double)bin;
		V_shell = 4.0/3.0*PI*( (dr+bin_width)*(dr+bin_width)*(dr+bin_width)-dr*dr*dr );
		fprintf (ofile, "%f %f\n",dr+0.5*bin_width,(double)bins[bin]/V_shell/(double)num_frames/(double)num_atoms/avg_density);
	}

	fclose(ofile);
}



/**
 * Write radial distribution function function
 */
void Traj::rhok_statistics(string ofilename,int num_subx){

	int num_nx = get_variable(variables,"num_nx",64);
	if(verbose>4) cout << "Number of k-vectors in x (y and z) direction: num_nx=" << num_nx << endl;

	//int num_subx = get_variable(variables,"num_subx",1);
	if(verbose>4) cout << "Number of subsystems in x (y and z) direction: num_subx=" << num_subx << endl;

	int num_sub=num_subx*num_subx*num_subx;
	if(verbose>4) cout << "Number of subsystems (num_subx^3): num_sub=" << num_sub << endl;

	double k1=2.0*PI*(double)num_subx/bbox[0]; // Smallest k vector
	if(verbose>4) cout << "Length of shortest k vector (2*PI*num_subx/L): k1=" << k1 << endl;



	double *Sk = new double[num_nx];
	for(int i=0;i<num_nx;i++) Sk[i]=0.0;

	double *Fourth = new double[num_nx];		// Fourth moment of distribution
	for(int i=0;i<num_nx;i++) Fourth[i]=0.0;

	double *RE = new double[num_sub*num_nx*NUM_DIM];
	double *IM = new double[num_sub*num_nx*NUM_DIM];



	for (int frame=0;frame<num_frames;frame++){
		for(int i=0;i<num_sub*num_nx*NUM_DIM;i++) RE[i]=0.0;
		for(int i=0;i<num_sub*num_nx*NUM_DIM;i++) IM[i]=0.0;


		for(int atom=0;atom<num_atoms;atom++){

			int sub=0; // subsystem index (assume that all particles are in primary image)
			sub+=                  (int)floor(get_x_img(0,atom,frame)/bbox[0]*(double)num_subx);
			sub+=num_subx*         (int)floor(get_x_img(1,atom,frame)/bbox[1]*(double)num_subx);
			sub+=num_subx*num_subx*(int)floor(get_x_img(2,atom,frame)/bbox[2]*(double)num_subx);
			if(verbose>8) cout << "atom=" << atom << " x=" << coords[atom*NUM_DIM+0] << " y=" << coords[atom*NUM_DIM+1] << "z=" << coords[atom*NUM_DIM+2]  << " sub=" << sub << endl;
			if(verbose>0 && (sub<0 || sub>num_sub-1)) cout << "Error: for atom " << atom << " subvolume index " << sub << "is out of range" << endl ;

			for(int nx=0;nx<num_nx;nx++){
				for(int d=0;d<NUM_DIM;d++){
					RE[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]+=cos((double)nx*k1*coords[atom*NUM_DIM+d]);
					IM[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]-=sin((double)nx*k1*coords[atom*NUM_DIM+d]);
				}
			}
		}

		// Print real and imaginary part of rho_k in subvolumes
		if(verbose>7){
			for(int sub=0;sub<num_sub;sub++){
				cout << endl << "rho_k,Real(subvol=" << sub << ",|k|,k direction): ";
				for(int i=sub*num_nx*NUM_DIM;i<(sub+1)*num_nx*NUM_DIM;i++) cout << " " << RE[i];
			}
			for(int sub=0;sub<num_sub;sub++){
				cout << endl << "rho_k,Imag(subvol=" << sub << ",|k|,k direction): ";
				for(int i=sub*num_nx*NUM_DIM;i<(sub+1)*num_nx*NUM_DIM;i++) cout << " " << IM[i];
			}
		}

		// print S(k)=<rhok^2>/N of this frame
		if(verbose>6) cout << endl << "S(k) : ";

		// Do special calculation for k=0
		double Sk_tmp=0.0;
		double Fourth_tmp = 0.0;
		for(int sub=0;sub<num_sub;sub++){
			int d=0;
			//for(int d=0;d<NUM_DIM;d++){
				int nx=0;
				Sk_tmp+= (RE[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]-(double)num_atoms/(double)num_sub)
						*(RE[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]-(double)num_atoms/(double)num_sub);
				Fourth_tmp+= (RE[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]-(double)num_atoms/(double)num_sub)
								*(RE[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]-(double)num_atoms/(double)num_sub)
								*(RE[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]-(double)num_atoms/(double)num_sub)
								*(RE[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]-(double)num_atoms/(double)num_sub);
				//Sk_tmp/=RE[sub*num_nx*NUM_DIM+0*NUM_DIM+d];
				//Sk_tmp+=IM[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]*IM[sub*num_nx*NUM_DIM+nx*NUM_DIM+d];
			//}
		}
		//Sk_tmp/=(double)num_sub;			// TODO Bessel correction ?
		//Sk_tmp/=(double)num_atoms/(double)num_sub;
		Sk_tmp/=(double)num_atoms;
		if(verbose>6) cout << " " << Sk_tmp << " " << Fourth_tmp;
		Sk[0]+=Sk_tmp;
		Fourth[0]+=Fourth_tmp;

		// do k>0
		for(int nx=1;nx<num_nx;nx++){
			Sk_tmp=0.0;
			Fourth_tmp=0.0;
			for(int sub=0;sub<num_sub;sub++){
				for(int d=0;d<NUM_DIM;d++){
					Sk_tmp+=RE[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]*RE[sub*num_nx*NUM_DIM+nx*NUM_DIM+d];
					Sk_tmp+=IM[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]*IM[sub*num_nx*NUM_DIM+nx*NUM_DIM+d];
					Fourth_tmp+=RE[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]*RE[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]*RE[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]*RE[sub*num_nx*NUM_DIM+nx*NUM_DIM+d];
					Fourth_tmp+=IM[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]*IM[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]*IM[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]*IM[sub*num_nx*NUM_DIM+nx*NUM_DIM+d];
					Fourth_tmp+=2.0*RE[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]*RE[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]*IM[sub*num_nx*NUM_DIM+nx*NUM_DIM+d]*IM[sub*num_nx*NUM_DIM+nx*NUM_DIM+d];
				}
			}
			Sk_tmp=Sk_tmp/(3.0*(double)num_atoms);  // TODO Bessel correction ?
			Fourth_tmp=Fourth_tmp/(3.0*(double)(num_atoms*num_atoms));
			if(verbose>6) cout << " " << Sk_tmp ;
			Sk[nx]+=Sk_tmp;
			Fourth[nx]+=Fourth_tmp;
		}
		if(verbose>6) cout << endl;

	}

	// Print the averaged Sk
	FILE * ofile = fopen(ofilename.c_str(),"w");
	fprintf (ofile, "# [k, S=<|r^2|> , <|r^4|> , <|r^4|>/(2*<|r^2|>^2)-1 , nx ] \n");

	for(int nx=0;nx<num_nx;nx++){
		//cout << k1*nx << " " << Sk[nx]/(double)num_frames << " " << Fourth[nx]/(double)num_frames << endl;
		double m2 = Sk[nx]/(double)num_frames;
		double m4 = Fourth[nx]/(double)num_frames;
		fprintf (ofile, "%f %f %f %f %i\n", k1*(double)nx , m2 , m4  , m4/(2*m2*m2)-1.0 , nx );
	}

	fclose(ofile);
}



/**
 * Write van Hoove correlation function
 */
void Traj::incoherent_self_scattering(string ofilename,double q_vector){

	if(verbose>4) cout << endl << " ..:: Incoherent self-scattering function, F_s(t,q=" << q_vector << ") (" << ofilename << ") ::.." << endl;

	double *Fs = new double[num_frames];

	get_self_intermediate_scattering_function(
		num_frames,
		num_atoms,
		bbox,
		img_coords,
		coords,
		Fs,
		q_vector,
		get_variable(variables,"trestart_Fs",10)
	);

	// Print characteristic times where
	double incoherent_self_scattering_y0_shift = get_variable(variables,"incoherent_self_scattering_y0_shift",0.1);
	if(verbose>4){
		for ( double y0 = 1.0 ; y0 > 0.0 ; y0-=incoherent_self_scattering_y0_shift ){
			cout << "Fs( t = " << frame_dt*first_x_where_y_equals_y0(Fs,num_frames-1,y0) << " , q = " << q_vector << " ) = " << y0 << endl ;
		}
	}

	double time_one_over_e = first_x_where_y_equals_y0(Fs,num_frames-1,exp(-1));
	time_one_over_e *= frame_dt;
	if(verbose>4){
		cout << "Fs( t = " << time_one_over_e << " , q = " << q_vector << " ) = " << exp(-1) << endl;
	}

	FILE * ofile = fopen(ofilename.c_str(),"w");
	fprintf (ofile, "# Incoherent self-scattering function, F_s(t,q=%f). [ time , F_s ] \n",q_vector);

	for (int frame=0 ; frame < num_frames-1 ; frame++ ){
		fprintf ( ofile , "%f %f\n",frame_dt*(double)frame,Fs[frame]);
	}

	fclose(ofile);
}


/**
 * Write van Hoove correlation function
 */
void Traj::incoherent_self_scattering_of_selected(string ofilename,double q_vector){

	if(verbose>4) cout << endl << " ..:: Incoherent self-scattering function, F_s(t,q=" << q_vector << ") (" << ofilename << ") ::.." << endl;

	double *Fs = new double[num_frames];

	get_self_intermediate_scattering_function_of_selected(
		num_frames,
		num_atoms,
		bbox,
		img_coords,
		coords,
		selected,
		Fs,
		q_vector,
		get_variable(variables,"trestart_Fs",10)
	);

	// Print characteristic times where
	double incoherent_self_scattering_y0_shift = get_variable(variables,"incoherent_self_scattering_y0_shift",0.1);
	if(verbose>4){
		for ( double y0 = 1.0 ; y0 > 0.0 ; y0-=incoherent_self_scattering_y0_shift ){
			cout << "Fs( t = " << frame_dt*first_x_where_y_equals_y0(Fs,num_frames-1,y0) << " , q = " << q_vector << " ) = " << y0 << endl ;
		}
	}

	double time_one_over_e = first_x_where_y_equals_y0(Fs,num_frames-1,exp(-1));
	time_one_over_e *= frame_dt;
	if(verbose>4){
		cout << "Fs( t = " << time_one_over_e << " , q = " << q_vector << " ) = " << exp(-1) << endl;
	}

	FILE * ofile = fopen(ofilename.c_str(),"w");
	fprintf (ofile, "# Incoherent self-scattering function, F_s(t,q=%f). [ time , F_s ] \n",q_vector);

	for (int frame=0 ; frame < num_frames-1 ; frame++ ){
		fprintf ( ofile , "%f %f\n",frame_dt*(double)frame,Fs[frame]);
	}

	fclose(ofile);
}


/**
 * Return distance
 */
double Traj::get_dr(int atom0 ,int atom1, int frame0 ,int frame1){
	// TODO return error if atoms or frames are out of boundaries
	return dr(atom0,atom1,frame0,frame1,num_frames,num_atoms,bbox,coords);
}

/**
 * Return distance with the minimal image convention
 */
double Traj::get_dr_img(int atom0 ,int atom1, int frame){
	// TODO return error if atoms or frames are out of boundaries
	return dr_img(atom0,atom1,frame,num_frames,  num_atoms,  bbox,  coords);
}

/**
 * Return coordinate in primary image
 */
double Traj::get_x_img(int dim, int atom ,int frame){
	// TODO return error if atoms or frames are out of boundaries
	double x = coords[frame*num_atoms*NUM_DIM+atom*NUM_DIM+dim];
	x-= bbox[dim]*floor(x/bbox[dim]);
	return x;
}

/**
 * Return distance in x-dimension using minimum image convention. ASSUME fixed box size
 */
double Traj::get_dx_img(int dim, int atom0, int atom1 ,int frame){
	// TODO return error if atoms or frames are out of boundaries
	double x0 = get_x_img(dim,atom0,frame);
	double x1 = get_x_img(dim,atom1,frame);
	double dx = x1-x0;
	dx -= bbox[dim]*floor( (dx/bbox[dim]) + 0.5 );
	return dx;
}

/**
 * Return distance in x-dimension using minimum image convention. ASSUME fixed box size
 */
double Traj::get_dx_img(int dim, int atom0, int atom1 ,int frame0,int frame1){
	// TODO return error if atoms or frames are out of boundaries
	double x0 = get_x_img( dim, atom0, frame0 );
	double x1 = get_x_img( dim, atom1, frame1 );
	double dx = x1-x0;
	dx -= bbox[dim]*floor( (dx/bbox[dim]) + 0.5 );
	return dx;
}

double Traj::get_dr_img(int atom0, int atom1 ,int frame0,int frame1){
	double dx= 0.0;
	double dr = 0.0;
	for (int dim=0;dim<NUM_DIM;dim++){
		dx = get_dx_img(dim,atom0,atom1,frame0,frame1);
		dr += dx*dx;
	}
	return pow(dr,0.5);
}

/**
 * Return time of this trajectory
 */
double Traj::get_time(int f){
	return frame_dt*(double)f;
}


/**
 * Deselect all atoms. TODO and bonds.
 */
void Traj::deselect_all(){
	for ( int frame=0 ; frame < num_frames ; frame++) {
		for ( int atom=0 ; atom < num_atoms ; atom++) {
			selected[num_atoms*frame + atom] = false;
		}
	}
}

/**
 * Convert a (full) trajectory of water, into one where only the molecules exits.
 * Return vector with information about network.
 *
 * 	n_type = 0 when for the doner
 *    and
 *  n_type = 1 for the accepter
 *
 * Elements in nl_selected is true for bonded elements
 *
 */
void Traj::construct_hbond_network() {

	//nl_selected

	if(verbose>4) cout << endl << " ..:: Construct Hydrogen-bond network ::.." << endl;

	int atoms_per_mol = 3;
	int num_mol = num_atoms/atoms_per_mol;
	if ( verbose>4 & !(num_atoms%atoms_per_mol==0) ) cout << "Warning: !(num_atoms%atoms_per_mol==0), this may not be a (pure) water trajectory" << endl;

	// Definition of hydrogen bond
	if(verbose>4) {
		cout << "   If h_bond_definition < 1 then no H-bonds." << endl;
		cout << "   If h_bond_definition = 1 then a H-bond is when angle_OOH < max_ooh_angle and distance_OO < max_oo_distance." << endl;
		cout << "   If h_bond_definition = 2 then a H-bond is when beta_h = cosine(angle_OOH) - h_bond_alpha * OO_distance > h_bond_beta_cut." << endl;
		cout << "   Note: the *_bound variables are for a more restrict definition of H-Bonds in the bounded state." << endl;
	}

	int h_bond_definition = get_variable(variables,"h_bond_definition",1);
	//int h_bond_definition_ts = get_variable(variables,"h_bond_definition_ts",h_bond_definition);

	double max_oo_distance = get_variable(variables,"max_oo_distance",3.25);
	double max_ooh_angle = get_variable(variables,"max_ooh_angle",30.0);
	double max_oo_distance_bound = get_variable(variables,"max_oo_distance_bound",max_oo_distance);
	double max_ooh_angle_bound = get_variable(variables,"max_ooh_angle_bound",max_ooh_angle);
	double min_ooh_cos_angle = cos(max_ooh_angle/180.0*PI);
	double min_ooh_cos_angle_bound = cos(max_ooh_angle_bound/180.0*PI);
	if(verbose>4) cout << "cos(max_ooh_angle/180.0*PI) = min_ooh_cos_angle=" <<  min_ooh_cos_angle << endl;
	if(verbose>4) cout << "cos(max_ooh_angle_ts/180.0*PI) = min_ooh_cos_angle_bound=" <<  min_ooh_cos_angle_bound << endl;

	double h_bond_alpha = get_variable( variables , "h_bond_alpha" , 0.25 );
	double h_bond_beta_cut = get_variable( variables, "h_bond_beta_cut" , -0.125 );
	double h_bond_beta_cut_bound = get_variable( variables, "h_bond_beta_cut_bound" , 0.075 );

	// Make file with hbond angles and lengths histogram
	int num_bins_cos_ooh_angle = get_variable(variables,"num_bins_cos_ooh_angle",200);
	double max_oo_distance_in_histogram = get_variable(variables,"max_oo_distance_in_histogram",8.0); // TODO  make sure that
	int num_bins_oo_distance = get_variable(variables,"num_bins_oo_distance",(int)floor(max_oo_distance_in_histogram*20.0));
	double free_energy_output_for_inf = get_variable(variables,"free_energy_output_for_inf",20.0);

	if(verbose > 4 && max_oo_distance_in_histogram<max_oo_distance ){
		cout << "Warning: max_oo_distance_in_histogram < max_oo_distance is not allowed, max_oo_distance_in_histogram have been reset to 2.0*max_oo_distance" << endl;
		max_oo_distance_in_histogram = 2.0*max_oo_distance;
	}

	int* hbond_geometry_histogram = new int[num_bins_cos_ooh_angle*num_bins_oo_distance];
	for ( int i = 0 ; i < num_bins_cos_ooh_angle*num_bins_oo_distance ; i++ ) hbond_geometry_histogram[i]=0;

	// File with histogram of h_bond_beta's
	double min_h_bond_beta_in_histogram = get_variable ( variables , "min_h_bond_beta_in_histogram" , -0.5 );
	double max_h_bond_beta_in_histogram = get_variable(variables , "max_h_bond_beta_in_histogram" , 0.5 );
	double width_of_h_bond_beta_in_histogram = max_h_bond_beta_in_histogram - min_h_bond_beta_in_histogram;
	int num_bins_h_bond_beta = get_variable ( variables , "num_bins_h_bond_beta" , (int)(width_of_h_bond_beta_in_histogram*200.0) );
	int * h_bond_beta_histogram = new int[num_bins_h_bond_beta];
	for ( int i = 0 ; i < num_bins_h_bond_beta ; i++ ) h_bond_beta_histogram[i]=0;

	// Arrays with 'hydrogen bond' vectors
	double* oo_vector = new double [NUM_DIM];		// For doner bonds
	double* ohA_vector= new double [NUM_DIM];
	double* ohB_vector= new double [NUM_DIM];

	double* oo_vector_reverse = new double [NUM_DIM];	// For accepter bonds
	double* ohC_vector= new double [NUM_DIM];
	double* ohD_vector= new double [NUM_DIM];


	// Clear neighbor list
	nl.clear();
	nl_first.clear();
	nl_last.clear();
	nl_num.clear();
	n_type.clear();
	nl_selected.clear();
	int max_number_of_hbonds=0;

	// Make neighbor list
	// N^2 search algorithm TODO make cell list for finding hydrogen-bonded neighbors if this is to slow
	if(h_bond_definition==0){
		// Make empty neighbor list
		nl.resize(1,0);
		nl_first.resize(num_mol*num_frames,0);
		nl_last.resize(num_mol*num_frames,0);
		nl_num.resize(num_mol*num_frames,0);
		n_type.resize(num_mol*num_frames,0);
		nl_selected.resize(num_mol*num_frames,false);
	}else{
		for ( int frame=0 ; frame < num_frames ; frame++) {
			for ( int mol0 = 0 ; mol0 < num_mol ; mol0++ ) {
				nl_first.push_back(nl.size());
				for (int mol1 = 0 ; mol1 < num_mol ; mol1++ ) {

					if(!(mol0==mol1)){

						// Calculate O-O distance
						double oo_length = get_dr_img ( mol0*atoms_per_mol , mol1*atoms_per_mol , frame ) ;

						if ( oo_length < max_oo_distance_in_histogram ) {

							// Calculate OOH angle
							for(int d=0;d<NUM_DIM;d++) oo_vector[d]=get_dx_img ( d, mol0*atoms_per_mol , mol1*atoms_per_mol , frame);
							for(int d=0;d<NUM_DIM;d++) ohA_vector[d]=get_dx_img( d, mol0*atoms_per_mol , mol0*atoms_per_mol+1 , frame);
							for(int d=0;d<NUM_DIM;d++) ohB_vector[d]=get_dx_img( d, mol0*atoms_per_mol , mol0*atoms_per_mol+2 , frame);

							for(int d=0;d<NUM_DIM;d++) oo_vector_reverse[d]=get_dx_img ( d, mol1*atoms_per_mol , mol0*atoms_per_mol , frame);
							for(int d=0;d<NUM_DIM;d++) ohC_vector[d]=get_dx_img( d, mol1*atoms_per_mol , mol1*atoms_per_mol+1 , frame);
							for(int d=0;d<NUM_DIM;d++) ohD_vector[d]=get_dx_img( d, mol1*atoms_per_mol , mol1*atoms_per_mol+2 , frame);

							double cosOOHA = get_cos_theta( oo_vector         , ohA_vector );
							double cosOOHB = get_cos_theta( oo_vector         , ohB_vector );
							double cosOOHC = get_cos_theta( oo_vector_reverse , ohC_vector );
							double cosOOHD = get_cos_theta( oo_vector_reverse , ohD_vector );

							// Add configuration to 'geometry' histogram (C and D hydrogens will be added as some point)
							hbond_geometry_histogram[	(int)floor((cosOOHA+1.0)/2.0*(double)num_bins_cos_ooh_angle) * num_bins_oo_distance
													 +	(int)floor( oo_length/max_oo_distance_in_histogram*(double)num_bins_oo_distance )  ]++;
							hbond_geometry_histogram[	(int)floor((cosOOHB+1.0)/2.0*(double)num_bins_cos_ooh_angle) * num_bins_oo_distance
													 +	(int)floor( oo_length/max_oo_distance_in_histogram*(double)num_bins_oo_distance )  ]++;


							// Calculate beta_H, i.e. the 'r_OO intersection's ... and add to histogram
							double h_bond_beta_A = cosOOHA - h_bond_alpha * oo_length;
							double h_bond_beta_B = cosOOHB - h_bond_alpha * oo_length;
							double h_bond_beta_C = cosOOHC - h_bond_alpha * oo_length;
							double h_bond_beta_D = cosOOHD - h_bond_alpha * oo_length;
							int h_bond_beta_A_index =  (int)((double)num_bins_h_bond_beta*( h_bond_beta_A - min_h_bond_beta_in_histogram ) / width_of_h_bond_beta_in_histogram);
							int h_bond_beta_B_index =  (int)((double)num_bins_h_bond_beta*( h_bond_beta_B - min_h_bond_beta_in_histogram ) / width_of_h_bond_beta_in_histogram);
							if ( h_bond_beta_A_index>0 && h_bond_beta_A_index < num_bins_h_bond_beta ) h_bond_beta_histogram[h_bond_beta_A_index]++;
							if ( h_bond_beta_B_index>0 && h_bond_beta_B_index < num_bins_h_bond_beta ) h_bond_beta_histogram[h_bond_beta_B_index]++;


							// Add H-Bond to neighbor list, dependent on the used Hbonds definition
							if(h_bond_definition==0){
								// No H-bonds selected
							}else if(h_bond_definition==1){

								if( ( cosOOHA > min_ooh_cos_angle || cosOOHB > min_ooh_cos_angle ) && oo_length < max_oo_distance ){	// Donor H-Bond (type = 0)
									nl.push_back(mol1);
									n_type.push_back(0);
									if( (      cosOOHA > min_ooh_cos_angle_bound || cosOOHB > min_ooh_cos_angle_bound )
											&& oo_length < max_oo_distance_bound ){nl_selected.push_back(true);}else{nl_selected.push_back(false);}
								}

								if( ( cosOOHC > min_ooh_cos_angle || cosOOHD > min_ooh_cos_angle ) && oo_length < max_oo_distance ){	// Accepter H-Bond (type = 1)
									//if(output_hbond_geometry==1) fprintf (ofile_hbond_geometry, "%f %f %f\n",oo_length,cosOOHA,cosOOHB);
									nl.push_back(mol1);
									n_type.push_back(1);
									if( (      cosOOHC > min_ooh_cos_angle_bound || cosOOHD > min_ooh_cos_angle_bound )
											&& oo_length < max_oo_distance_bound  ){nl_selected.push_back(true);}else{nl_selected.push_back(false);}
								}

							}else if(h_bond_definition==2){

								if ( (h_bond_beta_A > h_bond_beta_cut) || (h_bond_beta_B > h_bond_beta_cut) ){
									nl.push_back(mol1);
									n_type.push_back(0);
									if ( (h_bond_beta_A > h_bond_beta_cut_bound) || (h_bond_beta_B > h_bond_beta_cut_bound) ){nl_selected.push_back(true);}else{nl_selected.push_back(false);}
								}

								if ( (h_bond_beta_C > h_bond_beta_cut) || (h_bond_beta_D > h_bond_beta_cut) ){
									nl.push_back(mol1);
									n_type.push_back(1);
									if ( (h_bond_beta_C > h_bond_beta_cut_bound) || (h_bond_beta_D > h_bond_beta_cut_bound) ){nl_selected.push_back(true);}else{nl_selected.push_back(false);}
								}

							}else{
								// We don't expect to program to go here
							}
						}
					}
				}

				nl_last.push_back(nl.size());
				nl_num.push_back( nl_last.back() - nl_first.back() );
				if( max_number_of_hbonds<nl_num.back() ) max_number_of_hbonds=nl_num.back();
			}
		} // END for ( int frame=0 ; frame < num_frames ; frame++)
	}

	// Print H-Bond geometry histogram to file
	FILE * ofile_hbond_distribution = fopen("hbond_distribution.dat","w");
	fprintf (ofile_hbond_distribution, "# [OO distance, OOH angle, number of bonds ] \n");

	FILE * ofile_hbond_free_energy = fopen("hbond_free_energy.dat","w");
	fprintf (ofile_hbond_free_energy, "# [OO distance, OOH angle, Free energy relative to infinite oxygen separation in units where kT = %f ] \n",kT);

	double rho_avg = (double)num_mol/get_average_box_volume();
	for (int i = 0 ; i < num_bins_cos_ooh_angle ; i++ ){

		fprintf (ofile_hbond_distribution , "\n");	// Blank line for gnu-plot
		fprintf (ofile_hbond_free_energy , "\n");

		for ( int j = 0 ; j < num_bins_oo_distance ; j++ ){

			fprintf (ofile_hbond_distribution , "%f %f %f\n" ,
					(double)j/(double)num_bins_oo_distance*max_oo_distance_in_histogram ,
					((double)i/(double)num_bins_cos_ooh_angle*2.0-1.0) ,
					(double)hbond_geometry_histogram[	i * num_bins_oo_distance + j ]/(double)num_frames/(double)num_mol );


			double r0 = (double)(j+0)/(double)num_bins_oo_distance*max_oo_distance_in_histogram;
			double r1 = (double)(j+1)/(double)num_bins_oo_distance*max_oo_distance_in_histogram;
			double bin_volume = 4.0/3.0*PI*( r1*r1*r1 - r0*r0*r0 );
			bin_volume *= 2.0/(double)num_bins_cos_ooh_angle;

			double free_energy_of_bin = (double)hbond_geometry_histogram[	i * num_bins_oo_distance + j ];
			free_energy_of_bin /= (double)num_frames*(double)num_mol;
			free_energy_of_bin /= bin_volume;
			free_energy_of_bin /= rho_avg;
			free_energy_of_bin = -log(free_energy_of_bin);
			free_energy_of_bin *= kT;

			if(hbond_geometry_histogram[	i * num_bins_oo_distance + j ]>0){
				fprintf (ofile_hbond_free_energy , "%f %f %f\n" ,
						(double)j/(double)num_bins_oo_distance*max_oo_distance_in_histogram ,
						((double)i/(double)num_bins_cos_ooh_angle*2.0-1.0) ,
						free_energy_of_bin );
			}else{
				fprintf (ofile_hbond_free_energy , "%f %f %f\n" ,
						(double)j/(double)num_bins_oo_distance*max_oo_distance_in_histogram ,
						((double)i/(double)num_bins_cos_ooh_angle*2.0-1.0) ,
						free_energy_output_for_inf );
			}




		}

	}
	fclose(ofile_hbond_distribution);
	fclose(ofile_hbond_free_energy);


	// Print h_bond_beta histogram
	FILE * ofile_h_bond_beta = fopen("hbond_beta_parameter.dat","w");
	fprintf (ofile_h_bond_beta, "# [h_bond_beta=cosine(angle_OOH)-%f*distance_OO , frequency per frame per particle ]\n",h_bond_alpha);

	for ( int i = 0 ; i < num_bins_h_bond_beta ; i++){
		double volume_element = (double)num_frames*(double)num_atoms/width_of_h_bond_beta_in_histogram;
		fprintf (ofile_h_bond_beta, "%f %f\n"
				, width_of_h_bond_beta_in_histogram*(double)i/(double)num_bins_h_bond_beta + min_h_bond_beta_in_histogram
				, (double)h_bond_beta_histogram[i]/volume_element );
	}
	fclose(ofile_h_bond_beta);


	// Print some statistics on H-Bonds
	if ( verbose > 4 ) {
		cout << endl << "Found " << nl.size()/2 << " hydrogen bonds in total, i.e. " << (double)nl.size()/(double(num_frames))/2.0 << " per frame or " << (double)nl.size()/((double)num_frames)/((double)num_mol)/2.0 << " per molecule." << endl;
		int *count_all = new int[max_number_of_hbonds+1];
		for ( int i = 0 ; i < max_number_of_hbonds+1 ; i++ ) count_all[i]=0;
		for ( unsigned i = 0 ; i<nl_num.size() ; i++ ) count_all[nl_num.at(i)]++;
		for ( int i = 0 ; i<max_number_of_hbonds+1 ; i++ ) if(count_all[i]>0) cout << count_all[i] << " waters with " << i << " h-bonds, i.e a fraction of  " << (double)count_all[i]/(double)(num_mol*num_frames) << endl;

		int *count_ts = new int[max_number_of_hbonds+1];
		for ( int i = 0 ; i < max_number_of_hbonds+1 ; i++ ) count_ts[i]=0;
		int *count_bound = new int[max_number_of_hbonds+1];
		for ( int i = 0 ; i < max_number_of_hbonds+1 ; i++ ) count_bound[i]=0;
		int total_num_bound = 0;
		int total_num_ts = 0;
		for ( int fp = 0 ; fp < num_frames*num_mol ; fp++ ){
			int this_num_ts = 0;
			int this_num_bound = 0;
			for ( int i = nl_first.at(fp) ; i < nl_last.at(fp) ; i++ ){
				if(nl_selected.at(i)){this_num_bound++;}else{this_num_ts++;}
			}
			count_ts[this_num_ts]++;
			count_bound[this_num_bound]++;
			total_num_bound += this_num_bound;
			total_num_ts += this_num_ts;
		}

		cout << endl << "Found " << total_num_bound/2 << " bounded hydrogen bonds in total, i.e. " << (double)total_num_bound/(double(num_frames))/2.0 << " per frame or " << (double)total_num_bound/((double)num_frames)/((double)num_mol)/2.0 << " per molecule." << endl;
		for ( int i = 0 ; i<max_number_of_hbonds+1 ; i++ ) if(count_bound[i]>0) cout << count_bound[i] << " waters with " << i << " bound state h-bonds, i.e a fraction of  " << (double)count_bound[i]/(double)(num_mol*num_frames) << endl;

		cout << endl << "Found " << total_num_ts/2 << " transition state hydrogen bonds in total, i.e. " << (double)total_num_ts/(double(num_frames))/2.0 << " per frame or " << (double)total_num_ts/((double)num_frames)/((double)num_mol)/2.0 << " per molecule." << endl;
		for ( int i = 0 ; i<max_number_of_hbonds+1 ; i++ ) if(count_ts[i]>0) cout << count_ts[i] << " waters with " << i << " transition state h-bonds, i.e a fraction of  " << (double)count_ts[i]/(double)(num_mol*num_frames) << endl;

	}


	// Re-make position arrays so they only contains oxygen atoms
	{
		int *new_type;							// Atom type
		new_type = new int[num_mol];
		for (int i = 0 ; i<num_mol ; i++) new_type[i]=type[i*atoms_per_mol];

		int *new_img_coords;					// image coordinates
		new_img_coords = new int[num_frames*num_mol*NUM_DIM];

		double *new_coords;						// Coordinates
		new_coords = new double[num_frames*num_mol*NUM_DIM];

		for ( int f = 0 ; f < num_frames ; f ++){
			for ( int mol = 0 ; mol < num_mol ; mol++ ) {
				for ( int d = 0 ; d<NUM_DIM ; d++){

					int i_in_old = f*num_mol*atoms_per_mol*NUM_DIM + mol*atoms_per_mol*NUM_DIM + d;
					int i_in_new = f*num_mol*NUM_DIM + mol*NUM_DIM + d;

					new_img_coords[ i_in_new ] = img_coords [ i_in_old ] ;

					new_coords[ i_in_new ] = coords [ i_in_old ] ;
				}
			}
		}

		delete[] type;
		type = new int[num_mol];
		for ( int i = 0 ; i < num_mol ; i++) type[i]=new_type[i];

		delete[] img_coords;
		img_coords = new int[num_frames*num_mol*NUM_DIM];
		for ( int i = 0 ; i < num_frames*num_mol*NUM_DIM ; i++) img_coords[i]=new_img_coords[i];

		delete[] coords;
		coords = new double[num_frames*num_mol*NUM_DIM];
		for ( int i = 0 ; i < num_frames*num_mol*NUM_DIM ; i++) coords[i]=new_coords[i];

		if (verbose > 4) cout << "Reset num_atoms from " << num_atoms << " to " << num_mol << endl;
		num_atoms = num_mol;

		// num_types = ### TODO num_types should be set
	}



}


/**
 * Enumerate bonds in neighbor list.
 *
 *   index = index_of_doner_atom * number_of_atoms + index_of_accepter_atom
 *
 *	The negative value is given to accepters oxygens.
 *
 *  When calling this function n_type should be contain 0 or 1, for H-bond donors and accepters respectively.
 */
void Traj::enumerate_nl_hbonds() {
	if(verbose>4) cout << endl << " ..:: Enumerate hydrogen bonds with i = index_of_doner_atom * number_of_atoms + index_of_accepter_atom  ::.." << endl;

	// Enumerate
	for (int f = 0 ; f < num_frames ; f++ ){
		for (int a = 0 ; a < num_atoms ; a++ ){
			for ( int n = nl_first.at(f*num_atoms+a) ; n < nl_last.at(f*num_atoms+a) ; n++ ) {
				if ( n_type.at(n) == 0 ){	// this is a doner water
					n_type.at(n) = a*num_atoms+nl.at(n);
				} else if ( n_type.at(n) == 1 ) { // this is a accepter water
					n_type.at(n) = - nl.at(n)*num_atoms+a;
				}
			}
		}
	}

	// Write the n_type to a file (hbond_index.dat)
	if(verbose>4) cout << "Write hbond_index.dat if write_hbond_index=1 (note that this file might get large)" << endl;
	int write_hbond_index = get_variable(variables,"write_hbond_index",1);
	if (write_hbond_index == 1){
		FILE * ofile = fopen("hbond_index.dat","w");
		fprintf (ofile, "# [ frame ; hydrogen bond index ; atom ]\n");

		for ( int f = 0 ; f < num_frames ; f++ ){
			for (int a = 0 ; a < num_atoms ; a++ ){
				for ( int n = nl_first.at(f*num_atoms+a) ; n < nl_last.at(f*num_atoms+a) ; n++ ) {
					fprintf (ofile, "%d %d %d \n", f , n_type.at(n) , a );
				}
			}
		}
		fclose(ofile);
	}

}



/**
 * Enumerate hydrogen bonds in neighbor list.
 *
 * This functions redefine the array stating the neighbor type:
 * 		When calling this function n_type should be contain 0 or 1, for H-bond donors and accepters respectively. .
 * 		After, the donors have n_types > 0, and accepters have n_type < 0.
*/
void Traj::enumerate_nl_hbonds_by_occurrence() {

	if(verbose>4) cout << endl << " ..:: Enumerate hydrogen bonds by occurrence !!! TODO this function is not written !!!!  ::.." << endl;

	int counter = 0;	// Counter to enumerate H-bonds.

	// Enumerate hydrogen bonds in first frame.
	for (int a = 0 ; a < num_atoms ; a++ ){
		for ( int n = nl_first.at(a) ; n < nl_last.at(a) ; n++ ) {

			// If doner, bond type is set to (positive) counter value and accepter is set no negative counter value.
			if( n_type.at(n) == 0 ) {
				counter++;
				n_type.at(n)=counter;

				// Find the accepter in neighbor list.
				int this_counter_should_become_one = 0;
				for ( int nn=nl_first.at( nl.at(n) ) ; nn < nl_last.at( nl.at(n) ) ; nn++  ) {
					if ( nl.at(nn) == a ){
						n_type.at(nn)=-counter;
						this_counter_should_become_one++;
					}
				}
				if (this_counter_should_become_one!=1 && verbose > 3 ){
					cout << "Warning: Inconsistent found when looking for accepter hydrogen to atom " << a << endl;
					exit(0);
				}
			}
		}
	}

	// Enumerate hydrogen bonds in remaining frames.  TODO


	// Write the n_type to a file (hbond_index.dat)
	if(verbose>4) cout << "Write hbond_index.dat if write_hbond_index=1" << endl;
	int write_hbond_index = get_variable(variables,"write_hbond_index",1);
	if (write_hbond_index == 1){
		FILE * ofile = fopen("hbond_index.dat","w");
		fprintf (ofile, "# [ time ; hydrogen bond index ]\n");

		for ( int f = 0 ; f < num_frames ; f++ ){
			for ( int n = nl_first.at(f*num_atoms) ; n < nl_last.at((f+1)*num_atoms-1) ; n++ ) {
				fprintf (ofile, "%f %d\n", frame_dt*(double)f , n_type.at(n) );
			}
		}
		fclose(ofile);
	}
}


/**
 * Print <h_i(0)h_i(t)>_i where h_i(t)=1 if bond i is present at time t and h_i(t)=0 othervise.
 */
void Traj::nl_bond_time_autocorrelation(string ofilename) {

	if(verbose>4) cout << endl << " ..:: Hydrogen-bond time autocorrelation ::.." << endl;

	FILE * ofile = fopen(ofilename.c_str(),"w");
	fprintf (ofile, "# HBond time autocorrelation [t; < h_i(0) h_i(t) > ]\n");

	int nl_bond_time_autocorrelation_frame_stride = get_variable(variables,"nl_bond_time_autocorrelation_frame_stride",1);

	for ( int df = 0 ; df < num_frames ; df++ ) {

		int num_hbonds = 0;
		int num_persistent_hbonds = 0;

		for ( int f0 = 0 ; f0 < num_frames-df ; f0 += nl_bond_time_autocorrelation_frame_stride ) {
			for ( int a = 0 ; a < num_atoms ; a++ ){
				for ( int n = nl_first.at(f0*num_atoms+a) ; n < nl_last.at(f0*num_atoms+a) ; n++) {
					num_hbonds++;

					// See if this bond is present in f1=f0+df frame
					bool this_bonds_persist = false;
					for ( int n_future = nl_first.at((f0+df)*num_atoms+a) ; n_future < nl_last.at((f0+df)*num_atoms+a) ; n_future++) {
						if( (nl.at(n) == nl.at(n_future)) && (n_type.at(n) == n_type.at(n_future)) ){
							if(this_bonds_persist){
								cout << "Warning: Persisting hbond found twice," << endl;
							}else{
								num_persistent_hbonds ++;
								this_bonds_persist = true;
							}
						}
					}
				}
			}
		}

		fprintf (ofile, "%f %f\n", frame_dt*(double)df,(double)num_persistent_hbonds/(double)num_hbonds);
	}

	fclose(ofile);
}

/**
 * Find Hydrogen-bond braking and formations,
 * e.i. excitations, by finding trajectories that travel though transitions states
 * between bound and unbound.
 *
 * The neighbor bool vector nl_selected elements should be true for the neighbor list elements
 *   that are in the bounded state. The remaining H-Bonds are assumed to be in (a) transition state.
 *
 * After running this function, 'true' elements in the n_is_bound is
 * referes to H-Bonds that are moving though the transition region.
 *
 *    Change type list.
 * type = 0 : not in transition
 * type = 1 : in transition
 *
 */
void Traj::hbond_dynamic_exitations_ts() {
	if(verbose>4) cout << endl << " ..:: Find Hydrogen-bond braking and formations, e.i. excitations, by finding trajectories that travel though transitions states between bound and unbound ::.." << endl;

	// Swap the content of nl_selected into n_is_bound. At the end of this function, the content of is_in_transit is swapped into nl_selected.
	vector<bool> is_in_transit;
	vector<bool> n_is_bound;
	n_is_bound.swap(nl_selected);

	// File with all h-bonds that visit the transition region
	int write_exitations_to_file = get_variable(variables,"write_exitations_to_file",1);
	FILE * ofile_exitations = fopen("exitations.dat","w");
	if ( write_exitations_to_file == 1){
		fprintf (ofile_exitations, "# [ frame ; hydrogen bond index of excitation ]\n");
	} else {
		fprintf (ofile_exitations, "# This file is blank since write_exitations_to_file != 1");
	}

	// File with all h-bonds that transit the transition region region
	int write_exitations_transit_to_file = get_variable(variables,"write_exitations_transit_to_file",1);
	FILE * ofile_exitations_transit = fopen("exitations_transit.dat","w");
	if ( write_exitations_transit_to_file == 1){
		fprintf (ofile_exitations_transit, "# [ frame ; hydrogen bond index of excitation ]\n");
	} else {
		fprintf (ofile_exitations_transit, "# This file is blank since write_exitations_transit_to_file != 1");
	}

	// Variables for collecting statistics
	int num_transit_hbonds = 0;
	int instantaneous_hbond_transitions = 0;
	int * length_of_transit_region = new int[num_frames];
	for (int i = 0 ; i < num_frames ; i++ ) length_of_transit_region[i]=0;

	// Find paths that goes though the transition region from the bounded state to the un-bounded (or visor versa)
	for (int f = 0 ; f < num_frames ; f++ ){
		for (int a = 0 ; a < num_atoms ; a++){
			bool have_a_bond_in_transit = false;
			int fa = f*num_atoms+a;
			for (int n = nl_first.at(fa) ; n < nl_last.at(fa) ; n ++ ){

				// See if this is an instantaneous transition
				if(        f > 0
						&& f<num_frames-1
						&& n_is_bound.at(n)
						&& (
							   !do_n_type_exist_in_frame (n_type.at(n),f-1,a)
							|| !do_n_type_exist_in_frame (n_type.at(n),f+1,a)
							)
						){
					instantaneous_hbond_transitions++;
					length_of_transit_region[0]++;
				}


				// Try to find the time where this trajectory entered the transition region
				bool tmp_is_in_transit = false;
				if(!n_is_bound.at(n)){
					int n_index = n_type.at(n);

					int f_past;
					bool was_bound_in_past = false;
					for ( f_past = f ; f_past >= 0 ; f_past-- ){
						int n_ref = n_type_in_frame ( n_index , f_past , a );
						if ( n_ref == -1 ){ break ; }
						if ( n_ref > -1 && n_is_bound.at( n_ref ) ){ was_bound_in_past=true ; break ; }
					}

					int f_future;
					bool was_bound_in_future = false;
					for ( f_future = f ; f_future < num_frames ; f_future++ ){
						int n_ref = n_type_in_frame ( n_index , f_future , a );
						if ( n_ref == -1 ){ break ; }
						if ( n_ref > -1 && n_is_bound.at( n_ref ) ){ was_bound_in_future=true ; break ; }
					}

					if ( f_past >= 0 &&
						 f_future<num_frames &&
						 was_bound_in_past^was_bound_in_future ){
						tmp_is_in_transit = true;
						//cout << "Found transition path: n_type = " << n_type.at(n) << " from frame " << f_past << " to " <<  f_future << ". df = " << f_future - f_past - 1 << endl;
						num_transit_hbonds++;
						length_of_transit_region[ f_future - f_past - 1 ]++ ;
						if ( write_exitations_transit_to_file==1 ) fprintf ( ofile_exitations_transit , "%d %d\n" , f , n_type.at(n) );
					}
				}
				if( write_exitations_to_file == 1 && !n_is_bound.at(n) ){
					fprintf ( ofile_exitations , "%d %d\n" , f , n_type.at(n) );
				}

				is_in_transit.push_back(tmp_is_in_transit);

				have_a_bond_in_transit = ( have_a_bond_in_transit || tmp_is_in_transit);

				// Put index to one for oxygens with a bond in transit
				/*if(have_a_bond_in_transit){
					type[num_atoms*f+a] = 1;
					type[num_atoms*f+nl[n]] = 1;
				}*/
			} // end n
		} // end a
	} // end f
	fclose(ofile_exitations);
	fclose(ofile_exitations_transit);

	// Print some statistics
	if( verbose > 4 ) {
		cout << "num_transit_hbonds = " << num_transit_hbonds/2 << " or (estimated) " << (double)num_transit_hbonds/2/((double)num_frames-2.0)/frame_dt/(double)num_atoms << " transit H-Bonds per atom per time_unit " << endl;
		cout << "instantaneous_hbond_transitions = " << instantaneous_hbond_transitions << " ( N. B. instantaneous_hbond_transitions/num_transit_hbonds << 1 )" << endl;
	}

	FILE * ofile_length_of_transit_region = fopen("length_of_transit_region.dat","w");
	fprintf (ofile_length_of_transit_region, "# [ t ; unnormalized probability of being in a transition trajectory of length t ]\n");
	for (int f = 0 ; f < num_frames ; f++ )
		fprintf (ofile_length_of_transit_region, "%f %d\n",frame_dt*(double)f,length_of_transit_region[f] );
	;
	fclose(ofile_length_of_transit_region);

	// Let the nl_selected have true elements for neighbors that are in transit
	is_in_transit.swap(nl_selected);

	// make selected oxygens the ones that have at least one bonds in transit
	unsigned int number_of_selected_oxygen = 0;
	for (int f = 0 ; f < num_frames ; f++ ){
		for (int a = 0 ; a < num_atoms ; a++){
			selected[f*num_atoms+a] = false;
			for (int n = nl_first.at(f*num_atoms+a) ; n < nl_last.at(f*num_atoms+a) ; n ++ ){
				if(nl_selected.at(n)) selected[f*num_atoms+a] = true;
			}
			if(selected[f*num_atoms+a]) number_of_selected_oxygen++;
		}
	}
	if( verbose > 4 ) cout << number_of_selected_oxygen << " oxygens have one or more h-bonds in transit. These are now selected particles." << endl;
}



/**
 * Find and analyze excitations in the hydrogen network.
 */
void Traj::hbond_dynamic_exitations_sojourn() {
	if(verbose>4) cout << endl << " ..:: Find Hydrogen-bond braking and formations, e.i. excitations, using a \"sojourn\" criteria ::.." << endl;

	// File with output
	int write_exitations_to_file = get_variable(variables,"write_exitations_to_file",1);
	FILE * ofile_exitations = fopen("exitations.dat","w");
	if ( write_exitations_to_file == 1){
		fprintf (ofile_exitations, "# [ frame ; hydrogen bond index of excitation ]\n");
	} else {
		fprintf (ofile_exitations, "# This file is blank since write_exitations_to_file!=1");
	}

	int write_exitations_frame_to_file = get_variable(variables,"write_exitations_frame_to_file",1);
	FILE * ofile_exitations_frame = fopen("exitations_frame.dat","w");
	if ( write_exitations_to_file == 1){
		fprintf (ofile_exitations_frame, "# [ frame ; number of excitation in this frame ]\n");
	} else {
		fprintf (ofile_exitations_frame, "# This file is blank since write_exitations_frame_to_file!=1");
	}


	// Initialize variables
	if(verbose > 4) cout << "TODO if transition variables are not default values, transitions may be found multiple times. Fix This." << endl;
	int t_transition = get_variable(variables,"t_transition",1);
	int t_sojourn = get_variable(variables,"t_sojourn",4);

	int low_sojourn_integral_cut = get_variable(variables,"low_sojourn_integral_cut",4);
	int upper_sojourn_integral_cut = get_variable(variables,"upper_sojourn_integral_cut",0);

	int num_excitations_total = 0;

	// Find dynamic excitations loop
	for ( int f = t_transition+t_sojourn+t_transition ; f < num_frames-t_transition-t_sojourn ; f++ ){

		int num_excitations_frame = 0;

		for ( int a = 0 ; a < num_atoms ; a++ ){
			for ( int n = nl_first.at(f*num_atoms+a) ; n < nl_last.at(f*num_atoms+a) ; n++ ) {

				// Integrate back path
				int sum_back = 0 ;
				for ( int ff = f - t_transition ; ff > f - t_transition - t_sojourn ; ff-- ){
					if ( do_n_type_exist_in_frame ( n_type.at(n) , ff , a ) ) sum_back++;
				}

				// Integrate forward path
				int sum_forward = 0 ;
				for ( int ff = f + t_transition ; ff < f + t_transition + t_sojourn ; ff++ ){
					if ( do_n_type_exist_in_frame ( n_type.at(n) , ff , a ) ) sum_forward++;
				}

				//cout << sum_back << " " << sum_forward <<  endl;

				//
				if (	( sum_back >= low_sojourn_integral_cut && sum_forward <= upper_sojourn_integral_cut )
						||
						( sum_back <= upper_sojourn_integral_cut && sum_forward >= low_sojourn_integral_cut )
					) {
					num_excitations_frame++;
					fprintf ( ofile_exitations , "%d %d\n" , f , n_type.at(n) );
				}
			}
		}
		fprintf ( ofile_exitations_frame , "%d %d\n" , f , num_excitations_frame );
		if(verbose > 6 && write_exitations_frame_to_file==1 ) cout << num_excitations_frame << " dynamic exitations in frame " << f << endl;
		num_excitations_total += num_excitations_frame;
	}

	// Print statistics   TODO correct how number of frames is counted
	if( verbose > 4) {
		int num_frames_effectiv = num_frames-2*(t_sojourn+t_transition);
		cout << "Found " << num_excitations_total << " excitations, e.i. " << (double)num_excitations_total/(double)(num_frames_effectiv) << " per frame " << " or " << (double)num_excitations_total/(double)num_frames_effectiv/(double)num_atoms << " per atom " << endl;
		double excitation_rate = (double)num_excitations_total/(double)num_frames_effectiv/(double)num_atoms/frame_dt;
		cout << "Excitation rate (intensive): " << excitation_rate << endl;
		cout << "Excitation time (intensive): " << 1.0/excitation_rate << endl;
	}

	fclose(ofile_exitations);
	fclose(ofile_exitations_frame);
}



/**
 * Return true if 'type_index' exits in frame 'frame' in nl_type list
 */
bool Traj::do_n_type_exist_in_frame ( int type_index_to_finde , int frame , int atom ) {

	bool answer = false;

	for ( int n = nl_first.at(frame*num_atoms + atom) ; n < nl_last.at(frame*num_atoms+atom) ; n++ ) {
		if ( n_type.at(n) == type_index_to_finde )	answer = true;
	}

	return answer;
}



/**
 * Return the neighbor index if 'type_index' exits in frame 'frame' in nl_type list. Return -1 if neighbor type does not exist.
 */
int Traj::n_type_in_frame ( int type_index_to_finde , int frame , int atom ) {

	int answer = -1;

	for ( int n = nl_first.at(frame*num_atoms + atom) ; n < nl_last.at(frame*num_atoms+atom) ; n++ ) {
		if ( n_type.at(n) == type_index_to_finde )	answer = n;
	}

	return answer;
}



/**
 * Calculate the rotational auto-correlation function of the water dipole
 */
void Traj::water_dipole_rotational_correlation (  ) {

	if(verbose>4) cout << endl << " ..:: Calculate the rotational auto-correlation function of the water dipole (rotacf.dat) ::.." << endl;

	int atoms_per_mol = 3;
	int num_mol = num_atoms/atoms_per_mol;

	// Make dipole vectors. Assume that the first atom is oxygen
	double *dipole_vectors = new double[num_frames*num_mol*NUM_DIM];
	for (int frame = 0 ; frame < num_frames ; frame++ ) {
		for (int mol = 0 ; mol < num_mol ; mol++) {
			for (int dim = 0 ; dim < NUM_DIM ; dim++ ) {
				double x_HA = get_dx_img( dim, mol*atoms_per_mol+0 , mol*atoms_per_mol+1 , frame );
				double x_HB = get_dx_img( dim, mol*atoms_per_mol+0 , mol*atoms_per_mol+2 , frame );
				dipole_vectors[frame*num_mol*NUM_DIM+mol*NUM_DIM+dim] = x_HA + x_HB;
			}
		}
	}
	if(verbose>4) cout << "The first dipole vector is v_1( t = 0 ) = ( " << dipole_vectors[0] << " , " << dipole_vectors[1] << " , " <<  dipole_vectors[2] << " ) " <<  endl;

	// Calculate rotational auto correlation functions
	double *rotacf = new double[num_frames];
	double *rotacf2 = new double[num_frames];
	get_rotational_auto_correlation(
		num_frames,
		num_mol,
		dipole_vectors,
		rotacf,
		rotacf2,
		get_variable(variables,"rotacf_trestart",10)
	);

	// Print characteristic times
	double water_dipole_rotational_correlation_y0_shift = get_variable(variables,"water_dipole_rotational_correlation_y0_shift",0.1);
	if(verbose>4){

		for ( double y0 = 1.0 ; y0 > 0.0 ; y0 -= water_dipole_rotational_correlation_y0_shift ){
			cout << "C_1( t = " << frame_dt*first_x_where_y_equals_y0(rotacf,num_frames-1,y0) << " ) = " << y0 << endl ;
		}
		cout << endl;
	}
	if(verbose>4){
		for ( double y0 = 1.0 ; y0 > 0.0 ; y0-=0.1 ){
			cout << "C_2( t = " << frame_dt*first_x_where_y_equals_y0(rotacf2,num_frames-1,y0) << " ) = " << y0 << endl ;
		}
		cout << endl;
	}

	double C1_time_one_over_e = first_x_where_y_equals_y0(rotacf,num_frames-1,exp(-1));
	C1_time_one_over_e *= frame_dt;
	double C2_time_one_over_e = first_x_where_y_equals_y0(rotacf2,num_frames-1,exp(-1));
	C2_time_one_over_e *= frame_dt;
	if(verbose>4){
		cout << "C_1( t = " << C1_time_one_over_e << " ) = " << exp(-1) << endl;
		cout << "C_2( t = " << C2_time_one_over_e << " ) = " << exp(-1) << endl;
	}

	// Write C1 and C2 to file
	FILE * ofile = fopen("rotacf.dat","w");
	fprintf ( ofile , "# Rotational auto-correlation function of water dipole vector [ time , C_1 , C_2 ]\n" );

	for (int frame = 0 ; frame < num_frames ; frame++ ) {
		fprintf ( ofile , "%f %f %f\n" , frame_dt*(double)frame , rotacf[frame] , rotacf2[frame] );
	}

	fclose(ofile);
}

/**
 * Find kinks as displacements that stick
 */
void Traj::select_displacement_kinks ( double min_kink_displacement , int t_transition, int t_sojourn, int stride ) {

	int t_transition_half = t_transition/2;
	int t_a = t_transition + 2*t_sojourn;

	if(verbose>4){
		cout << endl << " ..:: Find kinks. A atom is 'kinking' when is displace  ::.." << endl;
		cout << "stride = " << stride << endl;
		cout << "min_kink_displacement = " << min_kink_displacement << endl;
		cout << "t_transition = " << t_transition << endl;
		cout << "t_transition*frame_dt = " << (double)t_transition*frame_dt << endl;
		cout << "t_sojourn = " << t_sojourn << endl;
		cout << "t_sojourn*frame_dt = " << (double)t_sojourn*frame_dt << endl;
		cout << "t_transition_half = " << t_transition_half << endl;
		cout << "t_a = t_transition + 2*t_sojourn " << t_a << endl;
		cout << "t_a*frame_dt = " << (double)t_a*frame_dt << endl;
	}

	if(t_transition%1==1){
		t_transition = 2*t_transition_half;
		if(verbose>2){
			cout << "Warning: input t_transition is odd, and is reset to 2*t_transition_half = " << 2*t_transition_half << endl;
			cout << "t_transition = " << t_transition << endl;
		}
	}


	// Deselect all particles
	deselect_all();

	// Write kinks to file
	if(verbose>4) cout << "Write kinks.dat with kink data." << endl;
	FILE * ofile = fopen("kinks.dat","w");
	fprintf ( ofile , "# Kinks[ frame , particle , average sojourn distance , shortest transition time , transition frame ]\n" );

	// Find kinks
	int num_kinks = 0;
	int frame_effective = 0; // (num_frames - t_transition - 2 * t_sojourn)/stride
	double disp_avg_all = 0.0;
	int average_shortest_transition_time = 0 ;
	for ( int atom=0 ; atom < num_atoms ; atom++) {
		for ( int frame=t_sojourn+t_transition_half ; frame < num_frames - t_sojourn - t_transition_half ; frame+=stride) {
			if(atom==0) frame_effective++;

			double is_kink = true;

			// Is this particle kinking?
			double disp_avg = 0.0;
			for (int delta_frame = t_transition_half ; delta_frame < t_transition_half+t_sojourn ; delta_frame++){

				double disp = get_dr_img(atom,atom,frame-delta_frame,frame+delta_frame);
				disp_avg += disp;
				//cout << frame + (double)atom/1000 << " " <<  disp  << " " << df << endl;

				if(disp<min_kink_displacement)
					is_kink = false;
			}

			if (is_kink){
				// Find shortest time where this kink could have been detected
				int smallest_dt=t_transition;
				int frame_transition = frame - t_transition_half ;
				for ( int frame0=frame-t_transition_half ; frame0<frame+t_transition_half-1 ; frame0++ ) {
					for ( int frame1=frame0+1 ; frame1<frame+t_transition_half ; frame1++ ) {
						if( frame1-frame0<smallest_dt && (get_dr_img(atom,atom,frame0,frame1)>min_kink_displacement) ){
							// Does is stik in t_sojourn
							bool it_stiks = true;
							for (int df = 0 ; df < t_sojourn ; df++ ){	// TODO This loop ilesss expensive
								if ( get_dr_img(atom,atom,frame0-df,frame1+df) < min_kink_displacement ) it_stiks=false;
							}
							if (it_stiks) {
								smallest_dt = frame1-frame0;
								frame_transition=frame0;
							}
						}
					}
				}
				average_shortest_transition_time += smallest_dt;

				disp_avg /= (double)t_sojourn;
				disp_avg_all += disp_avg;

				fprintf ( ofile , "%d %d %f %d %d \n", frame , atom , disp_avg , smallest_dt , frame_transition );
				selected[ (frame + t_transition/2)*num_atoms + atom ]=true;
				num_kinks++;
			}
		}
	}

	if(verbose>4) cout << "num_kinks = " << num_kinks << endl;
	//if(verbose>4) cout << "frame_effective = num_frames - t_transition - 2 * t_sojourn = " << frame_effective << endl;
	if(verbose>4) cout << "frame_effective = " << frame_effective << endl;
	if(verbose>4) cout << "Kink concentration: " << (double)num_kinks/(double)frame_effective << endl;
	if(verbose>4) cout << "Kink concentration per particle: " << (double)num_kinks/(double)frame_effective/(double)num_atoms << endl;
	if(verbose>4) cout << "Average sojourn displacement: "  << disp_avg_all/(double)num_kinks << endl;
	if(verbose>4) cout << "Average shortest transition time: " << (double)average_shortest_transition_time/(double)num_kinks << " ( " << frame_dt*(double)average_shortest_transition_time/(double)num_kinks << " time units )" << endl;
	if(verbose>2 && frame_effective<1) cout << "Warning: frame_effective<1" << endl;

	fclose(ofile);
}




#endif /* TRAJ_H_ */

