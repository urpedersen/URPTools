//============================================================================
//
//   Program for analyzing a system of acute triangular molecules like the LWOTP model:
//
//       A0     (non-central atom)
//      /
//     A1       (central atom)
//      \
//      A2     (non-central atom)
//
//
//   Version 0.1, October 20, 2010, by Ulf R. Pedersen (www.urp.dk).
//
//============================================================================

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

#include "trj.h"

using namespace std;

#ifndef PI
#define PI 3.14159265358979323846	// ratio of a circles circumference to its diameter
#endif

#ifndef NUM_DIM
#define NUM_DIM 3					// number of spatial dimensions
#endif

int verbose=9;						// Determines amount of print outs. 0 for silent and 9 for verbose.

int main(int narg, char **arg) {

	// Put variables from command line to variables string.
	string variables=set_variable_string(narg,arg);;


	verbose=get_variable(variables,"verbose",5);


	// Send greeting to dear user.
	if(verbose>4){
		cout << endl << "This program will analyze a trajectory of (acute) triangular molecules." << endl
			<< "Usage: " << arg[0] << " --load_traj_type=0" << endl
			<< "See documentation and output for variables that can and should be set with --variable=value flags." << endl
			<< "Command line: " << variables << endl
		;
	}



	// Load trajectory.
	if(verbose>4) cout << endl << " ..:: Load trajectory ::.." << endl;
	int num_atoms=0;											// Number of atoms.
	int num_frames=0;											// Number of frames.
	get_metadata_traj(variables,&num_frames,&num_atoms);
	double *bbox = new double[num_frames*NUM_DIM];				// Size of boundary box.
	int *type = new int[num_atoms];								// Atom type (0, 1 or 2).
	int *img_coords = new int[num_frames*num_atoms*NUM_DIM];    // Image coordinates.
	double *coords = new double[num_frames*num_atoms*NUM_DIM];  // Coordinates in primary image. TODO check that this is what is done, also when Gromacs files are loaded
	num_frames=load_traj(variables,num_frames,num_atoms,bbox,type,img_coords,coords);







	// Initialize and print variables.
	int atoms_per_molecule=3;
	if(verbose>4) cout << "atoms_per_molecules=" << atoms_per_molecule << endl;

	int num_molecules=num_atoms/atoms_per_molecule;
	if(verbose>4) cout << "num_molecules=" << num_molecules << endl;
	if(verbose>2 && num_atoms%atoms_per_molecule!=0)  cout << "Warning: Wrong number at atoms: num_atoms%atoms_per_molecules=" << num_atoms%atoms_per_molecule <<endl ;

	for(int a=0;a<num_atoms;a++) type[a]=a%3; // Set atom types. (TODO remove this variable. We assume the order of the atoms, and this type is not needed.)

	double frame_dt = get_variable(variables,"frame_dt",1.0);	// Time interval between frames
	if(verbose>4) cout << "frame_dt=" << frame_dt << endl;

	int num_log2_df = 0; // Number of possible frame differences on log2 scale (1,2,4,8 ...)
	for (int df = 1  ; df < num_frames ; df*=2) num_log2_df++;
	if(verbose>4) cout << "num_log2_df=" << num_log2_df << endl;








	// Set arrays with coordinates of A0, A1 and A2 atoms, and geometric center.
	double *coords_A0 = new double[num_frames*num_molecules*NUM_DIM];	// Atom 0
	int *img_coords_A0 = new int[num_frames*num_molecules*NUM_DIM];

	double *coords_A1 = new double[num_frames*num_molecules*NUM_DIM];	// Atom 1 (center atom)
	int *img_coords_A1 = new int[num_frames*num_molecules*NUM_DIM];

	double *coords_A2 = new double[num_frames*num_molecules*NUM_DIM];	// Atom 2
	int *img_coords_A2 = new int[num_frames*num_molecules*NUM_DIM];

	double *coords_GC = new double[num_frames*num_molecules*NUM_DIM];	// Geometric center
	int *img_coords_GC = new int[num_frames*num_molecules*NUM_DIM];

	for(int frame=0;frame<num_frames;frame++){
		for(int molecule=0;molecule<num_molecules;molecule++){
			for(int dim=0;dim<NUM_DIM;dim++){
				coords_A0[frame*num_molecules*NUM_DIM+molecule*NUM_DIM+dim]=coords[frame*num_atoms*NUM_DIM+molecule*atoms_per_molecule*NUM_DIM+0*NUM_DIM+dim];
				img_coords_A0[frame*num_molecules*NUM_DIM+molecule*NUM_DIM+dim]=img_coords[frame*num_atoms*NUM_DIM+molecule*atoms_per_molecule*NUM_DIM+0*NUM_DIM+dim];

				coords_A1[frame*num_molecules*NUM_DIM+molecule*NUM_DIM+dim]=coords[frame*num_atoms*NUM_DIM+molecule*atoms_per_molecule*NUM_DIM+1*NUM_DIM+dim];
				img_coords_A1[frame*num_molecules*NUM_DIM+molecule*NUM_DIM+dim]=img_coords[frame*num_atoms*NUM_DIM+molecule*atoms_per_molecule*NUM_DIM+1*NUM_DIM+dim];

				coords_A2[frame*num_molecules*NUM_DIM+molecule*NUM_DIM+dim]=coords[frame*num_atoms*NUM_DIM+molecule*atoms_per_molecule*NUM_DIM+2*NUM_DIM+dim];
				img_coords_A2[frame*num_molecules*NUM_DIM+molecule*NUM_DIM+dim]=img_coords[frame*num_atoms*NUM_DIM+molecule*atoms_per_molecule*NUM_DIM+2*NUM_DIM+dim];

				coords_GC[frame*num_molecules*NUM_DIM+molecule*NUM_DIM+dim]=get_geometric_center_img(molecule,dim,frame,num_frames,num_molecules,atoms_per_molecule,bbox,coords);
				img_coords_GC[frame*num_molecules*NUM_DIM+molecule*NUM_DIM+dim]=img_coords[frame*num_atoms*NUM_DIM+molecule*atoms_per_molecule*NUM_DIM+0*NUM_DIM+dim];
			}
		}
	}






	// Calculate Radial distribution function.

	if(verbose>4) cout << endl << " ..:: Radial distribution function (rdf.dat) ::.." << endl;
	int num_bins=get_variable(variables,"num_bins",5000);
	double r_max=get_variable(variables,"r_max",5);
	int frame_stride=get_variable(variables,"frame_stride",10);

	{
		double *rdf = new double[num_bins];
		get_radial_distribution_function(num_frames,num_atoms,bbox,coords,rdf,r_max,num_bins,frame_stride);

		double *rdf_A0 = new double[num_bins];
		get_radial_distribution_function(num_frames,num_molecules,bbox,coords_A0,rdf_A0,r_max,num_bins,frame_stride);

		double *rdf_A1 = new double[num_bins];
		get_radial_distribution_function(num_frames,num_molecules,bbox,coords_A1,rdf_A1,r_max,num_bins,frame_stride);

		double *rdf_GC = new double[num_bins];
		get_radial_distribution_function(num_frames,num_molecules,bbox,coords_GC,rdf_GC,r_max,num_bins,frame_stride);

		FILE * ofile_rdf = fopen("rdf.dat","w");
		fprintf (ofile_rdf, "# Radial distribution function: [pair distance] [all atoms] [A0 (bottom atom)] [A1 (center atom)] [geometric center] \n");
		for(int bin=0;bin<num_bins-1;bin++){
			fprintf (ofile_rdf, "%f %f %f %f %f\n",(double)bin*r_max/(double)num_bins,rdf[bin],rdf_A0[bin],rdf_A1[bin],rdf_GC[bin]);
		}
		fclose(ofile_rdf);
	}





	// Calculate Van Hove self-correlation function (Gs_LJ.dat)
	{
		if(verbose>4) cout << endl << " ..:: Self part of van Hove correlation function (Gs_LJ.dat and Gs_GC.dat) ::.." << endl;
		int frame_stride_Gs=get_variable(variables,"frame_stride_Gs",1);

		double *Gs = new double[num_bins*num_log2_df];

		// Lennard-Jones particles
		get_van_Hove_self_log2(num_frames,num_atoms,bbox,coords,Gs,r_max,num_bins,frame_stride_Gs);


		FILE * ofile_Gs = fopen("Gs_LJ.dat","w");

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

		// Geometric centers particles
		get_van_Hove_self_log2(num_frames,num_molecules,bbox,coords_GC,Gs,r_max,num_bins,frame_stride_Gs);

		ofile_Gs = fopen("Gs_GC.dat","w");

		fprintf (ofile_Gs, "# Self part of van Hove correlation function of molecular geometric center (4*pi*r^2*Gs): [Jump distance]");
		for ( int df = 1 ; df < num_frames ; df*=2 ) fprintf(ofile_Gs," [t=%f]",(double)df*frame_dt);
		fprintf (ofile_Gs, "\n");

		for(int bin=0;bin<num_bins-1;bin++){
			fprintf (ofile_Gs, "%f" , (double)bin*r_max/(double)num_bins );
			for (int df_index = 0  ; df_index < num_log2_df ; df_index++){
				fprintf ( ofile_Gs, " %f" , Gs [df_index*num_bins+bin] ) ;
			}
			fprintf (ofile_Gs, "\n");
		}

		fclose(ofile_Gs);
}







	// Calculate mean square displacements.
	if(verbose>4) cout << endl << " ..:: Mean square displacements (msd.dat) ::.." << endl;
	int trestart=get_variable(variables,"trestart",1);		// Interval between 'time zeros
	{
		double *msd=new double[num_frames];
		get_mean_square_displacement(num_frames,num_atoms,bbox,img_coords,coords,msd,trestart);

		double *msd_A0=new double[num_frames];
		get_mean_square_displacement(num_frames,num_molecules,bbox,img_coords_A0,coords_A0,msd_A0,trestart);

		double *msd_A1=new double[num_frames];
		get_mean_square_displacement(num_frames,num_molecules,bbox,img_coords_A1,coords_A1,msd_A1,trestart);

		double *msd_GC=new double[num_frames];
		get_mean_square_displacement(num_frames,num_molecules,bbox,img_coords_GC,coords_GC,msd_GC,trestart);

		FILE * ofile_msd = fopen("msd.dat","w");

		fprintf (ofile_msd, "# Mean square displacement: [Time] [all atoms] [A0 (base atom)] [A1 (top particle)] [geometric center] \n");
		for(int f=0;f<num_frames-1;f++){
			fprintf (ofile_msd, "%f %f %f %f %f\n",(double)f*frame_dt,msd[f],msd_A0[f],msd_A1[f],msd_GC[f]);
		}
		fclose(ofile_msd);
	}




	// ..:: Calculate self-intermediate (incoherent) scattering function ::..
	{
		if(verbose>4) cout << endl << " ..:: Self-intermediate scattering function (Fs.dat) ::.." << endl;

		double q_vector=get_variable(variables,"q_vector",2*PI);	// Length of q vector, |q|.

		double *Fs=new double[num_frames];
		get_self_intermediate_scattering_function(num_frames,num_atoms,bbox,img_coords,coords,Fs,q_vector,trestart);

		double *Fs_A0=new double[num_frames];
		get_self_intermediate_scattering_function(num_frames,num_molecules,bbox,img_coords_A0,coords_A0,Fs_A0,q_vector,trestart);

		double q_vector_GC=get_variable(variables,"q_vector_GC",q_vector);

		double *Fs_GC=new double[num_frames];
		get_self_intermediate_scattering_function(num_frames,num_molecules,bbox,img_coords_GC,coords_GC,Fs_GC,q_vector_GC,trestart);

		FILE * ofile_Fs = fopen("Fs.dat","w");
		fprintf (ofile_Fs, "# Self-intermediate scattering function: [Time] [All atoms] [A0 (base atom)] [Geometric center] \n");
		for(int f=0;f<num_frames-1;f++){
			fprintf (ofile_Fs, "%f %f %f %f\n",(double)f*frame_dt,Fs[f],Fs_A0[f],Fs_GC[f]);
		}
		fclose(ofile_Fs);
	}





	// ..:: Calculate rotational auto-correlation ::..
	if(verbose>4) cout << endl << " ..:: Rotational auto-correlation functions (rotacf.dat) ::.." << endl;


	// Base vectors
	double *base_vectors = new double[num_frames*num_molecules*NUM_DIM]; // Base vectors: vector pointing from atom A0 to A2.
	get_vectors ( num_frames*num_molecules , coords_A0 , coords_A2 , base_vectors );	// Note, we assume that molecules are connected.

	double avg_length_of_base_vector;	// Calculate and print the average length.
	for(int i=0;i<num_molecules*num_frames;i++) avg_length_of_base_vector += sqrt( base_vectors[i*NUM_DIM+0]*base_vectors[i*NUM_DIM+0]+base_vectors[i*NUM_DIM+1]*base_vectors[i*NUM_DIM+1]+base_vectors[i*NUM_DIM+2]*base_vectors[i*NUM_DIM+2] );
	avg_length_of_base_vector/=(double)(num_molecules*num_frames);
	if(verbose>4) cout << "Average length of base vector: avg_length_of_base_vector=" << avg_length_of_base_vector << endl;

	double *rotacf1_base_vector=new double[num_frames];
	double *rotacf2_base_vector=new double[num_frames];
	get_rotational_auto_correlation(num_frames,num_molecules,base_vectors,rotacf1_base_vector,rotacf2_base_vector,trestart);


	// Height vectors
	double *height_vectors = new double[num_frames*num_molecules*NUM_DIM]; // Height vectors: vector pointing from atom geometric center to A1.
	get_vectors ( num_frames*num_molecules , coords_GC , coords_A1 , height_vectors );	// Note, we assume that molecules are connected.

	double avg_length_of_height_vector;
	for(int i=0;i<num_molecules*num_frames;i++) avg_length_of_height_vector += sqrt( height_vectors[i*NUM_DIM+0]*height_vectors[i*NUM_DIM+0]+height_vectors[i*NUM_DIM+1]*height_vectors[i*NUM_DIM+1]+height_vectors[i*NUM_DIM+2]*height_vectors[i*NUM_DIM+2] );
	avg_length_of_height_vector/=(double)(num_molecules*num_frames);
	if(verbose>4) cout << "Average length of height vector: avg_length_of_height_vector=" << avg_length_of_height_vector << endl;

	double *rotacf1_height_vector=new double[num_frames];
	double *rotacf2_height_vector=new double[num_frames];
	get_rotational_auto_correlation(num_frames,num_molecules,height_vectors,rotacf1_height_vector,rotacf2_height_vector,trestart);


	// Plane vector
	double *plane_vectors = new double[num_frames*num_molecules*NUM_DIM]; // Plane vectors: cross product of base vector and height vector.
	get_cross_products ( num_frames*num_molecules , base_vectors , height_vectors , plane_vectors );

	double avg_length_of_plane_vector;
	for(int i=0;i<num_molecules*num_frames;i++) avg_length_of_plane_vector += sqrt( plane_vectors[i*NUM_DIM+0]*plane_vectors[i*NUM_DIM+0]+plane_vectors[i*NUM_DIM+1]*plane_vectors[i*NUM_DIM+1]+plane_vectors[i*NUM_DIM+2]*plane_vectors[i*NUM_DIM+2] );
	avg_length_of_plane_vector/=(double)(num_molecules*num_frames);
	if(verbose>4) cout << "Average length of plane vector: avg_length_of_plane_vector=" << avg_length_of_plane_vector << endl;

	double *rotacf1_plane_vector=new double[num_frames];
	double *rotacf2_plane_vector=new double[num_frames];
	get_rotational_auto_correlation(num_frames,num_molecules,plane_vectors,rotacf1_plane_vector,rotacf2_plane_vector,trestart);

	{
		// Print (numerical) Normalization factors.
		if(verbose>4) cout << "Normalization factors:"<< endl;
		if(verbose>4) cout << "rotacf1_base_vector[0]=" << rotacf1_base_vector[0] << endl;
		if(verbose>4) cout << "rotacf2_base_vector[0]=" << rotacf2_base_vector[0] << endl;
		if(verbose>4) cout << "rotacf1_height_vector[0]=" << rotacf1_height_vector[0] << endl;
		if(verbose>4) cout << "rotacf2_height_vector[0]=" << rotacf2_height_vector[0] << endl;
		if(verbose>4) cout << "rotacf1_plane_vector[0]=" << rotacf1_plane_vector[0] << endl;
		if(verbose>4) cout << "rotacf2_plane_vector[0]=" << rotacf2_plane_vector[0] << endl;

		// Print to file
		FILE * ofile_rotacf = fopen("rotacf.dat","w");
		fprintf (ofile_rotacf, "# Rotational auto-correlation, C_l(t) = < L_l(0)*L_l(t) > where L_1=cos(theta) and L_2=1/2(3*cos^2(theta)-1/2):\n # [Time] [C_1 of base vec.] [C_2 of base vec.] [C_1 of height vec.] [C_2 of height vec.] [C_1 of plane vec.] [C_2 of plane vec.] \n");
		for(int f=0;f<num_frames-1;f++){
			fprintf (ofile_rotacf, "%f %f %f %f %f %f %f\n",
					(double)f*frame_dt,
					rotacf1_base_vector[f]/rotacf1_base_vector[0],
					rotacf2_base_vector[f]/rotacf2_base_vector[0],
					rotacf1_height_vector[f]/rotacf1_height_vector[0],
					rotacf2_height_vector[f]/rotacf2_height_vector[0],
					rotacf1_plane_vector[f]/rotacf1_plane_vector[0],
					rotacf2_plane_vector[f]/rotacf2_plane_vector[0]);
		}
		fclose(ofile_rotacf);
	}






	// ..:: Calculate collective rotational auto-correlation ::..
	if(verbose>4) cout << endl << " ..:: Collective rotational auto-correlation functions (rotacf_col.dat) ::.." << endl;
	get_collective_rotational_auto_correlation(num_frames,num_molecules,plane_vectors,rotacf1_plane_vector,rotacf2_plane_vector,trestart);

	{
		// Print (numerical) Normalization factors.
		if(verbose>4) cout << "Normalization factors:"<< endl;
		if(verbose>4) cout << "rotacf1_base_vector[0]=" << rotacf1_base_vector[0] << endl;
		if(verbose>4) cout << "rotacf2_base_vector[0]=" << rotacf2_base_vector[0] << endl;
		if(verbose>4) cout << "rotacf1_height_vector[0]=" << rotacf1_height_vector[0] << endl;
		if(verbose>4) cout << "rotacf2_height_vector[0]=" << rotacf2_height_vector[0] << endl;
		if(verbose>4) cout << "rotacf1_plane_vector[0]=" << rotacf1_plane_vector[0] << endl;
		if(verbose>4) cout << "rotacf2_plane_vector[0]=" << rotacf2_plane_vector[0] << endl;

		// Print to file
		FILE * ofile_rotacf = fopen("rotacf_col.dat","w");
		fprintf (ofile_rotacf, "# Collective rotational auto-correlation, C_l(t) = < L_l(0)*L_l(t) > where L_1=cos(theta) and L_2=1/2(3*cos^2(theta)-1/2):\n # [Time] [C_1 of base vec.] [C_2 of base vec.] [C_1 of height vec.] [C_2 of height vec.] [C_1 of plane vec.] [C_2 of plane vec.] \n");
		for(int f=0;f<num_frames-1;f++){
			fprintf (ofile_rotacf, "%f %f %f %f %f %f %f\n",
					(double)f*frame_dt,
					rotacf1_base_vector[f]/rotacf1_base_vector[0],
					rotacf2_base_vector[f]/rotacf2_base_vector[0],
					rotacf1_height_vector[f]/rotacf1_height_vector[0],
					rotacf2_height_vector[f]/rotacf2_height_vector[0],
					rotacf1_plane_vector[f]/rotacf1_plane_vector[0],
					rotacf2_plane_vector[f]/rotacf2_plane_vector[0]);
		}
		fclose(ofile_rotacf);
	}









	// Qubatic order-parameter
	{
		if(verbose>4) cout << endl << " ..:: Qubatic orderparameter (qubatic.dat) ::.." << endl;

		// Get single molecule order-parameters
		double *Qi = new double[num_frames*num_molecules];
		get_qubatic_orderparameter(num_frames,num_molecules,base_vectors,Qi);

		// Count fraction of negative Q_i's
		{
			int counter = 0;
			for (int i = 0 ; i < num_frames*num_molecules ; i++){if(Qi[i]<0.0) counter++;}
			if(verbose>4) cout << "Fraction of negative Q_i's: " << (double)counter/(double)(num_frames*num_molecules) << endl;

		}

		// Write cubatic orderparameters to file, Q_i's
		FILE * ofile_qubatic = fopen("qubatic.dat","w");
		fprintf (ofile_qubatic, "# Qubatic orderparameters: [time] [Qi] [molecule]\n");
		for ( int f = 0 ; f < num_frames ; f++ ) {
			for ( int m = 0 ; m < num_molecules ; m++ ) {
				fprintf(ofile_qubatic,"%f %f %d\n",(double)f*frame_dt,Qi[f*num_molecules+m],m);
			}
		}
		fclose(ofile_qubatic);

		// Write Q = <Q_i>
		ofile_qubatic = fopen("qubatic_average.dat","w");
		fprintf (ofile_qubatic, "# Qubatic orderparameters: [time] [Q] [fraction of negative Q_i's] [number of defects]\n");
		for ( int f = 0 ; f < num_frames ; f++ ){
			int negative_Qi = 0;
			double average = 0.0;
			for ( int m = 0 ; m < num_molecules ; m++ ){
				average += Qi[f*num_molecules+m]/(double)num_molecules;
				if(Qi[f*num_molecules+m] < 0.0) negative_Qi++;
			}
			fprintf(ofile_qubatic,"%f %f %f %d\n",(double)f*frame_dt,average,negative_Qi/(double)num_molecules,negative_Qi);
		}
		fclose(ofile_qubatic);
	}

		// ..:: Jump angle statistics ::..
	{
		if(verbose>4) cout << endl << " ..:: Jump angle statistics (jump_angles.dat) ::.." << endl;

		int num_bins_angles=get_variable(variables,"num_bins_angles",180);

		int *histogram_base = new int[num_log2_df*num_bins_angles];
		int *histogram_height = new int[num_log2_df*num_bins_angles];

		int *df_counter = new int[num_log2_df];
		for ( int i = 0 ; i < num_log2_df ; i++ ) df_counter[i] = 0;

		double *vecA, *vecB;	// Pointers to the first element in vectors
		double *cos_thetas = new double[num_molecules];

		int df_index = 0;

		for ( int df = 1 ; df < num_frames ; df*=2 ){
			for ( int f0 = 0 ; f0 < num_frames-df ; f0++ ){

				// Base vector
				vecA = &base_vectors[f0*num_molecules*NUM_DIM];
				vecB = &base_vectors[(f0+df)*num_molecules*NUM_DIM];
				get_cos_thetas(num_molecules,vecA,vecB,cos_thetas);
				for (int i = 0 ; i < num_molecules ; i++ ) {
					histogram_base[df_index*num_bins_angles + (int)( (double)num_bins_angles*0.5*(cos_thetas[i]+1.0) ) ]++;
				}

				// Height vector
				vecA = &height_vectors[f0*num_molecules*NUM_DIM];
				vecB = &height_vectors[(f0+df)*num_molecules*NUM_DIM];
				get_cos_thetas(num_molecules,vecA,vecB,cos_thetas);
				for (int i = 0 ; i < num_molecules ; i++ ) {
					histogram_height[df_index*num_bins_angles + (int)( (double)num_bins_angles*0.5*(cos_thetas[i]+1.0) ) ]++;
				}

				// Count number of molecules for later normalization
				df_counter[df_index]+=num_molecules;
			}
			df_index++;
		}

		// Print for base vector
		FILE * ofile_jump_angles_base = fopen("jump_angles_base.dat","w");
		fprintf (ofile_jump_angles_base, "# Jump angles distribution for base vector: [cos(Delta theta)]");
		for ( int df = 1 ; df < num_frames ; df*=2 ) fprintf(ofile_jump_angles_base," [t=%f]",(double)df*frame_dt);
		fprintf (ofile_jump_angles_base, "\n");

		double bin_width = 2.0/(double)num_bins_angles;
		for(int bin=0;bin<num_bins_angles-1;bin++){
			fprintf (ofile_jump_angles_base, "%f",2*((double)bin+0.5)/((double)num_bins_angles)-1.0);
			for (int i = 0 ; i < num_log2_df ; i++){
				fprintf (ofile_jump_angles_base, " %f" , (double)histogram_base[i*num_bins_angles+bin]/(double)df_counter[i]/bin_width );
			}
			fprintf (ofile_jump_angles_base, "\n");
		}
		fclose(ofile_jump_angles_base);

		// Print for height vector
		FILE * ofile_jump_angles_height = fopen("jump_angles_height.dat","w");
		fprintf (ofile_jump_angles_height, "# Jump angles distribution for base vector: [cos(Delta theta)]");
		for ( int df = 1 ; df < num_frames ; df*=2 ) fprintf(ofile_jump_angles_height," [t=%f]",(double)df*frame_dt);
		fprintf (ofile_jump_angles_height, "\n");

		for(int bin=0;bin<num_bins_angles-1;bin++){
			fprintf (ofile_jump_angles_height, "%f",2*((double)bin+0.5)/((double)num_bins_angles)-1.0);
			for (int i = 0 ; i < num_log2_df ; i++){
				fprintf (ofile_jump_angles_height, " %f" , (double)histogram_height[i*num_bins_angles+bin]/(double)df_counter[i]/bin_width );
			}
			fprintf (ofile_jump_angles_height, "\n");
		}
		fclose(ofile_jump_angles_height);
	}

	// Say goodbye to the kind user.
	cout << endl << "Happy ending of " << arg[0] << endl;



}

