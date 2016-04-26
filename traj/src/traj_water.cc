//============================================================================
// Name        : traj_water.cpp
// Author      : Ulf R. Pedrsen
//============================================================================

#include <iostream>

#include "traj.h"
#include "../../cluster_analysis/src/cluster_analysis.h"

using namespace std;

int main(int narg, char **arg) {

	// Initialize.
	string tmp=set_variable_string(narg,arg);
	Traj traj(tmp);

	// Welcome user.
	if(traj.verbose>4){
		cout << "Post analyzes of water trajectory." << endl
			<< "Usage: " << arg[0] << " --num_frames=1" << endl
			<< "See documentation and output for variables that can and should be set with --variable=value flags." << endl;
	}

	// Load trajectory.
	traj.load_data_from_disk();

	// Analysis done on full trajectory
	traj.radial_distribution("rdf_all.dat");
	traj.water_dipole_rotational_correlation();

	// Smooth trajectory.
	//traj.smooth_by_time_integral();

	// Find hydrogen bonds (after this, traj only have oxygens coordinates hereafter )
	traj.construct_hbond_network();
	traj.enumerate_nl_hbonds();

	// H-Bond analysis
	traj.nl_bond_time_autocorrelation("hbond_autocorr.dat");

	// Find hbond kinks
	traj.hbond_dynamic_exitations_ts();

	// Cluster analysis of kinks - TODO move below to traj, and make it cluana on selected particles
	if(traj.verbose>4) cout << "\n ..:: Cluster analysis of excitations (cluster_analysis_exitations.dat)" << endl;
	//double d_cut = 3.6;
	FILE * ofile = fopen("cluster_analysis_exitations.dat","w");
	fprintf (ofile, "# [frame, cluster size, msd ]\n");

	for(int f = 0 ; f < traj.num_frames ; f++){
	//for(int f = 5 ; f < 7 ; f++){
		Cluster_analysis clu_ana;
		clu_ana.create_nodes(traj.num_atoms);
		clu_ana.verbose = 5;

		// Add connections between exited oxygens
		for(int a = 0 ; a < traj.num_atoms ; a++){
			for ( int n = traj.nl_first.at(f*traj.num_atoms+a) ;
					n < traj.nl_last.at(f*traj.num_atoms+a) ;
					n++ ) {
				if(	traj.nl_selected.at(n) && a < traj.nl.at(n)){
					clu_ana.create_node_connection(clu_ana.nodes.at(a),clu_ana.nodes.at(traj.nl.at(n)));
				}
			}
		}

		// Remove non-excited oxygens TODO move to cluster_analysis.h (as a 'remove single nodes' function)
		unsigned int counter = 0;
		for(int a = 0 ; a < traj.num_atoms ; a++){
			//if(traj.type[f*traj.num_atoms+a]==1){
			if(clu_ana.nodes.at(counter)->neighbor_nodes.size()>0){
				counter++;
			}else{
				clu_ana.nodes.erase(clu_ana.nodes.begin() + counter);
			}
		}

		// Assign clusters
		clu_ana.assign_nodes_to_clusters();

		// Print result of cluster analysis
		/*if(traj.verbose>4){
			cout << "\nFrame: " << f << endl;
			clu_ana.print();
		}*/
		for(unsigned int i=0;i<clu_ana.clusters.size();i++){
			int cluster_size = clu_ana.clusters.at(i)->nodes.size();

			// Calculate MSD of clusters (using minimum convention)
			double msd = 0.0;
			for ( int j = 0 ; j < cluster_size-1 ; j++ ){
				int p0 = clu_ana.clusters.at(i)->nodes.at(j)->index;
				for ( int k = j+1 ; k < cluster_size ; k++ ){
					int p1 = clu_ana.clusters.at(i)->nodes.at(k)->index;
					double dist = traj.get_dr_img(p0,p1,f);
					msd += dist*dist;
				}
			}
			msd /= (double)cluster_size * (double)(cluster_size-1);

			// Print to file
			fprintf (ofile, "%f %u %f\n",
					traj.get_time(f),
					cluster_size ,
					msd );
		}

		// Destroy clu_ana
	}
	fclose(ofile);
	// END of cluster analysis of hbond kinks

	// Analysis on oxygen trajectories
	traj.van_Hove_self("Gs.dat");
	traj.radial_distribution("rdf_oxygen.dat");
	traj.mean_squared_displacement("msd_oxygen.dat");
	traj.incoherent_self_scattering("Fs_q2.dat",2.0);
	traj.incoherent_self_scattering_of_selected("Fs_q2_brokenOH.dat",2.0);
	traj.incoherent_self_scattering("Fs_q3.dat",3.0);
	traj.incoherent_self_scattering_of_selected("Fs_q3_brokenOH.dat",3.0);

	// Find kinks using the Keys definition
	double min_kink_displacement = get_variable(traj.variables,"min_kink_displacement",1.0);
	int t_transition = get_variable(traj.variables,"t_transition",250);
	int t_sojourn = get_variable(traj.variables,"t_sojourn",250);
	int kinks_stride = get_variable(traj.variables,"kinks_stride",t_transition);
	traj.select_displacement_kinks(
			min_kink_displacement,
			t_transition,
			t_sojourn,
			kinks_stride );

	// Say goodbye to the kind user.
	if(traj.verbose>4) cout << endl << "Variable string: " << traj.variables << endl;
	if(traj.verbose>4) cout << "Happy ending. " << arg[0] <<endl;
}

