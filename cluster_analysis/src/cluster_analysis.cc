//============================================================================
// Name        : cluster_analysis.cpp
// Author      : Ulf R. Pedersen
// Version     :
// Copyright   : See COPY
// Description : Preform cluser analysis.
//============================================================================

#include <iostream>
#include <assert.h>
#include "cluster_analysis.h"

using namespace std;

int main() {

	Cluster_analysis clu_ana;
	clu_ana.verbose=6;

	clu_ana.load_data_from_file("node_connections.dat");

	clu_ana.assign_nodes_to_clusters();

	clu_ana.print();

	return 0;

}

