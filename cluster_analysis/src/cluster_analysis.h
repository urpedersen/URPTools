/*
 * cluster_analysis.h
 *
 *  Created on: Apr 6, 2011
 *      Author: Ulf R. Pedersen.
 */

#ifndef CLUSTER_ANAYSIS_H_
#define CLUSTER_ANAYSIS_H_


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <stack>
#include <assert.h>

using namespace std;

class Cluster_analysis_node;
class Cluster_analysis_cluster;
class Cluster_analysis;

class Cluster_analysis_node {
public:

	int index;

	Cluster_analysis_cluster* cluster;
	vector<Cluster_analysis_node*> neighbor_nodes;

	Cluster_analysis_node();		// Constructor

	void connect_to(Cluster_analysis_node*);
	void print();
};

class Cluster_analysis_cluster {
public:

	int index;

	vector<Cluster_analysis_node*> nodes;

	Cluster_analysis_cluster();	// Constructor
	

	void print();
};

class Cluster_analysis {
public:
	int verbose;

	vector<Cluster_analysis_node*> nodes;
	vector<Cluster_analysis_cluster*> clusters;

	Cluster_analysis();	// Constructor
	~Cluster_analysis(); // Destructor

	void load_data_from_file(string);

	void create_nodes(unsigned int);
	void create_node_connection(Cluster_analysis_node *,Cluster_analysis_node *);

	void assign_nodes_to_clusters();
	void create_cluster();
	
	void get_nodes_in_cluster(unsigned cluster_index, vector<unsigned>& out);

	int get_mass_of_largest_cluster();
	int get_index_of_largest_cluster();

	void print();
};

/**
 * (Standard) Constructor
 */
Cluster_analysis_node::Cluster_analysis_node(){
	index = -1;
	cluster = NULL;
}

void Cluster_analysis_node::connect_to(Cluster_analysis_node* node){
	neighbor_nodes.push_back(node);
}

void Cluster_analysis_node::print(){
	cout << "Node " << index;
	if(cluster==NULL){
		cout << " is not assigned a to a cluster ";
	}else{
		cout << " is in cluster " << cluster -> index;
	}
	cout << " and have ";
	cout << neighbor_nodes.size() << " connection(s):";
	for (unsigned int i=0;i<neighbor_nodes.size();i++) cout << " " << neighbor_nodes.at(i) -> index;

	cout << endl;
}

/**
 * (Standard) Constructor
 */
Cluster_analysis_cluster::Cluster_analysis_cluster(){
	index = -1;
}

/**
 * Print information about cluster to cout
 */
void Cluster_analysis_cluster::print(){
	cout << "Cluster " << index << " contains " << nodes.size() << " nodes:";
	for (unsigned int i = 0 ; i < nodes.size() ; i++ ) cout << " " << nodes.at(i) -> index;
	cout << endl;
}

/**
 * (Standard) Constructor
 */
Cluster_analysis::Cluster_analysis(){
	verbose=5;
}

Cluster_analysis::~Cluster_analysis(){
	for (unsigned int i = 0 ; i < nodes.size() ; i++ ) delete nodes[i];
	for (unsigned int i = 0 ; i < clusters.size() ; i++ ) delete clusters[i];
}

/**
 * Load data from files and allocate (all) variables
 */
void Cluster_analysis::load_data_from_file(string filename){

	//string filename="node_connections.dat";
	ifstream file;
	file.open(filename.c_str());

	if (file.is_open()) {
		string line;
		size_t pos;

		// Load & set number of nodes
		getline(file,line);
		pos = line.find(" ");
		create_nodes( atoi(line.substr(0,pos).c_str()) );

		// Load & set connection between nodes
		unsigned int node0,node1;
		while(! file.eof() ){
			getline(file,line);

			pos = line.find(" ");
			node0 = atoi(line.substr(0,pos).c_str());
			line=line.substr(pos+1);

			pos = line.find(" ");
			node1=atoi(line.substr(0,pos).c_str());

			if(verbose>6) cout << "connection: node0=" << node0 << " node1=" <<  node1 << endl;

			create_node_connection(nodes[node0],nodes[node1]);
		}
	}else{
		cout << "Error: Unable to open " << filename << ". Exit." << endl;
		cout << "10  # Example with ten nodes, ( 0, 1, 2, ... , 9 )\n0 1 #\n1 2 #    8-7-6\n1 8 #    |   |\n3 5 #  0-1   3-5  4-9\n8 7 #    |\n7 6 #    2\n6 3\n4 9\n";
		exit(0);
	}

	file.close();
}

/**
 * Create (empty) note. After running this, nodes should be connected.
 */
void Cluster_analysis::create_nodes(unsigned int number_of_nodes_to_add){
	nodes.resize(nodes.size()+number_of_nodes_to_add);
	for(unsigned int i=nodes.size()-number_of_nodes_to_add ; i<nodes.size() ; i++){
		nodes.at(i) = new Cluster_analysis_node();
		nodes.at(i)->index=i;
	}
}

/**
 *
 */
void Cluster_analysis::create_node_connection(Cluster_analysis_node* node0,Cluster_analysis_node* node1){
		node0->connect_to(node1);
		node1->connect_to(node0);
}

/**
 *  Assign nodes to cluster. Loop stack, while putting connecting nodes into top of stack, so that a cluster is found.
 */
void Cluster_analysis::assign_nodes_to_clusters(){

	stack<Cluster_analysis_node*> node_stack;
	for(unsigned int i = 0 ; i < nodes.size() ; i ++ ) node_stack.push(nodes.at(i));

	while ( !(node_stack.empty()) ){

		Cluster_analysis_node *current_node = node_stack.top();
		node_stack.pop();

		// First node in cluster.
		if( current_node -> cluster == NULL ){
			create_cluster();

			current_node -> cluster = clusters.back();
			current_node -> cluster -> nodes.push_back(current_node);
		}

		// If connecting node is a new in the cluster,
		//   assign it to current cluster,
		//   and put it in top of stack.
		for(unsigned int i = 0 ; i < current_node->neighbor_nodes.size() ; i++ ){

			Cluster_analysis_node * neighbor_node = current_node -> neighbor_nodes.at(i);

			if(neighbor_node -> cluster != current_node -> cluster){
				assert( neighbor_node -> cluster == NULL );

				neighbor_node -> cluster = current_node -> cluster;
				neighbor_node -> cluster -> nodes.push_back(neighbor_node);

				node_stack.push( neighbor_node );
			}
		}
	} // END while ( !(node_stack.empty()) )
}

void Cluster_analysis::create_cluster(){
	clusters.push_back(new Cluster_analysis_cluster());
	clusters.back()->index = clusters.size()-1;
}

void Cluster_analysis::get_nodes_in_cluster(unsigned c,vector<unsigned>& out){
	out.clear();
	for( unsigned i = 0 ; i<clusters[c]->nodes.size() ; i++ )
		out.push_back(clusters[c]->nodes.at(i)->index);
}

/**
 * Return mass of largest cluster.
 */
int Cluster_analysis::get_mass_of_largest_cluster(){

	unsigned int out=0;

	for (unsigned int i = 0 ; i < clusters.size() ; i++ )
		if ( clusters[i]->nodes.size() > out )
			out = clusters[i]->nodes.size()
	;

	return out;
}

/**
 * Return the index of largest cluster.
 */
int Cluster_analysis::get_index_of_largest_cluster(){

	unsigned int out=0;
	unsigned int size=0;
	
	for (unsigned int i = 0 ; i < clusters.size() ; i++ )
		if ( clusters[i]->nodes.size() > size ){
			size = clusters[i]->nodes.size();
			out=i;
		}
	;

	return out;
}



void Cluster_analysis::print ( ) {
	if(verbose>4){

		cout << "Number of nodes: " << nodes.size() << endl;
		cout << "Number of clusters: " << clusters.size() << endl;
		cout << "Mass of largest cluster: " << get_mass_of_largest_cluster() << endl;
		cout << "Index of largest cluster: " << get_index_of_largest_cluster() << endl;

		if(verbose>5){
			for(unsigned int i=0; i<clusters.size() ; i++ )
				clusters.at(i)->print();
			for(unsigned int i=0; i<nodes.size() ; i++ )
				nodes.at(i)->print();
		}
	}
}

#endif /* CLUSTER_ANAYSIS_H_ */
