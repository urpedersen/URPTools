/**
 * helpers_io.cpp
 *
 *  Created on: Oct 6, 2011
 *      Author: Ulf R. Pedersen
 */

#include <iostream>
using namespace std;

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include "helpers_io.h"




/* return a ''command line string'' with variables. This is send to the get_variable functions
* ... let the user change values of variables
*/
string set_variable_string(int narg, char **arg){
	string variables="";

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
*      ... this could look some like " --var1=146 --var2=2e5 and a optional comment --var3=-200"
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
		if(verbose>4) cout << "Read variable: "<< variable_name << default_value << "\n" ;
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
		if(verbose>4) cout << "Warning: Missing variable. Use default"<< variable_name << default_value << " \n" ;
	}
	else{
		string tmp_str=variable_string.substr(pos + variable_name.size());
		pos = tmp_str.find(" ");
		tmp_str=tmp_str.substr(0,pos);
		default_value=atof(tmp_str.c_str());
		if(verbose>4) cout << "Read variable: "<< variable_name << default_value << " \n";
	}
	return default_value;
}






/**
*   Return variable value given in variable_string.
*      ... this could look some like "var1=146 var2=2e5 and a optional comment var3=-200"
*      If the variable is set more than ones, the first one is used.
*
* @param variable_string
* @param variable_name
* @param default_value
*/
int get_header_variable(string variable_string,string variable_name,int default_value){
	//int verbose=DEFAULT_VERBOSE;
	int verbose=5;

	size_t pos;
	variable_name.insert(0," ");   // make sure that variable name is surrounded by " " and "="
	variable_name.append("=");
	pos = variable_string.find(variable_name,0);

	if(pos==-1){ // No match. NB: comparing signed and unsigned integer
		if(verbose>4) cout << "Warning: Missing header variable. Use default"<< variable_name << default_value << "\n" ;
	}
	else{
		string tmp_str=variable_string.substr(pos + variable_name.size());
		pos = tmp_str.find(" ");
		tmp_str=tmp_str.substr(0,pos);
		default_value=atol(tmp_str.c_str());
		if(verbose>4) cout << "Read header variable: "<< variable_name << default_value << "\n" ;
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
double get_header_variable(string variable_string,string variable_name,double default_value){
	//int verbose=DEFAULT_VERBOSE;
	int verbose=5;

	size_t pos;
	variable_name.insert(0," ");   // make sure that variable name is surrounded by white space
	variable_name.append("=");     //   and "="
	pos = variable_string.find(variable_name,0);

	if(pos==-1){  // No match. NB: comparing signed and unsigned integer
		if(verbose>4) cout << "Warning: Missing header variable. Use default"<< variable_name << default_value << " \n" ;
	}
	else{
		string tmp_str=variable_string.substr(pos + variable_name.size());
		pos = tmp_str.find(" ");
		tmp_str=tmp_str.substr(0,pos);
		default_value=atof(tmp_str.c_str());
		if(verbose>4) cout << "Read header variable: "<< variable_name << default_value << " \n";
	}
	return default_value;
}



