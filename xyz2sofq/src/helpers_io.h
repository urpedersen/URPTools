/*
 * helpers.h
 *
 *  Created on: Oct 6, 2011
 *      Author: Ulf R. Pedersen
 */

#ifndef HELPERS_IO_H_
#define HELPERS_IO_H_

#include <iostream>
using namespace std;

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>


string set_variable_string(int, char **);
int get_variable(string,string,int);
double get_variable(string,string,double);

int get_header_variable(string,string,int);
double get_header_variable(string,string,double);

#endif /* HELPERS_IO_H_ */
