#include "split.h"

#include <string>
#include <sstream>
#include <vector>

using namespace std;

/**
 * Split a string into a vector string
 */
vector<string> split(string str, char delimiter) {
    vector<string> output;
    stringstream ss(str);
    string substr;
    while(getline(ss, substr, delimiter)) 
        output.push_back(substr);
    return output;
}
