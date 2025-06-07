#include "Parameters.h"
#include <fstream>
#include <iostream>
#include <sstream>

Parameters::Parameters(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open parameters file " << filename << std::endl;
        exit(1);
    }

    std::string line;
    std::string key;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        ss >> key;
        if (key == "NX") ss >> nx;
        else if (key == "NY") ss >> ny;
        else if (key == "LX") ss >> lx;
        else if (key == "LY") ss >> ly;
        else if (key == "H") ss >> h;
        else if (key == "PERM") ss >> perm;
        else if (key == "PORO") ss >> poro;
        else if (key == "CR") ss >> cr;
        else if (key == "DT") ss >> dt;
        else if (key == "TIME") ss >> total_time;
        else if (key == "INITIAL_PRESSURE") ss >> initial_pressure;
        else if (key == "INITIAL_SATURATION") ss >> initial_saturation;
    }
    std::cout << "Parameters loaded from " << filename << std::endl;
} 