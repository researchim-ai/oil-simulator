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
        if (key == "nx") ss >> nx;
        else if (key == "ny") ss >> ny;
        else if (key == "lx") ss >> lx;
        else if (key == "ly") ss >> ly;
        else if (key == "h") ss >> h;
        else if (key == "perm") ss >> perm;
        else if (key == "poro") ss >> poro;
        else if (key == "cr") ss >> cr;
        else if (key == "dt") ss >> dt;
        else if (key == "total_time") ss >> total_time;
        else if (key == "initial_pressure") ss >> initial_pressure;
        else if (key == "initial_saturation") ss >> initial_saturation;
        else if (key == "pc_entry") ss >> pc_entry;
    }
    std::cout << "Parameters loaded from " << filename << std::endl;
} 