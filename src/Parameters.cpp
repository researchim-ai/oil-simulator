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
        if (key == "WIDTH") ss >> width;
        else if (key == "HEIGHT") ss >> height;
        else if (key == "DT") ss >> dt;
        else if (key == "TIME") ss >> total_time;
        else if (key == "ALPHA") ss >> alpha;
        else if (key == "INITIAL_PRESSURE") ss >> initial_pressure;
    }
    std::cout << "Parameters loaded from " << filename << std::endl;
} 