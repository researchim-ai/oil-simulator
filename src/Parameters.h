#pragma once
#include <string>

class Parameters {
public:
    Parameters(const std::string& filename);

    int width;
    int height;
    double dt;
    double total_time;
    double alpha;
    double initial_pressure;
}; 