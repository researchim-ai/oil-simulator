#pragma once
#include "Grid.h"
#include "FluidProperties.h"

class Solver {
public:
    Solver(Grid& grid, double alpha);

    void step(double dt);

private:
    Grid& m_grid;
    double m_alpha;
    FluidProperties m_fluid;
}; 