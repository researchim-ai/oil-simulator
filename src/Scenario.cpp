#include "Scenario.h"
#include "Solver.h"
#include "Parameters.h"
#include "Well.h"
#include "Grid.h" // Need grid to set initial conditions

void Scenario::create_five_spot(Solver& solver, const Parameters& params) {
    // --- Set Initial Conditions ---
    Grid& grid = solver.get_grid();
    const double initial_pressure = params.initial_pressure;
    const double initial_sw = params.initial_saturation;
    for (int i = 0; i < grid.getNx(); ++i) {
        for (int j = 0; j < grid.getNy(); ++j) {
            grid.pressure(i, j) = initial_pressure;
            grid.saturation(i, j) = initial_sw;
        }
    }

    // --- Scenario Setup: Five-spot pattern ---
    // 1. Injector in the center
    int inj_i = params.nx / 2;
    int inj_j = params.ny / 2;
    double inj_radius = 0.1; // meters
    double inj_bhp = 3.0e7;  // Pascals (300 bar), higher than initial pressure
    Well injector(inj_i, inj_j, inj_radius, inj_bhp);
    solver.add_well(injector);

    // 2. Producers at the corners
    double prod_radius = 0.1; // meters
    double prod_bhp = 1.0e7;  // Pascals (100 bar), lower than initial pressure
    int margin = 1; // place wells not exactly at 0,0 but slightly inside
    
    Well producer1(margin, margin, prod_radius, prod_bhp);
    solver.add_well(producer1);
    
    Well producer2(params.nx - 1 - margin, margin, prod_radius, prod_bhp);
    solver.add_well(producer2);

    Well producer3(margin, params.ny - 1 - margin, prod_radius, prod_bhp);
    solver.add_well(producer3);

    Well producer4(params.nx - 1 - margin, params.ny - 1 - margin, prod_radius, prod_bhp);
    solver.add_well(producer4);
} 