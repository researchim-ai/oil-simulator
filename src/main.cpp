#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

#include "Parameters.h"
#include "Grid.h"
#include "Solver.h"
#include "Scenario.h"
#include "Well.h"

int main(int argc, char* argv[]) {
    std::cout << "Starting Oil Simulator..." << std::endl;

    if (argc < 2) {
        std::cerr << "Error: Please provide a parameters file." << std::endl;
        std::cout << "Usage: " << argv[0] << " <path_to_params.dat>" << std::endl;
        return 1;
    }
    
    std::string params_file = argv[1];

    // Load parameters
    Parameters params(params_file);

    // Create grid and solver
    Grid grid(params);
    Solver solver(params, grid);

    // Setup wells using the Scenario class
    Scenario::create_five_spot(solver, params);

    std::cout << "\nInitial state set up." << std::endl;

    // Run simulation
    solver.solve();

    std::cout << "\nSimulation finished." << std::endl;
    
    // Save results
    solver.save_results("results.txt");
    
    // Visualize results
    std::cout << "Calling visualization script..." << std::endl;
    std::string command = "python3 ../python/visualize.py results.txt";
    int ret = system(command.c_str());
    if (ret != 0) {
        std::cerr << "Error: Visualization script failed." << std::endl;
    }

    std::cout << "Oil Simulator finished." << std::endl;
    return 0;
} 