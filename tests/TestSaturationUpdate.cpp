#include <iostream>
#include <cassert>
#include <cmath>
#include "Solver.h"
#include "Grid.h"
#include "Parameters.h"
#include "Well.h"
#include "FluidProperties.h"

// Custom assert for floating point comparison
void assert_almost_equal(double a, double b, double epsilon = 1e-9) {
    if (std::abs(a - b) > epsilon) {
        std::cerr << "ASSERT FAILED: " << a << " != " << b << ". Difference: " << std::abs(a-b) << std::endl;
        assert(false);
    }
}

void test_saturation_update_isolated() {
    std::cout << "--- Running Test: Isolated Saturation Update ---" << std::endl;

    // 1. Setup a 1x1 grid
    Parameters params;
    params.nx = 1;
    params.ny = 1;
    params.lx = 10;
    params.ly = 10;
    params.h = 10;
    params.poro = 0.2;
    params.perm = 1e-13; // 100 mD. THIS WAS THE MISSING PIECE.
    params.dt = 1.0; // 1 second timestep
    Grid grid(params);
    Solver solver(params, grid);

    // 2. Set known initial state
    double p_initial = 2.0e7;
    double s_initial = 0.1;
    grid.pressure(0, 0) = p_initial;
    grid.saturation(0, 0) = s_initial;

    // 3. Manually set the result of the pressure solve
    // This is the key part of the isolated test
    double p_final = 2.2e7; // Let's assume pressure increased
    solver.p_new.resize(1);
    solver.p_new(0) = p_final;

    // 4. Add a well
    Well well(0, 0, 0.1, 3.0e7); // 300 bar injector
    well.calculateWellIndex(grid);
    solver.add_well(well);

    // 5. Calculate expected result
    // This logic must EXACTLY match the solver's logic for a well in a single cell
    double q_w = 0;
    bool is_injector = well.getBHP() > p_final;
    if (is_injector) {
        double mob_water_inj = solver.fluid.krw(1.0) / solver.fluid.waterViscosity(p_final);
        q_w = well.getWellIndex() * mob_water_inj * (well.getBHP() - p_final);
    } else {
        double mob_w = solver.fluid.krw(s_initial) / solver.fluid.waterViscosity(p_final);
        double mob_o = solver.fluid.kro(s_initial) / solver.fluid.oilViscosity(p_final);
        double total_mobility = mob_w + mob_o;
        if (total_mobility > 1e-9) {
            double q_total = well.getWellIndex() * total_mobility * (well.getBHP() - p_final);
            q_w = (mob_w / total_mobility) * q_total;
        }
    }
    
    // Divergence is zero for a 1x1 grid
    double divergence_term = 0.0;
    double cell_volume = params.lx * params.ly * params.h;
    // The main equation is: s_new = s_old - (dt/Vp) * (Divergence - Source)
    double s_expected = s_initial - (params.dt / (params.poro * cell_volume)) * (divergence_term - q_w);
    s_expected = std::max(solver.fluid.s_wc, std::min(1.0 - solver.fluid.s_or, s_expected));
    
    // --- DEBUGGING THE MANUAL CALCULATION ---
    std::cout << "\n--- DEBUG INFO FOR MANUAL CALCULATION ---\n";
    std::cout << "  Well Index: " << well.getWellIndex() << std::endl;
    std::cout << "  BHP: " << well.getBHP() << ", P_final: " << p_final << std::endl;
    std::cout << "  Is Injector?: " << is_injector << std::endl;
    std::cout << "  krw(1.0): " << solver.fluid.krw(1.0) << std::endl;
    std::cout << "  waterViscosity(p_final): " << solver.fluid.waterViscosity(p_final) << std::endl;
    std::cout << "  mob_water_inj: " << (solver.fluid.krw(1.0) / solver.fluid.waterViscosity(p_final)) << std::endl;
    std::cout << "----------------------------------------\n" << std::endl;

    std::cout << "p_initial: " << p_initial << ", p_final: " << p_final << std::endl;
    std::cout << "q_w (manually calculated): " << q_w << std::endl;
    std::cout << "Expected s_new: " << s_expected << std::endl;

    // 6. Run the function to be tested
    solver.update_saturation(params.dt);

    // 7. Check the result
    double s_simulated = grid.saturation(0, 0);
    std::cout << "Simulated s_new: " << s_simulated << std::endl;
    assert_almost_equal(s_expected, s_simulated);

    std::cout << "--- Test Passed ---" << std::endl;
}


int main() {
    test_saturation_update_isolated();
    return 0;
} 