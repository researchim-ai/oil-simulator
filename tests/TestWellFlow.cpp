#include <iostream>
#include <cassert>
#include <cmath>
#include "Well.h"
#include "Grid.h"
#include "Parameters.h"
#include "FluidProperties.h"

// Custom assert for floating point comparison
void assert_almost_equal(double a, double b, double epsilon = 1e-9) {
    if (std::abs(a - b) > epsilon) {
        std::cerr << "ASSERT FAILED: " << a << " != " << b << ". Difference: " << std::abs(a-b) << std::endl;
        assert(false);
    }
}

void test_well_source_term() {
    std::cout << "--- Running Test: Well Source Term Calculation ---" << std::endl;

    // 1. Setup
    Parameters params;
    params.nx = 1;
    params.ny = 1;
    params.lx = 10;
    params.ly = 10;
    params.h = 10;
    params.perm = 1e-13; // 100 mD
    Grid grid(params);
    FluidProperties fluid;

    double block_pressure = 2.0e7; // 200 bar
    double s_w = 0.1; // connate water
    double bhp = 3.0e7; // 300 bar (injector)

    // 2. Create well and calculate its properties
    Well well(0, 0, 0.1, bhp);
    well.calculateWellIndex(grid);
    assert(well.getWellIndex() > 0);

    // 3. Calculate mobility and expected source term
    double mu_o = fluid.oilViscosity(block_pressure);
    double kro = fluid.kro(s_w); // At s_wc, krw is 0.
    double total_mobility = kro / mu_o;
    assert(total_mobility > 0);
    
    // This is the term q = WI * lambda * (P_bh - P_block)
    double expected_q = well.getWellIndex() * total_mobility * (bhp - block_pressure);
    assert(expected_q > 0);

    // 4. Use the class method to calculate the source term
    double calculated_q = well.calculateSourceTerm(block_pressure, total_mobility);

    // 5. Assert they are equal
    std::cout << "Expected Source Term q: " << expected_q << std::endl;
    std::cout << "Calculated Source Term q: " << calculated_q << std::endl;
    assert_almost_equal(expected_q, calculated_q);

    std::cout << "--- Test Passed ---" << std::endl;
}

int main() {
    test_well_source_term();
    return 0;
} 