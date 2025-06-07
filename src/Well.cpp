#include "Well.h"
#include "Grid.h"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Well::Well(int i, int j, double radius, double bhp)
    : m_i(i), m_j(j), m_radius(radius), m_bhp(bhp), m_wellIndex(0.0) {}

void Well::calculateWellIndex(const Grid& grid) {
    double dx = grid.getDx();
    double dy = grid.getDy();
    // Using cell-center permeability
    double k = grid.getPermeability(m_i, m_j); 
    double h = 1.0; // Assuming 2D reservoir, so thickness is 1.0 m.

    // Peaceman's equivalent radius for an isotropic permeability field in a grid block
    double re = 0.14 * std::sqrt(dx*dx + dy*dy);

    if (re <= m_radius) {
        // This is a configuration error (e.g., well radius is too large for the cell).
        // To prevent division by zero or log of a non-positive number, we set WI to zero.
        m_wellIndex = 0;
        return;
    }

    // Standard Peaceman well index formula
    m_wellIndex = (2.0 * M_PI * k * h) / std::log(re / m_radius);
}

double Well::calculateSourceTerm(double block_pressure, double total_mobility) const {
    // Source/sink term q = WI * lambda_total * (P_bh - P_block)
    // The result has units of m^3/s. The solver will handle division by cell volume.
    return m_wellIndex * total_mobility * (m_bhp - block_pressure);
} 