#pragma once

class Grid; // Forward declaration

class Well {
public:
    // A well is defined by its location (i, j), its radius, and its bottom-hole pressure (BHP)
    Well(int i, int j, double radius, double bhp);

    // Peaceman well index calculation. Must be called after grid is created.
    void calculateWellIndex(const Grid& grid);

    // Calculate the source/sink term for the flow equation
    double calculateSourceTerm(double block_pressure, double total_mobility) const;

    int getI() const { return m_i; }
    int getJ() const { return m_j; }
    double getBHP() const { return m_bhp; }
    double getWellIndex() const { return m_wellIndex; }

private:
    int m_i;
    int m_j;
    double m_radius;      // Wellbore radius, rw
    double m_bhp;         // Bottom-hole pressure, P_bh
    double m_wellIndex;   // Well Index, WI
}; 