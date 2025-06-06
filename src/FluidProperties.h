#pragma once

// Простая модель относительных фазовых проницаемостей (ОФП)
// k_rw = S_w
// k_ro = 1 - S_w
class FluidProperties {
public:
    FluidProperties();

    // Относительная проницаемость для воды (relative permeability of water)
    double krw(double saturation) const;

    // Относительная проницаемость для нефти (relative permeability of oil)
    double kro(double saturation) const;
}; 