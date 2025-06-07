#include "FluidProperties.h"
#include <cmath>
#include <algorithm>

FluidProperties::FluidProperties() {
    // Viscosity properties
    mu_ref_w = 1.0;    // cP
    mu_ref_o = 10.0;   // cP
    c_w = 1e-6;        // 1/psi for water viscosity
    c_o = 1e-5;        // 1/psi for oil viscosity
    p_ref_w = 0.0;     // Reference pressure, Pa
    p_ref_o = 0.0;     // Reference pressure, Pa

    // Corey model relative permeability properties
    s_wc = 0.1;        // Connate water saturation
    s_or = 0.1;        // Residual oil saturation
    krw_max = 0.8;     // Max water relative permeability
    kro_max = 1.0;     // Max oil relative permeability
    n_w = 2.0;         // Corey exponent for water
    n_o = 2.0;         // Corey exponent for oil

    // Capillary pressure
    pc_entry = 2000.0; // Pa (default value, ~0.3 psi)
}

// Модель ОФП Кори
// Swc - связанная вода, Sor - остаточная нефть
// k_rw = krw_end * ((S_w - Swc) / (1 - Swc - Sor))^exp_w
// k_ro = kro_end * ((1 - S_w - Sor) / (1 - Swc - Sor))^exp_o

double FluidProperties::krw(double s_w) const {
    if (s_w <= s_wc) {
        return 0.0;
    }
    if (s_w >= 1.0 - s_or) {
        return krw_max;
    }
    double s_norm = (s_w - s_wc) / (1.0 - s_wc - s_or);
    return krw_max * std::pow(s_norm, n_w);
}

double FluidProperties::kro(double s_w) const {
    if (s_w >= 1.0 - s_or) {
        return 0.0;
    }
     if (s_w <= s_wc) {
        return kro_max;
    }
    double s_norm = (1.0 - s_w - s_or) / (1.0 - s_wc - s_or);
    return kro_max * std::pow(s_norm, n_o);
}

double FluidProperties::capillaryPressure(double s_w) const {
    if (s_w >= 1.0 - s_or) {
        return 0.0;
    }
    // Effective saturation for capillary pressure
    double s_eff = (1.0 - s_w - s_or) / (1.0 - s_wc - s_or);
    if (s_eff <= 0) return 1e7; // A large pressure to prevent s_eff <= 0
    
    // Simple Brooks-Corey-like model
    // Using the same exponent 'n_o' as for oil rel-perm for simplicity
    return pc_entry * std::pow(s_eff, -1.0/n_o);
}

double FluidProperties::waterViscosity(double pressure) const {
    const double cP_to_PaS = 0.001;
    // Конвертируем Паскали в psi для формулы
    double pressure_psi = (pressure - p_ref_w) * 0.000145038;
    // Простая экспоненциальная модель вязкости
    return mu_ref_w * std::exp(c_w * pressure_psi) * cP_to_PaS;
}

double FluidProperties::oilViscosity(double pressure) const {
    const double cP_to_PaS = 0.001;
    double pressure_psi = (pressure - p_ref_o) * 0.000145038;
    return mu_ref_o * std::exp(c_o * pressure_psi) * cP_to_PaS;
} 