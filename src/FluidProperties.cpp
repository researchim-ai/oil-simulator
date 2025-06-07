#include "FluidProperties.h"
#include <cmath>
#include <algorithm>

FluidProperties::FluidProperties() {
    // Задаем стандартные значения
    m_viscosity_water_ref = 1.0;   // cP
    m_viscosity_oil_ref = 10.0;  // cP
    m_compressibility_viscosity_water = 1e-6; // 1/psi
    m_compressibility_viscosity_oil = 1e-5;   // 1/psi
}

// Модель ОФП Кори
// Swc - связанная вода, Sor - остаточная нефть
// k_rw = krw_end * ((S_w - Swc) / (1 - Swc - Sor))^exp_w
// k_ro = kro_end * ((1 - S_w - Sor) / (1 - Swc - Sor))^exp_o

double FluidProperties::krw(double saturation) const {
    const double swc = 0.1;      // Связанная вода
    const double sor = 0.1;      // Остаточная нефтенасыщенность
    const double krw_end = 0.8;  // ОФП воды при остаточной нефте-ти
    const double exp_w = 2.0;    // Экспонента Кори для воды

    if (saturation <= swc) {
        return 0.0;
    }
    if (saturation >= 1.0 - sor) {
        return krw_end;
    }
    double s_norm = (saturation - swc) / (1.0 - swc - sor);
    return krw_end * std::pow(s_norm, exp_w);
}

double FluidProperties::kro(double saturation) const {
    const double swc = 0.1;
    const double sor = 0.1;
    const double kro_end = 1.0;  // ОФП нефти при связанной воде
    const double exp_o = 2.0;    // Экспонента Кори для нефти

    if (saturation >= 1.0 - sor) {
        return 0.0;
    }
     if (saturation <= swc) {
        return kro_end;
    }
    double s_norm = (1.0 - saturation - sor) / (1.0 - swc - sor);
    return kro_end * std::pow(s_norm, exp_o);
}

double FluidProperties::waterViscosity(double pressure) const {
    const double cP_to_PaS = 0.001;
    // Конвертируем Паскали в psi для формулы
    double pressure_psi = pressure * 0.000145038;
    // Простая экспоненциальная модель вязкости
    return m_viscosity_water_ref * std::exp(m_compressibility_viscosity_water * pressure_psi) * cP_to_PaS;
}

double FluidProperties::oilViscosity(double pressure) const {
    const double cP_to_PaS = 0.001;
    double pressure_psi = pressure * 0.000145038;
    return m_viscosity_oil_ref * std::exp(m_compressibility_viscosity_oil * pressure_psi) * cP_to_PaS;
} 