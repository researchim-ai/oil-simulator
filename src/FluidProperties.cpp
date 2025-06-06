#include "FluidProperties.h"
#include <algorithm>

FluidProperties::FluidProperties() {}

double FluidProperties::krw(double saturation) const {
    // Ограничиваем насыщенность диапазоном [0, 1] для стабильности
    double s_norm = std::max(0.0, std::min(1.0, saturation));
    return s_norm;
}

double FluidProperties::kro(double saturation) const {
    // Ограничиваем насыщенность диапазоном [0, 1] для стабильности
    double s_norm = std::max(0.0, std::min(1.0, saturation));
    return 1.0 - s_norm;
} 