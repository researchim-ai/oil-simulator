#pragma once

// EN: Manages the properties of oil and water, and their relative permeabilities.
// RU: Управляет свойствами нефти и воды и их относительными проницаемостями.
class FluidProperties {
public:
    // EN: Constructor. Sets default fluid and rock-fluid properties.
    // RU: Конструктор. Устанавливает стандартные свойства флюидов и породы-флюида.
    FluidProperties();

    // EN: Calculates oil viscosity as a function of pressure.
    // RU: Рассчитывает вязкость нефти как функцию давления.
    double oilViscosity(double p) const;

    // EN: Calculates water viscosity as a function of pressure.
    // RU: Рассчитывает вязкость воды как функцию давления.
    double waterViscosity(double p) const;

    // EN: Calculates oil relative permeability based on the Corey model.
    // RU: Рассчитывает относительную проницаемость по нефти по модели Кори.
    double kro(double s_w) const;

    // EN: Calculates water relative permeability based on the Corey model.
    // RU: Рассчитывает относительную проницаемость по воде по модели Кори.
    double krw(double s_w) const;

    // EN: Connate water saturation (irreducible).
    // RU: Связанная (остаточная) водонасыщенность.
    double s_wc;

    // EN: Residual oil saturation.
    // RU: Остаточная нефтенасыщенность.
    double s_or;

private:
    // EN: Reference pressure for oil viscosity calculation.
    // RU: Референсное давление для расчета вязкости нефти.
    double p_ref_o;

    // EN: Oil viscosity at reference pressure.
    // RU: Вязкость нефти при референсном давлении.
    double mu_ref_o;

    // EN: Oil compressibility factor for viscosity.
    // RU: Коэффициент сжимаемости нефти для вязкости.
    double c_o;

    // EN: Reference pressure for water viscosity calculation.
    // RU: Референсное давление для расчета вязкости воды.
    double p_ref_w;

    // EN: Water viscosity at reference pressure.
    // RU: Вязкость воды при референсном давлении.
    double mu_ref_w;

    // EN: Water compressibility factor for viscosity.
    // RU: Коэффициент сжимаемости воды для вязкости.
    double c_w;

    // EN: Corey model exponent for oil.
    // RU: Показатель степени Кори для нефти.
    double n_o;

    // EN: Corey model exponent for water.
    // RU: Показатель степени Кори для воды.
    double n_w;

    // EN: Maximum relative permeability for water (endpoint).
    // RU: Максимальная относительная проницаемость для воды.
    double krw_max;

    // EN: Maximum relative permeability for oil (endpoint).
    // RU: Максимальная относительная проницаемость для нефти.
    double kro_max;
}; 