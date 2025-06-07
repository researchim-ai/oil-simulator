#pragma once
#include <vector>
#include <string>
#include "Parameters.h"
#include <Eigen/Dense>

// EN: Represents the 2D Cartesian grid of the reservoir.
// RU: Представляет 2D Декартову сетку пласта.
class Grid {
public:
    // EN: Constructor. Creates the grid based on simulation parameters.
    // RU: Конструктор. Создает сетку на основе параметров симуляции.
    Grid(const Parameters& params);

    // EN: Returns the linear index for a cell at (i, j).
    // RU: Возвращает линейный индекс для ячейки с координатами (i, j).
    int idx(int i, int j) const;

    // EN: Calculates the (i, j) coordinates from a linear index k.
    // RU: Вычисляет координаты (i, j) по линейному индексу k.
    void get_ij(int k, int& i, int& j) const;

    // EN: Gets and sets the pressure value at cell (i, j).
    // RU: Получает и устанавливает значение давления в ячейке (i, j).
    double& pressure(int i, int j);
    double pressure(int i, int j) const;

    // EN: Gets and sets the water saturation value at cell (i, j).
    // RU: Получает и устанавливает значение насыщенности водой в ячейке (i, j).
    double& saturation(int i, int j);
    double saturation(int i, int j) const;

    // EN: Gets the porosity at cell (i, j). Currently uniform.
    // RU: Получает пористость в ячейке (i, j). В данный момент однородная.
    double getPorosity(int i, int j) const;

    // EN: Gets the permeability at cell (i, j). Currently uniform.
    // RU: Получает проницаемость в ячейке (i, j). В данный момент однородная.
    double getPermeability(int i, int j) const;

    // EN: Returns the grid dimension in the x-direction (number of cells).
    // RU: Возвращает размер сетки по оси X (количество ячеек).
    int getNx() const;

    // EN: Returns the grid dimension in the y-direction (number of cells).
    // RU: Возвращает размер сетки по оси Y (количество ячеек).
    int getNy() const;

    // EN: Returns the physical size of a cell in the x-direction.
    // RU: Возвращает физический размер ячейки по оси X.
    double getDx() const;
    
    // EN: Returns the physical size of a cell in the y-direction.
    // RU: Возвращает физический размер ячейки по оси Y.
    double getDy() const;

    // EN: Returns the physical thickness of the cells (z-direction).
    // RU: Возвращает физическую толщину ячеек (по оси Z).
    double getDz() const;

    // EN: Returns the entire water saturation vector.
    // RU: Возвращает весь вектор насыщенности водой.
    Eigen::VectorXd get_sw_vector() const;

    // EN: Sets the water saturation for all cells from a vector.
    // RU: Устанавливает насыщенность водой для всех ячеек из вектора.
    void set_sw_vector(const Eigen::VectorXd& sw_vector);

    // EN: Sets the pressure for all cells from a vector.
    // RU: Устанавливает давление для всех ячеек из вектора.
    void set_p_vector(const Eigen::VectorXd& p_vector);

    void print() const;

private:
    int m_nx, m_ny;         // EN: Number of cells in X and Y directions. / RU: Количество ячеек по X и Y.
    double m_dx, m_dy, m_h; // EN: Cell dimensions in meters. / RU: Размеры ячеек в метрах.

    std::vector<double> m_pressure;    // EN: Pressure values for each cell. / RU: Значения давления для каждой ячейки.
    std::vector<double> m_saturation;  // EN: Saturation values for each cell. / RU: Значения насыщенности для каждой ячейки.
    std::vector<double> m_porosity;    // EN: Porosity values for each cell. / RU: Значения пористости для каждой ячейки.
    std::vector<double> m_permeability;// EN: Permeability values for each cell. / RU: Значения проницаемости для каждой ячейки.
}; 