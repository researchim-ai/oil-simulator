#include "Grid.h"
#include <iostream>
#include <vector>
#include <Eigen/Dense>

Grid::Grid(const Parameters& params)
    : m_nx(params.nx),
      m_ny(params.ny),
      m_dx(params.lx / params.nx),
      m_dy(params.ly / params.ny),
      m_perm(params.perm),
      m_poro(params.poro),
      m_pressure(params.nx * params.ny, 0.0),
      m_saturation(params.nx * params.ny, 0.0)
{
    std::cout << "Grid created with size " << m_nx << "x" << m_ny << std::endl;
    std::cout << "Cell dimensions: dx = " << m_dx << "m, dy = " << m_dy << "m" << std::endl;
}

void Grid::print() const {
    std::cout << "Grid size: " << m_nx << "x" << m_ny << std::endl;
    std::cout << "Sample values (P - pressure, S - saturation):" << std::endl;
    for (int i = 0; i < std::min(m_nx, 5); ++i) {
        for (int j = 0; j < std::min(m_ny, 5); ++j) {
            std::cout << "P(" << i << "," << j << ") = " << pressure(i, j)
                      << ", S(" << i << "," << j << ") = " << saturation(i, j) << std::endl;
        }
    }
}

int Grid::getNx() const {
    return m_nx;
}

int Grid::getNy() const {
    return m_ny;
}

double Grid::getDx() const {
    return m_dx;
}

double Grid::getDy() const {
    return m_dy;
}

double Grid::getPermeability(int i, int j) const {
    // For now, permeability is homogeneous
    return m_perm;
}

double Grid::getPorosity(int i, int j) const {
    // For now, porosity is homogeneous
    return m_poro;
}

double& Grid::pressure(int i, int j) {
    return m_pressure[i * m_ny + j];
}

const double& Grid::pressure(int i, int j) const {
    return m_pressure[i * m_ny + j];
}

double& Grid::saturation(int i, int j) {
    return m_saturation[i * m_ny + j];
}

const double& Grid::saturation(int i, int j) const {
    return m_saturation[i * m_ny + j];
}

// --- Implementation of new methods ---

int Grid::idx(int i, int j) const {
    return i * m_ny + j;
}

void Grid::get_ij(int k, int& i, int& j) const {
    i = k / m_ny;
    j = k % m_ny;
}

Eigen::VectorXd Grid::get_sw_vector() const {
    return Eigen::Map<const Eigen::VectorXd>(m_saturation.data(), m_saturation.size());
}

void Grid::set_p_vector(const Eigen::VectorXd& p) {
    Eigen::VectorXd::Map(&m_pressure[0], m_pressure.size()) = p;
}

void Grid::set_sw_vector(const Eigen::VectorXd& s) {
    Eigen::VectorXd::Map(&m_saturation[0], m_saturation.size()) = s;
} 