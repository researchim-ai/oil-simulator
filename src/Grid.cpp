#include "Grid.h"
#include <iostream>
#include <vector>

Grid::Grid(int width, int height) 
    : m_width(width), m_height(height), m_pressure(width * height, 0.0), m_saturation(width * height, 0.0)
{
    std::cout << "Grid created with size " << m_width << "x" << m_height << std::endl;
}

void Grid::print() const {
    std::cout << "Grid size: " << m_width << "x" << m_height << std::endl;
    std::cout << "Sample values (P - pressure, S - saturation):" << std::endl;
    for (int i = 0; i < std::min(m_width, 5); ++i) {
        for (int j = 0; j < std::min(m_height, 5); ++j) {
            std::cout << "P(" << i << "," << j << ") = " << pressure(i, j)
                      << ", S(" << i << "," << j << ") = " << saturation(i, j) << std::endl;
        }
    }
}

int Grid::getWidth() const {
    return m_width;
}

int Grid::getHeight() const {
    return m_height;
}

double& Grid::pressure(int i, int j) {
    return m_pressure[i * m_height + j];
}

const double& Grid::pressure(int i, int j) const {
    return m_pressure[i * m_height + j];
}

double& Grid::saturation(int i, int j) {
    return m_saturation[i * m_height + j];
}

const double& Grid::saturation(int i, int j) const {
    return m_saturation[i * m_height + j];
} 