#pragma once
#include <vector>

class Grid {
public:
    Grid(int width, int height);

    int getWidth() const;
    int getHeight() const;

    void print() const;

    double& pressure(int i, int j);
    const double& pressure(int i, int j) const;

    double& saturation(int i, int j);
    const double& saturation(int i, int j) const;

private:
    int m_width;
    int m_height;
    std::vector<double> m_pressure;
    std::vector<double> m_saturation;
}; 