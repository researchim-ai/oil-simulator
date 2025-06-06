#include "Solver.h"
#include <iostream>
#include <vector>

Solver::Solver(Grid& grid, double alpha) : m_grid(grid), m_alpha(alpha) {
    std::cout << "Solver created with alpha = " << m_alpha << std::endl;
}

void Solver::step(double dt) {
    const int width = m_grid.getWidth();
    const int height = m_grid.getHeight();
    const double dx = 1.0;
    const double dy = 1.0;

    Grid new_grid = m_grid; // Копия для новых значений

    // --- 1. Решение для давления ---
    // В двухфазной модели мы должны учитывать подвижность каждой фазы.
    // Коэффициент в уравнении теперь зависит от насыщенности.
    for (int i = 1; i < width - 1; ++i) {
        for (int j = 1; j < height - 1; ++j) {
            const double s_center = m_grid.saturation(i, j);
            const double total_mobility = m_fluid.krw(s_center) + m_fluid.kro(s_center);

            const double p_center = m_grid.pressure(i, j);
            const double p_left = m_grid.pressure(i - 1, j);
            const double p_right = m_grid.pressure(i + 1, j);
            const double p_bottom = m_grid.pressure(i, j - 1);
            const double p_top = m_grid.pressure(i, j + 1);

            const double laplacian = ((p_left - 2 * p_center + p_right) / (dx * dx)) +
                                     ((p_bottom - 2 * p_center + p_top) / (dy * dy));
            
            // Давление обновляется с учетом общей подвижности.
            // Для упрощения мы считаем подвижность константой на этом шаге,
            // что является большим упрощением, но подходит для начала.
            new_grid.pressure(i, j) = p_center + total_mobility * m_alpha * dt * laplacian;
        }
    }

    // --- 2. Решение для насыщенности ---
    // Используем новое поле давления (из new_grid) для расчета скоростей.
    for (int i = 1; i < width - 1; ++i) {
        for (int j = 1; j < height - 1; ++j) {
            // --- Upstream-схема для потока по X ---
            // Поток между (i-1, j) и (i, j)
            double s_upstream_x_left = (new_grid.pressure(i - 1, j) > new_grid.pressure(i, j)) ? m_grid.saturation(i - 1, j) : m_grid.saturation(i, j);
            double vx_left = -m_fluid.krw(s_upstream_x_left) * (new_grid.pressure(i, j) - new_grid.pressure(i - 1, j)) / dx;
            
            // Поток между (i, j) и (i+1, j)
            double s_upstream_x_right = (new_grid.pressure(i, j) > new_grid.pressure(i + 1, j)) ? m_grid.saturation(i, j) : m_grid.saturation(i + 1, j);
            double vx_right = -m_fluid.krw(s_upstream_x_right) * (new_grid.pressure(i + 1, j) - new_grid.pressure(i, j)) / dx;
            
            // --- Upstream-схема для потока по Y ---
            // Поток между (i, j-1) и (i, j)
            double s_upstream_y_bottom = (new_grid.pressure(i, j-1) > new_grid.pressure(i, j)) ? m_grid.saturation(i, j - 1) : m_grid.saturation(i, j);
            double vy_bottom = -m_fluid.krw(s_upstream_y_bottom) * (new_grid.pressure(i, j) - new_grid.pressure(i, j - 1)) / dy;
            
            // Поток между (i, j) и (i, j+1)
            double s_upstream_y_top = (new_grid.pressure(i, j) > new_grid.pressure(i, j+1)) ? m_grid.saturation(i, j) : m_grid.saturation(i, j+1);
            double vy_top = -m_fluid.krw(s_upstream_y_top) * (new_grid.pressure(i, j + 1) - new_grid.pressure(i, j)) / dy;

            // Дивергенция потока воды
            double div_v = (vx_right - vx_left) / dx + (vy_top - vy_bottom) / dy;
            
            // Явная схема для переноса насыщенности с ограничением
            double old_saturation = m_grid.saturation(i, j);
            double new_saturation = old_saturation - dt * div_v;
            new_grid.saturation(i, j) = std::max(0.0, std::min(1.0, new_saturation));
        }
    }

    // Обновляем состояние сетки
    m_grid = new_grid;
} 