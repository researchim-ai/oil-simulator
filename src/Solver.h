#pragma once
#include "Grid.h"
#include "Parameters.h"
#include "FluidProperties.h"
#include "Well.h"
#include <vector>
#include <Eigen/Sparse>

// EN: Main class for solving the reservoir simulation problem.
// RU: Основной класс для решения задачи гидродинамического моделирования.
class Solver {
public:
    // EN: Constructor. Initializes the solver with parameters and a grid.
    // RU: Конструктор. Инициализирует решатель параметрами и сеткой.
    Solver(const Parameters& params, Grid& grid);

    // EN: Adds a well to the simulation.
    // RU: Добавляет скважину в симуляцию.
    void add_well(const Well& well);

    // EN: Runs the full simulation for the total time specified in parameters.
    // RU: Запускает полную симуляцию на общее время, указанное в параметрах.
    void solve();

    // EN: Runs a single timestep of the simulation.
    // RU: Выполняет один шаг симуляции по времени.
    void solve_step(double dt);
    
    // EN: Saves the final pressure and saturation maps to a file.
    // RU: Сохраняет итоговые карты давления и насыщенности в файл.
    void save_results(const std::string& filename) const;

    // EN: Returns a reference to the grid object.
    // RU: Возвращает ссылку на объект сетки.
    Grid& get_grid() { return grid; }

// EN: This section is temporarily public for testing purposes.
// RU: Эта секция временно сделана публичной для целей тестирования.
// private: 
    // EN: Assembles the pressure matrix (A) and the right-hand side vector (b).
    // RU: Собирает матрицу давления (A) и вектор правой части (b).
    void assemble_system(double dt);

    // EN: Solves the linear system Ax=b for the new pressure distribution.
    // RU: Решает линейную систему Ax=b для нового распределения давления.
    void solve_pressure();

    // EN: Explicitly updates the saturation based on the new pressure field.
    // RU: Явно обновляет насыщенность на основе нового поля давления.
    void update_saturation(double dt);

    // EN: Constant reference to the simulation parameters.
    // RU: Константная ссылка на параметры симуляции.
    const Parameters& params;

    // EN: Reference to the simulation grid.
    // RU: Ссылка на сетку симуляции.
    Grid& grid;

    // EN: Object containing fluid properties.
    // RU: Объект, содержащий свойства флюидов.
    FluidProperties fluid;

    // EN: Vector of all wells in the simulation.
    // RU: Вектор всех скважин в симуляции.
    std::vector<Well> wells;

    // EN: System matrix for the pressure equation (implicit part).
    // RU: Матрица системы для уравнения давления (неявная часть).
    Eigen::SparseMatrix<double> A;

    // EN: Right-hand side vector for the pressure equation.
    // RU: Вектор правой части для уравнения давления.
    Eigen::VectorXd b;

    // EN: Vector to store the new pressure values after solving.
    // RU: Вектор для хранения новых значений давления после решения.
    Eigen::VectorXd p_new;
}; 