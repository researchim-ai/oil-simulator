#include <iostream>
#include "Grid.h"
#include "Solver.h"
#include "Parameters.h"

int main() {
    std::cout << "Starting Oil Simulator..." << std::endl;

    // Загружаем параметры из файла
    Parameters params("data/simple.dat");

    // Создаем сетку и задаем начальные условия
    Grid grid(params.width, params.height);
    
    // Начальные условия: нагнетательная скважина на левой границе
    for (int j = 0; j < grid.getHeight(); ++j) {
        grid.pressure(0, j) = params.initial_pressure; // Высокое давление
        grid.saturation(0, j) = 1.0; // 100% воды
    }
    
    // Создаем решатель
    Solver solver(grid, params.alpha);

    std::cout << "\nInitial state:" << std::endl;
    grid.print();

    // Запускаем симуляцию
    double time = 0.0;
    int step_count = 0;
    while (time < params.total_time) {
        solver.step(params.dt);
        time += params.dt;
        step_count++;
    }

    std::cout << "\nFinal state after " << step_count << " steps:" << std::endl;
    grid.print();
    
    std::cout << "Oil Simulator finished." << std::endl;
    return 0;
} 