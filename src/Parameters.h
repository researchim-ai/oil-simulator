#pragma once
#include <string>

// EN: Stores all simulation parameters loaded from a data file.
// RU: Хранит все параметры симуляции, загруженные из файла данных.
class Parameters {
public:
    // EN: Default constructor.
    // RU: Конструктор по умолчанию.
    Parameters() = default;

    // EN: Constructor that loads parameters from a file.
    // RU: Конструктор, загружающий параметры из файла.
    Parameters(const std::string& filename);

    int nx, ny;                 // EN: Number of blocks in X and Y / RU: Количество блоков по X и Y
    double lx, ly, h;           // EN: Reservoir dimensions in meters / RU: Размеры пласта в метрах
    double perm;                // EN: Permeability in m^2 / RU: Проницаемость в м^2
    double poro;                // EN: Porosity (fraction) / RU: Пористость (доли единицы)
    double cr;                  // EN: Rock compressibility in 1/Pa / RU: Сжимаемость породы в 1/Па
    double dt;                  // EN: Timestep size in seconds / RU: Размер шага по времени в секундах
    double total_time;          // EN: Total simulation time in seconds / RU: Общее время симуляции в секундах
    double initial_pressure;
    double initial_saturation;
}; 