#pragma once
#include "Solver.h"
#include "Parameters.h"

// EN: A utility class to set up predefined simulation scenarios.
// RU: Вспомогательный класс для настройки предопределенных сценариев симуляции.
class Scenario {
public:
    // EN: Sets up the "five-spot" well pattern.
    // RU: Устанавливает "пятиточечную" схему расположения скважин.
    static void create_five_spot(Solver& solver, const Parameters& params);
}; 