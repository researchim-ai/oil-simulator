# Oil Simulator

This is a simple 2D, two-phase (oil and water) reservoir simulator written in C++. It uses an IMPES (Implicit in Pressure, Explicit in Saturation) formulation.

## Physics & Numerics

The simulator solves the standard black-oil model equations for two-phase flow.

### Mathematical Model

The simulator implements a two-phase (oil and water) model in a 2D Cartesian grid. The mathematical model consists of two main equations solved sequentially.

The pressure equation is derived from the conservation of mass for both phases and is solved implicitly. It now includes the effect of capillary pressure ($P_c = P_o - P_w$):
$$
\nabla \cdot \left( \mathbf{K} \lambda_w \nabla p_w \right) + q_w = \phi \frac{\partial S_w}{\partial t}
$$
$$
\nabla \cdot \left( \mathbf{K} \lambda_o \nabla (p_w + P_c) \right) + q_o = \phi \frac{\partial S_o}{\partial t}
$$
Summing these and assuming total velocity formulation leads to a pressure equation for one of the phases (e.g., water pressure $p_w$) and an explicit saturation equation.

The saturation equation is solved explicitly after the pressure field is known:
$$
\phi \frac{\partial S_w}{\partial t} + \nabla \cdot \left( \mathbf{v}_w \right) + q_w = 0
$$

Where:
- $\lambda_t = \lambda_w + \lambda_o$ is the total mobility, where $\lambda_p = k_{rp} / \mu_p$ is the mobility of phase $p$.
- $\mathbf{K}$ is the permeability tensor.
- $p$ is the pressure.
- $q_t$ is the total source term from wells.
- $\phi$ is the porosity.
- $c_t$ is the total compressibility (rock + fluids).
- $S_w$ is the water saturation.
- $\mathbf{v}_t$ is the total Darcy velocity.
- $f_w$ is the fractional flow function for water.
- $q_w$ is the source term for water.
- **Relative Permeability**: Corey model.
- **Capillary Pressure**: Brooks-Corey model.
- **Fluid Viscosity**: Modeled as a function of pressure.
- **Linear Solver**: Eigen's BiCGSTAB for the pressure equation.
- **Parallelism**: Saturation update is parallelized using OpenMP.

### Well Model

Wells are modeled using the standard **Peaceman well model**. The source/sink term for a well is calculated as:
$$
q = WI \cdot \lambda \cdot (p_{bh} - p_{block})
$$
Where $WI$ is the well index, calculated as:
$$
WI = \frac{2 \pi k h}{\ln(r_e / r_w)}
$$
with the effective radius $r_e = 0.14 \sqrt{\Delta x^2 + \Delta y^2}$.

## Dependencies

- A C++17 compliant compiler (e.g., GCC, Clang)
- CMake (version 3.11 or higher)
- Git
- Python 3 with `matplotlib` and `numpy` for visualization.
  ```
  pip install matplotlib numpy
  ```

## How to Build

The project uses CMake for building.

1.  **Clone the repository:**
    ```bash
    git clone <repository_url>
    cd oil-simulator
    ```

2.  **Create a build directory:**
    ```bash
    mkdir build
    cd build
    ```

3.  **Run CMake and build:**
    ```bash
    cmake ..
    make
    ```
    This will create two main executables inside the `build` directory: `simulator` and `test_runner`.

## How to Run

The main executable `simulator` takes one argument: the path to a parameter file.

From inside the `build` directory, run:
```bash
./simulator ../data/five_spot.dat
```

This starts the simulation. After it completes, it will generate output files in the **current directory** (i.e., `build/`).

## Simulation Process

The simulation follows the IMPES (Implicit in Pressure, Explicit in Saturation) scheme:

1.  **Initialization**:
    *   The program reads parameters from the specified `.dat` file, including grid dimensions, rock/fluid properties, and time control settings.
    *   A `Grid` object is created to store the geometry and properties (porosity, permeability, pressure, saturation) for each cell.
    *   Wells are placed on the grid according to the scenario (e.g., `create_five_spot` places one injector in the center and four producers at the corners).

2.  **Main Time Loop**:
    *   The simulator iterates through time with a fixed step `dt` until the `total_time` is reached.
    *   In each time step, two primary equations are solved sequentially:

    a.  **Pressure Equation (Implicit Solve)**:
        *   A system of linear equations `Ax = b` is assembled, where `x` is the vector of cell pressures for the new time step.
        *   The `A` matrix (transmissibility matrix) accounts for permeability, fluid viscosity, and cell geometry.
        *   The right-hand side vector `b` includes terms for fluid accumulation, well flows, and **capillary pressure gradients**.
        *   This system is solved using the Eigen library to find the new pressure field across the reservoir.

    b.  **Saturation Equation (Explicit Update)**:
        *   Using the newly calculated pressure field, inter-cell flow rates are computed based on Darcy's Law.
        *   These flows are a function of the pressure gradient, fluid mobilities, and permeability. The **capillary pressure** term is included here as well to make the flow more physically accurate.
        *   The water saturation of each cell is then explicitly updated based on the net flow into it.

3.  **Save Results**:
    *   After the final time step, the resulting pressure and saturation fields are saved.
    *   A Python script is called to visualize the data.

## Output Files Description

After a successful run, the following files are created in the `build/` directory:

*   `results.txt`:
    *   **Format**: A comma-separated values (CSV) text file.
    *   **Content**: Each row represents a single grid cell and contains four values:
        1.  `i`: The X-index of the cell.
        2.  `j`: The Y-index of the cell.
        3.  `pressure`: The final pressure in the cell (in Pascals).
        4.  `saturation`: The final water saturation in the cell (a fraction from 0 to 1).

*   `pressure_map.png`:
    *   **Format**: PNG image.
    *   **Content**: A 2D color map of the pressure distribution. The "viridis" colormap shows pressure variations, helping to identify high and low-pressure zones.

*   `saturation_map.png`:
    *   **Format**: PNG image.
    *   **Content**: A 2D color map of the water saturation distribution. The "jet" colormap visualizes the water front (red/orange tones) and remaining oil (blue tones), clearly showing how the water has moved from the injector towards the producers.

## Project Structure

- `src/`: Contains all C++ source code and headers.
- `data/`: Contains parameter files for simulation scenarios.
- `python/`: Contains the Python script for visualization.
- `tests/`: Contains C++ source code for unit tests.
- `build/`: Build directory (created by you).
- `results.txt`: Main simulation output file (generated after a run).
- `*.png`: Visualization output files (generated after a run).
- `CMakeLists.txt`: Main CMake build script.
- `README.md`: This file. 