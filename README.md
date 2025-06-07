# Oil Simulator

This is a simple 2D, two-phase (oil and water) reservoir simulator written in C++. It uses an IMPES (Implicit in Pressure, Explicit in Saturation) formulation.

## Physics & Numerics

The simulator solves the standard black-oil model equations for two-phase flow.

### Mathematical Model

The simulator implements a two-phase (oil and water) model in a 2D Cartesian grid. The mathematical model consists of two main equations solved sequentially.

The pressure equation is derived from the conservation of mass for both phases and is solved implicitly:
$$
\nabla \cdot \left( \lambda_t \mathbf{K} \nabla p \right) + q_t = \phi c_t \frac{\partial p}{\partial t}
$$

The saturation equation is solved explicitly after the pressure field is known:
$$
\phi \frac{\partial S_w}{\partial t} + \nabla \cdot \left( \mathbf{v}_t f_w \right) + q_w = 0
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

### Other Models

- **Relative Permeability**: Corey model.
- **Fluid Viscosity**: Modeled as a function of pressure.
- **Linear Solver**: Eigen's BiCGSTAB for the pressure equation.
- **Parallelism**: Saturation update is parallelized using OpenMP.

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

### Running a Simulation

The main executable `simulator` takes one argument: the path to a parameter file.

From inside the `build` directory, run:
```bash
./simulator ../data/five_spot.dat
```
This will:
1.  Run the simulation based on the parameters in `five_spot.dat`.
2.  Save the final pressure and saturation for each grid cell into `results.txt` in the project's root directory.
3.  Automatically call the Python visualization script.
4.  Generate `pressure_map.png` and `saturation_map.png` in the project's root directory.

### Running Tests

The `test_runner` executable runs all compiled unit tests.

From inside the `build` directory, run:
```bash
./test_runner
```
This will execute the tests and report if they pass or fail.

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