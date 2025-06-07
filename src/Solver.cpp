#include "Solver.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <Eigen/Dense> // For VectorXd
#include <Eigen/Sparse>
#include <cstdio> // For fprintf

using Triplet = Eigen::Triplet<double>;

Solver::Solver(const Parameters& params, Grid& grid)
    : params(params), grid(grid) {
    int num_cells = params.nx * params.ny;
    A.resize(num_cells, num_cells);
    b.resize(num_cells);
    p_new.resize(num_cells);
    std::cout << "IMPES Solver created." << std::endl;

    // This needs the grid to be finalized.
    for (auto& well : wells) {
        well.calculateWellIndex(grid);
    }
}

void Solver::add_well(const Well& well) {
    wells.push_back(well);
}

void Solver::solve() {
    // Before starting the simulation, calculate the well indices for all wells.
    // This needs the grid to be finalized.
    for (auto& well : wells) {
        well.calculateWellIndex(grid);
    }

    double time = 0;
    int step_count = 0;
    while (time < params.total_time) {
        // For now, use a fixed small timestep for stability.
        // Auto dt calculation (CFL) needs to be revisited for IMPES.
        double dt = params.dt;
        if (time + dt > params.total_time) {
            dt = params.total_time - time;
        }

        solve_step(dt);
        time += dt;
        step_count++;

        std::cout << "Step: " << step_count
                  << ", Time: " << time / 86400.0 << " days"
                  << ", dt: " << dt << " s" << std::endl;
    }
}

void Solver::solve_step(double dt) {
    assemble_system(dt);
    solve_pressure();
    update_saturation(dt);
}

void Solver::assemble_system(double dt) {
    const int nx = params.nx;
    const int ny = params.ny;
    const int num_cells = nx * ny;
    
    // Clear previous system
    A.setZero();
    b.setZero();
    std::vector<Triplet> triplets;
    triplets.reserve(num_cells * 5);
    
    const double dx = params.lx / nx;
    const double dy = params.ly / ny;
    const double h = params.h;
    const double rock_compressibility = params.cr; 
    const double fluid_compressibility = 3e-6 / 6894.76;
    const double total_compressibility = rock_compressibility + fluid_compressibility;

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            const int k = grid.idx(i, j);
            const double s_w_old = grid.saturation(i, j);

            // 1. Accumulation Term
            const double Bw = 1.0;
            const double accumulation_coeff = (grid.getPorosity(i,j) * total_compressibility) / (dt * Bw);
            double diagonal_element = accumulation_coeff;
            b(k) = accumulation_coeff * grid.pressure(i, j);

            // 2. Transmissibility & Capillary Pressure Terms
            auto get_trans_and_pc = [&](int ci, int cj, int ni, int nj, double cross_sectional_area, double distance) {
                double p_c = grid.pressure(ci, cj);
                double p_n = grid.pressure(ni, nj);
                double s_c = grid.saturation(ci, cj);
                double s_n = grid.saturation(ni, nj);

                // Upstream weighting for mobility
                double s_up, p_up;
                if (p_c > p_n) { // Flow from c to n
                    s_up = s_c;
                    p_up = p_c;
                } else { // Flow from n to c
                    s_up = s_n;
                    p_up = p_n;
                }
                
                double total_mobility = (fluid.krw(s_up) / fluid.waterViscosity(p_up)) + (fluid.kro(s_up) / fluid.oilViscosity(p_up));
                
                double perm_c = grid.getPermeability(ci, cj);
                double perm_n = grid.getPermeability(ni, nj);
                double avg_perm = (2.0 * perm_c * perm_n) / (perm_c + perm_n); // Harmonic mean for permeability
                
                double T_total = avg_perm * total_mobility * cross_sectional_area / distance;

                // Capillary term (added to RHS)
                double mob_w_up = fluid.krw(s_up) / fluid.waterViscosity(p_up);
                double T_w = avg_perm * mob_w_up * cross_sectional_area / distance;
                double pc_c = fluid.capillaryPressure(s_c);
                double pc_n = fluid.capillaryPressure(s_n);
                double cap_term = T_w * (pc_n - pc_c);

                return std::make_pair(T_total, cap_term);
            };

            double T_right, T_left, T_top, T_bottom;
            double Pc_term_right = 0, Pc_term_left = 0, Pc_term_top = 0, Pc_term_bottom = 0;

            if (i < nx - 1) { 
                auto [T, Pc_term] = get_trans_and_pc(i, j, i + 1, j, dy * h, dx);
                T_right = T;
                Pc_term_right = Pc_term;
                diagonal_element += T_right;
                triplets.emplace_back(k, grid.idx(i + 1, j), -T_right);
            }
            if (i > 0) {
                auto [T, Pc_term] = get_trans_and_pc(i, j, i - 1, j, dy * h, dx);
                T_left = T;
                Pc_term_left = -Pc_term; // Note the sign change
                diagonal_element += T_left;
                triplets.emplace_back(k, grid.idx(i - 1, j), -T_left);
            }
            if (j < ny - 1) {
                auto [T, Pc_term] = get_trans_and_pc(i, j, i, j + 1, dx * h, dy);
                T_top = T;
                Pc_term_top = Pc_term;
                diagonal_element += T_top;
                triplets.emplace_back(k, grid.idx(i, j + 1), -T_top);
            }
            if (j > 0) {
                auto [T, Pc_term] = get_trans_and_pc(i, j, i, j - 1, dx * h, dy);
                T_bottom = T;
                Pc_term_bottom = -Pc_term; // Note the sign change
                diagonal_element += T_bottom;
                triplets.emplace_back(k, grid.idx(i, j - 1), -T_bottom);
            }
            
            b(k) += (Pc_term_right + Pc_term_left + Pc_term_top + Pc_term_bottom);

            // 3. Well Terms
            for (const auto& well : wells) {
                if (well.getI() == i && well.getJ() == j) {
                    double s_w = grid.saturation(i, j);
                    double p = grid.pressure(i,j);
                    double total_mobility = (fluid.krw(s_w) / fluid.waterViscosity(p)) + (fluid.kro(s_w) / fluid.oilViscosity(p));
                    double T_w = well.getWellIndex() * total_mobility;
                    
                    diagonal_element += T_w;
                    b(k) += T_w * well.getBHP();
                }
            }

            // 4. Add diagonal entry
            triplets.emplace_back(k, k, diagonal_element);
        }
    }

    A.setFromTriplets(triplets.begin(), triplets.end());
}

void Solver::solve_pressure() {
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    p_new = solver.solve(b);

    if(solver.info() != Eigen::Success) {
        std::cerr << "Pressure solver failed to converge. Info: " << solver.info() << std::endl;
        // Optionally, handle the failure, e.g., by returning or throwing an exception
    }
}

void Solver::update_saturation(double dt) {
    int nx = params.nx;
    int ny = params.ny;
    const double dx = params.lx / nx;
    const double dy = params.ly / ny;
    const double cell_volume = dx * dy * params.h;

    Eigen::VectorXd s_old = grid.get_sw_vector();
    Eigen::VectorXd s_new = s_old;

    auto calculate_flux = [&](int i_from, int j_from, int i_to, int j_to, double L, double area) {
        // Harmonic average of permeability
        double perm_from = grid.getPermeability(i_from, j_from);
        double perm_to = grid.getPermeability(i_to, j_to);
        double avg_perm = (2.0 * perm_from * perm_to) / (perm_from + perm_to);

        double p_from = p_new(grid.idx(i_from, j_from));
        double p_to = p_new(grid.idx(i_to, j_to));
        
        // Upstream weighting for saturation and mobility
        int i_up, j_up;
        double p_up;
        if (p_from > p_to) {
            i_up = i_from; j_up = j_from; p_up = p_from;
        } else {
            i_up = i_to; j_up = j_to; p_up = p_to;
        }
        
        double s_up = s_old(grid.idx(i_up, j_up));
        double mob_w = fluid.krw(s_up) / fluid.waterViscosity(p_up);
        
        // Adding capillary pressure effect to the water flux
        double s_from = s_old(grid.idx(i_from, j_from));
        double s_to = s_old(grid.idx(i_to, j_to));
        double pc_from = fluid.capillaryPressure(s_from);
        double pc_to = fluid.capillaryPressure(s_to);

        return -mob_w * avg_perm * area * ( (p_to - p_from) - (pc_to - pc_from) ) / L;
    };

    #pragma omp parallel for
    for (int k = 0; k < nx * ny; ++k) {
        int i, j;
        grid.get_ij(k, i, j);

        const double poro = grid.getPorosity(i,j);
        
        // Calculate fluxes for all cell faces, assuming no-flow at boundaries
        double flux_right = (i < nx - 1) ? calculate_flux(i, j, i + 1, j, dx, dy * params.h) : 0.0;
        double flux_left  = (i > 0)      ? calculate_flux(i - 1, j, i, j, dx, dy * params.h) : 0.0;
        double flux_top   = (j < ny - 1) ? calculate_flux(i, j, i, j + 1, dy, dx * params.h) : 0.0;
        double flux_bottom= (j > 0)      ? calculate_flux(i, j - 1, i, j, dy, dx * params.h) : 0.0;
        
        double divergence_term = (flux_right - flux_left) + (flux_top - flux_bottom);
        
        double q_w = 0.0; // water source term from wells in m^3/s
        for (const auto& well : wells) {
            if (well.getI() == i && well.getJ() == j) {
                double p_res = p_new(k);
                double s_w = s_old(k);

                // Determine if it's an injector or producer based on potential
                bool is_injector = well.getBHP() > p_res;

                if (is_injector) {
                    // For an injector, flow rate is determined by the mobility of the injected fluid (pure water)
                    double mob_water_inj = fluid.krw(1.0) / fluid.waterViscosity(p_res);
                    double q_total = well.getWellIndex() * mob_water_inj * (well.getBHP() - p_res);
                    q_w = q_total; // The entire flow is water
                } else {
                    // For a producer, flow rate is based on the combined mobility of fluids in the cell
                    double mob_w = fluid.krw(s_w) / fluid.waterViscosity(p_res);
                    double mob_o = fluid.kro(s_w) / fluid.oilViscosity(p_res);
                    double total_mobility = mob_w + mob_o;
                    if (total_mobility > 1e-9) { // Avoid division by zero
                        double q_total = well.getWellIndex() * total_mobility * (well.getBHP() - p_res);
                        q_w = (mob_w / total_mobility) * q_total;
                    }
                }
            }
        }
        
        s_new(k) = s_old(k) - (dt / (poro * cell_volume)) * (divergence_term - q_w);
        s_new(k) = std::max(fluid.s_wc, std::min(1.0 - fluid.s_or, s_new(k)));
    }
    
    // Update pressure and saturation in the grid
    grid.set_p_vector(p_new);
    grid.set_sw_vector(s_new);
}


void Solver::save_results(const std::string& filename) const {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open output file " << filename << std::endl;
        return;
    }

    outfile << "i,j,pressure,saturation\n"; // Header
    for (int i = 0; i < params.nx; ++i) {
        for (int j = 0; j < params.ny; ++j) {
            outfile << i << "," << j << "," 
                    << grid.pressure(i, j) << "," 
                    << grid.saturation(i, j) << "\n";
        }
    }
    std::cout << "Results saved to " << filename << std::endl;
}