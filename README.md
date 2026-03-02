# Finite-element-method-solver
C++ implementation of the Finite Element Method (MES/FEM) for thermal analysis, including grid generation and simulation parameters.

# Finite Element Method (MES) Solver

## Project Overview
This project is a C++ implementation of the **Finite Element Method (MES)** designed for thermal conduction simulations. It was developed as part of the "Technologie Informatyczne" or "Metody Numeryczne" coursework at **AGH University of Science and Technology**.

## Features
* **Grid Processing:** Parses simulation data including nodes, elements, and boundary conditions.
* **Thermal Analysis:** Calculates temperature distribution based on conductivity, density, and specific heat.
* **Simulation Control:** Supports custom simulation time, step time, and initial temperatures.
* **Input Handling:** Reads data from structured text files (e.g., MixGrid format).

## Technical Details
* **Language:** C++
* **Development Environment:** Microsoft Visual Studio
* **Key Parameters handled:** Conductivity, Alfa, Density, Specific Heat, Initial Temperature.

## Example Data
The solver processes input files containing node coordinates and element definitions, such as:
* Nodes number: 16
* Elements number: 9
* Initial Temp: 100
* Simulation Time: 500s

## How to Run
1. Open `MES.sln` in Visual Studio.
2. Build the project in Release/Debug mode.
3. Provide an input file (e.g., from the `data/` folder) to run the simulation.
