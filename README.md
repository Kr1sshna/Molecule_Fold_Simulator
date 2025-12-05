# Molecular Fold Simulator 🧬

![Language](https://img.shields.io/badge/language-C-blue?style=for-the-badge&logo=c)
![License](https://img.shields.io/badge/license-MIT-green?style=for-the-badge)
![Field](https://img.shields.io/badge/Field-Biophysics-purple?style=for-the-badge)

> A simplified computational biophysics engine to simulate molecular folding and energy minimization, written in C.

## 🔬 About The Project

This repository contains a C-based simulation of molecular folding. The project explores how linear chains of atoms (representing amino acids) collapse into stable structures based on simplified interaction rules.

This simulation is a study in **Theoretical Biology** and **Computational Physics**, focusing on:
* **Conformational Search:** How a molecule explores space to find stable states.
* **Energy Landscapes:** Visualizing the minimization of potential energy.
* **Stochastic Processes:** Utilizing random walk or Monte Carlo methods to drive the folding process.

It serves as a tool to understand the computational complexity behind **Levinthal's Paradox** and the protein folding problem.

## 📊 Features

* **Dynamic Chain Generation:** Simulates atomic chains of varying lengths.
* **Energy Minimization Logic:** Calculates interaction energies (e.g., hydrophobic interactions or simple distance constraints).
* **Terminal Visualization:** Outputs the coordinates  in the trajectory.xyz file and radius of gyration in the analysis.csv.
* **Lightweight C Code:** Zero dependencies (standard libraries only) for high-performance execution.

## 🚀 Getting Started

Follow these instructions to compile and run the simulator on your local machine.

### Prerequisites
* A C Compiler (GCC, Clang, or MSVC) and VMD or any molecular dynamics viewing platform.

### Installation & Compilation

1.  **Clone the repository**
    ```bash
    git clone https://github.com/kr1sshna/Molecular_Fold_Simulator.git
    ```

2.  **Navigate to the folder**
    ```bash
    cd Molecular_Fold_Simulator
    ```

3.  **Compile the simulation**
    *Note: The `-lm` flag is often required to link the math library for physics calculations in unix/linux.*
    ```bash
    gcc molecule_fold_sim.c -o fold_sim -lm
    ```

## 🧪 Usage

Run the compiled executable:

```bash
./fold_sim
```
## 🔮 Future Roadmap (Bioinformatics)

    [ ] Scaling-Up: Scaling up to measure big proteins.

    [ ] Better Force Approximation: Implement Better force approximation for a protein in solution.

    [ ] Free Energy Landscape: Calculate the Free Energy Landscape of the protein conformational changes.

## 🤝 Contributing
This is an open-source project. If you are interested in Computational Biology and want to improve the energy functions or optimization algorithms, feel free to fork and PR!

## 📜 License
Distributed under the MIT License. See LICENSE for more information.

## 👤 Author
Krishna

    Interests: Proteomics, Biophysics, Biomathematics, Simulation.

    GitHub: @kr1sshna
