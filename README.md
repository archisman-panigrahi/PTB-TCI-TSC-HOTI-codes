# Codes for "Projected branes as platforms for crystalline, superconducting, and higher-order topological phases"

Preprint: https://arxiv.org/pdf/2507.23783

Published version: https://doi.org/10.1103/vlk8-59x3

This project uses Jupyter notebooks with the [Julia](https://julialang.org/) programming language.

## Building block matrices

The building blocks for the tight binding model are available at `generate_matrices_1D.jl`, `generate_matrices_2D.jl` and `generate_matrices_dislocation2D.jl`. The code for generating spin matrices for arbitrary angular momentum is available at `angmom.jl`. The file `Hermitian_Check.jl` contains a function to check whether a matrix is Hermitian.

## Jupyter notebooks:

- `PTB_brane_demo.ipynb` demonstrates how the points in the brane are takes, with a simple square lattice.

### Topological Crystalline Insulator
- `PTB_hscti.ipynb`: Computes the eigensystem of the topological crystalline insulator in a projected brane.
- `PTB_hscti_dislocation.ipynb`: Computes the eigensystem of the topological crystalline insulator when a dislocation core falls in a projected brane.
- `hscti_Chern_phase_diagram.ipynb`: Generates the phase diagram of topological crystalline insulator using the Fukui-Hatsugai-Suzuki method.

### Topological superconductors
- `PTB_topo_super.ipynb`: Computes the eigensystem of the topolgoical superconductor in a projected brane.
- `PTB_topo_super_dislocation.ipynb`: Computes the eigensystem of the topolgoical superconductor in a projected brane with a dislocaiton core.
- `PTB_tipo_super_commute_projection.ipynb`: Computes the strong topological invariant where the normal phase Hamiltonian is first projected to the brane, and then superconductivity is induced by manually adding the on-site pairing terms.
- `PTB_ssh2d_SC_phase_weak_invariant.ipynb`: Computes the weak invariant
- `PTB_ssh2d_SC_phase_weak_invariant_commute_order.ipynb`: Computes the weak invariant where the superconductor is formed by first projecting a normal topological phase to the brane, and then superconductivity is induced by manually adding the on-site paring terms.
- `BHZ_model_superconductivity_phase_diagram.ipynb`: Generates the phase diagram of topological superconductor by computing the Chern number using the Fukui-Hatsugai-Suzuki method.

### Higher Order (second order) Topological Insulators
- `2D_higher_order_TI.ipynb`: Generates HOTI on a parent square lattice.
- `PTB_2D_higher_order_TI.ipynb`: Code for generating the higher order topological insulators on the projected brane.
- `PTB_2D_higher_order_TI_multi_runs.ipynb` - Computes the localizer index as a function of $m_0$.

## Codes for H11

These are the coeds describing the systems where we only consider $H_{11}$, i.e., we don't renormalize the hoppings.

- `PTB_hscti_H11.ipynb`: Topological crystalline insulator
- `PTB_topo_super_H11.ipynb`: The strong topological superconductor
- `PTB_ssh2d_SC_phase_weak_invariant_H11_only.ipynb`: The weak topological superconductor 

## Data

Data for the parent lattice of Higher order topological insulator (HOTI) is available in the directory `HOTI_parent_data`. All other data for the projected brane are available in the folder `data`. The data for topological crystalline insulators, topological superconductor, and higher order topological insulator are available at `data/hstci`, ``data/topo_super`, and `data/HOTI`, respectively.
The data for scaling of localizer index with system size is present at `data/HOTI/localizer_index/scaling_with_system_size/`.

## Mathematica codes for plotting

While the Julia code is utilized to generate data, we use Mathematica to plot these data. These are available in the directory `mathematica plot files`. Some Mathematica notebooks are also available inside the directory containing data files.
