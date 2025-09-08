# Codes for "Projected branes as platforms for crystalline, superconducting, and higher-order topological phases"
https://arxiv.org/pdf/2507.23783

This project uses Jupyter notebooks with the Julia programming language.

## Building block matrices

The building blocks for the tight binding model are available at `generate_matrices_1D.jl`, `generate_matrices_2D.jl` and `generate_matrices_dislocation2D.jl`. The code for generating spin matrices for arbitrary angular momentum is available at `angmom.jl`.

## Jupyter notebooks:

### Topological Crystalline Insulator
The code for generating the eigensystem is available at `PTB_hscti.ipynb`.

### Topological superconductors
- 
- PTB_ssh2d_SC_phase_weak_invariant`: Computes the weak invariant

### Higher Order (second order) Topological Insulators
- `2D_higher_order_TI.ipynb`: Generates HOTI on a parent square lattice.
- `PTB_2D_higher_order_TI`: Code for generating the higher order topological insulators on the projected brane.
- `PTB_2D_higher_order_TI_multi_runs.ipynb` - Computes the localizer index as a function of $m_0$.

## Codes for H11
- Topological crystalline insulator: `PTB_hscti_H11.ipynb`
- Topological superconductor: `PTB_ssh2d_SC_phase_weak_invariant_H11_only.ipynb` and `PTB_topo_super_H11.ipynb`

Data for the parent lattice of Higher order topological insulator (HOTI) is available in the directory `HOTI_parent_data`. All other data for the projected brane are available in the folder `HOTI`.

## Mathematica codes for plotting

While the Julia code is utilized to generate data, we use Mathematica to plot these data. These are available in the directory `mathematica plot files`.