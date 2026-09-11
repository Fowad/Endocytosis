# Code and Publications

This document explains how the public repository relates to the two papers produced from the Ph.D. work.

## 1. Thermodynamic foundation

**Paper:** F. Motahari and A. E. Carlsson, *Physical Review E* **100**, 042409 (2019)  
DOI: https://doi.org/10.1103/PhysRevE.100.042409

The paper studies a growing biopolymer interacting with a smoothly varying obstacle potential \(U(r)\). Its central methodological result is the detailed-balance relation

\[
\frac{k_{\mathrm{on}}(r)}
{k_{\mathrm{off}}(r-\delta)}
=
\exp\left(
-\frac{U(r-\delta)-U(r)}{k_B T}
\right)
\frac{k^0_{\mathrm{on}}}{k^0_{\mathrm{off}}}.
\]

This ensures that the correct thermodynamic stall force is obtained. The paper then adopts a minimal rate prescription in which the instantaneous rates do not exceed their free-filament values.

The current repository does **not** appear to contain a complete standalone snapshot of the single-filament parameter-sweep program used to generate every figure in this paper. Instead, the thermodynamic framework from this work is carried into the later multifilament source.

### Where it appears in `2Dmain.cpp`

The historical C++ program evaluates the interaction potential at the current gap and at gaps shifted by one projected monomer step, then updates `k_on` and `k_off` using exponential energy differences. It caps each rate at the free value when the corresponding energy change is favorable. This is the same rate construction used in the publications.

Because the code uses shifted potentials and stores force-like coefficients in its `A_*` and `B_*` variables, its parameter symbols should not be assumed to be numerically identical to the paper’s \(A\) and \(B\) energy coefficients without tracing the definitions.

## 2. Multifilament pulling-force model

**Paper:** F. Motahari and A. E. Carlsson, *Physical Biology* **17**, 016005 (2020)  
DOI: https://doi.org/10.1088/1478-3975/ab59bd

The public `2Dmain.cpp` source maps directly onto the architecture described in this paper:

- 144 filaments arranged as a 12 × 12 square array;
- 36 central filaments in a 6 × 6 puller region;
- 108 surrounding pusher filaments;
- filament-obstacle potentials that differ between the two populations;
- stochastic polymerization/depolymerization;
- obstacle Brownian motion;
- filament-tip bending fluctuations;
- elastic base displacement coupled to neighboring filaments;
- per-filament forces and force distributions;
- MPI-based parameter sweeps.

The paper’s baseline parameter table includes an actin step size of 2.7 nm before projection, 5.3 µM bulk actin, free on/off constants of 11.6 s⁻¹ µM⁻¹ and 1.4 s⁻¹, obstacle diffusion coefficient \(10^4\) nm²/s, tip and base diffusion coefficients \(10^5\) nm²/s, tip bending spring constant 4.17 pN/nm, medium gel spring constant 0.53 pN/nm, and timestep \(5\times10^{-10}\) s.

The source uses `delta = 2.21` nm, which is the approximately projected 2.7 nm actin step for a 35° filament angle.

## 3. Output-to-analysis map

The following is a best-effort map from the checked-in source to the published analyses. It is intentionally conservative: exact one-command reproduction of each journal figure is **not** claimed.

| Publication analysis | Relevant source/output | Status |
|---|---|---|
| Pulling-force buildup versus time (paper Fig. 4) | `F_pulling_vs_time.txt`, `F_total_vs_time.txt` | Directly compatible with the analysis; exact plotting script not present. |
| Gel deformation versus time (Fig. 5) | `Gel_deformation_avg_in_time.txt`, `dz_elas_*_vs_time.txt` | Directly compatible; exact plotting script not present. |
| Effect of pusher attraction on force (Fig. 6) | MPI parameter table varying `A_push`, `B_push`; summary in `Data.txt` | Parameter sweep is represented; aggregation/plotting script absent. |
| Spatial force distribution (Figs. 7–8) | `Force_Dist.txt`, `Force_Dist_Symmetrized.txt`, `Force_Dist_Row_3.txt`, `Force_Dist_Row_5.txt` | Direct output. Original `3D_Force_Dist.py` is a diagnostic visualizer; modern helper also supplied. |
| Time development of force distribution (Fig. 9) | per-filament `Force_vs_time_<row>_<col>.txt` | Underlying time-resolved output is present; exact figure script absent. |
| Mean-force F–V analysis (Fig. 10) | external-load parameter `F_ext`; `F_ext_vs_v_memb.txt` | Source contains the required external-load machinery; exact sweep/aggregation workflow not fully archived. |
| Pulling force vs puller growth rate (Fig. 11) | pulling-force summaries and calculated puller polymerization speed in `Data.txt` | Quantities are produced; exact aggregation script absent. |
| Polymerized subunit count (Fig. 12) | per-filament `Num_subunits` summaries in `Data.txt` | Underlying counts are produced; exact plotting script absent. |
| Detachment analyses (Figs. 13–14) | `t_break`, per-filament gap files `r_<row>_<col>.txt`, force time series | Relevant state variables are archived by the simulation; the journal-specific heatmap workflow is not present. |

## 4. Historical plotting and visualization files

### `3D_Force_Dist.py`

This script reads `Force_Dist.txt` and plots the spatial force distribution. It uses `matplotlib.mlab.griddata`, an API removed from modern Matplotlib. It is retained unchanged as historical code.

`plot_force_distribution.py` is a modern replacement for this one plotting task. It does not change the original simulation.

### `OpenGL_Visualization.cpp`

This program reads time-resolved per-filament force files and draws an OpenGL visualization of the 12 × 12 array. It depends on GLEW, GLFW, a local `Shader.h`, shader source files under `res/shaders/`, and run data under `res/`. Those support files are not in the current public repository, so this source should be regarded as an archival visualization component rather than a standalone build target.

## 5. Scientific results represented by the code

The ensemble paper found that the faster-growing outer filaments push while the more strongly bound/slower-growing central filaments pull. Pulling force increases as central filament growth is suppressed, reaching its maximum when the central pullers effectively do not grow. The steady-state force is comparatively insensitive to gel rigidity, but softer gels equilibrate more slowly and can detach from the obstacle. The force distribution becomes approximately flat in the pushing and pulling regions at steady state.

These are model results under the assumptions described in the papers; this repository should not be read as a general-purpose cellular endocytosis simulator.
