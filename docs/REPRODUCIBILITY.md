# Reproducibility and Build Notes

## Status

This repository is an archival snapshot of Ph.D.-era research code. It is being documented rather than rewritten so that the provenance of the original calculations remains visible.

The source is scientifically meaningful, but the repository is **not** currently a complete turnkey reproduction package for every figure in the associated papers.

## Main C++ simulation

`2Dmain.cpp` requires an MPI-capable C++ compiler and C++11 support.

Typical compile command:

```bash
mpic++ -O2 -std=c++11 2Dmain.cpp -o endocytosis_sim
```

The source currently defines parameter sets for MPI ranks **0 through 13**. Although the internal `MAXPROC` constant is larger, ranks beyond the explicitly populated parameter sets should not be launched without first extending the parameter table.

### Computational scale

The checked-in source uses:

```text
Ntstep = 20 × 10^9
dt     = 5 × 10^-10 s
```

so the default production trajectory corresponds to 10 s of simulated time and is computationally expensive.

This is not a recommended “hello world” run. A scientifically valid reduced test would require coordinated changes to timestep-count, equilibration time, averaging windows, and related bookkeeping. Those changes are not supplied here because they would alter the historical production configuration.

## Main output structure

Each MPI rank creates a directory named `Run00`, `Run01`, etc. Outputs include:

- `Parameters.txt`
- `Time_Course.txt`
- `F_pulling_vs_time.txt`
- `F_pushing_vs_time.txt`
- `F_total_vs_time.txt`
- `Gel_deformation_avg_in_time.txt`
- `dz_elas_*_vs_time.txt`
- `k_on.txt`
- `k_off.txt`
- `Force_Dist.txt`
- `Force_Dist_Symmetrized.txt`
- `Force_Dist_Row_3.txt`
- `Force_Dist_Row_5.txt`
- `Gel_Deformation_Dist.txt`
- `F_ext_vs_v_memb.txt`
- `Data.txt`
- per-filament time series such as `Force_vs_time_<row>_<col>.txt` and `r_<row>_<col>.txt`

`Force_Dist.txt` contains four columns:

```text
row  column  mean_force  estimated_uncertainty
```

## Historical parameter conventions

The publication writes the interaction potential using energy coefficients such as \(A\) and \(B\), with force obtained by differentiating the potential.

The historical C++ source stores coefficients in a shifted form and its `A_push`, `B_push`, `A_pull`, and `B_pull` variables appear directly in the `Force()` expression. Therefore repository variables should not be equated one-for-one with same-named paper symbols without following the code’s definitions.

The code’s `delta = 2.21` nm is approximately the projection of a 2.7 nm actin step along the obstacle-normal direction at the 35° filament angle used in the ensemble paper.

## Python plotting

The original `3D_Force_Dist.py` uses the removed `matplotlib.mlab.griddata` API.

A modern helper is included:

```bash
python -m pip install -r requirements.txt
python plot_force_distribution.py --input Run00/Force_Dist.txt
```

This helper reads the regular 12 × 12 force grid directly and does not modify simulation results.

## OpenGL visualization

`OpenGL_Visualization.cpp` requires at least:

- OpenGL
- GLEW
- GLFW
- a historical local `Shader.h`
- vertex/fragment shader files under `res/shaders/`
- simulation outputs copied under `res/`

The support header and shader files are not currently present in the repository. The visualization source is therefore retained for provenance, but no claim is made that it builds standalone from this repository.

## What is and is not reproducible from this snapshot

### Present

- central many-filament simulation logic;
- thermodynamically constrained position-dependent on/off rates;
- Brownian obstacle, tip, and base dynamics;
- elastic nearest-neighbor coupling;
- MPI parameter sweep;
- force and deformation distributions;
- time-resolved per-filament quantities;
- historical plotting and visualization source.

### Not fully archived

- complete single-filament production code for all figures in the 2019 *Physical Review E* paper;
- exact aggregation/plotting scripts for every published figure;
- original large output files used to produce the journal plots;
- support files for the OpenGL visualization;
- a frozen compiler/library environment from 2018–2019.

## Recommended use of this repository

Use it as:

1. an archival record of the computational methods underlying the publications;
2. a readable example of the stochastic multifilament model;
3. a starting point for understanding the implementation of the position-dependent thermodynamic rate prescription.

Do not treat it as a maintained end-user software library or as a one-command replication capsule for the papers.
