# Endocytosis — Actin Filament–Membrane Interaction Simulations

Historical research code for stochastic simulations of actin polymerization, filament–membrane interactions, and pulling-force generation during endocytosis.

This repository originated in my Ph.D. research in computational biophysics at Washington University in St. Louis. The code accompanies two first-author papers with A. E. Carlsson:

1. **F. Motahari and A. E. Carlsson**, “Thermodynamically Consistent Treatment of the Growth of a Biopolymer in the Presence of a Smooth Obstacle Interaction Potential,” *Physical Review E* **100**, 042409 (2019).  
   https://doi.org/10.1103/PhysRevE.100.042409

2. **F. Motahari and A. E. Carlsson**, “Pulling-Force Generation by Ensembles of Polymerizing Actin Filaments,” *Physical Biology* **17**, 016005 (2020).  
   https://doi.org/10.1088/1478-3975/ab59bd

## Scientific overview

Actin polymerization can generate both pushing and pulling forces on cellular membranes. The work represented here asks how those forces depend on the microscopic interaction between a growing filament tip and a membrane/obstacle, and how spatially distinct populations of filaments can cooperate to generate a net pulling region during endocytosis.

The central multifilament model uses:

- a **12 × 12 square array** of actin filaments;
- a **6 × 6 central population** of more strongly membrane-bound, slower-growing “puller” filaments;
- **108 surrounding “pusher” filaments** that grow more rapidly;
- stochastic polymerization and depolymerization;
- smooth filament–obstacle interaction potentials with repulsive and/or attractive components;
- thermodynamically constrained position-dependent polymerization and depolymerization rates;
- biased Brownian dynamics for obstacle motion and filament-tip fluctuations;
- stochastic filament-base motion with nearest-neighbor elastic coupling to represent actin-gel deformation.

In the ensemble simulations, the outer filaments push while the central filaments pull. Stronger central binding slows central growth and can increase the total pulling force; softer gels take longer to reach steady force and can become more susceptible to filament–membrane detachment.

## Relationship between the two papers

The *Physical Review E* paper establishes the thermodynamic constraint linking the smooth interaction potential \(U(r)\) to the instantaneous polymerization and depolymerization rates. In compact form,

\[
\frac{k_{\mathrm{on}}(r)}
{k_{\mathrm{off}}(r-\delta)}
=
\exp\left(
-\frac{U(r-\delta)-U(r)}{k_B T}
\right)
\frac{k^0_{\mathrm{on}}}{k^0_{\mathrm{off}}}.
\]

The later *Physical Biology* study uses this framework in a many-filament model of pulling-force generation. The public `2Dmain.cpp` snapshot is primarily associated with that ensemble study and implements the corresponding position-dependent rate prescription.

See [`docs/CODE_AND_PAPERS.md`](docs/CODE_AND_PAPERS.md) for a more detailed paper-to-code map.

## Repository contents

| File | Purpose |
|---|---|
| `2Dmain.cpp` | Main historical MPI/C++ stochastic multifilament simulation. |
| `3D_Force_Dist.py` | Original plotting helper for a force-distribution output file; preserved as historical code. |
| `OpenGL_Visualization.cpp` | Historical OpenGL visualization program for time-resolved filament-force output. |
| `plot_force_distribution.py` | Modern lightweight plotting helper added during archival documentation. |
| `docs/CODE_AND_PAPERS.md` | Maps the repository to the two publications and their analyses. |
| `docs/REPRODUCIBILITY.md` | Build notes, outputs, limitations, and reproducibility status. |
| `CITATION.cff` | Software citation metadata. |
| `requirements.txt` | Python requirements for the modern plotting helper. |

## Building the main simulation

The main source uses MPI and C++11 features. A typical compile command is:

```bash
mpic++ -O2 -std=c++11 2Dmain.cpp -o endocytosis_sim
```

**Important:** the checked-in source is a historical production-scale research program, not a short example. Its default configuration uses \(2\times10^{10}\) timesteps, corresponding to a 10 s simulated trajectory at \(\Delta t=5\times10^{-10}\) s. Do not start a full run casually.

The current parameter table explicitly defines parameter sets for MPI ranks **0 through 13**. See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) before attempting a run.

## Plotting a force distribution

`2Dmain.cpp` writes `Force_Dist.txt` in each run directory. The original `3D_Force_Dist.py` depends on an old Matplotlib API, so a modern helper is included:

```bash
python -m pip install -r requirements.txt

python plot_force_distribution.py \
    --input Run00/Force_Dist.txt \
    --output Run00/force_distribution.png
```

The script accepts the four-column `Force_Dist.txt` format produced by the C++ program:

```text
row  column  mean_force  uncertainty
```

## Reproducibility status

This is an **archival research-code repository**. The original simulation files are preserved rather than silently rewritten.

The public snapshot captures the principal multifilament simulation and historical visualization tools, but it is not a frozen turnkey reproduction package for every figure in the two papers. Some publication plots were generated from parameter sweeps and analysis steps whose exact plotting/aggregation scripts are not present in this repository. The OpenGL program also refers to historical shader-support files that are not currently included.

These limitations are documented explicitly in [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).

## Citation

If you use or discuss this code, please cite the repository and the relevant publication(s). Machine-readable citation metadata are provided in [`CITATION.cff`](CITATION.cff).

## Author

**Fowad Motahari**  
Ph.D. research performed in the Department of Physics and Center for Engineering Mechanobiology, Washington University in St. Louis.

## Archival note

Original scientific source files in this repository date from the Ph.D. project. Documentation and the modern plotting helper were added in 2026 to make the historical research code easier to understand and inspect without changing the scientific provenance of the original simulation.
