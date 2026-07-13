# Paramo

> "Vine a Comala porque me dijeron que acá vivía mi padre, un tal Pedro Páramo." - Juan Rulfo

`paramo` stands for: PArticles and RAdiation MOnitor. `paramo` is a Fortran 95 simulation library. It solves the Fokker-Planck equation to model the time evolution of relativistic particle distributions in relativistic plasma, coupled to synchrotron and inverse Compton radiative transfer. It is also parallelized with OpenMP, and uses HDF5 for output.

The Fortran core is organized into modular source files separating physical processes, numerical methods, and I/O routines.

[`Arriero.py`](./Arriero.py) is the user Python layer that provides a class-based interface for parameter management, compilation orchestration across GCC and Intel compilers, and HPC job submission supporting both PBS and SLURM schedulers. This separation between the performance-critical Fortran core and the Python workflow layer is a design pattern I carried forward into Tleco, its Rust and Python successor described in [Davis, Rueda-Becerril & Giannios (2022)](https://academic.oup.com/mnras/article/513/4/5766/6583293).

## Prerequisites

- Fortran (GNU/Intel)
- HDF5

## How to run it

1. **Clone the repository**

2. **Customize `Makefile`**. The `Makefile` has a series of variables that allows to configure compilation.

    - `COMPILER`: 0 for GCC and 1 for Intel compilers.
    - `DEBUGGING`: Compile in debugging or optimized mode.
    - `USEHDF5`: Compile with HDF5 libraries and save data in that format.

3. **Create parameters file**

4. **Run `make` in terminal**

5. **Run executable**

The alternative is calling [`Arriero.py`](./Arriero.py).

## Output

Output files can be analyzed using [SAPytho](https://github.com/altjerue/SAPytho), a companion Python toolkit for computing spectra, light curves, and SED fits from simulation data.

## References

These are some of the references used to build this code:

- [PVP14] Pannanen, Vurm, Poutanen, 2014, A&A, 564, A77
- [PM09]  Petropoulou, Mastichiadis, 2009, A&A, 507, 599
- [RM92]  Rees & Meszaros, 1992, MNRAS, 258, 41P
- [DM09]  Dermer & Menon, 2009, "High Energy Radiation from Black Holes: Gamma Rays, Cosmic Rays, and Neutrinos", Princeton Series in Astrophysics
- [SPN98] Sari, Piran, Narayan, 1998, ApJ, 497, L17
- [PK00]  Panaitescu & Kumar, 2000, ApJ, 543, 66
- [CL00]  Chavalier & Li, 2000, ApJ, 536, 195
- [CD99]  Chiang & Dermer, 1999, ApJ, 512, 699

Publications using `paramo`:

- [Rueda-Becerril, Harrison & Giannios (2021)](https://academic.oup.com/mnras/article/501/3/4092/6043215)
- [Davis, Rueda-Becerril & Giannios (2022)](https://academic.oup.com/mnras/article/513/4/5766/6583293)
- [Combi & Siegel (2023)](https://iopscience.iop.org/article/10.3847/1538-4357/acac29/meta)
