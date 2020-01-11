# Paramo

"Vine a Comala porque me dijeron que acá vivía mi padre, un tal Pedro Páramo."

## What it does

Paramo stands for: PArticles and RAdiation MOnitor. In a few words: this code solves que Fokker-Planck equation and calculates synchrotron and inverse Compton emission.

## Prerequisites

- Fortran (either GNU or Intel)
- HDF5
- OpenMP (optional)
- OpenACC (in progress)

### Python interface

- time
- [SAPyto](https://github.com/altjerue/SAPyto)
- [extractor](https://github.com/altjerue/extractor)

## Parameters

Name | Default Value | Concept
-----|---------------|--------
`R` | 1E16 (cm) | Radius of the emitting region
`R0` | 1E15 (cm)  | Distance from explosion
`dLum` | 4.0793E26 (cm)  | luminosity distance (Mrk 421)
`z` | 0.03 | redshift (Mrk 421)
`theta_obs` | 5.0 (deg) | observer viewing angle
`gamma_bulk` | 100 | bulk Lorentz factor of the emitting region
`mu_mag` | 1.0 | baryon loading ($\mu$)
`sigma` | 1.0 | magnetization
`f_rec` | 1.0 | magnetic reconectiondissipative efficiency
`b_index` | 0.0 | magnetic field decay index
`Bfield` | 1.0 (G) | Magnetic field strength
`eps_B` | 0.03 | $epsilon_B$
`eps_e` | 0.1 | $epsilon_e$
`theta_e` | 10.0 | electrons temperature (k T / m_e c^2)
`zeta_e` | 0.99 | fraction of non-thermal particles
`tstep` | 0.01 | time-step factor
`tmax` | 1E5 (s) | maximum time
`tmin` | 0.0 (s) | minimum time
`tvar` | 2.0 | variability time scale
`L_j` | 1E45 (erg / s) | jet luminosity
`E0` | 1E50 (erg) | isotropic energy of the blast wave
`n_ext` | 1.0 | number density of the external medium
`g1` | 100 | minimum Lorentz factor of non-thermal particles
`g2` | 1E4  | maximum Lorentz factor of non-thermal particles
`gmin` | 1.01 | minimum Lorentz factor of particles distribution
`gmax` | 2E$ | maximum Lorentz factor of particles distribution
`qind` | 2.5 | power law index of non-thermal particles
`nu_ext` | 1E14 | frequency of external radiation field
`u_ext` | 1E-4 | energy density of external radiation field
`numin` | 1E7 | minimum frequency
`numax` | 1E15 | maximum frequency
`numbins` | 128 | number of bin in the particles distribution
`numdt` | 300 | number of time-steps
`numdf` | 256 | number of frequencies
`time_grid` | 1 | type of time discretization
`params_file` | input.par | name of the parameters file

## How to run it

### Compilation

The exectuables are generated through `Makefile`. To generate all executables available you have to run
```
$ make all
```

To generate the blazar program:
```
$ make xBlazMag
```

To generate the afterglow program:
```
$ make xAglow
```

### Compiler flags

### Run it

### Run it with `Python`

### Run it in a server (e.g., Brown)

## How to red output

# References

These are the most referenced works on which I based all the modeling

## Blazars module

## GRB afterglow module
- [PVP14] Pannanen, Vurm, Poutanen, 2014, A&A, 564, A77
- [PM04]  Petropoulou, Mastichiadis, 2009, A&A, 507, 599
- [RM92]  Rees & Meszaros, 1992, MNRAS, 258, 41P
- [DM09]  Dermer & Menon, 2009, "High Energy Radiation from Black Holes: Gamma Rays, Cosmic Rays, and Neutrinos", Princeton Series in Astrophysics
- [SPN98] Sari, Piran, Narayan, 1998, ApJ, 497, L17
- [PK00]  Panaitescu & Kumar, 2000, ApJ, 543, 66
- [CL00]  Chavalier & Li, 2000, ApJ, 536, 195
- [CD99]  Chiang & Dermer, 1999, ApJ, 512, 699
