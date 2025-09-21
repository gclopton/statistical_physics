


| Week | Kardar |
| ---- | ------ |
| 1    | 1      |
| 2    |        |
| 3    |        |
| 4    |        |
| 5    |        |
| 6    |        |
| 7    |        |
| 8    |        |
| 9    |        |
| 10   |        |
| 11   |        |
| 12   |        |
| 13   |        |
| 14   |        |




# Simulations

Absolutely—here’s a curated, **code‑first** set of Python simulations matched to every topic in your syllabus schedule. I focused on resources that include **complete, runnable Python scripts or notebooks** (many use open scientific Python stacks like **HOOMD‑blue, OpenMM, ASE, MDAnalysis, NumPy/SciPy**). I mapped rows to the schedule items in your syllabus for **MSE 485 – Atomic Scale Simulations**.

> **How to read this:** The middle column links directly to the code (repo or notebook). Where a topic is broad, I picked one strong exemplar that you can drop into class with minimal prep; I also note what the example demonstrates or which package it uses.

---

### Python simulation examples aligned to your syllabus

|Syllabus topic|Python simulation with full solution|Why this fits / what it shows|
|---|---|---|
|**1. Overview**|**LJ Fluid Lab (Notebook + code)** — self‑contained Lennard‑Jones MD in Python (Colab + repo). ([Google Colab](https://colab.research.google.com/github/ggrazioli/LennardJonesFluidLab/blob/main/LennardJonesFluidLab_Colab_public.ipynb?utm_source=chatgpt.com "Lennard-Jones Fluid Simulation Lab - Colab"), [GitHub](https://github.com/ggrazioli/LennardJonesFluidLab?utm_source=chatgpt.com "ggrazioli/LennardJonesFluidLab: Lennard-Jones fluid simulation ..."))|A gentle but real MD workflow (init → integrate → analyze) suitable for a first live demo (positions, velocities, energy, RDF).|
|**2. Molecular Models**|||
|Periodic Boundary Conditions|**Minimal LJ MD (from scratch)** — PBC + minimum‑image + LJ forces + Verlet, pure Python. ([GitHub](https://github.com/m-rivera/ljmd?utm_source=chatgpt.com "m-rivera/ljmd: Lennard-Jones Molecular Dynamics for ..."))|Compact code where students can see and edit the PBC/minimum‑image logic.|
|Intramolecular Potentials|**OpenMM “first script”** — bonds/angles/dihedrals via force fields; full Python driver. ([docs.openmm.org](https://docs.openmm.org/latest/userguide/application/02_running_sims.html?utm_source=chatgpt.com "3. Running Simulations"))|Shows how bonded terms enter a real engine; easy to swap in simple molecules or polymers.|
|Short‑ranged Potentials|**HOOMD‑blue LJ tutorial (Python)** — define pair.LJ, integrate MD. ([hoomd-blue.readthedocs.io](https://hoomd-blue.readthedocs.io/en/v3.11.0/tutorial/01-Introducing-Molecular-Dynamics/01-Molecular-Dynamics-Simulations.html?utm_source=chatgpt.com "Molecular Dynamics Simulations - HOOMD-blue"))|Clean Python API for LJ (ε, σ, r_cut), neighbor lists, timestepping.|
|Electrostatics / Ewald Summation|**DScribe Ewald Sum Matrix (Python)** — compute Ewald‑summed energies for periodic charges. (Alt: `pymatgen.analysis.ewald.EwaldSummation`.) ([Singroup](https://singroup.github.io/dscribe/1.0.x/tutorials/descriptors/ewald_sum_matrix.html?utm_source=chatgpt.com "Ewald sum matrix — DScribe 1.0.x documentation"), [PyPI](https://pypi.org/project/pymatgen-analysis-defects/?utm_source=chatgpt.com "pymatgen-analysis-defects"))|Direct, reproducible Ewald implementation callable from Python; great for comparing real‑space vs reciprocal‑space terms on toy crystals.|
|**3. Statmech & Thermo**|||
|Statistics / Python|**Stat Mech course Python primer** — Jupyter basics + NumPy/Matplotlib workflow. ([nthu-yiping-huang.github.io](https://nthu-yiping-huang.github.io/Statistical_Mechanics_I_2022_Spring/Python_tutorial_0.html?utm_source=chatgpt.com "Basic introduction of python - 2022 Statistical Mechanics (I)"))|Fast on‑ramp for students new to scientific Python.|
|Statistical Mechanics|**Maxwell–Boltzmann in 2D (NumPy)** — simulate colliding particles → MB speed distribution. ([scipython.com](https://scipython.com/blog/the-maxwellboltzmann-distribution-in-two-dimensions/?utm_source=chatgpt.com "The Maxwell–Boltzmann distribution in two dimensions"))|Concrete link from microdynamics to equilibrium distributions.|
|Thermodynamics, Properties|**MD in Python for beginners** — velocity‑Verlet + virial pressure + thermodynamic properties (open code). ([Preprints](https://www.preprints.org/manuscript/202012.0179/v1/download?utm_source=chatgpt.com "[PDF] Practical Implementation of Molecular Dynamics code for beginners ..."))|A from‑scratch MD showing pressure/temperature/virial calculation.|
|Structural Correlations|**MDAnalysis: InterRDF** — compute g(r) from a trajectory in Python. ([MDAnalysis Documentation](https://docs.mdanalysis.org/stable/documentation_pages/analysis/rdf.html?utm_source=chatgpt.com "4.8.2.1. Radial Distribution Functions — MDAnalysis.analysis.rdf"))|Clear, tested API for RDF; drop‑in on HOOMD/OpenMM/LAMMPS trajectories.|
|Transport Coefficients|**MDAnalysis: MSD & diffusion** (Einstein relation) — full example + code. (Alt for viscosity: `visco.py` Green–Kubo script.) ([MDAnalysis Documentation](https://docs.mdanalysis.org/2.7.0/documentation_pages/analysis/msd.html?utm_source=chatgpt.com "4.7.2.2. Mean Squared Displacement"), [GitHub](https://github.com/omidshy/aMD?utm_source=chatgpt.com "omidshy/aMD: A collection of Python codes to calculate ..."))|One‑liner MSD + linear fit for D; optional Green–Kubo viscosity script for stress ACF.|
|**4. Molecular Dynamics**|||
|Classical Mechanics|**NVE MD notebook** — forces, integration, energies; a clean “first MD” in Python. ([GitHub](https://github.com/JayLau123/Molecular-dynamics-simulation-NVE-?utm_source=chatgpt.com "JayLau123/Molecular-dynamics-simulation-NVE-: Tutorial"))|Uses dimensionless units; great for linking Newton → MD.|
|Verlet Integration|**LJ‑Argon MD (Python)** — velocity‑Verlet on LJ liquid; run as `python simulate.py`. ([GitHub](https://github.com/KenNewcomb/LJ-Argon?utm_source=chatgpt.com "KenNewcomb/LJ-Argon: A molecular dynamics simulation ..."))|Full script: init → neighbor calc → force → VV integrate → observables.|
|Thermostats|**OpenMM Langevin/NVT in Python** — set temperature, friction; run canonical MD. ([docs.openmm.org](https://docs.openmm.org/latest/userguide/application/02_running_sims.html?utm_source=chatgpt.com "3. Running Simulations"))|Shows thermostatting in a few lines; switch between NVE/NVT easily.|
|Langevin Dynamics|**OU / Langevin notebook** — simulate the Langevin equation; analyze stats. ([IPython Books](https://ipython-books.github.io/134-simulating-a-stochastic-differential-equation/?utm_source=chatgpt.com "13.4. Simulating a stochastic differential equation"))|Pedagogical SDE implementation connecting stochastic force ↔ diffusion.|
|Constraints for MD|**OpenMM constraints (Python API)** — `System.addConstraint`, tolerance, etc. ([docs.openmm.org](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.System.html?utm_source=chatgpt.com "System — OpenMM Python API 8.3.0. ..."))|Demonstrates SHAKE/SETTLE‑style constraints at the Python level.|
|MD Packages|**HOOMD‑blue MD tutorial** — end‑to‑end MD entirely from Python. (Alt: ASE calculators overview.) ([hoomd-blue.readthedocs.io](https://hoomd-blue.readthedocs.io/en/v3.11.0/tutorial/01-Introducing-Molecular-Dynamics/01-Molecular-Dynamics-Simulations.html?utm_source=chatgpt.com "Molecular Dynamics Simulations - HOOMD-blue"), [wiki.fysik.dtu.dk](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html?utm_source=chatgpt.com "Calculators — ASE documentation"))|Modern, fast, Python‑native engine; simple to teach and to scale up.|
|Initialization / Neighbor search|**Neighbor lists in Python** — build cell/Verlet lists, complexity & updates. (Alt: HOOMD nlist docs.) ([Python in Plain English](https://python.plainenglish.io/molecular-dynamics-neighbor-lists-in-python-a8351e6fa3a8?utm_source=chatgpt.com "Molecular Dynamics: Neighbor Lists in Python"), [hoomd-blue.readthedocs.io](https://hoomd-blue.readthedocs.io/en/v2.9.6/nlist.html?utm_source=chatgpt.com "Neighbor lists — HOOMD-blue 2.9.6 documentation"))|Code you can annotate in lecture to show O(N)→O(N·const) scaling.|
|HPC, Discussion of MD|**HOOMD: MPI domain decomposition** — parallel MD from the same Python script. ([hoomd-blue.readthedocs.io](https://hoomd-blue.readthedocs.io/en/v3.1.0/tutorial/03-Parallel-Simulations-With-MPI/02-Domain-Decomposition.html?utm_source=chatgpt.com "Domain Decomposition — HOOMD-blue 3.1.0 documentation"))|Demonstrates scaling out with `mpirun` while keeping Python driver.|
|**5. Monte Carlo**|||
|Random Numbers|**Monte Carlo π (Jupyter, NumPy)** — RNG seeding, convergence, error ~1/√N. ([Risk Engineering](https://risk-engineering.org/notebook/monte-carlo-pi.html?utm_source=chatgpt.com "Monte Carlo methods: basic introduction - Risk Engineering"))|Compact RNG showcase; good to introduce reproducibility & variance.|
|Monte Carlo|**LJ Monte Carlo (NVT) in Python** — MC for LJ fluid; reproduces NIST benchmarks. ([GitHub](https://github.com/KenNewcomb/LJ-MC-NVT?utm_source=chatgpt.com "KenNewcomb/LJ-MC-NVT: A Monte Carlo simulation of the Lennard ..."))|Straight Metropolis sampler in continuous space; clear, short code.|
|Metropolis|**Ising model (Metropolis) notebook** — 1D/2D Ising, E/M/C_V vs T. ([hockygroup.com](https://hockygroup.com/exercise/ising-1d.html?utm_source=chatgpt.com "Computational HW 6, Monte Carlo, 1-d ising model"))|Classic Metropolis acceptance & detailed balance; easy to extend.|
|NVT / GC Monte Carlo|**Grand‑canonical lattice MC (Python)** — μVT for lattice gas/Ising. ([GitHub](https://github.com/bjmorgan/lattice_mc?utm_source=chatgpt.com "bjmorgan/lattice_mc: Lattice gas Monte Carlo simulation ..."))|A minimal GC sampler you can adapt to chemical potential demos.|
|Phase Coexistence|**Gibbs ensemble MC (notebooks; Colab)** — LJ liquid–vapor coexistence via GEMC. ([summer.ccp5.ac.uk](https://summer.ccp5.ac.uk/WORKSHOP/Day_5/Phase_Equilibria/MC-Gibbs.html?utm_source=chatgpt.com "Gibbs ensemble simulations - CCP5 Summer School"))|Directly targets coexistence curve with particle/volume exchanges.|
|Finite Size Effects|**Diffusion finite‑size correction (Python)** — D vs 1/L using multiple supercells. ([SCM](https://www.scm.com/doc/Tutorials/MolecularDynamicsAndMonteCarlo/MDintroPython/intro.html?utm_source=chatgpt.com "Molecular Dynamics with Python - SCM"))|Crisp demonstration that transport properties depend on box size.|

> **Schedule source:** The rows correspond 1:1 to the schedule in your syllabus (Overview; Molecular Models → PBC, Intramolecular, Short‑range, Ewald; Statmech & Thermo → Statistics/Python, Stat Mech, Thermodynamics/Properties, Structural Correlations, Transport; Molecular Dynamics → Classical Mechanics, Verlet, Thermostats, Langevin, Constraints, MD Packages, Initialization/Neighbor search, HPC; Monte Carlo → Random Numbers, Monte Carlo, Metropolis, NVT/GC, Phase Coexistence, Finite Size Effects).

---

## Notes & quick tips

- **All links include full Python code** (scripts or notebooks). Some examples call established engines (e.g., **HOOMD‑blue, OpenMM, ASE**), but the **driver and analysis are pure Python**, which is ideal for teaching and reproducibility. ([hoomd-blue.readthedocs.io](https://hoomd-blue.readthedocs.io/en/v3.11.0/tutorial/01-Introducing-Molecular-Dynamics/01-Molecular-Dynamics-Simulations.html?utm_source=chatgpt.com "Molecular Dynamics Simulations - HOOMD-blue"), [docs.openmm.org](https://docs.openmm.org/latest/userguide/application/02_running_sims.html?utm_source=chatgpt.com "3. Running Simulations"), [wiki.fysik.dtu.dk](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html?utm_source=chatgpt.com "Calculators — ASE documentation"))
    
- If you want a **single stack** for the MD track, HOOMD‑blue + MDAnalysis covers almost everything above (LJ MD, thermostats including Langevin/NVT, neighbor lists, parallel runs, RDF/MSD analysis). ([hoomd-blue.readthedocs.io](https://hoomd-blue.readthedocs.io/en/v3.11.0/tutorial/01-Introducing-Molecular-Dynamics/01-Molecular-Dynamics-Simulations.html?utm_source=chatgpt.com "Molecular Dynamics Simulations - HOOMD-blue"), [MDAnalysis Documentation](https://docs.mdanalysis.org/stable/documentation_pages/analysis/rdf.html?utm_source=chatgpt.com "4.8.2.1. Radial Distribution Functions — MDAnalysis.analysis.rdf"))
    
- For **electrostatics**, using **DScribe** or **pymatgen** gives you a compact, inspectable Ewald reference you can compare against MD package energies. ([Singroup](https://singroup.github.io/dscribe/1.0.x/tutorials/descriptors/ewald_sum_matrix.html?utm_source=chatgpt.com "Ewald sum matrix — DScribe 1.0.x documentation"), [PyPI](https://pypi.org/project/pymatgen-analysis-defects/?utm_source=chatgpt.com "pymatgen-analysis-defects"))
    
- For **Monte Carlo**, a neat progression is π‑estimation (random numbers) → Ising (Metropolis) → LJ MC (continuous) → **Gibbs ensemble** (coexistence). ([Risk Engineering](https://risk-engineering.org/notebook/monte-carlo-pi.html?utm_source=chatgpt.com "Monte Carlo methods: basic introduction - Risk Engineering"), [hockygroup.com](https://hockygroup.com/exercise/ising-1d.html?utm_source=chatgpt.com "Computational HW 6, Monte Carlo, 1-d ising model"), [GitHub](https://github.com/KenNewcomb/LJ-MC-NVT?utm_source=chatgpt.com "KenNewcomb/LJ-MC-NVT: A Monte Carlo simulation of the Lennard ..."), [summer.ccp5.ac.uk](https://summer.ccp5.ac.uk/WORKSHOP/Day_5/Phase_Equilibria/MC-Gibbs.html?utm_source=chatgpt.com "Gibbs ensemble simulations - CCP5 Summer School"))
    

---

### Want this as a living class resource?

I can export this into a CSV/Google Sheet or mirror the exact notebooks into a **single GitHub classroom repo** with a uniform environment file (`conda`/`pip`) so students can run each topic locally or on Colab. If you prefer one toolkit (e.g., only HOOMD + MDAnalysis), I can also replace mixed examples with equivalents in that stack.