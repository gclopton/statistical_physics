

Ceria ( $\mathrm{CeO}_2$ ) is a standard surrogate for fluorite nuclear fuels and a common target in swift heavy-ion (SHI) studies. In fuels, the lattice thermal conductivity $k$ sets centerline temperatures at a given power and feeds fission-gas, melting-margin, and pellet-clad mechanics models. In SHI physics, $k$ controls thermal-spike dissipation and, conversely, is altered by ion-track damage. Oxygen vacancies are ubiquitous in reduced ceria and scatter phonons, so $k$ should depend on temperature $T$ and vacancy fraction $x$. Our goal is to compute a small, uncertainty-controlled map $k(T, x)$ for $\mathrm{CeO}_{2-x}$ via classical MD that is immediately usable in fuel-performance and irradiation models.

We will implement a reproducible LAMMPS workflow. For each ( $T, x$ ), we build a fluorite supercell, remove $x N_{\mathrm{Ce}}$ oxygen sites, and for each vacancy reduce two neighboring $\mathrm{Ce}^{4+}$ to $\mathrm{Ce}^{3+}$ (explicit species and charges) to maintain neutrality, then relax (NPT-NVT) and run NVE production. Vacancy sites are chosen with a minimum-distance criterion to avoid artificial clustering; $\mathrm{Ce}^{3+}$ assignments follow nearestneighbor logic. We target $10 \times 10 \times 10$ and larger supercells so $x N_{\mathrm{Ce}}$ is integral and GK size effects are small, and we run multiple seeds to control statistical error.

Thermal conductivity will be computed primarily by Green-Kubo. We will record the microscopic heat flux $\mathbf{J}(t)$ using compute heat/flux with per-atom energies and per-atom stress including the PPPM kspace term, then build the heat-current autocorrelation and integrate to a plateau with a fixed analysis protocol (blocking or Sokal window) to report $\kappa \pm$ confidence intervals. Because our model is purely pairwise (Buckingham + Coulomb) with no angles/dihedrals, the recommended stress/atom route is consistent with long-range electrostatics; we will not switch to screened Coulombics unless strictly necessary for stability. As a cross-check, we will perform a single Müller-Plathe rNEMD calculation on the validated stoichiometric case and confirm consistent $k$.

We will validate on perfect $\mathrm{CeO}_2$ before introducing vacancies: reproduce a reference-scale $k$ at one $T$, and match the potential's lattice parameter and density. Convergence tests will cover timestep stability, thermostat/barostat settings that preserve fluctuations, GK plateau stability with increasing trajectory length and origin count, and a finite-size check by extending one box dimension. Literature ranges for bulk $\mathrm{CeO}_2$ (defect-free theory in the mid-teens $\mathrm{W} \mathrm{m}^{-1} \mathrm{~K}^{-1}$ at 300 K ; GK EMD high-teens as size grows; porous pellets lower) will be used to anchor expectations. 

Upon completion of the project, we will have: 
1.) a supercell + vacancy/ $\mathrm{Ce}^{3+}$ generator with reproducible seeds
2.) a driver that stages NPT $\rightarrow$ NVT $\rightarrow$ NVE and logs the quantities needed for GK/rNEMD and 
3.) a post-processor that outputs $k$ with documented CIs from a fixed protocol

We will produce a $k(T, x)$ data set showing the reduction of $k$ with increasing $x$ and its temperature dependence, suitable for centerline-temperature and thermal-spike models, and extensible to added $x$, dopants, or defect structures without changing the method.


# Sources

1.) [Nanometric Transformation of the Matter by Short and Intense Excitation: Experimental Data Versus Inelastic Thermal Spike Model](https://www.sciencedirect.com/science/article/pii/S0168583X11011645?utm_source=chatgpt.com)
2.) [Compute Stress/Atom Command–LAMMPS Documentation](https://docs.lammps.org/compute_stress_atom.html)
3.) [Atomistic and Experimental Study on Thermal Conductivity of Bulk and Porous Cerium Dioxide](https://pmc.ncbi.nlm.nih.gov/articles/PMC6474893/)
