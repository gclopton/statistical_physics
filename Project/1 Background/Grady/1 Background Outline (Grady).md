

# 1 Transport target and modeling scope
- **Lattice thermal conductivity $k$ in solids (Fourier + linear response):** In crystals, $\mathrm{q}=-k \nabla T$ emerges from linear response; MD estimates $k$ either from equilibrium fluctuations (Green-Kubo) or steady non-equilibrium gradients (NEMD). $\rightarrow$ Frames why GK/NEMD are both valid and equivalent in the linear regime.
- **Phonon picture vs. atomistic MD:** $k$ arises from phonon velocities, heat capacities, and scattering (defects, anharmonicity). MD resolves these implicitly without solving BTE. $\rightarrow$ Justifies using classical MD for defect-controlled trends $k(T, x)$.



# 2 Material physics you're encoding
- **Ceria redox chemistry and charge neutrality:** Removing $\mathrm{O}^{2-}$ creates two $\mathrm{Ce}^{3+}$ to keep neutrality; vacancies bind locally and scatter phonons. $\rightarrow$ Explains vacancy counting, $\mathrm{Ce}^{3+}$ assignment, and the expectation that $k \downarrow$ with $x$.
- **Interatomic potentials and transferability:** Buckingham + Coulomb (and possible polarizability) fit elastic constants, lattice parameter, defect energetics. $\rightarrow$ Why you must validate lattice/density and check that the chosen potential was fit for reduced ceria, not just stoichiometric.




# 3 Long-range electrostatics and energy partition
- **Ewald/PPPM as the correct Coulomb summation under PBC:** Ensures converged forces/energies for ionic solids. $\rightarrow$ Why we keep PPPM rather than "screening away" Coulomb; ad-hoc screening changes the physics.
- **Per-atom energies/stresses for heat flux:** The Green-Kubo heat flux $\mathbf{J}=\frac{d}{d t} \sum_i \mathbf{r}_i e_i$ depends on a consistent partition of potential energy and virial, including k -space contributions. $\rightarrow$ Dictates the specific LAMMPS computes we use for GK.



# 4 State preparation: NPT $\rightarrow$ NVT $\rightarrow$ NVE
- **Why NPT first?** At target $T, P$, volume (and thus density, sound speeds) must be right; $k$ is volumesensitive. Barostats (e.g., Nosé-Hoover) sample the isothermal-isobaric ensemble. $\rightarrow$ We relax the supercell in NPT to get the correct density before measuring transport.
- **Why NVT second?** Thermostat equilibrates kinetic/potential modes without changing volume; damps residual barostat oscillations. $\rightarrow$ Produces a stable microstructure and temperature distribution prior to measurement.
- **Why NVE for GK production?** Fluctuation formulas assume unperturbed Hamiltonian dynamics; thermostats distort time correlations. $\rightarrow$ We switch off thermostats to record $\mathbf{J}(t)$ faithfully.



# 5 Green-Kubo (equilibrium) methodology
- **Formal result:** $k=\frac{1}{3 V k_B T^2} \int_0^{\infty}\langle\mathbf{J}(0) \cdot \mathbf{J}(t)\rangle d t$. Stationarity and ergodicity let us replace ensemble with time averages. $\rightarrow$ Justifies computing the HCACF and integrating to a plateau.
- **Correlation statistics:** HCACF is noisy with long-time tails; effective sample size is set by the integrated autocorrelation time. $\rightarrow$ Motivates multiple seeds, long windows, and conservative uncertainty estimation (blocking/Sokal window).
- **Finite-size and finite-time effects:** Small cells cap long-wavelength phonons; short runs bias plateaus. $\rightarrow$ We choose large supercells and check plateau stability vs. run length; one finite-size check along a box axis.


# 6 rNEMD (Müller-Plathe) as cross-check/backup
- **Principle:** Impose a microscopic heat flux by swapping particle momenta; measure steady $\nabla T$. In linear response, $k=-q / \nabla T$. $\rightarrow$ Offers an independent estimate sharing the same potential but different systematics.
- **Key caveats:** Need a clear linear gradient away from exchange slabs; avoid thermostats in the measurement region; check that imposed flux is small enough (no nonlinear heating). $\rightarrow$ Explains gradient quality checks and rate selection.


# 7 Vacancy generation and sampling theory
- **Ensemble of defect configurations:** Real materials have correlated vacancy/ $\mathrm{Ce}^{3+}$ arrangements; one configuration is a sample from a distribution. $\rightarrow$ We use reproducible random placement with minimum-distance constraints, then let dynamics relax; optionally average over seeds to reduce configuration bias.
- **$\mathrm{Ce}^{3+}$ assignment rule:** Nearest-neighbor reduction minimizes local electrostatic energy and reproduces known $\mathrm{Ce}^{3+}-\mathrm{V}_{\mathrm{O}}$ association trends. $\rightarrow$ Justifies the deterministic "two nearest Ce become $\mathrm{Ce}^{3+}$ " rule at insertion.

# 8 Integrators, thermostats, and barostats
- **Velocity-Verlet and timestep:** Stability tied to highest vibrational frequencies (O-Ce stretch); too large $\Delta t$ biases energy and transport. $\rightarrow$ Sets the $0.5-1.0$ fs guidance and the need to re-check virial/temperature drifts.
- **Nosé-Hoover (chains) and relaxation times:** Proper ensemble sampling requires coupling times matched to phonon periods; over-tight coupling distorts fluctuations, over-loose slows equilibration. $\rightarrow$ Explains chosen damping constants and chain lengths.
- **Parrinello-Rahman vs. isotropic barostat:** Box-shape fluctuations can matter in anisotropic phases; fluorite is cubic, so isotropic volume change suffices. $\rightarrow$ Justifies a simpler NPT for speed/stability.

# 9 Data analysis and uncertainty quantification
- **Building the HCACF:** Multiple time origins reduce variance; window functions balance bias/variance. $\rightarrow$ Why we use overlapping origins and a fixed cutoff policy.
- **Plateau picking:** Bias-variance tradeoff: integrate until the noise floor, then stop to avoid random walk of partial sums. $\rightarrow$ Motivates Sokal-style window or blocking analysis to report $k \pm \mathrm{Cl}$ transparently.
- **Autocorrelation and effective samples:** $N_{\text {eff }}=N / \kappa$ with $\kappa$ the integrated autocorrelation time; CIs should scale with $1 / \sqrt{N_{\text {eff }}}$. $\rightarrow$ Justifies run lengths and the "seeds $\times$ length" budgeting.

# 10 Validation targets and pass/fail criteria
- **Equation of state and structure:** Lattice parameter, density, and bulk modulus must match the potential's training/validation set at $T, P . \rightarrow$ Guards against using an ill-posed parameterization.
- **Reference $k$ for stoichiometric $\mathrm{CeO}_2$:** Theory/experiment give a bracket; MD value should live in the expected window when size and time are adequate. $\rightarrow$ Provides a concrete gate before exploring $x>0$.

# 11 Simulation design and grid selection
- **Commensurability and integer vacancies:** For $x$ to be exact, $x N_{\text {Ce }} \in \mathbb{Z}$; choose supercell sizes accordingly. $\rightarrow$ Explains the $10 \times 10 \times 10$ style boxes.
- **Design of experiments:** Balance temperature points and vacancy fractions against total correlation time budget; prioritize conditions with highest sensitivity for models that will consume $k(T, x) . \rightarrow$ Justifies a "small, high-quality map" over a broad, shallow sweep.



# 12 Reproducibility engineering
- **Random seeds and provenance:** In stochastic MD workflows, seeds are part of the scientific input. We fix/log seeds, code versions, input hashes; outputs include raw flux time series for re-analysis.
- **Deterministic post-processing:** Analysis scripts should be pure functions of inputs, with configuration captured in the output metadata. $\rightarrow$ Ensures anyone can recompute $k$ and CIs bit-forbit.


# 13 Known limitations and interpretive boundaries
- **Classical nuclei:** Quantum heat capacity and high-frequency modes are misrepresented at low $T$; trends vs. $x$ remain meaningful at moderate/high $T . \rightarrow$ Sets expectation on absolute values at, say, 300 K .
- **Potential physics:** No explicit electron-phonon or polaron effects; defect energetics limited by the force field. $\rightarrow$ Explains why we validate and avoid over-claiming beyond defect-scattering trends.