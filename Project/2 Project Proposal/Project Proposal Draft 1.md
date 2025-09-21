



# Project Implementation Summary


Ceria ( $\mathrm{CeO}_2$ ) is a standard surrogate for fluorite nuclear fuels and is widely used in swift heavy-ion (SHI) studies. In fuels, the lattice thermal conductivity $k$ sets centerline temperatures at a given power and enters models of fission-gas release, melting margins, and pellet-clad mechanics. In SHI physics, $k$ controls thermal-spike dissipation and is itself altered by ion-track defects. Oxygen vacancies are common in reduced ceria and act as phonon scatterers, so $k$ should depend on temperature $T$ and vacancy fraction $x$. The aim of this project is to compute a small, uncertainty-controlled map $k(T, x)$ for $\mathrm{CeO}_{2-x}$ using classical molecular dynamics, suitable for direct use in fuel-performance and irradiation models.

We will implement a reproducible thermal-conductivity workflow in LAMMPS. For each state ( $T, x$ ), we build a fluorite supercell, introduce vacancies to reach the target $x$ with uniform spatial distribution, equilibrate with NPT followed by NVT, and then run production dynamics. The primary route will be the Green-Kubo method: record the microscopic heat flux, compute the heat-current autocorrelation function, and integrate to a plateau. If long-range electrostatics complicate Green-Kubo flux accounting, we will use a screened Coulomb scheme (Wolf/DSF). If that is not stable on our timeline, we will switch to a Müller-Plathe non-equilibrium calculation with PPPM, impose a steady heat flux, measure the temperature gradient away from exchange regions, and compute $k=-q / \nabla T$. One route will be chosen and held fixed for the main results.

We will start with validation of a perfect ceria supercell before introducing vacancies. We will reproduce a reference-scale $k$ for stoichiometric $\mathrm{CeO}_2$ at one temperature within reported uncertainty for the chosen potential, and check basic properties (lattice parameter, density) against the potential's literature values. Convergence will be shown explicitly: stable time step, thermostat/barostat settings that do not distort fluctuations, plateau stability (or steady gradient) with increasing trajectory length, and a one-time finite-size test by extending one box dimension.

We will provide input generators for supercells and vacancy placement, a driver that executes the NPT $\rightarrow$ NVT $\rightarrow$ production sequence and logs required quantities, and a post-processor that outputs $k$ with confidence intervals using a fixed analysis protocol (blocking or Sokal-window). We expect to obtain a small set of $k(T, x)$ values with documented uncertainties, showing the reduction of $k$ with increasing $x$ and its temperature dependence. These numbers can be inserted into fuel centerline-temperature calculations or into thermal-spike models where $k$ is an input, and they provide a baseline for future extensions (additional $x$, dopants, or defect structures) without changing the method.


# Project Schedule


Got it—here’s a clear, step-by-step plan. I’ll use phases with concrete exit criteria so you can track progress and keep risk low.

### Phase 0 — Repo, environment, and paper trail (Day 0–2)

1.) **Set up GitHub + structure:** `README.md`, `env/` (conda YAML), `lammps/` (inputs), `scripts/` (orchestration, post-proc), `data/` (raw), `results/` (processed), `figures/`, `paper/` (Overleaf link).
2.) **Pin software:** LAMMPS version; Python 3.11; packages (numpy, scipy, matplotlib, pint, pyyaml). Save a `conda-env.yml`.
 3.)   **Choose the electrostatics path:**  
    _Option A (GK): Buckingham + Wolf/DSF screened Coulomb:_ Advantage is that `compute heat/flux` is valid out-of-the-box.  
    _Option B (NEMD fallback):_ Buckingham + PPPM and Müller–Plathe heat swap. Pick one as primary now.
4.) **Select the potential:** One well-cited fixed-charge Ce–O Buckingham model. Collect its paper + parameters in `lammps/potentials/` and cite in `README.md`.
5.) **Create Overleaf project:** Add a 1-page “living abstract” and a figure checklist.
    
**Exit criteria:** repo compiles; LAMMPS runs a tiny CeO₂ cell; Overleaf linked in README.

---

### Phase 1 — Stoichiometric CeO₂ sanity & EOS check (Day 2–5)

1.) **Build fluorite supercell.** Start 4×4×4 (≈192 atoms) for quick checks; target 6×6×6 to 8×8×8 (≈2.6–6.1k atoms) for production.
2.) **Relax lattice at 300 K (or target T).** NPT, gentle damping (τ_T≈0.5–1 ps, τ_P≈5–10 ps), Δt≈0.5–1.0 fs depending on stiffness.
3.) **Equation of state (optional but quick):** 5–7 volumes, NPT relax to get a, ρ; compare to potential’s reference values.

**Exit criteria:** equilibrium lattice constant and density match published potential within a few percent; stable MD without energy drift.

---

### Phase 2 — Core workflow build (GK or NEMD) (Day 5–10)

- **GK path (Option A).**
    
    1. NPT → NVT equilibration to the target T.
    2. NVE production; log per-step heat flux `J(t)` via `compute heat/flux` (pair, bond, kspace terms—kspace terms only if using DSF/Wolf, not PPPM).
    3. Post-proc: compute HCACF CJJ(t)C_{JJ}(t); cumulative integral ∫0tCJJ(t′)dt′\int_0^t C_{JJ}(t') dt'; identify plateau.
        
- **NEMD path (Option B).**
    
    1. NPT → NVT;
    2. Müller–Plathe heat swap to impose constant heat flux qq;
    3. Accumulate steady-state temperature profile; fit ∇T away from swap regions; compute κ=−q/∇T\kappa = -q/\nabla T.
        
- **Freeze analysis recipe.** Pick one HCACF windowing rule (e.g., Sokal/autocorrelation time or first zero-crossing); or one NEMD fitting window. Don’t change it later.
    

**Exit criteria:** one clean κ\kappa estimate for stoichiometric CeO₂ at one T with a visible plateau (GK) or linear gradient (NEMD).

---

### Phase 3 — Uncertainty + convergence hardening (Day 10–14)

1.) **Statistical error bars:** Blocking or batch means; report $95 \% \mathrm{Cl}$. For GK, define $t_{\text {cut }}$ and show plateau stability vs. run length.
2.) **Time step & thermostat sanity:** Sweep $\Delta t$ small $\pm 20 \%$; show $\kappa$ invariance within error. Verify thermostat/barostat do not bias production (NVE for GK).
3.) **Finite size:** One check: double the length in one direction; confirm $\kappa$ change within Cl or document trend.

**Exit criteria:** a figure with $\kappa+\mathrm{Cl}$, a plateau/gradient plot, and a one-slide convergence panel you could show tomorrow.

---

### Phase 4 — Introduce a single vacancy fraction (Day 14–20)

1.) **Build $x>0$ supercell:** Remove 0 atoms to reach $x \approx 0.03-0.05$; distribute vacancies evenly; relax at target T .
2.) **Run the same workflow:** Identical equilibration, production, and analysis settings as Phase 2.
3.) **Compare to $x=0$:** Expect $\kappa$ to drop; keep the same CI machinery.

**Exit criteria:** one $\kappa$ point at $x>0$ with Cl at the same temperature; a plot showing the drop vs $x$.



---

### Phase 5 — Fill a tiny (T,x)(T,x) grid (Day 20–32)

1.) **Temperatures:** Add 1–2 more T values (e.g., ~600, 1000, 1400 K) for both x=0x=0 and x>0x>0.
2.) **Automate:** Write a launcher script to generate folders, submit jobs, and register results in a CSV/JSON with metadata.
3.) **QC loop:** Auto-flag runs with drift, no plateau, or noisy gradients; rerun longer where needed.
    

**Exit criteria:** a $2 \times 3$ or $2 \times 2$ panel of $\kappa(T, x)$ with Cl ; reproducible scripts that regenerate the figures.


---

### Phase 6 — Optional sensitivity (only if ahead) (Day 32–38)

1.) **One comparator potential:** Run one state point under a second literature potential; show difference bars.
2.) **Robustness appendix:** If your primary was GK, do **one** NEMD cross-check at the baseline state point (or vice-versa).
    
**Exit criteria:** a single sensitivity figure; clearly labeled as “supplementary.”

---

### Phase 7 — Paper-ready polish (Day 38–End)

1.) **Figures.** Finalize: (a) κ\kappa vs TT at fixed xx; (b) κ\kappa vs xx at fixed TT; (c) plateau/gradient; (d) convergence panel.
2.) **Reproducibility.** `make results` style target; seed control; YAML config dump with every figure.
3.) **Overleaf.** Tight 6–8 page report: intro → method → validation → results → limits; include parameter table in SI.
4.) **Talk.** 10–12 slides: motivation, method schematic, validation, main panel, caveats, “future work”.
    

**Exit criteria:** repo clone → `make` → figures regenerate; Overleaf compiles; dry-run talk is <15 minutes including questions.

---

### Guardrails & defaults (use unless you must change them)

- **Time step:** 0.5-1.0 fs (pick one; document).
- **Thermostat/barostat: Nosé-Hoover;** $\tau \_T \approx 0.5-1 \mathrm{ps} ; \tau \_P \approx 5-10 \tau \_T$; NVE for GK production.
- **Supercell:** Start $6 \times 6 \times 6$; scale only if plateau/gradient demands it.
- **Production length:** Aim for multi-ns per state point; extend until plateau/CI stable.
- **Electrostatics:** DSF/Wolf for GK; PPPM for NEMD. Don't mix GK + PPPM unless you've implemented flux accounting for $k$-space.
- **Vacancies:** Fixed-charge model; distribute uniformly; relax before production; don't chase $\mathrm{Ce}^{3+}$ chemistry this semester.
    

---

### What to show at the status report (mid-semester)

- One validated $\kappa$ at $x=0$ with Cl and a plateau/gradient plot.
- The scripted pipeline (a single command that reproduces that figure).
- A roadmap slide with the $2 \times 3$ target grid and what's already done.

