

# 1 Transport target and modeling scope

Awesome—here’s a tight, curated reading pack for Section 1. I grouped items by idea and, where useful, pointed to concrete “how-to” pages or example scripts.

### Continuum energy balance, linear irreversible thermodynamics, and Fourier’s law

For the macro-to-micro bridge (continuity of energy, linear response, Onsager reciprocity), a classic, readable route is de Groot & Mazur’s text; it frames Fourier’s law as the first term in a linear phenomenology and explains when reciprocity holds. For a modern lecture treatment keyed to Green–Kubo, see the noted chapter below.

- S. R. de Groot & P. Mazur, _Non-Equilibrium Thermodynamics_ (Dover). Especially Chs. on linear laws and reciprocal relations. ([Dover Publications](https://store.doverpublications.com/products/9780486647418?srsltid=AfmBOor6I8e4RYrvmJMVbl84TpTBRJN7PMWDxsCSfPoInuh8gcw8s6BR&utm_source=chatgpt.com "Non-Equilibrium Thermodynamics"))
    
- J. Vinals, “Chapter 4: The Green–Kubo Relations” (course notes). A compact derivation from linear response with clear assumptions (stationarity/ergodicity). ([College of Science and Engineering](https://www-users.cse.umn.edu/~vinals/tspot_files/Chapter_4.pdf?utm_source=chatgpt.com "Chapter 4. The Green Kubo Relations"))
    

### Microscopic foundations: Irving–Kirkwood fluxes → Green–Kubo

These give you the microscopic definition of the energy/heat current used in MD and how the current–current correlation yields kk.

- J. H. Irving & J. G. Kirkwood, “The Statistical Mechanical Theory of Transport Processes. IV. The Equations of Hydrodynamics,” _J. Chem. Phys._ **18**, 817 (1950). Original continuum-from-microscopic derivation. ([AIP Publishing](https://pubs.aip.org/aip/jcp/article/18/6/817/201367/The-Statistical-Mechanical-Theory-of-Transport?utm_source=chatgpt.com "The Statistical Mechanical Theory of Transport Processes. IV ..."))
    
- D. J. Evans & G. P. Morriss, _Statistical Mechanics of Nonequilibrium Liquids_, 2e (Cambridge). See chapters “The Green–Kubo relations” and “Linear response.” ([Cambridge University Press & Assessment](https://www.cambridge.org/core/books/statistical-mechanics-of-nonequilibrium-liquids/E97F04B0D7D423540AA3F9B91F155E7D?utm_source=chatgpt.com "Statistical Mechanics of Nonequilibrium Liquids"))
    
- Y. Chen, “Physical foundation and consistent formulation of atomic-level fluxes,” _Phys. Rev. E_ **98**, 052113 (2018). Helpful for subtleties of flux definitions with many-body forces. ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevE.98.052113?utm_source=chatgpt.com "Physical foundation and consistent formulation of atomic-level ..."))
    

### Green–Kubo for thermal conductivity (formalism + modern implementations)

These show the exact GK expression, practical issues, and first-principles variants (good for context even if you’ll use classical MD).

- J. Kang _et al._, “First-principles Green–Kubo method for thermal conductivity,” _Phys. Rev. B_ **96**, 020302 (2017), plus preprint with details. ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevB.96.020302?utm_source=chatgpt.com "First-principles Green-Kubo method for thermal conductivity ..."))
    
- D. B. Zhang _et al._, “Green–Kubo formalism…,” _Phys. Rev. B_ **108**, 104307 (2023). Recent, clear GK exposition (tight-binding context but the GK part is general). ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevB.108.104307?utm_source=chatgpt.com "Green-Kubo formalism for thermal conductivity with Slater ..."))
    
- A. Pereverzev & E. Zivkovic, “Heat-current filtering for Green–Kubo…,” _Int. J. Heat Mass Transf._ **188** (2022). Practical discussion of noise/long-time tails and windowing. ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0017931022001296?utm_source=chatgpt.com "Heat-current filtering for Green–Kubo and Helfand-moment ..."))
    

### Equivalence and contrasts: GK vs. NEMD

These papers make precise when EMD/GK and rNEMD/direct methods agree, plus finite-size/time systematics—exactly what you summarized.

- P. K. Schelling, S. R. Phillpot, P. Keblinski, “Comparison of atomic-level simulation methods for computing thermal conductivity,” _Phys. Rev. B_ **65**, 144306 (2002). Benchmarks GK vs NEMD for Si; shows convergence/equivalence when each is done right. ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevB.65.144306?utm_source=chatgpt.com "Comparison of atomic-level simulation methods for computing thermal conductivity | Phys. Rev. B - Physical Review Link Manager"))
    
- H. Dong _et al._, “Equivalence of EMD and NEMD… from bulk to nanowire silicon,” _Phys. Rev. B_ **97**, 094305 (2018) (+ arXiv). Maps GK “running kk” to an effective length to compare with NEMD. ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevB.97.094305?utm_source=chatgpt.com "From bulk to nanowire silicon | Phys. Rev. B"))
    
- H. Dong _et al._, “Exactly equivalent thermal conductivity in finite systems…” _J. Therm. Anal. Calorim._ (2022). Finite-domain equivalence with periodic BCs. ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S1386947722002405?utm_source=chatgpt.com "Exactly equivalent thermal conductivity in finite systems ..."))
    

### rNEMD (Müller–Plathe) original + reviews

The canonical source for the velocity-swap method, plus follow-ups applying it to liquids.

- F. Müller-Plathe, “A simple nonequilibrium MD method for calculating thermal conductivity,” _J. Chem. Phys._ **106**, 6082 (1997). The original rNEMD algorithm. ([Astrophysics Data System](https://ui.adsabs.harvard.edu/abs/1997JChPh.106.6082M/abstract?utm_source=chatgpt.com "A simple nonequilibrium molecular dynamics method for ..."))
    
- M. Zhang _et al._, “Thermal conductivities of molecular liquids by reverse NEMD,” _J. Chem. Phys._ **123**, 114501 (2005). Clear worked example in fluids. ([PubMed](https://pubmed.ncbi.nlm.nih.gov/16852906/?utm_source=chatgpt.com "Thermal conductivities of molecular liquids by reverse ..."))
    

### Phonon-gas (BTE) perspective for context

For the modal formula k=13∑λCλvλ2τλk=\tfrac{1}{3}\sum_\lambda C_\lambda v_\lambda^2\tau_\lambda and how defects/anharmonicity enter τλ\tau_\lambda.

- J.-P. Hansen & I. R. McDonald, _Theory of Simple Liquids_ (Elsevier). Classic reference that also links Green–Kubo and transport coefficients. ([ScienceDirect](https://www.sciencedirect.com/book/9780123870322/theory-of-simple-liquids?utm_source=chatgpt.com "Theory of Simple Liquids - ScienceDirect.com"))
    
- W. Lv & A. Henry, “Examining the validity of the phonon gas model in amorphous solids,” _Sci. Rep._ **6**, 37675 (2016). Nice review-style discussion of when the phonon-gas picture holds. ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC5137137/?utm_source=chatgpt.com "Examining the Validity of the Phonon Gas Model in ..."))
    

### LAMMPS: concrete, reproducible implementations you can adapt

These are the exact docs and example inputs you can mine for your GK and rNEMD runs (syntax, computes, what to average, pitfalls).

- **How-to κ page**: four ways to compute kk in LAMMPS; points to `examples/KAPPA` and gives the GK and rNEMD workflows side-by-side. ([LAMMPS Documentation](https://docs.lammps.org/Howto_kappa.html?utm_source=chatgpt.com "10.3.6. Calculate thermal conductivity"))
    
- **`compute heat/flux`** (GK): how to assemble the heat-flux vector from per-atom KE/PE/stress; the doc also flags scope/limitations and cross-links to κ how-to. ([LAMMPS Documentation](https://docs.lammps.org/compute_heat_flux.html?utm_source=chatgpt.com "compute heat/flux command"))
    
- **`fix thermal/conductivity`** (Müller–Plathe): the swap mechanics, binning, periodicity factor-of-2 caveat, and how to turn tallied energy into qq. ([LAMMPS Documentation](https://docs.lammps.org/fix_thermal_conductivity.html?utm_source=chatgpt.com "fix thermal/conductivity command"))
    
- **Example inputs**: `examples/KAPPA/` with parsing notes on extracting kk for each method. Good as a skeleton before you swap in your ceria model. ([C4Science](https://c4science.ch/source/lammps/browse/master/examples/KAPPA/%3Ba3cbd21b8467b49d75e7d9fac6e3f1bd0493e686?utm_source=chatgpt.com "KAPPA · lammps"))
    
- Tutorial slides (Sandia/Plimpton): minimal working rNEMD scripts + GK outline; handy to see the end-to-end calculation on one page. ([LAMMPS](https://www.lammps.org/tutorials/italy14/italy_kappa_viscosity_Mar14.pdf?utm_source=chatgpt.com "Modeling Thermal Transport and Viscosity with ..."))
    

### Bonus: extra perspectives and course-notes if you want more derivations

- A. Kundu & A. Dhar, “Kubo formula for heat conduction in open systems,” _J. Stat. Mech._ (2009). Useful to understand boundary-driven vs. bulk GK in low-D anomalies. ([ICTS Home](https://home.icts.res.in/~abhi/Papers/46-greenkubo.pdf?utm_source=chatgpt.com "Kubo formula for heat conduction in open systems"))
    
- J.-S. Wang & L. Zhang, “Phonon Hall thermal conductivity from Green–Kubo,” _Phys. Rev. B_ **80**, 012301 (2009). An application that walks through GK tensor structure. ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevB.80.012301?utm_source=chatgpt.com "Phonon Hall thermal conductivity from the Green-Kubo formula"))
    

If you’d like, I can turn this into a little “reading map” PDF with one-line takeaways and the exact LAMMPS commands to try under each source, but the links above already give you the definitive derivations and the working input decks for GK and rNEMD.


# 2 Material physics you're encoding


Perfect—here’s a compact, high-signal source pack keyed to the exact claims in your Section 2, with a mix of reviews, primary papers, and “how to implement it” examples.

### Redox bookkeeping in ceria (Kröger–Vink, two Ce³⁺ per VO_\mathrm{O}, small-polaron association)

- **DFT(+U) evidence for two Ce³⁺ near each O vacancy.** Sun _et al._ show the excess electrons localize on two Ce ions adjacent to VO_\mathrm{O} (i.e., Ce⁴⁺→Ce³⁺) and discuss the well-known sensitivity to UU, but the _localization and pairing near the vacancy_ is robust. ([Physical Review](https://link.aps.org/accepted/10.1103/PhysRevB.95.245101?utm_source=chatgpt.com "Disentangling the role of small polarons and oxygen ..."))
    
- **Vacancy–polaron complex picture.** A focused study disentangles the binding between small polarons and oxygen vacancies in CeO2_2 (oxygen vacancies act as double donors that can bind up to two small polarons). Useful conceptual groundwork for your “two nearest Ce become Ce³⁺” assignment. ([ResearchGate](https://www.researchgate.net/publication/317311668_Disentangling_the_role_of_small_polarons_and_oxygen_vacancies_in_Ce_O_2?utm_source=chatgpt.com "Disentangling the role of small polarons and oxygen ..."))
    
- **Vacancy/charge ordering trends.** Broader ordering arguments and surface/bulk comparisons (vacancy formation induces two Ce³⁺; ordering along ⟨111⟩ discussed). Good for thinking about _correlated_ arrangements beyond a single defect. ([Physical Review](https://link.aps.org/accepted/10.1103/PhysRevB.92.144105?utm_source=chatgpt.com "Origins and implications of the ordering of oxygen ..."))
    

### Ionic radii (Ce³⁺ larger than Ce⁴⁺) → strain relief argument

- **Shannon radii tables** (8-coordination typical for fluorite Ce): Ce³⁺(VIII) ≈ 128 pm vs Ce⁴⁺(VIII) ≈ 111 pm. This underpins your statement that larger Ce³⁺ helps relax the tensile field around VO_\mathrm{O}. ([WebElements](https://www.webelements.com/cerium/atom_sizes.html?utm_source=chatgpt.com "radii of atoms and ions - Cerium"))
    

### Point-defect phonon scattering (Rayleigh ω4\omega^4 term, “disorder parameter”)

- **Modern review of Klemens-type theory.** Gurunathan _et al._ give a clear, equation-forward treatment connecting mass/size/force-constant contrasts to a Rayleigh ∝ω4\propto\omega^4 rate and validate against first-principles. Excellent for quantifying why VO_\mathrm{O}+Ce³⁺ complexes are strong scatterers. ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevApplied.13.034011?utm_source=chatgpt.com "Analytical Models of Phonon--Point-Defect Scattering"))
    
- **Classic experimental benchmark** on point-defect scattering reducing kk. Historical but still pedagogical. ([Physical Review](https://link.aps.org/doi/10.1103/PhysRev.131.1433?utm_source=chatgpt.com "Phonon Scattering by Point Defects | Phys. Rev."))
    
- **First-principles case study (actinide fluorite, methodological analog).** Malakkal _et al._ quantify kk degradation from point defects in ThO2_2 (same fluorite topology), helpful as a modern ab-initio anchor for the “τλ\tau_\lambda↓ ⇒ kk↓” logic. ([Chris Marianetti](https://www.chrismarianetti.org/static/publications/prm_8_025401_2024_malakkal.pdf?utm_source=chatgpt.com "First-principles determination of the phonon-point defect ..."))
    

### Interatomic models you can actually use (fixed-charge Buckingham+Q vs polarizable DIPPIM)

- **Critical assessment of ceria potentials (what’s been fit, what works).** Survey comparing six popular Buckingham+Q parameterizations (Grimes, Götte, etc.) against lattice/elastic/defect properties—useful to justify your validation gates. ([ResearchGate](https://www.researchgate.net/publication/220037543_A_critical_assessment_of_interatomic_potentials_for_ceria_with_application_to_its_elastic_properties?utm_source=chatgpt.com "A critical assessment of interatomic potentials for ceria with ..."))
    
- **Polarizable DIPPIM for stoichiometric, reduced, and doped ceria.** Burbano _et al._ introduce a dipole-polarizable model (induced dipoles on O, sometimes cations) fitted to DFT; good discussion of why polarizability improves transferability across Ce³⁺/Ce⁴⁺ mixes. ([PubMed](https://pubmed.ncbi.nlm.nih.gov/21654047/?utm_source=chatgpt.com "A dipole polarizable potential for reduced and doped CeO ..."))
    
- **Recent perspective on getting defect chemistry “right.”** Zhang _et al._ (Chem. Mater. 2022) compare shell-model/Buckingham forms and emphasize which properties are captured well (phonons, VO_\mathrm{O} barriers) and where caveats remain—exactly the mindset you want before trusting trends in k(T,x)k(T,x). ([ACS Publications](https://pubs.acs.org/doi/10.1021/acs.chemmater.2c03019?utm_source=chatgpt.com "Toward a Consistent Prediction of Defect Chemistry in CeO2"))
    

### Concrete implementations with explicit Ce³⁺/Ce⁴⁺ species (what your scripts will mirror)

- **Large-cell MD with _randomly distributed_ O vacancies and Ce³⁺ ions.** Classic oxygen self-diffusion study uses boxes up to ~33 000 ions with explicit Ce³⁺/Ce⁴⁺ and random VO_\mathrm{O} placement—this is the workflow pattern you’ll emulate for generating mixed-valence lattices. ([ResearchGate](https://www.researchgate.net/publication/229260556_Molecular_Dynamics_Study_of_Oxygen_Self-Diffusion_in_Reduced_CeO2?utm_source=chatgpt.com "Molecular Dynamics Study of Oxygen Self-Diffusion in ..."))
    
- **Ce³⁺/Ce⁴⁺ ordering from interionic potentials + DFT.** Shows how explicit species and their short-range parameters control local structure; helpful when choosing which Ce³⁺–O and Ce⁴⁺–O parameter sets to adopt. ([PubMed](https://pubmed.ncbi.nlm.nih.gov/19691397/?utm_source=chatgpt.com "Exploring Ce3+/Ce4+ cation ordering in reduced ceria ..."))
    
- **Background on the Götte/Grimes family (fluorite focus).** Thesis and comparisons that include the fluorite topology and Buckingham parameter choices—useful when you want to verify your parameter files against literature values. ([DIVA Portal](https://www.diva-portal.org/smash/get/diva2%3A169372/FULLTEXT01.pdf?utm_source=chatgpt.com "Dynamics in Ceria and Related Materials from Molecular ..."))
    

### Why you must keep proper long-range electrostatics (PPPM/Ewald) for ionic ceria

Even in a “chemistry” section, it’s worth grounding your PPPM/Ewald claim with primary documentation since VO_\mathrm{O}–Ce³⁺ interactions are long-ranged in an ionic solid.

- **LAMMPS k-space docs (Ewald/PPPM)**—what’s being summed and to what tolerance; pairs well with your argument against ad-hoc screening. ([LAMMPS Documentation](https://docs.lammps.org/kspace_style.html?utm_source=chatgpt.com "kspace_style command"))
    

---

## How to translate these into your workflow (quick pointers)

- **Ce³⁺ assignment rule (implementation):** Use the “oxygen self-diffusion in reduced CeO2_2” setup as your template: after sampling $VO_\mathrm{O}$ sites, convert the _two nearest_ Ce to Ce³⁺ (distinct species types with their own Buckingham parameters and charges), then equilibrate—this exactly mirrors the practice in large-cell MD with explicit mixed valence. ([ResearchGate](https://www.researchgate.net/publication/229260556_Molecular_Dynamics_Study_of_Oxygen_Self-Diffusion_in_Reduced_CeO2?utm_source=chatgpt.com "Molecular Dynamics Study of Oxygen Self-Diffusion in ..."))
    
- **Choosing parameters:** Start from a vetted Buckingham+Q set compared in the assessment paper (e.g., Götte/Grimes family) and verify your equilibrium a0a_0, density, and bulk modulus against the reported values before introducing VO_\mathrm{O}. If you need better dielectric/phonon fidelity, consider the DIPPIM set (costlier but more transferable). ([ResearchGate](https://www.researchgate.net/publication/220037543_A_critical_assessment_of_interatomic_potentials_for_ceria_with_application_to_its_elastic_properties?utm_source=chatgpt.com "A critical assessment of interatomic potentials for ceria with ..."))
    
- **Defect-scattering intuition you can quantify:** When you later analyze k(T,x)k(T,x), you can quote the Klemens-type expression for the ω4\omega^4 Rayleigh term and discuss how mass removal (O), force-constant changes (lost Ce–O bonds), and size mismatch (Ce³⁺>Ce⁴⁺) all feed the “disorder parameter,” hence τλ↓\tau_\lambda\downarrow and k↓k\downarrow. ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevApplied.13.034011?utm_source=chatgpt.com "Analytical Models of Phonon--Point-Defect Scattering"))
    

# 3. Long-Range Electrostatics and Energy Partition


Great—here’s a compact, _implementation-ready_ reading kit for your Section 3. I grouped items by the exact claims you made and added the most useful LAMMPS doc pages and notes so you can wire this straight into your scripts.

### Ewald/PPPM: why you need it (conditional convergence, parameter-independence once converged)

- **Textbook/notes (clear derivation):** Allen & Tildesley’s treatment + short course notes: Ewald turns the conditionally convergent Coulomb lattice sum into two absolutely convergent pieces (real + reciprocal) with a Gaussian split; neutrality is required. ([levich.ccny.cuny.edu](https://levich.ccny.cuny.edu/koplik/molecular_simulation/AT2.pdf?utm_source=chatgpt.com "Computer Simulation of Liquids"))
    
- **“Why Ewald” one-pager:** Ceperley’s lecture notes explicitly state the conditional convergence and the Ewald cure; good for a crisp quote. ([mcc.uiuc.edu](https://www.mcc.uiuc.edu/summerschool/2001/David%20Ceperley/dmc_lec2.htm?utm_source=chatgpt.com "dmc_lec2"))
    
- **Classic primary sources on boundary treatment:** de Leeuw–Perram–Smith I/II (lattice sums, dielectric boundary conditions) are the canonical references that set the “tin-foil vs vacuum” context for Ewald sums. ([Royal Society Publishing](https://royalsocietypublishing.org/doi/10.1098/rspa.1980.0135?utm_source=chatgpt.com "Simulation of electrostatic systems in periodic boundary ..."))
    
- **Parameter independence in practice:** Demonstrations that once the real-space cutoff and k-space mesh/accuracy are converged, energies/forces are independent of the splitting parameter. ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevE.63.056703?utm_source=chatgpt.com "Ewald sum for electronic bilayer systems | Phys. Rev. E"))
    

### “Tin-foil” (conducting) boundary condition in Ewald

- **What “tin-foil” means:** In Ewald, dropping the surface dipole term corresponds to surrounding the infinite lattice by a perfect conductor; many MD codes implement this by default. Useful review-style statements and formulas. ([Deep Blue](https://deepblue.lib.umich.edu/bitstream/handle/2027.42/70827/JCPSA6-101-6-5024-1.pdf%3Bsequence%3D2?utm_source=chatgpt.com "How the unit cell surface charge distribution affects the ..."))
    

### Why Wolf/DSF cutoffs are **not** appropriate for ionic crystals when you care about phonons/defects

- **DSF method (what it is, when it works):** Fennell & Gezelter show damped-shifted-force electrostatics can approximate Ewald in some homogeneous fluids; this is the best-known “pairwise alternative.” Use it as the authoritative description—then explain you’re _not_ using it for an ionic crystal’s long-wavelength physics. ([PubMed](https://pubmed.ncbi.nlm.nih.gov/16821904/?utm_source=chatgpt.com "Is the Ewald summation still necessary? Pairwise ..."))
    
- **Limits of cutoffs:** Comparative studies emphasize these methods trade physical fidelity for speed; accuracy depends strongly on environment and can distort dielectric response—precisely what matters for long-wavelength phonons and defect fields in ionic solids. ([glass.ruc.dk](https://glass.ruc.dk/pdf/articles/2012_JPhysChemB_116_5738.pdf?utm_source=chatgpt.com "Simplistic Coulomb Forces in Molecular Dynamics"))
    

### Microscopic heat/energy flux and “gauge consistency” (Irving–Kirkwood/Hardy → GK)

- **Foundational derivation:** Irving & Kirkwood derive the hydrodynamic equations and identify the microscopic energy flux—the starting point for GK. ([AIP Publishing](https://pubs.aip.org/aip/jcp/article/18/6/817/201367/The-Statistical-Mechanical-Theory-of-Transport?utm_source=chatgpt.com "The Statistical Mechanical Theory of Transport Processes. IV ..."))
    
- **Modern consistency clarifications:** Chen (2018) discusses physically consistent atomic fluxes, gauge freedom, and many-body potentials—a good, careful read before you finalize your flux assembly. ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevE.98.052113?utm_source=chatgpt.com "Physical foundation and consistent formulation of atomic-level ..."))
    
- **Many-body caveat (important for flux in polymers/complex FFs):** Surblys _et al._ show naive “atomic stress” approximations fail for many-body interactions and give the corrected form—useful background even if your ceria FF is pairwise. ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevE.99.051301?utm_source=chatgpt.com "Application of atomic stress to compute heat flux via molecular ..."))
    
- **Further reading on Hardy/IK variants and near-equivalence:** Recent analyses comparing Hardy vs. Wigner/other formulations and why, when implemented consistently, transport coefficients are gauge-invariant. ([arXiv](https://arxiv.org/pdf/2202.02246?utm_source=chatgpt.com "arXiv:2202.02246v2 [cond-mat.mtrl-sci] 6 Jul 2022"))
    

### LAMMPS: the exact docs you’ll use to build a **GK-correct** flux with PPPM

- **`compute heat/flux` (GK workhorse):** How to assemble J\mathbf J from per-atom KE, per-atom PE, and per-atom _stress_; notes on using the right stress flavor. ([LAMMPS Documentation](https://docs.lammps.org/compute_heat_flux.html?utm_source=chatgpt.com "compute heat/flux command"))
    
- **Per-atom stress with **k-space** included:** `compute stress/atom` / `compute centroid/stress/atom` explain how Ewald/PPPM contributions to the per-atom virial are computed (Heyes/Sirk methods). This is the key link proving your flux _includes k-space_. ([LAMMPS Documentation](https://docs.lammps.org/compute_stress_atom.html?utm_source=chatgpt.com "compute centroid/stress/atom command"))
    
- **Per-atom potential energy with **k-space** included:** `compute pe/atom` doc states that KSpace contributions are tallied (with extra FFTs). That’s the “convective” term’s PE piece handled correctly. ([LAMMPS Documentation](https://docs.lammps.org/compute_pe_atom.html?utm_source=chatgpt.com "compute pe/atom command"))
    
- **K-space controls you’ll tighten for accuracy:** `kspace_style pppm` and `kspace_modify` (order, mesh, pressure/virial notes). These are the knobs you converge so mesh noise is well below HCACF noise. ([LAMMPS Documentation](https://docs.lammps.org/kspace_style.html?utm_source=chatgpt.com "kspace_style command"))
    

#### Minimal LAMMPS stanza you can copy (pairwise ionic FF + PPPM)

```
kspace_style pppm 1.0e-5            # tighten until GK tails are stable
kspace_modify order 5               # typical; converge vs 4/6/7
compute ke all ke/atom
compute pe all pe/atom              # includes k-space contribution
compute s all centroid/stress/atom NULL virial   # per-atom virial incl. PPPM
compute j all heat/flux ke pe s
fix samp all ave/correlate 1 ${ns} ${nblk} c_j[1] c_j[2] c_j[3] type auto file J.out ave running
# (Then integrate the HCACF offline with your windowing protocol.)
```

Docs backing each line are above. The `centroid/stress/atom` variant is recommended when bonded/many-body terms exist; for purely pairwise Buckingham+Coulomb, `stress/atom` is also acceptable, but the centroid form is future-proof. ([LAMMPS Documentation](https://docs.lammps.org/compute_stress_atom.html?utm_source=chatgpt.com "compute centroid/stress/atom command"))

### PPPM performance/diagnostics (so you can _prove_ k-space noise is sub-dominant)

- **Accuracy/performance notes and “why extra FFTs appear” when you ask for per-atom virial/energy:** handy when you see a time-step cost jump. ([LAMMPS Documentation](https://docs.lammps.org/compute_stress_atom.html?utm_source=chatgpt.com "compute centroid/stress/atom command"))
    
- **Tuning PPPM meshes/orders:** practical guidance on `kspace_modify order` and expected scaling. ([LAMMPS Documentation](https://docs.lammps.org/kspace_modify.html?utm_source=chatgpt.com "kspace_modify command"))
    

---

If you want, next I can do the same treatment for Section 4 (NPT→NVT→NVE) and Section 5 (GK statistics), with quick-start input fragments and links to exemplar scripts that already produce clean HCACFs.



# 4. NPT → NVT → NVE (Set the state →settle the microstate → measure without touching )


Fantastic—here’s a focused source pack for your **Section 4 (NPT → NVT → NVE)** with theory anchors plus practical, LAMMPS-ready references. I’ve grouped them by the specific claims in your text, and where it helps I’ve pointed to the exact knobs to use.

### Why NPT first: get the right density/elastic background at (T,P)

- **Foundations of variable-cell/pressure MD.** Parrinello–Rahman’s extended-Lagrangian (cell degrees of freedom) is the classic reference for relaxing volume/shape so phonon spectra live on the correct elastic background. ([Astrophysics Data System](https://ui.adsabs.harvard.edu/abs/1981JAP....52.7182P/abstract?utm_source=chatgpt.com "Polymorphic transitions in single crystals: A new molecular ..."))
    
- **Modern constant-pressure Nosé–Hoover family (MTK).** Martyna–Tobias–Klein derive stable, modularly invariant NPTNPT equations—the backbone of most “Nosé–Hoover barostat” implementations you’ll use. ([Astrophysics Data System](https://ui.adsabs.harvard.edu/abs/1994JChPh.101.4177M/abstract?utm_source=chatgpt.com "Constant pressure molecular dynamics algorithms"))
    
- **Practical barostat guidance.** Short, pedagogical notes on when to prefer isotropic NPT vs. Parrinello–Rahman (use PR for anisotropic solids/transformations; isotropic is fine—and stabler—for cubic fluorites like ceria). ([computecanada.github.io](https://computecanada.github.io/molmodsim-md-theory-lesson-novice/08-barostats/index.html?utm_source=chatgpt.com "Controlling Pressure – Practical considerations for Molecular ..."))
    
- **LAMMPS implementation knobs.** `fix npt` uses Nosé–Hoover thermostat+barostat; doc shows chain lengths, coupling constants, and restart behavior (so your NPT state cleanly hands off to NVT/NVE). ([LAMMPS Documentation](https://docs.lammps.org/fix_npt_body.html?utm_source=chatgpt.com "fix npt/body command"))
    

### Why NVT second: settle modes without changing volume

- **Nosé–Hoover chains (NHC) for canonical sampling.** The canonical reference showing why chain thermostats equilibrate efficiently without over-damping short-time dynamics—exactly what you want to kill barostat ring-down while preserving local fluctuations. ([AIP Publishing](https://pubs.aip.org/aip/jcp/article/97/4/2635/927962/Nose-Hoover-chains-The-canonical-ensemble-via?utm_source=chatgpt.com "Nosé–Hoover chains: The canonical ensemble via continuous ..."))
    
- **LAMMPS thermostat controls.** `fix nvt` / `fix nh` docs: set `tchain` (thermostat chain length) and coupling times; also `pchain` if you ever use PR/MTK. Good page to justify “longer than the fastest optical period” guidance. ([LAMMPS Documentation](https://docs.lammps.org/fix_nh.html?utm_source=chatgpt.com "fix nvt command"))
    

### Why NVE for GK production: measure without touching

- **Linear response requires unperturbed Hamiltonian dynamics.** Textbook treatments (Evans–Morriss; Tuckerman) lay out GK derivation and the caution about external couplings; ideal citations for your “thermostats imprint their own correlation times” paragraph. ([mini.ourphysics.org](https://mini.ourphysics.org/wiki/images/0/0b/Denis_J._Evans%2C-Statistical_Mechanics_of_Nonequilibrium_Liquids_%28Cambridge_2008%29.pdf?utm_source=chatgpt.com "STATISTICAL MECHANICS OF NONEQUILIBRIUM LIQUIDS"))
    
- **LAMMPS κ how-to.** The official “Calculate thermal conductivity” page spells out that GK uses equilibrium fluctuations of per-atom KE/PE/stress from an equilibrated state—i.e., run the _measurement_ under pure dynamics. Pair this with your flux assembly notes. ([LAMMPS Documentation](https://docs.lammps.org/Howto_kappa.html?utm_source=chatgpt.com "10.3.6. Calculate thermal conductivity"))
    

### rNEMD edge-case note (if you switch methods)

- **Original Müller–Plathe scheme.** Source for “swap hottest/coldest” and the practice of keeping the bulk region unthermostatted while pinning only the buffer ends; helpful when you contrast GK vs rNEMD in your write-up. ([ResearchGate](https://www.researchgate.net/profile/Mark-Tuckerman/publication/234921675_Nos-Hoover_chains_The_canonical_ensemble_via_continuous_dynamics/links/0c96052d4117cbd47b000000/Nos-Hoover-chains-The-canonical-ensemble-via-continuous-dynamics.pdf?utm_source=chatgpt.com "Nosé–Hoover chains: The canonical ensemble via ..."))
    

### Concrete, copy-able implementation references (LAMMPS)

- **GK input skeleton and pipeline.** `Howto_kappa` and `compute heat/flux` docs (what fields enter J\mathbf J, how to average with `fix ave/correlate`, and why to sample after full equilibration). Use these as your exact “what to type” anchors. ([LAMMPS Documentation](https://docs.lammps.org/Howto_kappa.html?utm_source=chatgpt.com "10.3.6. Calculate thermal conductivity"))
    
- **Thermostat/Barostat parameterization in LAMMPS.** The `fix nh`/`fix nvt` page lists `tchain`/`pchain` and relaxation-time settings—you can cite this when you justify choices like “τT∼0.1−0.5\tau_T \sim 0.1{-}0.5 ps; τP∼1−5\tau_P \sim 1{-}5 ps.” ([LAMMPS Documentation](https://docs.lammps.org/fix_nh.html?utm_source=chatgpt.com "fix nvt command"))
    

If you want, I can also drop a minimal two-file example: one NPT→NVT prep script that writes a restart, and one NVE GK script that reads it, sets `compute heat/flux`, and emits a clean HCACF ready for your windowed integral.



# 5. Green-Kubo (equilibrium) methodology


Awesome—here’s a _reader-friendly_ source pack for your Section 5 text, grouped by the exact claims you make and with concrete, copy-able implementation refs where helpful.

### Green–Kubo formalism, assumptions (stationarity, ergodicity), and the GK integral

- Derivations and linear-response background (GK from FDT): Evans & Morriss, chs. “Green–Kubo relations” / “Linear response.” ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevB.110.184203?utm_source=chatgpt.com "Diffusion, long-time tails, and localization in classical and ..."))
    
- Modern GK exposition (concise, materials-focused): Zhang et al., _Phys. Rev. B_ 108 (2023). ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevB.108.104307?utm_source=chatgpt.com "Green-Kubo formalism for thermal conductivity with Slater ..."))
    
- First-principles GK (nice companion read even if you’ll use classical MD): Kang et al., _Phys. Rev. B_ 96 (2017). ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevB.96.020302?utm_source=chatgpt.com "First-principles Green-Kubo method for thermal conductivity ..."))
    

### Practical GK in LAMMPS (how to implement exactly what you wrote)

- “How to compute κ” page (GK workflow + example input for solid Ar): LAMMPS Howto_kappa. ([LAMMPS Documentation](https://docs.lammps.org/Howto_kappa.html?utm_source=chatgpt.com "10.3.6. Calculate thermal conductivity"))
    
- Heat-flux assembly for GK: `compute heat/flux` (what goes into **J**, pitfalls). ([LAMMPS Documentation](https://docs.lammps.org/compute_heat_flux.html?utm_source=chatgpt.com "compute heat/flux command"))
    
- Time-correlation accumulation for HCACF: `fix ave/correlate` (overlapping time origins, integration). ([LAMMPS Documentation](https://docs.lammps.org/fix_ave_correlate.html?utm_source=chatgpt.com "fix ave/correlate command"))
    
- Older mirrored docs & example notes (handy when googling): legacy `compute_heat_flux` pages. ([ENEA](https://www.afs.enea.it/software/lammps/doc17/html/compute_heat_flux.html?utm_source=chatgpt.com "compute heat/flux command — LAMMPS documentation"))
    
- Worked GK tutorials / code to parse **J0Jt** and integrate (good for a first run): UIUC argon GK project; example notebook (Python) for GK post-processing. ([Physics Courses at Illinois](https://courses.physics.illinois.edu/phys466/sp2011/projects/2004/Team1/index.html?utm_source=chatgpt.com "Thermal Conductivity of Argon from the Green-Kubo Method"))
    

### Why you run the GK measurement in NVE (no thermostat during correlations)

- Textbook cautions about external couplings in linear-response GK (see same Evans & Morriss and Tuckerman-style lecture notes used above). ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevB.110.184203?utm_source=chatgpt.com "Diffusion, long-time tails, and localization in classical and ..."))
    
- LAMMPS κ page explicitly shows equilibrate → switch to pure dynamics for GK. ([LAMMPS Documentation](https://docs.lammps.org/Howto_kappa.html?utm_source=chatgpt.com "10.3.6. Calculate thermal conductivity"))
    

### Statistics of the HCACF: overlapping origins, IAT, blocking, Sokal window

- Blocking for correlated data (gold standard): Flyvbjerg & Petersen, _J. Chem. Phys._ 91, 461–466 (1989). ([Astrophysics Data System](https://ui.adsabs.harvard.edu/abs/1989JChPh..91..461F/abstract?utm_source=chatgpt.com "Error estimates on averages of correlated data - ADS"))
    
- Integrated autocorrelation time & effective samples (Sokal’s classic notes + accessible summaries/tools). ([SpringerLink](https://link.springer.com/chapter/10.1007/978-1-4899-0319-8_6?utm_source=chatgpt.com "Monte Carlo Methods in Statistical Mechanics: Foundations ..."))
    
- Application-oriented overview of GK uncertainty handling in MD (good “what can go wrong” read). ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0017931016341643?utm_source=chatgpt.com "Uncertainty quantification of thermal conductivities from ..."))
    

### Long-time tails and plateau picking (why you window the integral)

- Hydrodynamic long-time tails (original VACF results + modern reviews; establishes power-law tails that motivate windowed integrals). ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevA.1.18?utm_source=chatgpt.com "Decay of the Velocity Autocorrelation Function | Phys. Rev. A"))
    

### Finite-size / finite-time effects; GK vs NEMD equivalence checks

- Canonical comparison and finite-size analysis: Schelling, Phillpot & Keblinski, _Phys. Rev. B_ 65, 144306 (2002). ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevB.65.144306?utm_source=chatgpt.com "Comparison of atomic-level simulation methods for computing ..."))
    
- Quantitative EMD↔NEMD equivalence and length/time mapping (bulk → nanowires): Dong et al., _Phys. Rev. B_ 97, 094305 (2018), plus follow-ups on finite systems. ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevB.97.094305?utm_source=chatgpt.com "From bulk to nanowire silicon | Phys. Rev. B"))
    

---

#### Minimal LAMMPS→Python pathway you can mirror

1. In LAMMPS, after NVE starts:
    

```
compute ke all ke/atom
compute pe all pe/atom
compute s  all stress/atom NULL virial
compute J  all heat/flux ke pe s
fix hcacf all ave/correlate 1 ${Ns} ${Nb} c_J[1] c_J[2] c_J[3] type auto file J0Jt.dat ave running
```

(then disable all fixes except pure integration; write `J0Jt.dat`). See Howto_kappa and `compute heat/flux`. ([LAMMPS Documentation](https://docs.lammps.org/Howto_kappa.html?utm_source=chatgpt.com "10.3.6. Calculate thermal conductivity"))

2. Offline, read `J0Jt.dat`, form CJJ(t)C_{JJ}(t), integrate with a **window** t⋆=c τintt^\star=c\,\tau_\text{int} (e.g., c ⁣∼ ⁣5c\!\sim\!5), and get CIs via **blocking** or shard-wise variance. For blocking theory and IAT estimation recipes, use Flyvbjerg–Petersen and Sokal. ([Astrophysics Data System](https://ui.adsabs.harvard.edu/abs/1989JChPh..91..461F/abstract?utm_source=chatgpt.com "Error estimates on averages of correlated data - ADS"))
    

If you want, I can produce a tiny, plug-in Python script that (i) estimates τint\tau_\text{int}, (ii) picks t⋆t^\star via a Sokal-style rule, and (iii) returns k±k\pmCI from `J0Jt.dat`.


# 6. rNEMD (Müller-Plathe) as cross-check/backup


Here’s a tight reading/implementation pack for your rNEMD section, matched to each claim and with concrete LAMMPS hooks so you can reproduce everything you describe.

The original method and its logic are Müller-Plathe’s 1997 JCP paper. It introduces the velocity-swap scheme, defines the imposed heat flux qq from the cumulative exchanged kinetic energy, and shows how to extract k=−q/∇Tk=-q/\nabla T from the bulk-linear segment of the temperature profile. This is the canonical citation for the “reverse” idea (impose flux, measure gradient). ([AIP Publishing](https://pubs.aip.org/aip/jcp/article/106/14/6082/181799/A-simple-nonequilibrium-molecular-dynamics-method?utm_source=chatgpt.com "A simple nonequilibrium molecular dynamics method for calculating ..."))

For a modern, code-level implementation, the LAMMPS command `fix thermal/conductivity` implements the Müller-Plathe algorithm directly, with documentation that spells out swap frequency, slab geometry, periodicity caveats (the “two directions” factor in periodic boxes), and the formula to convert tallied energy into a heat flux. This is the page you’ll cite when you explain exactly how you computed qq and set the swap cadence. ([LAMMPS Documentation](https://docs.lammps.org/fix_thermal_conductivity.html?utm_source=chatgpt.com "fix thermal/conductivity command"))

Your requirement that the driving remain in the linear regime—and your “halve/double the swap rate” invariance check—follows best practice established in comparative studies that show EMD/GK and NEMD converge when nonlinearity and size effects are controlled. Schelling, Phillpot & Keblinski’s classic PRB comparison is the go-to reference for demonstrating equivalence and for discussing finite-size/time issues and rate dependence; it’s frequently cited when justifying that rNEMD and GK should agree within error once both are in their proper limits. ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevB.65.144306?utm_source=chatgpt.com "Comparison of atomic-level simulation methods for computing ..."))

Everything about how you build and read the temperature profile has direct, documented support in LAMMPS. To _measure_ ∇T\nabla T away from the exchange slabs, define spatial slabs with `compute chunk/atom bin/1d` along the transport axis and time-average bin temperatures with `fix ave/chunk`; this is the supported replacement for the older, deprecated `fix ave/spatial`. These pages give you the exact knobs for bin width, bounds, and averaging, and they’re what you cite for “exclude a few bins near the exchange regions, fit a straight line over the central window, and report R2R^2/residuals.” ([LAMMPS Documentation](https://docs.lammps.org/compute_chunk_atom.html?utm_source=chatgpt.com "compute chunk/atom command"))

If you need a temperature definition that removes any slab-wise streaming component before computing kinetic temperature, `compute temp/profile` is the sanctioned tool and its docs explain the subtraction of the spatially averaged velocity field per bin. That’s the right citation for your “one-dimensional on average; subtract any flow before computing TT” sentence. ([LAMMPS Documentation](https://docs.lammps.org/compute_temp_profile.html?utm_source=chatgpt.com "compute temp/profile command"))

Your note about keeping the _bulk_ microcanonical while optionally pinning only buffer zones is squarely in line with both the original rNEMD paper and the LAMMPS documentation; the latter emphasizes that the algorithm itself imposes the flux, so thermostats are not needed (and can be counterproductive) in the measurement region. Use the `fix thermal/conductivity` page for the practical warning and Müller-Plathe for the theoretical rationale. ([LAMMPS Documentation](https://docs.lammps.org/fix_thermal_conductivity.html?utm_source=chatgpt.com "fix thermal/conductivity command"))

Removing spurious center-of-mass momentum is also by the book. In LAMMPS, `fix momentum` (or `fix momentum/chunk` if you ever need per-slab corrections) zeroes linear momentum at a chosen interval; the manual and core developer notes make the distinction with `fix recenter`, which _doesn’t_ remove momentum. Cite these when you mention de-drifting the system so convective motion doesn’t contaminate the heat current. ([LAMMPS Documentation](https://docs.lammps.org/fix_momentum.html?utm_source=chatgpt.com "fix momentum/chunk command"))

If you want additional worked examples beyond the LAMMPS manual, there are tutorial-style walks that use rNEMD to build temperature gradients and read interfacial conductance or bulk kk; these are handy to see the full pipeline—slab binning, linear fit windows, rate scans—in one place before porting the pattern to ceria. ([QuantumATK Documentation](https://docs.quantumatk.com/tutorials/interfacial_thermal_conductance/interfacial_thermal_conductance.html?utm_source=chatgpt.com "Interfacial thermal conductance - QuantumATK"))

Finally, if you need to anchor that rNEMD is widely used outside simple fluids and retains the same “swap → flux → gradient” mechanics, there are many application papers adapting the method to polymers, composites, and nanoconfinement; they’re useful as practical demonstrations of the exact diagnostics you listed (slope symmetry upon folding the profile, linearity checks, superposed profiles over seeds/blocks). ([PubMed](https://pubmed.ncbi.nlm.nih.gov/16852906/?utm_source=chatgpt.com "Thermal conductivities of molecular liquids by reverse ... - PubMed"))

If it’s helpful, I can also drop a minimal two-file example: one script that equilibrates, then runs `fix thermal/conductivity` with a conservative swap rate, and a second that bins the profile with `compute chunk/atom`/`fix ave/chunk`, writes the per-slab temperatures, and fits the central window to report kk with a bootstrap CI over time blocks.


# 7 Vacancy generation and sampling theory


Great—here’s a compact source + implementation pack matched to your vacancy-generation/sampling section. I’ll flag one tiny typo first: it should be xNCe∈ZxN_{\mathrm{Ce}}\in\mathbb{Z} (not NGeN_{\mathrm{Ge}}).

### Why two Ce³⁺ per O-vacancy; local association and strain logic

- **Vacancy leaves two electrons that localize on neighboring Ce (small polarons)**—clear DFT(+U) evidence and discussion of vacancy–polaron complexes. These establish the “VO∙∙+2Ce3+V_{\mathrm O}^{\bullet\bullet}+2\mathrm{Ce}^{3+}” bookkeeping you implement. ([Physical Review](https://link.aps.org/accepted/10.1103/PhysRevB.95.245101?utm_source=chatgpt.com "Disentangling the role of small polarons and oxygen ..."))
    
- **Charge disproportionation around VOV_{\mathrm O}** (two Ce become Ce³⁺, two remain Ce⁴⁺) and anisotropic strain in cubic ceria—nice link to your “Ce³⁺ larger radius relaxes tensile field” argument. ([National Science Foundation](https://par.nsf.gov/servlets/purl/10198111?utm_source=chatgpt.com "Anisotropic Chemical Strain in Cubic Ceria due to Oxygen- ..."))
    
- **Explicit MD workflows in reduced ceria** that do exactly what you describe: remove O’s, then _replace double the number of Ce⁴⁺ by Ce³⁺_, typically the nearest neighbours. Use these as concrete precedents for your generator. ([Northwestern Civil Engineering](https://www.civil.northwestern.edu/people/qu/Qu%20Papers/Cui%20-%20SS%20Ionics%20-%20Molecular%20dynamics%20simulation%20of%20reduced%20CeO2.pdf?utm_source=chatgpt.com "Molecular dynamics simulation of reduced CeO2"))
    

### How people place defects in practice (random vs. correlated) and why “blue-noise” is sensible

- **Large-cell MD with randomly distributed VOV_{\mathrm O} and Ce³⁺** (4k–33k ions): a direct precedent for “draw, relax, then measure”—good model of your _draw-relax-measure_ pipeline. ([ResearchGate](https://www.researchgate.net/publication/229260556_Molecular_Dynamics_Study_of_Oxygen_Self-Diffusion_in_Reduced_CeO2?utm_source=chatgpt.com "Molecular Dynamics Study of Oxygen Self-Diffusion in ..."))
    
- **Cation/charge ordering in reduced ceria**: shows that correlations _do_ emerge from the force field, so it’s better not to hard-wire patterns at insertion time; let NPT→NVT express them. ([AIP Publishing](https://pubs.aip.org/aip/jcp/article/131/6/064701/188580/Exploring-Ce3-Ce4-cation-ordering-in-reduced-ceria?utm_source=chatgpt.com "cation ordering in reduced ceria nanoparticles using interionic ..."))
    
- **Blue-noise / Poisson-disk sampling with a minimum-distance constraint**: practical algorithms and the conflict/coverage metrics you can quote when you explain your rmin⁡r_{\min} rule. (Computer-graphics literature, but ideal for your “forbid tiny separations, avoid pathological clustering” rationale.) ([Visual Computing Lab](https://vcg.isti.cnr.it/Publications/2012/CCS12/TVCG-2011-07-0217.pdf?utm_source=chatgpt.com "Efficient and Flexible Sampling with Blue Noise Properties of ..."))
    

### Configuration averaging, self-averaging, and finite-size images

- **Random-placement sensitivity and “randomly distributed vacancies” improving transport** (older but widely cited): helps justify reporting vacancy–vacancy separation stats and avoiding adjacent-image pathologies. ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0167273800006755?utm_source=chatgpt.com "Molecular dynamics calculations on ceria-based solid ..."))
    
- General MD best-practice: under PBC each defect feels its images; mitigate by large cells and by reporting minimum image separations (you already propose this). A short practical reminder page with multi-run scripting patterns (seed control) from LAMMPS workshops is handy to cite for reproducible seeds/shards. ([LAMMPS](https://www.lammps.org/workshops/Aug17/pdf/tutorial_beginner_quick_tour.pdf?utm_source=chatgpt.com "A Quick Tour of LAMMPS"))
    

### Implementation pointers you can copy/adapt

- **Create vacancies and re-label nearest Ce as Ce³⁺ (script pattern):**  
    – Use a pre-processing script (Python) to: (i) list all O-site indices, (ii) sample nvacn_\mathrm{vac} under a Poisson-disk/rmin⁡r_{\min} constraint, (iii) for each chosen O, delete it and find the _two nearest Ce_ by distance under PBC, and (iv) change their species/type to Ce³⁺ (with Ce³⁺–O parameters and charge). The Cui _Solid State Ionics_ paper and Lawrence et al. explicitly describe the “remove O → switch two nearest Ce⁴⁺ to Ce³⁺” step; cite them in your methods. ([Northwestern Civil Engineering](https://www.civil.northwestern.edu/people/qu/Qu%20Papers/Cui%20-%20SS%20Ionics%20-%20Molecular%20dynamics%20simulation%20of%20reduced%20CeO2.pdf?utm_source=chatgpt.com "Molecular dynamics simulation of reduced CeO2"))  
    – If you prefer an _in-LAMMPS_ prototype for the mechanics of making point defects, this FCC vacancy tutorial shows the pattern (delete atoms by ID, then continue)—you’ll generalize it to your fluorite lattice and species retyping. ([cavs.msstate.edu](https://www.cavs.msstate.edu/icme/code/lammps/tutorials/lammps/vacancy.php?utm_source=chatgpt.com "LAMMPS Vacancy | Code Repository | ICME - cavs.msstate.edu"))
    
- **Nearest-neighbor lookup / coordination under PBC:** in LAMMPS you can dump neighbor lists or use analysis tools; for pre-processing, most folks do it in Python with KD-trees + minimum-image convention. A community thread shows runtime neighbor info, useful if you ever sanity-check distances in-sim. ([Materials Science Community Discourse](https://matsci.org/t/nearest-neighbours-information-at-simulation-time/45978?utm_source=chatgpt.com "Nearest neighbours information at simulation time"))
    

### Extra context you can cite when discussing correlations that emerge after relaxation

- **Vacancy/polaron binding, migration, and ordering trends** (recent reviews/perspectives): use these when you explain that short-range associations should re-appear naturally after NPT→NVT even when you start from a blue-noise placement. ([UCL Discovery](https://discovery.ucl.ac.uk/id/eprint/10162438/1/acs.chemmater.2c03019.pdf?utm_source=chatgpt.com "Toward a Consistent Prediction of Defect Chemistry in CeO2"))
    

If you want, I can also hand you a small, ready-to-run Python pre-processor that (1) reads a fluorite POSCAR/LAMMPS data file, (2) performs Poisson-disk sampling on O sites with your chosen rmin⁡r_{\min}, (3) re-labels the two nearest Ce as Ce³⁺ per vacancy, and (4) writes out the modified structure plus CSV diagnostics (min-image separations; histograms of VOV_{\mathrm O}–VOV_{\mathrm O} and Ce³⁺–VOV_{\mathrm O} distances) so you can drop those plots straight into your methods.


# 8 Integrators, thermostats, and barostats


Here’s a compact reading + implementation pack for your Section 8, matched to the exact claims you make and with concrete LAMMPS hooks.

Velocity–Verlet is the standard because it is time-reversible and symplectic, i.e., it approximately conserves a “shadow” Hamiltonian over long times; that’s why it controls secular energy drift better than non-symplectic schemes. Good derivations and diagnostics live in Allen–Tildesley’s text and in papers that explicitly analyze shadow energy for symplectic integrators (you can quote these when you justify your timestep choice and NVE drift checks). ([PagePlace](https://api.pageplace.de/preview/DT0400.9780192524706_A35505741/preview-9780192524706_A35505741.pdf?utm_source=chatgpt.com "Computer Simulation of Liquids"))

The practical timestep bound is set by the highest vibrational frequency (here, stiff Ce–O modes). You don’t need a closed-form ωmax; the standard practice is to test ∆t by measuring NVE energy drift over the same horizon used for GK and verifying it’s orders of magnitude below the stochastic spread of the running integral. Allen–Tildesley’s treatment of integrator stability, plus modern notes on reversible/symplectic schemes, are the right citations for this “choose ∆t, then verify by NVE drift” workflow. ([PagePlace](https://api.pageplace.de/preview/DT0400.9780192524706_A35505741/preview-9780192524706_A35505741.pdf?utm_source=chatgpt.com "Computer Simulation of Liquids"))

Nosé–Hoover chains (Martyna–Klein–Tuckerman, MKT) are the canonical route for clean canonical sampling in solids. The original papers derive (i) Nosé–Hoover chains for NVT and (ii) the modularly invariant constant-pressure Nosé–Hoover equations that generate NPT, together with explicit reversible integrators; these are exactly what LAMMPS implements in `fix nvt`/`fix npt`. Use them to defend your guidance on coupling times (τT for NVT; τP for NPT) and short chain lengths. ([Duke Statistical Science](https://www2.stat.duke.edu/~scs/Projects/REMD/NoseHooverChains1992.pdf?utm_source=chatgpt.com "Nose--Hoover chains"))

LAMMPS documentation then gives you the precise knobs: `fix nvt`/`fix nh` for thermostatting (chain length via `tchain`, coupling time via the third argument to `temp`), and `fix npt` for isotropic NPT with MTK barostatting. The “How-to: barostats” page summarizes when to use isotropic vs. anisotropic control; for cubic fluorite ceria, isotropic is both sufficient and stabler than full Parrinello–Rahman. ([LAMMPS Documentation](https://docs.lammps.org/fix_nh.html?utm_source=chatgpt.com "fix nvt command"))

Two small mechanics matter for the numerical cleanliness you’re after. First, remove residual center-of-mass momentum during prep (and during rNEMD, if needed) so a tiny drift doesn’t masquerade as convective heat current: `fix momentum` (or its per-chunk variant) is the sanctioned tool; don’t confuse it with `fix recenter`, which only shifts coordinates and does not remove momentum. ([LAMMPS Documentation](https://docs.lammps.org/fix_momentum.html?utm_source=chatgpt.com "fix momentum/chunk command")) Second, mind neighbor-list and k-space tolerances: loose PPPM or rare rebuilds introduce force/virial “texture” that no integrator can cure. The neighbor/Verlet-buffer docs explain skin-size vs. rebuild tradeoffs (and “dangerous builds”), while the k-space developer notes summarize PPPM’s real/reciprocal split; together, these are the dials you tighten until NVE drift and virial noise fall below your HCACF noise floor. ([LAMMPS Documentation](https://docs.lammps.org/neighbor.html?utm_source=chatgpt.com "neighbor command"))

If you want a minimal, citation-backed stanza to drop into your inputs:

```lmp
# --- NPT: set the (T,P) state on a cubic fluorite cell ---
variable T equal 300
variable P equal 1.0
timestep 0.001                     # 1.0 fs; verify by NVE drift (see refs)
fix prep all npt temp ${T} ${T} 0.2 iso ${P} ${P} 2.0   # NH chain; MTK barostat
#                       ↑ τT ~0.2 ps (optical<τT<acoustic few periods)  ↑ τP ~2 ps
run 200000

unfix prep
fix settle all nvt temp ${T} ${T} 0.2 tchain 3          # NVT to kill barostat ring-down
run 200000
unfix settle

# --- NVE: measure without touching (GK stage) ---
reset_timestep 0
fix removeDrift all momentum 100 linear 1 1 1           # de-drift COM periodically
timestep 0.0005                                         # 0.5 fs if needed; verify
# PPPM must be tight enough that GK tails are stable
kspace_style pppm 1.0e-5
kspace_modify order 5

compute ke all ke/atom
compute pe all pe/atom
compute s  all centroid/stress/atom NULL virial
compute J  all heat/flux ke pe s                         # GK flux assembly

fix hcacf all ave/correlate 1 10000 10000 c_J[1] c_J[2] c_J[3] type auto file J0Jt.dat ave running
run 2000000
```

`fix nvt`/`fix npt` (Nosé–Hoover chains; MTK), `fix momentum`, and the neighbor/k-space guidance are straight from the LAMMPS manual; the integrator/ensemble theory is from MKT and Allen–Tildesley, and the “shadow Hamiltonian → drift check” logic is from the energy-conservation literature on symplectic integrators. ([LAMMPS Documentation](https://docs.lammps.org/fix_nh.html?utm_source=chatgpt.com "fix nvt command"))

For completeness, if you ever need anisotropic cell relaxation (not the case for cubic ceria), the Parrinello–Rahman family and modern variants are your theoretical anchors; otherwise, stick with isotropic NPT for speed and stability. ([MDPI](https://www.mdpi.com/2076-3417/12/3/1139?utm_source=chatgpt.com "Molecular Dynamics of Solids at Constant Pressure and ..."))

That set of references and the snippet above should let you justify every number you pick (∆t, τT, τP, chain lengths) and show exactly how you verified “no hidden numerical heating or artificial damping” before you turn on GK.



# 9 Data analysis and uncertainty quantification


Here’s a focused reading kit for your Section 9, mapped to each idea in your text and with concrete implementation hooks where they matter.

Start with the linear-response foundation and the exact Green–Kubo object you’re estimating. Evans & Morriss lay out the GK relations, their assumptions (stationarity, ergodicity), and why the transport coefficient is an integral of a current–current correlation; it’s the cleanest “why your recipe is valid” reference for the paragraph that opens your section.

For building the HCACF from a discrete MD record and doing it **the LAMMPS way**, read the two correlation fixes. `fix ave/correlate` is the standard overlapping-origins accumulator your procedure describes; the doc explains how samples are formed and averaged, which is precisely why your estimator variance drops at fixed run length when you use overlaps. `fix ave/correlate/long` implements a multiple-τ scheme that extends the correlation horizon efficiently with progressively coarser resolution—handy if you want to check that the tail truly sits at the noise floor without exploding memory. Both pages give exact syntax you can mirror. ([LAMMPS Documentation](https://docs.lammps.org/fix_ave_correlate.html?utm_source=chatgpt.com "fix ave/correlate command"))

Long-time tails justify **windowed** integrals instead of chasing the asymptote forever. The classic results of Alder & Wainwright (and modern confirmations) show algebraic decay of velocity autocorrelations from hydrodynamic mode coupling (e.g., t−3/2t^{-3/2} in 3D), which is the canonical reason plateaus can meander if you integrate too far into noise. These are ideal citations for your “random walk of partial sums” warning and the need to stop at a physically motivated cutoff. ([Physical Review](https://link.aps.org/doi/10.1103/PhysRevA.1.18?utm_source=chatgpt.com "Decay of the Velocity Autocorrelation Function | Phys. Rev. A"))

On **uncertainty quantification for correlated data**, the gold standard is Flyvbjerg–Petersen’s blocking analysis: recursively coarse-grain until adjacent blocks decorrelate, then use the variance of block means. It’s exactly the tool you point to when you say “blocking (a.k.a. batch means).” If you want an updated, automated take (with proofs and code), the more recent “automated blocking” paper is a useful modern companion. Pair those with Sokal’s classic notes on integrated autocorrelation time (IAT) and effective sample size; Sokal is the canonical source for Neff≈T/τintN_{\text{eff}}\approx T/\tau_{\text{int}} and the 1/Neff1/\sqrt{N_{\text{eff}}} scaling that your section emphasizes. ([Academia](https://www.academia.edu/63287314/Error_estimates_on_averages_of_correlated_data?utm_source=chatgpt.com "Error estimates on averages of correlated data"))

Your **plateau-picking** rule (“Sokal-style window” t⋆=c τintt^\star=c\,\tau_{\text{int}} with conservative cc) is standard practice in MD/MC error analysis; it is essentially Sokal’s prescription (estimate τint\tau_{\text{int}}, choose a cutoff proportional to it, report the windowed integral and its standard error), transplanted to the GK integral. Use Sokal here again for the derivation and language around bias–variance tradeoffs. ([SpringerLink](https://link.springer.com/chapter/10.1007/978-1-4899-0319-8_6?utm_source=chatgpt.com "Monte Carlo Methods in Statistical Mechanics: Foundations ..."))

For **finite-size/time systematics** and the sanity of cross-checking against rNEMD when needed, Schelling–Phillpot–Keblinski is still the most-cited comparative study. It documents exactly the phenomena you describe: missing long-wavelength content in small EMD boxes, length dependence in NEMD, and the conditions under which both approaches converge to the same kk. It’s the right anchor for your “large supercell + one explicit size check + multiple long shards with independent seeds” protocol. ([Materials Science Community Discourse](https://matsci.org/uploads/short-url/gUoU92bPJUEkDeqLHGRHglhEx1g.pdf?utm_source=chatgpt.com "Comparison of atomic-level simulation methods for ..."))

Finally, to tie the **implementation** straight into your scripts, keep the LAMMPS docs for the GK pipeline at hand: use `fix ave/correlate` (or `/long`) to dump Jα(0)Jα(t)J_\alpha(0)J_\alpha(t), then integrate offline with your window t⋆t^\star and error bars from blocking or shard-wise variance; that’s exactly what the documentation intends. If you want extra color on why overlapping-origin estimators behave the way you note (the “triangular availability” with lag), the `fix ave/correlate` page explicitly explains the sampling cadence and averaging windows. ([LAMMPS Documentation](https://docs.lammps.org/fix_ave_correlate.html?utm_source=chatgpt.com "fix ave/correlate command"))

If you’d like, I can also drop a small, cited Python helper that (i) estimates τint\tau_{\text{int}} from a noisy HCACF, (ii) picks t⋆=c τintt^\star=c\,\tau_{\text{int}}, and (iii) returns k±k\pmCI via Flyvbjerg–Petersen blocking—so your analysis steps match these sources one-for-one.



# 10 Validation targets and pass/fail criteria


Awesome—here’s a tight, high-signal source pack for **Section 10 (Validation targets & pass/fail)**. I’ve grouped items by the specific gates you described and included practical docs where you’ll actually touch knobs.

### A. Equation of state & structure at (T,P)(T,P)

- **Birch–Murnaghan EOS (what you fit, what B, V0B,\,V_0 mean).** Birch’s original Phys. Rev. paper is the canonical citation for the finite-strain EOS you’ll use to extract BB and V0V_0. ([Physical Review Links](https://link.aps.org/doi/10.1103/PhysRev.71.809?utm_source=chatgpt.com "Finite Elastic Strain of Cubic Crystals | Phys. Rev."))
    
- **Murnaghan EOS (background/contrast).** Handy as a second reference if you do a quick linear-in-PP bulk-modulus check before Birch–Murnaghan. ([Wikipedia](https://en.wikipedia.org/wiki/Murnaghan_equation_of_state?utm_source=chatgpt.com "Murnaghan equation of state"))
    
- **Constant-pressure MD with Nosé–Hoover (MTK).** The derivation that underpins modern NPT (what LAMMPS implements in `fix npt`), i.e., the equations you rely on while scanning pressures for the EOS fit. ([physics.ujep.cz](https://physics.ujep.cz/~mlisal/md/martyna-tobias-klein_md.pdf?utm_source=chatgpt.com "Constant pressure molecular dynamics algorithms"))
    
- **LAMMPS barostat/thermostat docs (the “how”).** `Howto_barostat` overview and `fix nvt`/`fix npt` details; shows where to set isotropic control for a cubic fluorite and how to pick relaxation times. ([LAMMPS Documentation](https://docs.lammps.org/Howto_barostat.html?utm_source=chatgpt.com "10.2.5. Barostats"))
    

### B. Long-range electrostatics: convergence & “quiet numerics”

- **PPPM accuracy & pressure/virial cautions.** The `kspace_style pppm` doc explicitly warns that energy/pressure accuracy can lag force accuracy and should be converged when using a barostat (exactly your gate about halving the tolerance & enlarging cutoffs without shifting the mean pressure beyond equilibrium noise). ([LAMMPS Documentation](https://docs.lammps.org/kspace_style.html?utm_source=chatgpt.com "kspace_style command"))
    
- **`kspace_modify` tuning.** How to set absolute/relative accuracy and mesh/order; cite this when you define the PPPM tightening you’ll lock across the grid. ([LAMMPS Documentation](https://docs.lammps.org/kspace_modify.html?utm_source=chatgpt.com "kspace_modify command"))
    

### C. Stoichiometric kk gate: GK setup, plateau behavior, and cross-check

- **LAMMPS “How to compute κ\kappa” (GK + NEMD).** The official playbook that your workflow mirrors: equilibrate, switch to pure dynamics, build HCACF from per-atom KE/PE/stress, integrate to a plateau; also links the rNEMD alternative you’ll use for a sanity cross-check. ([LAMMPS Documentation](https://docs.lammps.org/Howto_kappa.html?utm_source=chatgpt.com "10.3.6. Calculate thermal conductivity"))
    
- **`compute heat/flux` (what actually goes into J\mathbf J).** Exact inputs (per-atom KE, PE, stress) and cross-links to the GK page—good to cite right where you claim to assemble J\mathbf J consistently. ([LAMMPS Documentation](https://docs.lammps.org/compute_heat_flux.html?utm_source=chatgpt.com "compute heat/flux command"))
    
- **GK ↔ NEMD equivalence & finite-size/time systematics.** Schelling–Phillpot–Keblinski’s classic comparison paper is the standard reference for your procedural gates: steady plateau insensitivity, shard scaling ∼1/M\sim1/\sqrt{M}, and “elongate one axis” tests; it also motivates the rNEMD cross-check. ([Physical Review Links](https://link.aps.org/doi/10.1103/PhysRevB.65.144306?utm_source=chatgpt.com "Comparison of atomic-level simulation methods for computing ..."))
    

### D. “Bracket” for defect-free CeO2\mathrm{CeO_2} at 300 K (to sanity-check your absolute)

- **Ab-initio/BTE side (clean crystal).** Recent first-principles work reports room-temperature lattice kk for stoichiometric CeO2_2 in the mid-teens W m−1^{-1} K−1^{-1} (e.g., ∼17\sim17 W m−1^{-1} K−1^{-1} at 0 GPa), furnishing the upper, defect-free “theory” bracket you referenced. ([ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0375960121006332?utm_source=chatgpt.com "Lattice thermodynamic behavior in CeO 2 from first principles"))
    
- **Dense-pellet experiments (microstructure lowers kk).** Experimental datasets on bulk CeO2_2 pellets (laser-flash diffusivity → kk) show significant spread with density/porosity and impurities—use these as your lower/realistic bracket and to justify why you compare to high-density specimens only. ([OSTI](https://www.osti.gov/servlets/purl/1469554?utm_source=chatgpt.com "Thermal and mechanical properties of CeO2"))
    

### E. What to quote in your methods (concise, defensible gates)

- **EOS gates.** “Lattice parameter & density within ~1–2% of published potential values; BB within 10–15% from a Birch–Murnaghan fit” → anchored by Birch (EOS theory), MTK (proper NPT sampling), and LAMMPS barostat docs. ([Physical Review Links](https://link.aps.org/pdf/10.1103/PhysRev.71.809?utm_source=chatgpt.com "Finite Elastic Strain of Cubic Crystals"))
    
- **PPPM gate.** “Halving the k-space tolerance and extending real-space cutoffs does not shift mean pressure beyond equilibrium fluctuations” → directly motivated by `kspace_style` accuracy notes. ([LAMMPS Documentation](https://docs.lammps.org/kspace_style.html?utm_source=chatgpt.com "kspace_style command"))
    
- **GK plateau gates.** “Windowed integral shows a broad flat region; shard means collapse as 1/M1/\sqrt{M}; doubling total correlation time changes kk by <1σ; elongating one axis leaves kk within error” → backed by LAMMPS κ how-to and the Schelling PRB comparison. ([LAMMPS Documentation](https://docs.lammps.org/Howto_kappa.html?utm_source=chatgpt.com "10.3.6. Calculate thermal conductivity"))
    
- **Absolute sanity bracket.** “Final kk sits between BTE clean-crystal estimates (mid-teens W m−1^{-1} K−1^{-1} at 300 K) and high-density pellet measurements” → use Bao et al. (first-principles) plus high-density experimental compilations (e.g., Suzuki/OSTI; Morrison review). ([ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0375960121006332?utm_source=chatgpt.com "Lattice thermodynamic behavior in CeO 2 from first principles"))
    

If you want, I can also tack on a one-page checklist with the exact LAMMPS commands you’ll run for each gate (pressure sweep for EOS, PPPM tightening deltas, GK/rNEMD convergence plots) so the validation becomes a literal “run this → check that” protocol.


# 11 Simulation design and grid selection


Awesome—here’s a tight reading/implementation pack for **Section 11 (Simulation design & grid selection)**, mapped directly to the claims in your text.

### Commensurability in fluorite CeO₂ (why Z=4Z=4, how NCe=4n3N_{\text{Ce}}=4n^3)

- **Crystal facts you need** (structure, ZZ, space group): multiple crystallographic entries list **fluorite (Fm3ˉ\bar3m) with Z=4Z=4** for CeO₂—this is the basis of NCe=4n3N_{\text{Ce}}=4n^3 and NO=8n3N_{\text{O}}=8n^3 in an n×n×nn\times n\times n supercell. ([SpringerMaterials](https://materials.springer.com/isp/crystallographic/docs/sd_1923902?utm_source=chatgpt.com "CeO2 Crystal Structure"))
    
- **Lattice parameter near 5.40–5.41 Å** (so a 10×10×1010\times10\times10 box is ≈5.4 nm): recent experimental and review sources report a≈5.40–5.41a\approx 5.40\text{–}5.41 Å for stoichiometric CeO₂. ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9565251/?utm_source=chatgpt.com "Structure and Surface Relaxation of CeO2 Nanoparticles ..."))
    

### Why composition must be _exact_ under PBC (and why “round” boxes help)

- **Image interactions & separations**: vacancy images sit one or more unit cells away under PBC; explicit examples in fluorites use a≈5.41a\approx 5.41 Å to quantify image distances—useful to sanity-check that your box puts images far enough apart at small xx. ([OSTI](https://www.osti.gov/servlets/purl/1185884?utm_source=chatgpt.com "Chemical expansion affected oxygen vacancy stability in ..."))
    

### Finite-size effects for kk (missing long-qq phonons) and the need for a size check

- **Canonical comparison & size analysis** (GK vs NEMD): Schelling–Phillpot–Keblinski analyze **finite-size and finite-time** effects and show how properly done EMD/NEMD converge—this underwrites your “one elongated cell” gate and the preference for large nn. ([Materials Science Community Discourse](https://matsci.org/uploads/short-url/gUoU92bPJUEkDeqLHGRHglhEx1g.pdf?utm_source=chatgpt.com "Comparison of atomic-level simulation methods for ..."))
    

### Correlation-time budgeting (why error bars scale like 1/Neff1/\sqrt{N_{\text{eff}}})

- **Effective sample size with autocorrelation**: clear, citable statement that uncertainties scale with **ESS** NeffN_{\text{eff}} rather than raw samples—exactly your Neff≈Trun/τintN_{\text{eff}}\approx T_{\text{run}}/\tau_{\text{int}} reasoning for allocating time per (T,x)(T,x). ([Stan](https://mc-stan.org/docs/2_19/reference-manual/effective-sample-size-section.html?utm_source=chatgpt.com "15.4 Effective Sample Size | Stan Reference Manual"))
    

### Practical tools you’ll actually use while building the grid

- **Overlapping-origin correlation accumulation** for GK (so one _long_ shard gives maximal information): LAMMPS `fix ave/correlate` docs describe the overlapping-origin estimator you reference and the sampling cadence—handy when you pin a common GK setup for every grid point. ([LAMMPS Documentation](https://docs.lammps.org/fix_ave_correlate.html?utm_source=chatgpt.com "fix ave/correlate command"))
    

---

## How these sources map to your specific numbers

- **“10×10×1010\times10\times10 gives NCe=4000N_{\text{Ce}}=4000”** follows directly from Z=4Z=4 fluorite data (SpringerMaterials entries); the linear size L≈10a≈5.4L\approx 10a\approx 5.4 nm uses the reported a≃5.40–5.41a\simeq 5.40\text{–}5.41 Å (PMC/Elsevier sources). This justifies your exact-integer counts at x=0.03x=0.03 and 0.060.06 (120 and 240 vacancies) and the box length you quote. ([SpringerMaterials](https://materials.springer.com/isp/crystallographic/docs/sd_1923902?utm_source=chatgpt.com "CeO2 Crystal Structure"))
    
- **“Do at least one size check”** (elongate one axis) is precisely the finite-size protocol advocated in Schelling et al.; cite it when you report that the plateau stays within error upon elongation. ([Materials Science Community Discourse](https://matsci.org/uploads/short-url/gUoU92bPJUEkDeqLHGRHglhEx1g.pdf?utm_source=chatgpt.com "Comparison of atomic-level simulation methods for ..."))
    
- **“Spend correlation time, not just steps”** uses the ESS logic in the Stan reference manual; it’s a standard, clear citation for your 1/Neff1/\sqrt{N_{\text{eff}}} scaling and the “few high-quality (T,x)(T,x) points” design. ([Stan](https://mc-stan.org/docs/2_19/reference-manual/effective-sample-size-section.html?utm_source=chatgpt.com "15.4 Effective Sample Size | Stan Reference Manual"))
    
- **“Use the same GK plumbing everywhere”** is supported by the LAMMPS correlation docs—pin the `fix ave/correlate` configuration once (stride, window, outputs) and keep it fixed across the entire (T,x)(T,x) grid. ([LAMMPS Documentation](https://docs.lammps.org/fix_ave_correlate.html?utm_source=chatgpt.com "fix ave/correlate command"))
    

If you’d like, I can bundle these into a one-page checklist with the exact arithmetic for commensurate nn at your chosen xx values, plus a tiny helper that, given xx and aa, proposes the smallest nn with x⋅4n3∈Zx\cdot4n^3\in\mathbb{Z} and estimates the linear size L=naL=na so you can pick between “adjust xx” vs “go to the next nn”.



# 12 Reproducibility engineering


Awesome—here’s a compact, implementation-ready source pack for **Section 12 (Reproducibility engineering)**. I’ve grouped items by the exact claims in your text and added the most useful primary docs for “how” details.

### Seeds, manifests, and provenance you can defend in a methods section

- **Ten Simple Rules for Reproducible Computational Research** — canonical guidance to _record exact versions, parameters, and workflows_; perfect to justify your run manifests and frozen protocol. ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3812051/?utm_source=chatgpt.com "Ten Simple Rules for Reproducible Computational Research"))
    
- **The Turing Way (reproducible research + version control)** — practical handbook on provenance, VCS, testing, CI; cite when you describe content-addressed manifests and code-hash checks in your folder layout. ([book.the-turing-way.org](https://book.the-turing-way.org/reproducible-research/reproducible-research?utm_source=chatgpt.com "Guide for Reproducible Research - The Turing Way"))
    
- **LiveCoMS “Best practices”** (Braun/Mobley/Shirts et al.) — simulation-specific reproducibility conventions (log seeds, software versions, parameter files; separate inputs/outputs; ship analysis scripts). ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6884151/?utm_source=chatgpt.com "Best Practices for Foundations in Molecular Simulations ..."))
    
- **FAIR Principles** — why you keep self-describing data + metadata (so others can find/understand/reuse), i.e., your `manifest.json`, `control.json`, and checksums. ([Nature](https://www.nature.com/articles/sdata201618?utm_source=chatgpt.com "The FAIR Guiding Principles for scientific data ..."))
    

### Content-addressed manifests & hashing (why you store code/data checksums)

- **Pro Git — “Git is a content-addressable filesystem”** — clean, citable explanation of object hashes; use when you say “content-addressed manifest with code/potential hashes.” ([Git](https://git-scm.com/book/en/v2/Git-Internals-Git-Objects?utm_source=chatgpt.com "10.2 Git Internals - Git Objects"))
    
- **Git’s SHA-256 transition** — official docs noting repositories can be SHA-256; good footnote for your “SHA-256 of flux arrays” and “code-hash drift” notes. ([Git](https://git-scm.com/docs/hash-function-transition?utm_source=chatgpt.com "hash-function-transition Documentation"))
    

### Binary formats for flux/HCACF (why NetCDF/HDF5 over text)

- **HDF5 User Guide** — self-describing, portable binary format, parallel-I/O capable; perfect to justify `flux.bin` as HDF5. ([The HDF Group Support Site](https://support.hdfgroup.org/documentation/hdf5/latest/_u_g.html?utm_source=chatgpt.com "HDF5 User Guide"))
    
- **NetCDF Users Guide** — classic + netCDF-4 (HDF5-backed) formats; concise data-model description you can cite for “structured, self-describing arrays + metadata.” ([docs.unidata.ucar.edu](https://docs.unidata.ucar.edu/nug/current/file_structure_and_performance.html?utm_source=chatgpt.com "NetCDF Users Guide: File Structure and Performance"))
    

### Why bitwise identity is hard in MD (parallel reductions, neighbor lists, RNG streams)

- **Floating-point non-associativity → non-deterministic reductions** — Demmel/Ahrens et al. on reproducible summation; excellent to justify “bitwise stability is tricky in MPI/threads.” ([EECS at UC Berkeley](https://www2.eecs.berkeley.edu/Pubs/TechRpts/2015/EECS-2015-229.pdf?utm_source=chatgpt.com "Efficient Reproducible Floating Point Summation and BLAS"))
    
- **MPI reduction/reordering reproducibility** — classical note from Argonne on why operation order changes answers; good one-line cite. ([mcs.anl.gov](https://www.mcs.anl.gov/papers/P4093-0713_1.pdf?utm_source=chatgpt.com "On the Reproducibility of MPI Reduction Operations"))
    
- **NVIDIA/HPC notes** on non-determinism from asynchronous parallel reductions; handy general reference (even if you’re on CPUs). ([NVIDIA Developer Forums](https://forums.developer.nvidia.com/t/numerical-reproducibility-randomness/285547?utm_source=chatgpt.com "Numerical Reproducibility & Randomness"))
    
- **LAMMPS-specific determinism** — INTEL package note: deterministic _given the same parallel config/libraries_; pairs with your “fix MPI layout” advice. ([LAMMPS Documentation](https://docs.lammps.org/Speed_intel.html?utm_source=chatgpt.com "9.4.2. INTEL package"))
    
- **LAMMPS RNG & processor-dependent behavior** — docs/forum threads describing per-processor RNG streams and why results can differ with core counts; cite where you say “don’t reuse seeds; keep shard MPI layouts constant.” ([LAMMPS Documentation](https://docs.lammps.org/fix_gld.html?utm_source=chatgpt.com "fix gld command"))
    
- **Neighbor-list details** — developer doc on Verlet neighbor lists (order/updates can differ across decompositions), supporting your note about “optimistic” reordering and rebuild settings. ([LAMMPS Documentation](https://docs.lammps.org/Developer_par_neigh.html?utm_source=chatgpt.com "4.4.3. Neighbor lists"))
    

### Practical templates that map 1-to-1 to your workflow

- **Stan/ESS style reference for NeffN_{\text{eff}}** — concise, standard citation for using effective sample size logic in your CI/QA (ties to your windowed-integral CIs). ([livecomsjournal.org](https://livecomsjournal.org/index.php/livecoms/catalog/category/bestpractices?utm_source=chatgpt.com "Best Practices"))
    
- **FAIR pages** — short, quotable bullets you can drop beside your folder tree explaining _why_ you store checksums and metadata with every artifact. ([GO FAIR](https://www.go-fair.org/fair-principles/?utm_source=chatgpt.com "FAIR Principles"))
    

---

## “Do it like this” (implementation notes tied to the sources)

- **Manifests**: include _code commit_, _LAMMPS version/build flags_, _potentials’ SHA-256_, _RNG seeds per phase/shard_, _PPPM/neighbor tolerances_, _MPI/OMP layout_, _platform info_. This mirrors Sandve Rules 1 & 3 and The Turing Way’s provenance chapters. ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3812051/?utm_source=chatgpt.com "Ten Simple Rules for Reproducible Computational Research"))
    
- **Bitwise-ish stability**: keep _identical_ MPI layout per shard; set a fixed PPPM mesh/order; avoid features that reorder neighbors unpredictably; if available to you, the **LAMMPS INTEL** package can enforce determinism on CPU for a given parallel config. ([LAMMPS Documentation](https://docs.lammps.org/Speed_intel.html?utm_source=chatgpt.com "9.4.2. INTEL package"))
    
- **Binary outputs**: write **HDF5/netCDF-4** for `flux` and `hcacf` arrays with attributes `{T,V,nsteps,dt,pppm_tol,neighbor_skin,seed,…}`; store a **SHA-256** next to each file and record it in `k_result.json`. ([The HDF Group Support Site](https://support.hdfgroup.org/documentation/hdf5/latest/_u_g.html?utm_source=chatgpt.com "HDF5 User Guide"))
    
- **Why results differ across core counts**: order of floating-point reductions changes; neighbor-list & domain decomposition differ; RNG stream partitioning may change — cite the Argonne/NVIDIA notes and LAMMPS RNG docs in your methods. ([mcs.anl.gov](https://www.mcs.anl.gov/papers/P4093-0713_1.pdf?utm_source=chatgpt.com "On the Reproducibility of MPI Reduction Operations"))
    

If you want, I can also hand you a tiny `manifest_writer.py` that computes SHA-256 for inputs, snapshots `sys.version`, `lmp -h` build flags, and your JSON configs, then writes a single `manifest.json` exactly in the shape you outlined.




# 13 Known limitations and interpretive boundaries


Here’s a tight, annotated source pack for **Section 13 (Known limitations & interpretive boundaries)**—each item maps to specific claims in your text and gives you something quotable or directly usable.

### Classical nuclei → low-T bias; “quantum corrections” caveats

- J. Liu _et al._, “Assessing the quantum effect in classical thermal transport simulations,” _J. Appl. Phys._ **129**, 235104 (2021). Clear demonstration that classical MD overpopulates high-ω modes; shows why post-hoc “quantum corrections” are not generally reliable for kk. ([AIP Publishing](https://pubs.aip.org/aip/jap/article/129/23/235104/286381/Assessing-the-quantum-effect-in-classical-thermal?utm_source=chatgpt.com "Assessing the quantum effect in classical thermal ..."))
    
- H. Zhou _et al._, “Impact of classical statistics on thermal conductivity…,” _Appl. Phys. Lett._ **125**, 172202 (2024) + OSTI open version. Quantifies specific-heat and lifetime errors from classical statistics; good language for your low-T caveat. ([AIP Publishing](https://pubs.aip.org/aip/apl/article/125/17/172202/3318140/Impact-of-classical-statistics-on-thermal?utm_source=chatgpt.com "Impact of classical statistics on thermal conductivity ..."))
    
- J. A. Thomas & A. J. H. McGaughey, “Assessing the applicability of quantum corrections…,” preprint classic often cited to warn against simple heat-capacity rescaling of kk. ([ResearchGate](https://www.researchgate.net/publication/228652410_Assessing_the_Applicability_of_Quantum_Corrections_to_Classical_Thermal_Conductivity_Predictions?utm_source=chatgpt.com "Assessing the Applicability of Quantum Corrections to ..."))
    

### Force-field physics boundaries (fixed-charge Buckingham+Coulomb; missing polaron/e–ph)

- Q. Mao _et al._, “Classical and reactive MD: principles and applications,” _Prog. Mater. Sci._ **142**, 101072 (2023). Broad, up-to-date review: limitations of classical FFs; transferability cautions—useful to justify “results are only as good as the fit.” ([ScienceDirect](https://www.sciencedirect.com/science/article/pii/S036012852300014X?utm_source=chatgpt.com "Classical and reactive molecular dynamics: Principles and ..."))
    
- L. Sun _et al._, “Disentangling the role of small polarons and oxygen vacancies in CeO2_2,” _Phys. Rev. B_ **95**, 245101 (2017). DFT(+U) study: oxygen vacancies are double donors; small polarons localize on Ce and bind near VOV_\mathrm O; anchors your “polaronic conduction absent in our model” point. ([Physical Review Links](https://link.aps.org/accepted/10.1103/PhysRevB.95.245101?utm_source=chatgpt.com "Disentangling the role of small polarons and oxygen ..."))
    
- J. S. Pelli Cresi _et al._, “Ultrafast formation of small polarons and the optical gap in CeO2_2,” _J. Phys. Chem. Lett._ **11**, 6211 (2020). Experimental/spectroscopic evidence of polaron physics in ceria—good external validation for the narrative. ([American Chemical Society Publications](https://pubs.acs.org/doi/10.1021/acs.jpclett.0c01590?utm_source=chatgpt.com "Ultrafast Formation of Small Polarons and the Optical Gap in ..."))
    

### Methodological limits: finite size/time; long-wavelength truncation; rNEMD nonlinearity

- P. K. Schelling, S. R. Phillpot & P. Keblinski, “Comparison of atomic-level simulation methods for computing thermal conductivity,” _Phys. Rev. B_ **65**, 144306 (2002). The canonical GK↔NEMD comparison; discusses finite-size/time effects and conditions for agreement—exactly supports your convergence gates. ([AIP Publishing](https://pubs.aip.org/aip/apl/article/125/17/172202/3318140/Impact-of-classical-statistics-on-thermal?utm_source=chatgpt.com "Impact of classical statistics on thermal conductivity ..."))
    
- S. Shenogina _et al._, “Finite-size effects in determination of thermal conductivities…,” _J. Heat Transfer_ **126**, 577 (2004). Focused analysis of size dependence and missing long-qq content; handy for your “GK can be suppressed in small boxes” line. ([ASME Digital Collection](https://asmedigitalcollection.asme.org/heattransfer/article/126/4/577/464249/Finite-Size-Effects-in-Determination-of-Thermal?utm_source=chatgpt.com "Finite Size Effects in Determination of Thermal Conductivities"))
    
- Z. Liang & P. Keblinski, “Finite-size effects on interfacial thermal resistance…,” _Phys. Rev. B_ **90**, 075411 (2014). While interfacial, the treatment of ballistic transport and domain-length artifacts is a good citation for rNEMD “too-strong driving/too-short box” warnings. ([Physical Review Links](https://link.aps.org/doi/10.1103/PhysRevB.90.075411?utm_source=chatgpt.com "Finite-size effects on molecular dynamics interfacial thermal ..."))
    

### Boundary conditions in electrostatics (PPPM/Ewald “tin-foil”)

- J. E. Roberts & J. Schnitker, _J. Chem. Phys._ **101**, 5024 (1994). Explicitly states that the usual Ewald sum corresponds to conducting (“tin-foil”) boundary conditions—use this when you note your bulk-limit assumption. ([Deep Blue Repositories](https://deepblue.lib.umich.edu/bitstream/handle/2027.42/70827/JCPSA6-101-6-5024-1.pdf%3Bsequence%3D2?utm_source=chatgpt.com "How the unit cell surface charge distribution affects the ..."))
    
- S. W. de Leeuw, J. W. Perram & E. R. Smith, _Proc. R. Soc. A_ **373**, 27 (1980). Classic series on lattice-sum electrostatics and boundary terms; authoritative background for your PPPM remark. ([Astrophysics Data System](https://ui.adsabs.harvard.edu/abs/1980RSPSA.373...27D/abstract?utm_source=chatgpt.com "Simulation of Electrostatic Systems in Periodic Boundary ..."))
    
- O. Andreussi _et al._, _Phys. Rev. B_ **90**, 245101 (2014). Modern discussion of electrostatic boundary conditions under PBC; useful one-liner about “tin-foil” vs. other choices. ([Physical Review Links](https://link.aps.org/doi/10.1103/PhysRevB.90.245101?utm_source=chatgpt.com "Electrostatics of solvated systems in periodic boundary ..."))
    

### Comparing to experiments: porosity/microstructure corrections (Maxwell–Eucken, etc.)

- J. Wang _et al._, “A unifying equation for effective thermal conductivity models,” _Int. J. Heat Mass Transf._ **49**, 3075 (2006). Summarizes Series/Parallel and both Maxwell–Eucken forms with limits/assumptions—good when you explain why pellet data sit below bulk and why direct mapping needs a microstructural model. ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0017931006001293?utm_source=chatgpt.com "A new approach to modelling the effective thermal ..."))
    
- M. Zhang _et al._, “Influence of pore distribution on equivalent thermal conductivity…,” _Ceram. Int._ **45**, 18408 (2019). Shows deviations from simple ME models when pore morphology matters; a nice cautionary anchor. ([ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0272884218319011?utm_source=chatgpt.com "Influence of pore distribution on the equivalent thermal ..."))
    
- A. Rai _et al._, “Conduction heat transfer through porous materials,” OSTI (2021). Review-style discussion of EMA limits; quotable for “EMA with only kk and porosity may be insufficient.” ([OSTI](https://www.osti.gov/servlets/purl/1813261?utm_source=chatgpt.com "Conduction Heat Transfer through Solid in Porous Materials"))
    

### What MD can (and cannot) decompose without extra analysis

- J. A. Thomas _et al._, “Predicting phonon dispersion relations and lifetimes from the spectral energy density,” _Phys. Rev. B_ **81**, 081411 (2010) (+ open repository). The go-to SED method you can cite when you say “modal breakdown requires separate tools.” ([Physical Review Links](https://link.aps.org/doi/10.1103/PhysRevB.81.081411?utm_source=chatgpt.com "Predicting phonon dispersion relations and lifetimes from the ..."))
    
- T. Feng & X. Ruan, “Prediction of spectral phonon mean free path and thermal conductivity…,” _J. Nanomater._ **2014**, 206370. Concise review of BTE vs. MD-based approaches; great for framing what extra assumptions SED/BTE bring. ([feng.mech.utah.edu](https://feng.mech.utah.edu/wp-content/uploads/sites/152/2025/07/Feng-Ruan-2014-Prediction-of-Spectral-Phonon-Mean-Free-Path-and-Thermal-Conductivity-with-Applications-to-Thermoelectrics-and-Therm.pdf?utm_source=chatgpt.com "Review Article Prediction of Spectral Phonon Mean Free ..."))
    

---

If you want, I can turn these into a one-page bibliography with pull-quotes keyed to your manuscript margins (e.g., a margin note next to “tin-foil boundary” with the exact line from Roberts & Schnitker), plus a tiny appendix explaining how to apply a Maxwell–Eucken porosity correction when you later compare your bulk kk to pellet data.


