
We first fix a single thermodynamic state-temperature, composition (EC:LiPF ${ }_6$ salt fraction), and density -representing bulk electrolyte conditions. A periodic liquid box is prepared by packing EC and LiPF ${ }_6$ to the target composition, then equilibrating to the target density and temperature. This preparation is performed once and archived as a common starting configuration. From that configuration we generate matched replicas by assigning velocity seeds; each replica index has the same initial coordinates and velocity seed across all model variants so that differences in outcomes can be attributed to the force-field choice rather than initial conditions.

ReaxFF is then exercised in four "paired" experiments: parameter set A with QEq, A with ACKS2, B with QEq, and B with ACKS2. In every case we use the same neighbor settings, integration scheme (velocityVerlet), a sub-femtosecond timestep appropriate for reactive bond-order dynamics, and the same thermostat to maintain the target temperature during the short reactive window. Prior to production, each replica is briefly relaxed under the chosen scheme so that charges and local structure settle; production trajectories then proceed without changing any settings other than the force-field lever being tested.

Because we care about the very first chemical change, sampling emphasizes many short, independent trajectories rather than a few long ones. Each model variant is run for $O\left(10^1-10^2\right)$ replicas over a window long enough to register first reactions with reasonable probability (tens to hundreds of picoseconds depending on temperature), with the understanding that some trajectories will remain unreacted within the window. To keep runs comparable, wall-clock length, writeout cadence, and analysis cadence (for bond order, charges, and coordinates) are identical across variants.

Event detection is based on bond-order time series. We define bond formation/breaking as threshold crossings with hysteresis and minimum-duration criteria to suppress flicker. The first elementary event in a trajectory is the earliest time at which a chemically interpretable change occurs-e.g., EC ring opening (C-O scission), solvent deprotonation, or $\mathrm{PF}_6^{-}$dissociation. Once an event is flagged, we freeze a short post-event window for product typing. Product identities are assigned from connectivity graphs constructed by binarizing the bond-order matrix; this yields counts for CO vs. $\mathrm{CO}_2$, small oligomers, and, in $\mathrm{LiPF}_6$ systems, the incidence of $\mathrm{Li}-\mathrm{F}$ contacts. $\mathrm{Li}-\mathrm{F}$ contact lifetimes are measured by monitoring minimum-image Li-F distances, with a cutoff taken from the first minimum of the Li-F radial distribution function computed from the same trajectories.

For each variant, the distribution of first-event times is estimated with a Kaplan-Meier curve that correctly handles right-censored runs (those with no event before the cutoff). Median times, hazard rates over the early window, and confidence intervals are computed by nonparametric bootstrap over replicas. Differences between variants are assessed with logrank tests on the survival curves and with bootstrap confidence intervals on pairwise median time ratios. Branching fractions among competing first events are estimated from the same replica pool; uncertainty is quantified with multinomial bootstrap, and cross-variant differences are checked with exact tests on contingency tables.

To situate reactivity in its structural environment, we compute radial distribution functions and running coordination numbers for $\mathrm{Li}-\mathrm{O}, \mathrm{Li}-\mathrm{F}$, and $\mathrm{C}-\mathrm{O}$ pairs over the pre-event portions of each trajectory, pooling replicas within a variant. Where helpful, we condition these statistics on whether a replica ultimately exhibits a given first event, to identify structural predictors of pathway choice. We also track charge statistics (e.g., Li and $\mathrm{PF}_6$ partial-charge distributions) to make the impact of QEq vs. ACKS2 explicit.

Throughout, controls and diagnostics ensure a like-for-like comparison. Timestep and neighbor settings are validated with short NVE checks for energy drift. Thermostat parameters are held fixed across variants, charge-equilibration tolerances are matched so the linear solves are of comparable accuracy for QEq and ACKS2, and replicate-level randomization (common seeds across variants) is used to reduce between-variant noise.

This design yields three families of outputs per variant: (i) survival curves and summary statistics for first-event times, (ii) branching fractions across mechanistic pathways with uncertainties, and (iii) structural context (RDFs, coordinations, Li-F contact lifetimes and frequencies). Comparing these across the four experiments directly answers whether EC decomposition in ReaxFF is robust to parameter-set choice and to the charge-equilibration scheme-or whether mechanistic conclusions are conditional and require qualification and higher-level validation.



# Sources

## Background on ReaxFF + Why Parameter/Charge-Scheme Choices Matter

1.) **ReaxFF formalism and transferability:** bond-order mechanics and typical analysis choices are covered in the original hydrocarbon ReaxFF paper and a broad review. ReaxFF encodes chemistry through bond order + fitted parameter sets. [ReaxFF:  A Reactive Force Field for Hydrocarbons](https://pubs.acs.org/doi/10.1021/jp004368u)
2.) QEq (electronegativity equalization) is the classic flexible-charge scheme used in ReaxFF. It is known to over-delocalize charge in extended systems ("metallicity problem"), motivating ACKS2/SQE-style fixes that penalize long-range charge transfer. This motivates testing QEq vs. ACKS2. [Charge Equilibration For Molecular Dynamics Schemes](https://pubs.acs.org/doi/10.1021/j100161a070)
3.) ACKS2 as a drop-in alternative in LAMMPS ReaxFF (and its expected impact on polarization/fragment neutrality). There's some practical notes on using `fix acks2/reaxff` that might be useful for us. We can use these to ground the implementation knob and to argue it can change early-time reactivity. [fix acks2/reaxff command](https://docs.lammps.org/fix_acks2_reaxff.html)
4.) Comparative studies that explicitly evaluate QEq vs. ACKS2 (and related SQE/EEM families) support the claim that charge-assignment method shifts dipoles, fields, and response-hence kinetics.  ReaxFF for Li-ion electrolyte chemistry (including EC reduction) and eReaxFF+ACKS2 usage in battery contexts: good precedent that our chemistry target is standard for ReaxFF and that ACKS2 is used in practice. [Reductive Decomposition Reactions of Ethylene Carbonate by Explicit Electron Transfer from Lithium: An eReaxFF Molecular Dynamics Study](https://pubs.acs.org/doi/10.1021/acs.jpcc.6b08688)


## EC/LiPF6 Early Chemistry (What "First Events" Might Be)
1.) Ab-initio and mechanism studies identifying EC ring-opening, one- vs two-electron routes, and surface effects give us canonical "first events" to look for. We can use these to motivate ring-opening, deprotonation, and salt-involved pathways as targets. [Two-electron reduction of ethylene carbonate: A quantum chemistry re-examination of mechanisms](https://www.sciencedirect.com/science/article/pii/S0009261412009384)
2.) For $LiPF_{6}$, modern DFT work clarifies chemical decomposition channels (PF5/POF3 chemistry, LiF emergence) without requiring water-useful to justify Li-F contact tracking as a proxy for LiF precursors. [Elementary Decomposition Mechanisms of Lithium Hexafluorophosphate in Battery Electrolytes and Interphases](https://pubs.acs.org/doi/10.1021/acsenergylett.2c02351?)


## Trajectory Design: Many Short Replicas to Capture First Events

General rare-events/first-passage-time treatments in simulation justify emphasizing many short, independent trajectories and analyzing first-event time distributions. Peters' text is the cleanest single citation; recent work on short-time kinetics (e.g., infrequent metadynamics papers) reinforces the "many short" logic.


## MD Plubming; Integration, Thermostatting, and Validation

Velocity-Verlet and Nosé-Hoover are the standard for reactive MD at sub-fs timesteps. Allen & Tildesley is a good textbook anchor for RDFs/coordination numbers. [A computer simulation method for the calculation of equilibrium constants for the formation of physical clusters of molecules: Application to small water clusters](https://pubs.aip.org/aip/jcp/article/76/1/637/397800/A-computer-simulation-method-for-the-calculation)



## Event Detection from Bond-Order Time Series; Connectivity Graphs; Bond-Order Cutoffs

1.) Methodologies that use ReaxFF bond order with thresholds/hysteresis to detect reactions and then assign products by graph connectivity are well-established. VARxMD (reaction parsing from bond orders) and documentation specifying standard bond-order cutoffs make our detection pipeline by-the-book. [Reaction analysis and visualization of ReaxFF molecular dynamics simulations](https://pubmed.ncbi.nlm.nih.gov/25064439/)
2.) Bond-order time-series-based event detection (explicitly discussed in the JCTC piece) lets us justify hysteresis/minimum-duration rules to suppress flicker. [Bond-Order Time Series Analysis for Detecting Reaction Events in _Ab Initio_ Molecular Dynamics Simulations](https://pubs.acs.org/doi/10.1021/acs.jctc.9b01039)
3.) LAMMPS community guidance and ReaxFF engine docs back the claim that cutoffs affect analysis only (not forces), which is useful when you state your post-processing thresholds. [ReaxFF Bond Order Cutoffs](https://matsci.org/t/reaxff-bond-order-cutoffs/16360)


## Structural Context: RDFs, Coordination Numbers, Li-F Contact Lifetimes

Standard practice is to define coordination by integrating $4 \pi \rho \int_0^{r_{\min }} r^2 g(r) d r$ to the first minimum of the RDF; cite a canonical text/tutorial and an example application. This supports your Li-O/Li-F/C-O coordination metrics and the "first-minimum as cutoff" rule used for Li-F contacts. [Computer Simulation of Liquids by M.P. Allen and D.J. Tildesley](https://datagrid.hu/boda/Boda-sajat/Rush/Books/Allen-Tildesley.pdf)


## Survival Analysis for time-to-first-event with Censoring

1.) Use Kaplan-Meier for right-censored first-event times and the Mantel (log-rank) test for cross-variant comparisons. A modern tutorial plus the originals cover both method and interpretation. [Nonparametric Estimation from Incomplete Observations](https://web.stanford.edu/~lutian/coursepdf/KMpaper.pdf)

2.) Bootstrap CIs for medians/hazard summaries (Efron-Tibshirani) and exact tests for pathway-branching contingency tables (Fisher) justify your uncertainty quantification and pathway-fraction comparisons. [An Introduction to Bootstrap](https://www.taylorfrancis.com/books/mono/10.1201/9780429246593/introduction-bootstrap-bradley-efron-tibshirani)


