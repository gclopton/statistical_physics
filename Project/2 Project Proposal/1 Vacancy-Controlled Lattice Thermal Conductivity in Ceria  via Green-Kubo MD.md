



**_How does lattice thermal conductivity $\kappa(T)$ change with oxygen-vacancy fraction $x$ and temperature, and how sensitive are the conclusions to the chosen interatomic potential?_**

**What we Expect to See from the Model:** More oxygen vacancies lower $\kappa$ by strengthening phonon-defect scattering; the relative reduction grows with temperature as the dominant mean free paths shrink.

We will calculate $\kappa$ for a perfect cerium oxide cell, as well as cerium oxide cells with different oxygen-vacancy fractions.

For each system, we will relax the volume in the NPT ensemble, then stabilize the temperature in NVT. Finally we'll run NVE for production. During the production run, we compute the microscopic heat flux in LAMMPS, accumulate the heat-current autocorrelation function (HCACF), and integrate it until the conductivity integral reaches a stable plateau.





1. [Atomistic and experimental study on thermal conductivity of bulk and porous cerium dioxide](https://pmc.ncbi.nlm.nih.gov/articles/PMC6474893/?)
2. [A critical assessment of interatomic potentials for ceria with application to its elastic properties](https://www.sciencedirect.com/science/article/pii/S0167273810000986?)
3. [Molecular Dynamics simulations of reduced CeO2: bulk and surface](https://www.diva-portal.org/smash/record.jsf?pid=diva2%3A95716&)
4. [Toward a Consistent Prediction of Defect Chemistry in CeO2](https://pmc.ncbi.nlm.nih.gov/articles/PMC9835833/?)



# Atomistic and Experimental Study on Thermal Conductivity of Bulk and Porous Cerium Dioxide



---

**Overview:**

- (i) a parameter-free ab-initio BTE places bulk single-crystal  near $17 \mathrm{~W} \mathrm{~m}^{-1} \mathrm{~K}^{-1}$ at 300 K ; 
- (ii) porosity and microstructure readily explain much lower experimental values; 
- (iii) optical modes are nonnegligible carriers; and 
- (iv) MFP spectra identify ~tens-of-nanometers as the key scale to engineer thermal transport in ceria. These results give you clean validation targets for stoichiometric $\mathrm{CeO}_2$ and concrete size/porosity regimes to test in GK MD.


---



The paper tackles the long-standing spread in reported room-temperature thermal conductivity $\kappa$ for ceria by combining first-principles lattice-dynamics/Boltzmann transport, classical MD, and laser-flash measurements. Using DFT/DFPT-derived force constants and a full iterative solution of the phonon Boltzmann transport equation (via ShengBTE), the authors predict a defect-free single-crystal value of $\kappa(300 \mathrm{~K})=16.71 \mathrm{~W} \mathrm{~m}^{-1} \mathrm{~K}^{-1}$ and **recover the expected decrease** of $\kappa$ with temperature as Umklapp scattering grows. They argue this provides an "upper bound," and note that past numbers ranged widely from $\sim 6$ to $19 \mathrm{~W} \mathrm{~m}^{-1} \mathrm{~K}^{-1}$-depending on density, microstructure, and methodology; their 300 K prediction aligns best with the higher-quality $\sim 14 \mathrm{~W} \mathrm{~m}^{-1} \mathrm{~K}^{-1}$ experiments on $\sim 95 \%$ dense pellets.



A central physics result is the mode-resolved breakdown: optical phonons contribute a substantial ~27\% of the total $\kappa$ at 300 K (with $\sim 13.1 \%$ from one transverse optical branch, $\sim 13.7 \%$ from the triply degenerate Raman-active branch, and ~1\% from the longitudinal optical branch). This challenges the common heuristic that optics are negligible in oxides. 


To bridge ab-initio "perfect crystal" predictions to real samples, the team measures porous, spark-plasma-sintered pellets by laser flash and models porosity explicitly with equilibrium Green-Kubo MD using a many-body EAM potential. Pure-crystal MD at 300 K gives $\kappa \approx 18-19 \mathrm{~W} \mathrm{~m}^{-1} \mathrm{~K}^{-1}$ (consistent across 6-96k atoms), but introducing ~ 5\% porosity suppresses $\kappa$ markedly, in line with selected experiments; the magnitude of the drop depends on pore spacing relative to phonon mean free paths (smaller cells with closer pores show a larger reduction). They also rationalize which porosity-correction formulas are appropriate, arguing that accurate theory helps choose corrections rather than fitting ad hoc.


Nanostructuring guidance comes from cumulative mean-free-path (MFP) analysis: phonons with MFPs below ~65 nm carry ~80\% of the heat, implying that feature sizes below this threshold will strongly depress $\kappa$. Thin-film calculations show both in-plane and cross-plane $\kappa$ fall with decreasing thickness and with increasing temperature, framing expectations for coatings and micro-/nano-devices that use ceria.







