

# 1. The XYZ Format

`ethane.xzy` is in a format called XZY format, and it carries the identities of the atoms and Cartesian coordinates.


``` title:ethane.xyz
8
Ethane, length in nm
C 0.000 0.000 0.000        
H 0.000 0.000 0.110  
H 0.110 0.000 0.000  
H -0.110 0.000 0.000 
C 0.000 0.150 0.000  
H 0.000 0.150 0.110  
H 0.000 0.150 -0.110 
H 0.110 0.150 0.000  
```


Here, the first line is the atom count. `8` tells the parser there are eight atoms in this frame. Many tools rely on this exact integer to know how many coordinate lines to read next.

The second line in the file is a free-form comment. It can say anything you want. Here, the line tells us that the units of the cartesian coordinates are in nanometers. This is important to know because most chemistry files default to ångstrom (Å). 


The subsequent lines provide the positions of the atoms in terms of Cartesian coordinates. There is one atom per line. Each line has four whitespace separated fields: `<element> x y z`, where `<element>` is the chemical symbol, and `xyz` are the Cartesian coordinates in nm.

The ordering defines the atom indices that will be used everywhere else, (e.g., bonds like (0,4) for C-C).


``` title:ethane.xyz
8
Ethane, length in nm
C 0.000 0.000 0.000        # atom 0
H 0.000 0.000 0.110        # atom 1
H 0.110 0.000 0.000        # atom 2
H -0.110 0.000 0.000       # atom 3
C 0.000 0.150 0.000        # atom 4
H 0.000 0.150 0.110        # atom 5
H 0.000 0.150 -0.110       # atom 6
H 0.110 0.150 0.000        # atom 7
```


The geometry looks like two carbons along +y separated by 1.50 nm (1.50 Å), with hydrogens placed roughly 0.110 nm (1.10 Å) from their respective carbons along orthogonal axes. Looking at the molecule in OVITO, we can see that it is in a high-energy eclipsed conformation. Normally, in an ethane molecule, the azimuthal angle separating the hydrogen atoms on a single carbon should be 120 degrees. The HCH bond angles should be 109.5 degrees. The C-H bond length should be 1.09 Å (0.109 nm).


![[Pasted image 20250919142454.png|500]]







# 2. Defining Molecular Topology and Force Fields with OPLS (Optimized Potentials for Liquid Simulations) Nomenclature
This is exactly the information the “OPLS nomenclature” provides: the names tie your topology patterns to the correct parameters and the bookkeeping rules tell you which interactions to include and how.

OPLS (Optimized Potentials for Liquid Simulations) Nomenclature defines the topology of a molecules by providing information about atom types (with non-bonded parameters and charges), connectivity (bonds and angles), and torsional types (the Fourier series that fixes the rotamer energetics). The OPLS notation ties the topology patterns to the correct parameters and the bookkeeping rules indicate which interactions to include and how. *==This information can be used to turn the XYZ into tuples the functions expect==*.


## 2.1 Atom Typing

In OPLS-AA, an aliphatic $\mathrm{sp}^3$ carbon is `CT` and the hydrogens bound to that carbon are `HC`. Those names are force-field types (not just elements): the type carries the Lennard-Jones size/energy, a partial charge, and which bonded/angle/dihedral parameters to use when this atom appears in a pattern. So for ethane we convert the eight element symbols to types as: the two carbons $\rightarrow$ CT, the six hydrogens attached to them $\rightarrow$ HC. (In the code, the `nb_params` dictionary keys will be `CT` and `HC`, and those entries hold $\sigma, \varepsilon$, and $q$ for each type.)


## 2.2 Representing Bonds, Angles, and Dihedrals in OPLS Notation


### 2.2.1 Bonds

XYZ files don't encode bonds, so we have to declare them. For ethane, there is one C-C bond and six C-H bonds. Using the XYZ file's line order as indices (C0 H1 H2 H3 C4 H5 H6 H7), the bonds are ( 0,4 ), then ( 0,1$)$, ( 0,2 ), ( 0,3 ), and (4,5), (4,6), (4,7). 


![[Pasted image 20250919174116.png|500]]


In OPLS, a bond type is named by the ordered pair of atom types on that bond, e.g., `CT-CT` for the central carbon-carbon bond and `CT-HC` for the six carbon-hydrogen bonds. In the code, each bond tuple generally carries the type label and the two indices, e.g., ( `CTCT`, 0, 4) and six ('CT-HC', $\mathrm{i}, \mathrm{j}$ ) entries. The bond type label is what is used to look up the harmonic parameters $r_0$ and $k_r$ in `bond_params`.


### 2.2.2 Angles

Angles are represented as triples of bonded atoms with the middle atom as the vertex. In ethane, each carbon is tetrahedral, so around C0 you have $\mathrm{H}-\mathrm{C}-\mathrm{H}$ angles for the three hydrogen pairs and $\mathrm{H}-\mathrm{C}-\mathrm{C}$ angles for each hydrogen with the neighboring carbon. Concretely, at C0 the HCH angles are ( $1,0,2$ ), ( $1,0,3$ ), ( $2,0,3$ ), and the HCC angles are ( $1,0,4$ ), ( $2,0,4$ ), ( $3,0,4$ ). At C4 the triples are: ( $5,4,7$ ), ( $5,4,6$ ), ( $7,4,6$ ) for HCH, and ( $5,4,0$ ), (6,4,0), (7,4,0) for HCC. In OPLS the angle type is a three-type string `HC-CT-HC` for HCH and `HC-CT-CT` (or similarly `CT-CT-HC`) for HCC. Each such type points to a **harmonic equilibrium angle** $\theta_0$ and **stiffness** $k_\theta$ in `angle_params` . An important thing to remember is that the ordering must stay consistent so that dictionary lookups are unambiguous. One approach is to store both permutations that are physically the same. Another is to normalize the type string (e.g., always put the central type in the middle and sort the flanking types when they're identical).


### 2.2.3 Dihedrals

A proper dihedral is a sequence of four atoms i-j-k-I where i-j-k and j-k-I are bonded. In ethane, the chemically meaningful torsions are all HC-CT-CT-HC connecting one hydrogen on C0 through the $\mathrm{C}-\mathrm{C}$ bond to one hydrogen on C4 . There are $3 \times 3=9$ such dihedrals: (1,0,4,5), (1,0,4,6), ..., (3,0,4,7). In OPLS, a dihedral type is written by its four atom types, e.g., `HC-CT-CT-HC`. That type is associated with a Fourier series with up to four coefficients $V_1, \ldots, V_4$ that produce the staggered minima and eclipsed energy barriers for ethane. In general $V_{n}$ represent the amplitudes of the first four cosine harmonics in a compact Fourier model of the torsional potential $U(\phi)$ about a bond. For OPLS, the standard functional form is the following truncated cosine series:


$$
U_{\text {dihedral }}(\phi)=\frac{1}{2} V_1(1+\cos \phi)+\frac{1}{2} V_2(1-\cos 2 \phi)+\frac{1}{2} V_3(1+\cos 3 \phi)+\frac{1}{2} V_4(1-\cos 4 \phi)
$$

The index on $V_{n}$ is the periodicity: $n=1$ varies once per $360^{\circ}$, $n=2$ twice, $n=3$ three times, and $n=4$ four times. OPLS caps the expansion at $n=4$ as a pragmatic fit; it's flexible enough to reproduce the known barrier shapes for organic torsions without overfitting or making force/torque calculations messy.

There is an important thing to note in the case of ethane (HC-CT-CT-HC). Here $V_{3}$ is the threefold term that captures the tetrahedral symmetry around a single C-C bond. It represents how much higher the eclipsed states sit than the staggered ones because of $H \cdots H$ repulsion and hyperconjugation loss when bonds eclipse. So, considering the special angles; (1) at eclipsed conformations (e.g., $\phi=0^{\circ}, 120^{\circ}, 240^{\circ}$ ), $\cos (3 \phi)=+1$, so the  contribution is $\frac{1}{2} V_3(1+1)=V_3$ – it raises the energy. At staggered conformations $\left(\phi=60^{\circ}, 180^{\circ}, 300^{\circ}\right), \cos (3 \phi)=-1$, so the $V_3$ contribution is $\frac{1}{2} V_3(1-1)=0-$ the term vanishes. So in ethane, when $V_{3} > 0$ , the three eclipsed states are penalized and the three staggered states become minima at $\phi= 60^{\circ}, 180^{\circ}, 300^{\circ}$. bonds eclipse. If $V_3=1$, the $n=3$ term contributes exactly +1 at eclipsed and 0 at staggered; if $V_3=2$, it contributes +2 at eclipsed, and so on. The total torsion energy at any $\phi$ is the sum of all the $V_n$ pieces evaluated there.


The other terms fine-tune details that the pure threefold symmetry can't handle. The $V_1$ (onefold term) can tilt the landscape so the barrier at $0^{\circ}$ differs from features near . The $V_2$ (twofold term) can split trans vs. gauche when substitution breaks symmetry. For ethane though, this term is usually small or zero because the three staggered states are equivalent. $V_4$ is rarely needed for simple alkanes, but can sharpen eclipse profiles.


> [!NOTE] **Force field family" and "unit system**
> Different families (OPLS-AA, AMBER, CHARMM) and different unit conventions tabulate numbers that aren't directly comparable unless the functional form and units are checked. Some families also prefer the Ryckaert-Bellemans (RB) polynomial form $\sum_{m=0}^5 C_m \cos ^m(\phi)$ rather than the OPLS cosine series. One can convert between OPLS $V_n$ and RB $C_m$ one-to-one; the converted parameters will still embed a positive "threefold" influence for alkanes.



The `dihedral_params[HC-CT-CT-HC]` entry in the code holds a dict of the form `{ 'V1': ..., 'V2': ..., 'V3': ..., 'V4': ...}` (where zeros stand in place of the omitted terms).




## 2.3 Non-Bonded Interactions in the OPLS System

OPLS-AA includes non-bonded interactions (Lennard-Jones and Coulomb) for atom pairs that are separated by more than an angle, but it applies special rules:

**Rule #1: 1-2 (bonded) and 1-3 (angle) pairs are excluded from non-bonded interactions entirely.**
For ethane, the 1-2 bonded pairs are the connectivity edges: (0,4),(0,1),(0,2),(0,3),(4,5),(4,6),(4,7). For the angles (1-3 pairs): take every angle $i-j-k$ and mark the endpoints $(i, k)$ as 1-3. Around each carbon we get the three $\mathrm{H}-\mathrm{H}$ pairs and the three $\mathrm{H}-\mathrm{C}$ (across the $\mathrm{C}-\mathrm{C}$ ) pairs.
- At $\mathrm{CO}(\mathrm{j}=0)$ : $\mathrm{H}-\mathrm{H}$ across $\mathrm{CO} \rightarrow(1,2),(1,3),(2,3) ; \mathrm{H}-\mathrm{C}$ across $\mathrm{COC} 4 \rightarrow (1,4),(2,4),(3,4)$.
- At C4 $(j=4)$ : H -H across C4 $\rightarrow(5,6),(5,7),(6,7) ; \mathrm{H}-\mathrm{C}$ across C4C0 $\rightarrow (5,0),(6,0),(7,0)$.

Collecting all the pairs, we have:

$$
(1,2),(1,3),(2,3),(1,4),(2,4),(3,4),(5,6),(5,7),(6,7),(5,0),(6,0),(7,0) .
$$

These twelve pairs are also excluded from non-bonded.


**Rule #2: 1-4 pairs (atoms at the ends of a single torsion i-j-k-l) are included but scaled by 0.5 for both LJ and electrostatics.**

For torsions (1-4 pairs, the endpoints of i-j-k-I chains for ethane), the chemically relevant dihedrals are HC-CT-CT-HC, i.e., hydrogens across the central C-C bond. Every H on C0 paired with every H on C4 is a 1-4 pair:

$$
(1,5),(1,6),(1,7), \quad(2,5),(2,6),(2,7), \quad(3,5),(3,6),(3,7) .
$$


These nine pairs are included with 0.5 scaling for both LJ and Coulomb in OPLS-AA.


**Rule #3: All other pairs (1-5 and beyond) are are included unscaled.**
For ethane, this rule is irrelevant.


In the code, the topology is separated from the physics. It starts by precomputing the three disjoint sets: `excluded_12`, `excluded_13`, and `scaled_14` from the bonds/angles/and dihedrals. Then it loops over all unordered pairs $i < j$: 
- if $(i, j) \in$ `excluded_12` $\cup$ `excluded_13`: skip
- else if $(i,j) \in$ `scaled_14`: pass `scale_lj=0.5, scale_coul=0.5` to the `nonbond_potential`.
- else: pass `scale_lj=1.0, scale_coul=1.0`.

*==Be careful! Do not double-apply scaling if you later add group-based 1-4 scaling elsewhere. Keep it centralized in the non-bonded pair evaluation==*.


> [!NOTE] **A useful sanity check**
> With eight atoms there are $\binom{8}{2}=28$ unique pairs. Subtract 7 (12) and 12 (1-3) leaves 9 , which are exactly the nine 1-4 pairs above. Ethane has no 1-5 or longer pairs, so in this molecule every non-bonded interaction that survives is a scaled 1-4.


### 2.3.1 Non-Bonded Parameters are Attached to Atom Types

OPLS assigns each atom an atom type (e.g., `CT` for an $\mathrm{sp}^3$ aliphatic carbon; `HC` for a hydrogen attached to `CT`). The type carries its own non-bonded parameters:
- a partial charge $q$ (in units of the elementary charge $e$ ),
- a Lennard-Jones size $\sigma$ (a distance scale), and
- a Lennard-Jones well depth $\varepsilon$ (an energy scale).

We store ( $q, \sigma, \varepsilon$ ) once per type and derive cross-type parameters on the fly using a combining rule (next section). Bonds/angles/dihedrals have their own bonded parameters (harmonic $r_0, k_r ; \theta_0, k_\theta ; V_n$ for torsions), but they do not define non-bonded behavior.


```Python
nb_params = {
    'CT': {'q': -0.18, 'sigma': 0.350, 'epsilon': 0.276, 'coulomb_factor': 138.935456},
    'HC': {'q':  0.06, 'sigma': 0.250, 'epsilon': 0.125, 'coulomb_factor': 138.935456},
}
```


### 2.3.2 Geometric Combining Rules in OPLS

When two atoms $i$ and $j$ interact non-bondedly, you must build a pair Lennard-Jones parameter set $\left(\sigma_{i j}, \varepsilon_{i j}\right)$ from the two type sets $\left(\sigma_i, \varepsilon_i\right)$ and $\left(\sigma_j, \varepsilon_j\right)$. OPLS uses the geometric-geometric rule:

$$
\sigma_{i j}=\sqrt{\sigma_i \sigma_j}, \quad \varepsilon_{i j}=\sqrt{\varepsilon_i \varepsilon_j} .
$$


That choice is part of the force field's design (other families sometimes use Lorentz-Berthelot: $\sigma_{i j}=\left(\sigma_i+\sigma_j\right) / 2, \varepsilon_{i j}=\sqrt{\varepsilon_i \varepsilon_j}$ ). Since the code is implementing OPLS, it uses geometric for both.


#### Example
Concretely, suppose:
- CT: $\sigma_{\mathrm{CT}}=0.350 \mathrm{~nm}, \varepsilon_{\mathrm{CT}}=0.276 \mathrm{~kJ} \mathrm{~mol}^{-1}(\approx 0.066 \mathrm{kcal} / \mathrm{mol})$,
- $\mathrm{HC}: \sigma_{\mathrm{HC}}=0.250 \mathrm{~nm}, \varepsilon_{\mathrm{HC}}=0.125 \mathrm{~kJ} \mathrm{~mol}^{-1}(\approx 0.030 \mathrm{kcal} / \mathrm{mol})$.

Then for a CT-HC pair we'll use

$$
\sigma_{\mathrm{CT}-\mathrm{HC}}=\sqrt{0.350 \times 0.250} \approx 0.296 \mathrm{~nm}, \quad \varepsilon_{\mathrm{CT}-\mathrm{HC}}=\sqrt{0.276 \times 0.125} \approx 0.186 \mathrm{~kJ} \mathrm{~mol}^{-1} .
$$

Note: the code precompute all unique pairs of types once and caches them.


### 2.3.3 The Leonard-Jones Interaction Potential

With $\sigma_{i j}$ and $\varepsilon_{i j}$ in hand, the Lennard-Jones pair potential is

$$
U_{i j}^{\mathrm{LJ}}(r)=4 e_{i j}\left[\left(\frac{\sigma_{i j}}{r}\right)^{12}-\left(\frac{\sigma_{i j}}{r}\right)^6\right] .
$$


The shape (repulsive core and attractive well) comes from the exponents 12/6; the scale and position of the minimum come entirely from $\varepsilon_{i j}$ and $\sigma_{i j}$, which is why the mixing rule is so important. Using the wrong rule shifts well depths and equilibrium separations for cross interactions and breaks the fit the force field was parameterized to.


### 2.3.4 The Coulomb Potential

The coulomb potential uses the partial charges on the two atoms:

$$
U_{i j}^{\text {Coul }}(r)=k_e \frac{q_i q_j}{\varepsilon_r r} .
$$


Condensed-phase MD almost always sets the relative dielectric $\varepsilon_r=1$ and accounts for screening explicitly with solvent and long-range summation (PME/PPPM). *==In a GROMACS-style unit system (distances in nm, energies in kJ/mol, charges in e), you'll take==*

$$
k_e \approx 138.935456 \frac{\mathrm{~kJ} \mathrm{~nm}}{\mathrm{~mol} \mathrm{e}^2} .
$$


The code can carries a single coulomb_factor equal to  and multiply it by $q_i q_j / r$ for each pair. Note: When adding OPLS 1-4 scaling, the code must multiply the LJ and Coulomb pair energies (and forces) by 0.5 for those 1-4 pairs.



The Leonard-Jones interaction is $4 \varepsilon_{i j}\left[\left(\sigma_{i j} / r\right)^{12}-\left(\sigma_{i j} / r\right)^6\right]$ and Coulomb is $k_e q_i q_j / r$, with the code carrying a `coulomb_factor` equal to $k_e / \varepsilon_r$ *==expressed in whatever unit system you choose (e.g., $\mathrm{kJ} \cdot \mathrm{nm} \cdot \mathrm{mol}^{-1} \cdot \mathrm{e}^{-2}$ in GROMACS-style units)==*. 



## 2.4 Units and Consistency

Because the XYZ file is in nanometers, but the code implements energies in $\mathrm{kJ} / \mathrm{mol}$, it's convenient to keep the whole force field in the same GROMACS-style unit system: distances in nm, energies in kJ/mol, charges in elementary charge, Coulomb factor $k_e \approx 138.935456 \mathrm{~kJ} \mathrm{~nm} \mathrm{~mol}^{-1} \mathrm{e}^{-2}$. 


> [!Warning]
> Published OPLS numbers are sometimes tabulated in $\mathrm{kcal} / \mathrm{mol}$ and Å, which need to be converted once they are loaded to keep everything internally consistent.( $1 Å=0.1 \mathrm{~nm}, 1 \mathrm{kcal} / \mathrm{mol} \approx 4.184 \mathrm{~kJ} / \mathrm{mol}$.




## 2.4 Defining the Molecular Topology: Python Implementation


The `define_molecule` function reads the `.xyz` file and assign type strings `['CT','HC','HC','HC','CT','HC','HC','HC']` in the file order, then populates the following lists:

1.) `bonds = [('CT-CT',0,4), ('CT-HC',0,1), ('CT-HC',0,2), ('CT-HC',0,3), ('CT-HC',4,5), ('CT-HC',4,6), ('CT-HC',4,7)]`
    
2.) `angles = [('HC-CT-HC',1,0,2), (1,0,3), (2,0,3), ('HC-CT-CT',1,0,4), (2,0,4), (3,0,4), ('HC-CT-HC',5,4,6), (5,4,7), (6,4,7), ('HC-CT-CT',5,4,0), (6,4,0), (7,4,0)]` (with the type string repeated on each tuple)
    
3.) `dihedrals = [('HC-CT-CT-HC',1,0,4,5), ..., ('HC-CT-CT-HC',3,0,4,7)]` covering the 9 H–C–C–H combinations.


From there, the parameter dictionaries key off the type strings: 
- `bond_params['CT-CT'] = {'r0': ..., 'k': ...}`
- `angle_params['HC-CT-HC'] = {'th0': ..., 'k': ...}`
- `dihedral_params['HC-CT-CT-HC'] = {'V1': ..., 'V2': ..., 'V3': ..., 'V4': ...}` and 
- `nb_params['CT'] = {'q': ..., 'sigma': ..., 'epsilon': ..., 'coulomb_factor': 138.935456}` (with a similar entry for `'HC'`). 

Finally, when the non-bonded pairs are assembled, the exclusion lists are built directly from the bonds/angles/dihedrals: skip 1–2 and 1–3; include 1–4 with 0.5 scaling; include all others unscaled. 




```Python
def define_molecule(filename: str):

    with open(filename) as f:
        N = int(f.readline().strip())   
        _ = f.readline()                

        positions = np.empty((N, 3), dtype=float)
        elements  = []
        for i in range(N):
            sym, x, y, z = f.readline().split()
            positions[i] = [float(x), float(y), float(z)]
            elements.append(sym)

    types = [('CT' if e == 'C' else 'HC') for e in elements]

    C0, H1, H2, H3, C4, H5, H6, H7 = range(8)

    bonds = (
        ('CT-CT', C0, C4),
        ('CT-HC', C0, H1), ('CT-HC', C0, H2), ('CT-HC', C0, H3),
        ('CT-HC', C4, H5), ('CT-HC', C4, H6), ('CT-HC', C4, H7),
    )

    angles = (
        ('HC-CT-HC', H1, C0, H2), ('HC-CT-HC', H1, C0, H3), ('HC-CT-HC', H2, C0, H3),
        ('HC-CT-CT', H1, C0, C4), ('HC-CT-CT', H2, C0, C4), ('HC-CT-CT', H3, C0, C4),
        ('HC-CT-HC', H5, C4, H6), ('HC-CT-HC', H5, C4, H7), ('HC-CT-HC', H6, C4, H7),
        ('HC-CT-CT', H5, C4, C0), ('HC-CT-CT', H6, C4, C0), ('HC-CT-CT', H7, C4, C0),
    )

    dihedrals = (
        ('HC-CT-CT-HC', H1, C0, C4, H5), ('HC-CT-CT-HC', H1, C0, C4, H6), ('HC-CT-CT-HC', H1, C0, C4, H7),
        ('HC-CT-CT-HC', H2, C0, C4, H5), ('HC-CT-CT-HC', H2, C0, C4, H6), ('HC-CT-CT-HC', H2, C0, C4, H7),
        ('HC-CT-CT-HC', H3, C0, C4, H5), ('HC-CT-CT-HC', H3, C0, C4, H6), ('HC-CT-CT-HC', H3, C0, C4, H7),
    )

    return positions, types, bonds, angles, dihedrals

```


> [!NOTE]
> A few expectations to keep in mind as you use these returns downstream
> - The `positions` array is immediately usable for distance/angle/dihedral calculations and for building forces. 
> - The `types` list is the lookup key into your `nb_params` dictionary, which stores each type’s `q`, `sigma`, `epsilon`, and a `coulomb_factor` constant. 
> - The `bonds`, `angles`, and `dihedrals` sequences are declarative: they don’t carry any numbers—just the type string that keys into `bond_params`, `angle_params`, or `dihedral_params`, plus the **indices** so your energy/force routines know which atoms to read and where to accumulate forces.



## 2.5 Defining the Force Field

Perfect-let's drop in a clean, self-consistent OPLS-AA set in the same units you're using elsewhere:
- distances: nm
- energies: $\mathrm{kJ}^{-} \mathrm{mol}^{-1}$
- angles: radians
- charges: e
- Coulomb factor: $k_e=138.935456 \mathrm{~kJ}$ cdotpnm cdotpmol ${ }^{-1} e^{-2}$

The numbers below are the classic OPLS-AA alkane values (converted from $\AA / \mathrm{kcal}$ where needed), and they match the OPLS cosine form you're coding (the $1 / 2$ factors are already in the equation). $\mathrm{CT}=\mathrm{sp}^3$ aliphatic carbon; $\mathrm{HC}=\mathrm{its}$ hydrogen.

> [!NOTE]
> Why these choices (quick sanity check)
> - CT/HC nonbonded: classic OPLS-AA alkanes ( $\sigma: 3.50 / 2.50 \AA \rightarrow 0.350 / 0.250 \mathrm{~nm}$; $\varepsilon$ : $0.066 / 0.030 \mathrm{kcal} \rightarrow 0.276 / 0.1255 \mathrm{~kJ}$ ). Charges sum to $\sim 0$ for each $\mathrm{CH}_3(-0.18+ 3 \times 0.06=0$ ).
> - Bonds: $r_0(\mathrm{C}-\mathrm{C}) \approx 1.529 \AA, r_0(\mathrm{C}-\mathrm{H}) \approx 1.090 \AA$; force constants converted to $\mathrm{kJ} \cdot \mathrm{mol}^{-1} \cdot \mathrm{~nm}^{-2}$.
> - Angles: tetrahedral region; HC-CT-HC $\approx 107.8^{\circ}$, HC-CT-CT $\approx 110.7^{\circ}$, with standard OPLS angle stiffnesses.
> - Dihedral: positive V3 ( $\approx 0.300 \mathrm{kcal}=1.2552 \mathrm{~kJ}$ ) penalizes eclipsed and favors staggered; V1, V2, V4 kept 0 for symmetric ethane.
> 
> From here, keep your nonbonded loop doing OPLS bookkeeping: skip 1-2 \& 1-3, include 1-4 with 0.5 scaling on both LJ and Coulomb, and use geometric mixing for $\sigma / \varepsilon$.




```Python
def define_force_field():

    k_e = 138.935456  # kJ·nm·mol⁻¹·e⁻²

    # Nonbonded (per atom type)
    nb_params = {
        'CT': {'q': -0.18, 'sigma': 0.350, 'epsilon': 0.276,  'coulomb_factor': k_e},
        'HC': {'q':  0.06, 'sigma': 0.250, 'epsilon': 0.1255, 'coulomb_factor': k_e},
    }

    # Bonds: U = 1/2 k (r - r0)^2
    bond_params = {
        'CT-CT': {'r0': 0.15290, 'k': 224262.4},  # kJ·mol⁻¹·nm⁻²
        'CT-HC': {'r0': 0.10900, 'k': 284512.0},  # kJ·mol⁻¹·nm⁻²
    }

    # Angles: U = 1/2 k (θ - θ0)^2
    angle_params = {
        'HC-CT-HC': {'th0': np.deg2rad(107.8), 'k': 276.144},  # kJ·mol⁻¹·rad⁻²
        'HC-CT-CT': {'th0': np.deg2rad(110.7), 'k': 313.800},  # kJ·mol⁻¹·rad⁻²
    }

    # Proper dihedral (OPLS Fourier, kJ·mol⁻¹)
    dihedral_params = {
        'HC-CT-CT-HC': {'V1': 0.0, 'V2': 0.0, 'V3': 1.2552, 'V4': 0.0}
    }

    return nb_params, bond_params, angle_params, dihedral_params

```




# 3 Calculating the non-bonded, bonded, angle, and dihedral potentials:**


## 3.1 The Bond Potential

At the bond level we're modeling a covalent connection as a 1D spring embedded in 3D space. In OPLS-AA the functional form is the simple harmonic:

$$
U_{\text {bond }}(R)=\frac{1}{2} k\left(R-R_0\right)^2,
$$

where $R=\left\|\mathbf{r}_i-\mathbf{r}_j\right\|$ is the instantaneous inter-atomic distance, $R_0$ is the equilibrium bond length, and $k$ is the bond force constant. In your unit system $R$ and $R_0$ are in nm and $k$ is in $\mathrm{kJ} \cdot^2 \mathrm{~mol}^{-1} \cdot \mathrm{~nm}^{-2}$, so $U$ is in $\mathrm{kJ}^2 \cdot \mathrm{~mol}^{-1}$.

The force is the negative gradient of the energy with respect to the atomic positions. Because the energy depends on positions only through $R$, the chain rule gives a central scalar factor and a direction:

$$
\frac{d U}{d R}=k\left(R-R_0\right), \quad \hat{\mathbf{u}}=\frac{\mathbf{r}_{i j}}{R}=\frac{\mathbf{r}_i-\mathbf{r}_j}{\left\|\mathbf{r}_i-\mathbf{r}_j\right\|} .
$$


The pairwise forces are then equal and opposite along the bond axis,

$$
\mathrm{F}_i=-k\left(R-R_0\right) \hat{\mathrm{u}}, \quad \mathrm{~F}_j=+k\left(R-R_0\right) \hat{\mathrm{u}} .
$$



Two implementation details matter. First, guard against $R→0$ to avoid division by zero when forming $\hat{\mathrm{u}}$ (an ε clamp is fine; physically you’ll never be at exactly zero in a sane configuration). Second, be clear about the function’s input and output shapes so it snaps into your skeleton: your `compute_energy_forces` loop does `u, f = ...` and then `forces[[i, j]] += f`, which strongly suggests `bond_potential` should return a scalar energy and a `(2,3)` array of forces ordered as $\left[F_i, F_j\right]$. It’s convenient to pass the two positions as a `(2,3)` array and let the bond routine compute RRR, the unit vector, and the two forces.


```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def bond_potential(r, r0, k):
    """
    Harmonic bond potential U = 1/2 k (R - r0)^2.

    Parameters
    ----------
    r_ij : (2,3) array
        Positions of the bonded atoms [ri, rj] in nm.
    r0 : float
        Equilibrium bond length in nm.
    k : float
        Bond force constant in kJ·mol^-1·nm^-2.

    Returns
    -------
    energy : float
        Bond energy in kJ·mol^-1.
    forces : (2,3) array
        Forces on atoms i and j in kJ·mol^-1·nm^-1 (same as “force units”).
        forces[0] is on atom i, forces[1] on atom j.
    """
    ri, rj = r_ij
    rij = ri - rj
    R = np.linalg.norm(rij)
    # Safe unit vector
    if R < 1e-12:
        # No defined direction; zero force but very large energy in practice.
        # You could also pick an arbitrary axis; zero is safer here.
        uhat = np.zeros(3)
    else:
        uhat = rij / R

    dU_dR = k * (R - r0)                # kJ·mol^-1·nm^-1
    Fi = -dU_dR * uhat                  # on i
    Fj = -Fi                            # Newton's third law
    energy = 0.5 * k * (R - r0) ** 2    # kJ·mol^-1
    forces = np.vstack((Fi, Fj))
    return energy, forces

```


## 3.2 The Angle Potential

The OPLS angle treats three atoms $i-j-k$ as a hinge with a harmonic penalty on the deviation of the internal angle $\theta$ at the vertex $j$ from its equilibrium value $\theta_0$. In the unit system we're using (nm, $\mathrm{kJ}^{-} \mathrm{mol}^{-1}$, radians, e), the potential is

$$
U_{\text {angle }}(\theta)=\frac{1}{2} k\left(\theta-\theta_0\right)^2,
$$

where $k$ is in $\mathrm{kJ}^{\cdot} \mathrm{mol}^{-1} \cdot \mathrm{rad}^{-2}$ and $\theta, \theta_0$ are in radians.
To evaluate $U$ and forces from Cartesian coordinates, form the two bond vectors meeting at the vertex:

$$
\mathbf{a}=\mathbf{r}_i-\mathbf{r}_j, \quad \mathbf{b}=\mathbf{r}_k-\mathbf{r}_j .
$$


Let $a=\|\mathrm{a}\|, b=\|\mathrm{b}\|$. The angle follows from

$$
\cos \theta=\frac{\mathbf{a} \cdot \mathbf{b}}{a b}, \quad \sin \theta=\frac{\|\mathbf{a} \times \mathbf{b}\|}{a b} .
$$


Using $\sin \theta$ via the cross product is numerically more stable near $\theta \approx 0$ or $\pi$ than computing $\sqrt{1-\cos ^2 \theta}$. You still clamp both $\cos \theta$ into $[-1,1]$ and $\sin \theta$ away from zero by a tiny $\varepsilon$ to avoid division blow-ups.

Forces are the negative gradient of $U$ with respect to positions. By the chain rule,

$$
\frac{\partial U}{\partial \theta}=k\left(\theta-\theta_0\right), \quad \frac{\partial \theta}{\partial \mathrm{r}}=-\frac{1}{\sin \theta} \frac{\partial \cos \theta}{\partial \mathrm{r}} .
$$


Therefore the Cartesian forces are

$$
\mathbf{F}=-\nabla U=\frac{k\left(\theta-\theta_0\right)}{\sin \theta} \nabla \cos \theta .
$$


Using $\mathbf{a}$ and $\mathbf{b}$ as the independent vectors, the gradients of $\cos \theta$ are standard:

$$
\begin{gathered}
\frac{\partial \cos \theta}{\partial \mathbf{a}}=\frac{\mathbf{b}}{a b}-\frac{\mathbf{a} \cdot \mathbf{b}}{a^2 b} \frac{\mathbf{a}}{a}=\frac{\mathbf{b}}{a b}-\cos \theta \frac{\mathbf{a}}{a^2}, \\
\frac{\partial \cos \theta}{\partial \mathbf{b}}=\frac{\mathbf{a}}{a b}-\cos \theta \frac{\mathbf{b}}{b^2} .
\end{gathered}
$$


Let

$$
g \equiv \frac{k\left(\theta-\theta_0\right)}{\sin \theta} .
$$


Then the three forces are

$$
\begin{gathered}
\mathbf{F}_i=g\left(\frac{\mathbf{b}}{a b}-\cos \theta \frac{\mathbf{a}}{a^2}\right), \quad \mathbf{F}_k=g\left(\frac{\mathbf{a}}{a b}-\cos \theta \frac{\mathbf{b}}{b^2}\right) \\
\mathbf{F}_j=-\left(\mathbf{F}_i+\mathbf{F}_k\right)
\end{gathered}
$$


which guarantees Newton's third law and total momentum conservation. Units work out: $g$ has units $\mathrm{kJ}^{-} \mathrm{mol}^{-1} \cdot \mathrm{rad}^{-1}$ divided by rad (dimensionless) $\rightarrow \mathrm{kJ}^{-1} \cdot \mathrm{~mol}^{-1}$, multiplied by inverse lengths in the brackets to give force units $\mathrm{kJ} \cdot \mathrm{mol}^{-1} \cdot \mathrm{~nm}^{-1}$.

Two numerical subtleties matter in practice. First, clamp $a$ and $b$ away from zero with a tiny floor $\left(10^{-12} \mathrm{~nm}\right)$ so the unit vectors are well defined even if two atoms momentarily coincide. Second, clamp $\sin \theta$ with a small $\varepsilon\left(\right.$ e.g., $\left.10^{-8}\right)$ when the angle is extremely straight or bent; the analytic force truly diverges at exactly 0 or $\pi$, but a stable integrator or minimizer should never let you sit exactly there-this safeguard prevents catastrophic steps.

Here is a drop-in implementation that matches your skeleton's calling convention and return shapes: a scalar energy and a $(3,3)$ array of forces ordered as $\left[F_i, F_j, F_k\right]$.



```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def angle_potential(r_ijk, th0, k):
    """
    Harmonic angle potential U = 1/2 k (theta - th0)^2 for atoms (i-j-k).

    Parameters
    ----------
    r_ijk : (3,3) array
        Positions [[ri],[rj],[rk]] in nm, with j the vertex.
    th0 : float
        Equilibrium angle in radians.
    k : float
        Angle force constant in kJ·mol^-1·rad^-2.

    Returns
    -------
    energy : float
        Angle energy in kJ·mol^-1.
    forces : (3,3) array
        Forces on atoms i, j, k in kJ·mol^-1·nm^-1, ordered [Fi, Fj, Fk].
    """
    ri, rj, rk = r_ijk
    a = ri - rj              # vector j->i
    b = rk - rj              # vector j->k
    a2 = np.dot(a, a)
    b2 = np.dot(b, b)
    a_len = np.sqrt(a2)
    b_len = np.sqrt(b2)

    # Guard near-zero bond lengths
    eps_len = 1e-12
    if a_len < eps_len or b_len < eps_len:
        return 0.0, np.zeros((3,3))

    inv_ab = 1.0 / (a_len * b_len)
    cos_th = np.dot(a, b) * inv_ab
    # Clamp for acos
    cos_th = max(-1.0, min(1.0, cos_th))

    # Use atan2-like stable sin(theta) via cross product
    cross = np.cross(a, b)
    sin_th = np.linalg.norm(cross) * inv_ab
    # Clamp sin to avoid blowups near 0 or pi
    sin_th = max(sin_th, 1e-8)

    theta = np.arccos(cos_th)
    dU_dtheta = k * (theta - th0)
    g = dU_dtheta / sin_th

    # Force components on i and k
    Fi = g * ( b * inv_ab - cos_th * a / a2 )
    Fk = g * ( a * inv_ab - cos_th * b / b2 )
    Fj = -(Fi + Fk)

    energy = 0.5 * k * (theta - th0)**2
    forces = np.vstack((Fi, Fj, Fk))
    return energy, forces

```


This routine preserves symmetry, conserves momentum, and is robust in pathological geometries. It also slots directly into your `compute_energy_forces` loop:


```Python
for t, i, j, k in angles:
    th0 = angle_params[t]['th0']
    ka  = angle_params[t]['k']
    u, f = angle_potential(positions[[i, j, k]], th0, ka)
    energy['angle'] += u
    forces[[i, j, k]] += f
```


As with bonds, keep everything in radians and kJ·mol⁻¹ throughout, and you’ll avoid silent unit bugs.


## 3.3 The Dihedral Potential

### 3.3.1 Computing the signed dihedral angle robustly
Define the bond vectors ( $I$ 'll use $r_a$ for atom $a$ ):
- $\mathbf{b}_1=\mathbf{r}_2-\mathbf{r}_1$
- $\mathbf{b}_2=\mathbf{r}_3-\mathbf{r}_2$ (the central bond)
- $\mathbf{b}_3=\mathbf{r}_4-\mathbf{r}_3$

Build the two plane normals:
- $\mathbf{n}_1=\mathbf{b}_1 \times \mathbf{b}_2$
- $\mathbf{n}_2=\mathbf{b}_3 \times \mathbf{b}_2$

Let $b_2=\left\|\mathrm{b}_2\right\|$. A numerically stable signed $\phi$ is:

$$
\phi=\operatorname{atan} 2\left(b_2\left(\mathbf{b}_1 \cdot \mathbf{n}_2\right), \mathbf{n}_1 \cdot \mathbf{n}_2\right) .
$$


This avoids catastrophic cancellation when the planes become nearly parallel.


### 3.3.2 Energy derivative and mapping to forces
We first need $\frac{d U}{d \phi}$. Differentiate the series term-by-term:

$$
\frac{d U}{d \phi}=-\frac{1}{2} V_1 \sin \phi+V_2 \sin 2 \phi-\frac{3}{2} V_3 \sin 3 \phi+2 V_4 \sin 4 \phi .
$$


For the Cartesian forces we use the chain rule:

$$
\mathrm{F}_a=-\frac{d U}{d \phi} \frac{\partial \phi}{\partial \mathrm{r}_a} \quad(a=1,2,3,4) .
$$


A standard, symmetric and well-behaved set of $\partial \phi / \partial \mathbf{r}_a$ (used in MD codes and derived in, e.g., Swope et al., Mayo et al.) is:

Let
- $\mathbf{n}_1=\mathbf{b}_1 \times \mathbf{b}_2, \quad \mathbf{n}_2=\mathbf{b}_3 \times \mathbf{b}_2$
- $n_1^2=\left\|\mathbf{n}_1\right\|^2, n_2^2=\left\|\mathbf{n}_2\right\|^2, b_2^2=\left\|\mathbf{b}_2\right\|^2, b_2=\sqrt{b_2^2}$.

Then define

$$
\begin{aligned}
\mathbf{g}_1=-\frac{b_2}{n_1^2} \mathbf{n}_1, & \mathbf{g}_4=\frac{b_2}{n_2^2} \mathbf{n}_2, \\
\alpha=\frac{\mathbf{b}_1 \cdot \mathbf{b}_2}{b_2^2}, & \beta=\frac{\mathbf{b}_3 \cdot \mathbf{b}_2}{b_2^2}, \\
\mathbf{g}_2=-\mathbf{g}_1+\alpha \mathbf{g}_1-\beta \mathbf{g}_4, & \mathbf{g}_3=-\mathbf{g}_4+\beta \mathbf{g}_4-\alpha \mathbf{g}_1 .
\end{aligned}
$$

These $\mathbf{g}_a$ are the gradients $\partial \phi / \partial \mathbf{r}_a$. Finally,

$$
\mathrm{F}_a=-\frac{d U}{d \phi} \mathrm{~g}_a \text { for } a=1,2,3,4 .
$$


This construction guarantees:
- forces are equal-and-opposite overall ( $\sum_a \Gamma_a=0$ ),
- torques are balanced about the central bond,
- smooth behavior except at genuinely singular geometries (collinear triplets).


### 3.3.3 Numerical Stability and Edge Cases


- If $b_2$ is extremely small (central bond collapsed) or $n_1$ or $n_2$ nearly vanishes (colinear triplets), we'd divide by $\sim 0$. In practice, clamp denominators with small epsilons. When the geometry is that pathological, returning energy with very small/zero forces for that term is acceptable for a minimizer step (the bonded/angle terms will push you back to sanity).
- Keep $\phi$ in $(-\pi, \pi]$ via `atan2`, and use radians consistently.
- The dihedral does not handle non-bonded scaling; 1-4 LJ/Coulomb scaling is separate (apply it only in your non-bonded loop, as we discussed).


### 3.3.4 Python Implementation

Here `r` is a `(4,3)` array `[[ri],[rj],[rk],[rl]]` and `V` can be a dict with `V1..V4` or a 4-tuple. The code below accepts either.

```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
import numpy as np

def dihedral_potential(r_ijkl, V):
    """
    OPLS proper dihedral for atoms (i-j-k-l).
    U(phi) = 1/2 V1(1+cos φ) + 1/2 V2(1−cos 2φ) + 1/2 V3(1+cos 3φ) + 1/2 V4(1−cos 4φ)

    Parameters
    ----------
    r_ijkl : (4,3) array
        Positions [[ri],[rj],[rk],[rl]] in nm.
    V : dict or sequence
        Fourier coefficients V1..V4 in kJ·mol^-1. If dict, keys 'V1','V2','V3','V4'.

    Returns
    -------
    energy : float
        Torsion energy in kJ·mol^-1.
    forces : (4,3) array
        Forces on atoms i, j, k, l in kJ·mol^-1·nm^-1, ordered [Fi, Fj, Fk, Fl].
    """
    ri, rj, rk, rl = r_ijkl
    b1 = rj - ri
    b2 = rk - rj
    b3 = rl - rk

    # Precompute cross products & norms
    n1 = np.cross(b1, b2)
    n2 = np.cross(b3, b2)
    n1_sq = np.dot(n1, n1)
    n2_sq = np.dot(n2, n2)
    b2_sq = np.dot(b2, b2)
    b2_len = np.sqrt(b2_sq)

    # Guard against singular geometries
    eps = 1e-12
    if b2_len < eps or n1_sq < eps or n2_sq < eps:
        # Degenerate geometry; return finite energy and zero forces
        # (bond/angle terms will regularize the structure).
        # Compute phi in the safest way possible; if still degenerate, set to 0.
        try:
            x = np.dot(n1, n2)
            y = b2_len * np.dot(b1, n2)
            phi = np.arctan2(y, x)
        except Exception:
            phi = 0.0
        V1 = V.get('V1', 0.0) if isinstance(V, dict) else (V[0] if len(V) > 0 else 0.0)
        V2 = V.get('V2', 0.0) if isinstance(V, dict) else (V[1] if len(V) > 1 else 0.0)
        V3 = V.get('V3', 0.0) if isinstance(V, dict) else (V[2] if len(V) > 2 else 0.0)
        V4 = V.get('V4', 0.0) if isinstance(V, dict) else (V[3] if len(V) > 3 else 0.0)
        energy = 0.5*V1*(1+np.cos(phi)) + 0.5*V2*(1-np.cos(2*phi)) \
               + 0.5*V3*(1+np.cos(3*phi)) + 0.5*V4*(1-np.cos(4*phi))
        return float(energy), np.zeros((4,3))

    # Signed dihedral
    x = np.dot(n1, n2)
    y = b2_len * np.dot(b1, n2)
    phi = np.arctan2(y, x)

    # Unpack V1..V4
    if isinstance(V, dict):
        V1 = V.get('V1', 0.0); V2 = V.get('V2', 0.0); V3 = V.get('V3', 0.0); V4 = V.get('V4', 0.0)
    else:
        V1 = V[0] if len(V) > 0 else 0.0
        V2 = V[1] if len(V) > 1 else 0.0
        V3 = V[2] if len(V) > 2 else 0.0
        V4 = V[3] if len(V) > 3 else 0.0

    # Energy
    energy = (0.5*V1*(1+np.cos(phi))
              + 0.5*V2*(1-np.cos(2*phi))
              + 0.5*V3*(1+np.cos(3*phi))
              + 0.5*V4*(1-np.cos(4*phi)))

    # dU/dphi
    dUdphi = (-0.5*V1*np.sin(phi)
              +     V2*np.sin(2*phi)
              - 1.5*V3*np.sin(3*phi)
              +   2.0*V4*np.sin(4*phi))

    # Gradients of phi w.r.t. Cartesian coords (stable formulas)
    g1 = -(b2_len / n1_sq) * n1
    g4 =  (b2_len / n2_sq) * n2
    alpha = np.dot(b1, b2) / b2_sq
    beta  = np.dot(b3, b2) / b2_sq
    g2 = -g1 + alpha * g1 - beta * g4
    g3 = -g4 + beta  * g4 - alpha * g1

    # Forces: F = - dU/dphi * grad(phi)
    Fi = -dUdphi * g1
    Fj = -dUdphi * g2
    Fk = -dUdphi * g3
    Fl = -dUdphi * g4

    forces = np.vstack((Fi, Fj, Fk, Fl))
    return float(energy), forces

```




## 3.4 The Non-Bonded Potential





```Python
def nonbond_potential(r, indexes, q, sigma, epsilon, coulomb_factor):

    ri, rj = r
    rij = ri - rj
    r2 = float(np.dot(rij, rij))
    rlen = np.sqrt(r2)
    # Optional scales from indexes (for OPLS 1-4 pairs)
    scale_lj = scale_coul = 1.0
    if len(indexes) >= 4:
        # indexes can be (i, j, scale_lj, scale_coul)
        scale_lj = float(indexes[2])
        scale_coul = float(indexes[3])

    # Guard against zero distance
    if rlen < 1e-12:
        return 0.0, 0.0, np.zeros((2,3))

    uhat = rij / rlen

    # Unpack charges
    qi, qj = (float(q[0]), float(q[1]))

    # Geometric mixing if per-atom sigma/epsilon provided
    try:
        # try as array-like
        si, sj = float(sigma[0]), float(sigma[1])
        ei, ej = float(epsilon[0]), float(epsilon[1])
        sij = np.sqrt(si * sj)
        eij = np.sqrt(ei * ej)
    except (TypeError, IndexError):
        # treat as scalars already mixed
        sij = float(sigma)
        eij = float(epsilon)

    # --- Lennard-Jones energy and force ---
    sr = sij / rlen
    sr6 = sr**6
    sr12 = sr6 * sr6
    # Energy
    u_lj = 4.0 * eij * (sr12 - sr6)
    # Force magnitude along rhat
    f_lj_mag = (24.0 * eij / rlen) * (2.0 * sr12 - sr6)
    F_lj_i = f_lj_mag * uhat
    F_lj_j = -F_lj_i

    # --- Coulomb energy and force ---
    u_el = coulomb_factor * qi * qj / rlen
    f_el_mag = coulomb_factor * qi * qj / (rlen * rlen)
    F_el_i = f_el_mag * uhat
    F_el_j = -F_el_i

    # Apply optional 1-4 scaling (caller can pass scales via indexes)
    u_lj *= scale_lj
    u_el *= scale_coul
    Fi = scale_lj * F_lj_i + scale_coul * F_el_i
    Fj = scale_lj * F_lj_j + scale_coul * F_el_j

    forces = np.vstack((Fi, Fj))
    return float(u_lj), float(u_el), forces

```


## 4.) We need to then write a function to compute the total energy and forces on the molecule:





```Python
def compute_energy_forces(positions, types, bonds, angles, dihedrals,
                          nb_params, bond_params, angle_params, dihedral_params):
    # running totals
    energy = {'bond':0.0, 'angle':0.0, 'dihedral':0.0, 'LJ':0.0, 'coulomb':0.0}
    N = len(positions)
    forces = np.zeros((N,3))

    # --- 1) bonded sums ---
    for t, i, j in bonds:
        r0 = bond_params[t]['r0']; kb = bond_params[t]['k']
        u, f = bond_potential(positions[[i, j]], r0, kb)
        energy['bond'] += u
        forces[[i, j]] += f

    for t, i, j, k in angles:
        th0 = angle_params[t]['th0']; ka = angle_params[t]['k']
        u, f = angle_potential(positions[[i, j, k]], th0, ka)
        energy['angle'] += u
        forces[[i, j, k]] += f

    for t, i, j, k, l in dihedrals:
        V = dihedral_params[t]   # dict with V1..V4 (kJ/mol)
        u, f = dihedral_potential(positions[[i, j, k, l]], V)
        energy['dihedral'] += u
        forces[[i, j, k, l]] += f

    # --- 2) build OPLS exclusion/scale sets from topology ---
    s12 = set()
    for _, i, j in bonds:
        s12.add((min(i,j), max(i,j)))

    s13 = set()
    for _, i, j, k in angles:
        s13.add((min(i,k), max(i,k)))

    s14 = set()
    for _, i, j, k, l in dihedrals:
        p = (min(i,l), max(i,l))
        if p not in s12 and p not in s13:
            s14.add(p)

    # --- 3) non-bonded sums (LJ + Coulomb, with OPLS 1–4 scaling) ---
    for i in range(N):
        for j in range(i+1, N):
            p = (i, j)
            if p in s12 or p in s13:
                continue
            scale = 0.5 if p in s14 else 1.0

            ti, tj = types[i], types[j]
            qi, qj = nb_params[ti]['q'], nb_params[tj]['q']
            si, sj = nb_params[ti]['sigma'], nb_params[tj]['sigma']
            ei, ej = nb_params[ti]['epsilon'], nb_params[tj]['epsilon']
            ke     = nb_params[ti]['coulomb_factor']  # same for all types

            u_lj, u_el, f = nonbond_potential(
                positions[[i, j]],
                (i, j, scale, scale),     # pass 1–4 scaling here if needed
                (qi, qj),
                (si, sj),                 # per-atom; function does geometric mixing
                (ei, ej),
                ke
            )
            energy['LJ']      += u_lj
            energy['coulomb'] += u_el
            forces[[i, j]]    += f

    return energy, forces

```


# Steepest Descent




We will calculate the initial energies of the molecule:


```Python
positions, types, bonds, angles, dihedrals =
nb_params, bond_params, angle_params, dihedral_params =
# forces and energies of initial configuration
U,F =
```


```Python
# steepest descent
gamma = 1.0e-6                       # nm^2·mol/kJ
tol = 10.0                           # kJ/mol/nm (max force criterion)
max_iterations = 1000

for i in range(max_iterations):
    U, F = compute_energy_forces(positions, types,
                                 bonds, angles, dihedrals,
                                 nb_params, bond_params, angle_params, dihedral_params)
    Fmax = np.max(np.linalg.norm(F, axis=1))
    if Fmax < tol:
        break
    positions += gamma * F

if i == max_iterations:
    print('warning: may not have converged')
else:
    print(f'converged in {i} iterations')


```


And then we use steepest descent to minimize the energy and get the .xyz file for the lowest energy conformation of ethane:


```Python
U_final,F_final =

print('Final energy and forces')
print('-----------------------')
print('total: {:.5g}'.format(sum(U.values())))
for k,v in U.items():
    print('{}: {:.5g}'.format(k,v))
print(F)
print('-----------------------')
```



And finally, we will write out the minimized structure as an .xyz file.

![[Pasted image 20250919230832.png|500]]
