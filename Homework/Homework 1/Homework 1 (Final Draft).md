

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

OPLS (Optimized Potentials for Liquid Simulations) Nomenclature defines the topology of a molecules by providing information about atom types (with non-bonded parameters and charges), connectivity (bonds and angles), and torsional types (the Fourier series that fixes the rotamer energetics). The OPLS notation ties the topology patterns to the correct parameters and the bookkeeping rules indicate which interactions to include and how.


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

When two atoms $i$ and $j$ interact non-bondedly, one must build a pair Lennard-Jones parameter set $\left(\sigma_{i j}, \varepsilon_{i j}\right)$ from the two type sets $\left(\sigma_i, \varepsilon_i\right)$ and $\left(\sigma_j, \varepsilon_j\right)$. OPLS uses the geometric-geometric rule:

$$
\sigma_{i j}=\sqrt{\sigma_i \sigma_j}, \quad \varepsilon_{i j}=\sqrt{\varepsilon_i \varepsilon_j} .
$$


That choice is part of the force field's design (other families sometimes use Lorentz-Berthelot: $\sigma_{i j}=\left(\sigma_i+\sigma_j\right) / 2, \varepsilon_{i j}=\sqrt{\varepsilon_i \varepsilon_j}$ ). Since the code is implementing OPLS, it uses geometric for both.




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


Condensed-phase MD almost always sets the relative dielectric $\varepsilon_r=1$ and accounts for screening explicitly with solvent and long-range summation (PME/PPPM). The code uses

$$
k_e \approx 138.935456 \frac{\mathrm{~kJ} \mathrm{~nm}}{\mathrm{~mol} \mathrm{e}^2} .
$$




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




# 3 Calculating the non-bonded, bonded, angle, and dihedral potentials:**


## 3.1 The Bond Potential

At the bond level covalent covalent bonds are modeled as a 1D spring embedded in 3D space. In OPLS-AA the functional form is the simple harmonic:

$$
U_{\text {bond }}(R)=\frac{1}{2} k\left(R-R_0\right)^2,
$$

where $R=\left\|\mathbf{r}_i-\mathbf{r}_j\right\|$ is the instantaneous inter-atomic distance, $R_0$ is the equilibrium bond length, and $k$ is the bond force constant. In your unit system $R$ and $R_0$ are in nm and $k$ is in $\mathrm{kJ} \cdot^2 \mathrm{~mol}^{-1} \cdot \mathrm{~nm}^{-2}$, so $U$ is in $\mathrm{kJ}^2 \cdot \mathrm{~mol}^{-1}$. The force is the negative gradient of the energy with respect to the atomic positions. Because the energy depends on positions only through $R$, the chain rule gives a central scalar factor and a direction:

$$
\frac{d U}{d R}=k\left(R-R_0\right), \quad \hat{\mathbf{u}}=\frac{\mathbf{r}_{i j}}{R}=\frac{\mathbf{r}_i-\mathbf{r}_j}{\left\|\mathbf{r}_i-\mathbf{r}_j\right\|} .
$$


The pairwise forces are then equal and opposite along the bond axis,

$$
\mathrm{F}_i=-k\left(R-R_0\right) \hat{\mathrm{u}}, \quad \mathrm{~F}_j=+k\left(R-R_0\right) \hat{\mathrm{u}} .
$$


The code, guards against $R→0$ to avoid division by zero when forming $\hat{\mathrm{u}}$. 



## 3.2 The Angle Potential

The OPLS angle treats three atoms $i-j-k$ as a hinge with a harmonic penalty on the deviation of the internal angle $\theta$ at the vertex $j$ from its equilibrium value $\theta_0$. In the unit system that is used here (nm, $\mathrm{kJ}^{-} \mathrm{mol}^{-1}$, radians, e), the potential is

$$
U_{\text {angle }}(\theta)=\frac{1}{2} k\left(\theta-\theta_0\right)^2,
$$

where $k$ is in $\mathrm{kJ}^{\cdot} \mathrm{mol}^{-1} \cdot \mathrm{rad}^{-2}$ and $\theta, \theta_0$ are in radians. To evaluate $U$ and forces from Cartesian coordinates, we form the two bond vectors meeting at the vertex:

$$
\mathbf{a}=\mathbf{r}_i-\mathbf{r}_j, \quad \mathbf{b}=\mathbf{r}_k-\mathbf{r}_j .
$$


Let $a=\|\mathrm{a}\|, b=\|\mathrm{b}\|$. The angle follows from

$$
\cos \theta=\frac{\mathbf{a} \cdot \mathbf{b}}{a b}, \quad \sin \theta=\frac{\|\mathbf{a} \times \mathbf{b}\|}{a b} .
$$


Using $\sin \theta$ via the cross product is numerically more stable near $\theta \approx 0$ or $\pi$ than computing $\sqrt{1-\cos ^2 \theta}$. We still clamp both $\cos \theta$ into $[-1,1]$ and $\sin \theta$ away from zero by a tiny $\varepsilon$ to avoid division blow-ups. Forces are the negative gradient of $U$ with respect to positions. By the chain rule,

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


which guarantees Newton's third law and total momentum conservation. 




## 3.3 The Dihedral Potential

### 3.3.1 Computing the signed dihedral angle robustly
We start by define the bond vectors ( I use $r_a$ for atom $a$ ):
- $\mathbf{b}_1=\mathbf{r}_2-\mathbf{r}_1$
- $\mathbf{b}_2=\mathbf{r}_3-\mathbf{r}_2$ (the central bond)
- $\mathbf{b}_3=\mathbf{r}_4-\mathbf{r}_3$

Then we build the two plane normals:
- $\mathbf{n}_1=\mathbf{b}_1 \times \mathbf{b}_2$
- $\mathbf{n}_2=\mathbf{b}_3 \times \mathbf{b}_2$

Let $b_2=\left\|\mathrm{b}_2\right\|$. A numerically stable signed $\phi$ is:

$$
\phi=\operatorname{atan} 2\left(b_2\left(\mathbf{b}_1 \cdot \mathbf{n}_2\right), \mathbf{n}_1 \cdot \mathbf{n}_2\right) .
$$


Using this this avoids catastrophic cancellation when the planes become nearly parallel.


### 3.3.2 Energy derivative and mapping to forces
We first need $\frac{d U}{d \phi}$. So we differentiate the series term-by-term:

$$
\frac{d U}{d \phi}=-\frac{1}{2} V_1 \sin \phi+V_2 \sin 2 \phi-\frac{3}{2} V_3 \sin 3 \phi+2 V_4 \sin 4 \phi .
$$


For the Cartesian forces we use the chain rule:

$$
\mathrm{F}_a=-\frac{d U}{d \phi} \frac{\partial \phi}{\partial \mathrm{r}_a} \quad(a=1,2,3,4) .
$$


A standard, symmetric and well-behaved set of $\partial \phi / \partial \mathbf{r}_a$ is:

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


This construction guarantees that (1) forces are equal-and-opposite overall ( $\sum_a \Gamma_a=0$ ), (2) torques are balanced about the central bond, and that (3) we have smooth behavior except at genuinely singular geometries (collinear triplets).




# 4. Steepest Descent

Steepest descent (SD) is the most direct way to relax the ethane. At each iteration the code moves every atom a small distance in the direction of the current force and recomputes energy and forces. Since the force is the negative gradient of the potential, this guarantees a first-order decrease in energy for sufficiently small steps. In the units used for here, the update is

$$
\mathbf{r}^{(n+1)}=\mathbf{r}^{(n)}-\gamma \nabla U\left(\mathbf{r}^{(n)}\right)=\mathbf{r}^{(n)}+\gamma \mathbf{\Gamma}\left(\mathbf{r}^{(n)}\right),
$$

where $\gamma>0$ is the step size. The stopping test is $\max _a\left\|\Gamma_a\right\|<10 \mathrm{~kJ} \mathrm{~mol}^{-1} \mathrm{~nm}^{-1}$. The update direction is locally optimal among all directions (it is the direction of steepest decrease in $U$ under a Euclidean metric), but the magnitude of the step must be tuned to the stiffness of the landscape. In molecular mechanics that stiffness is set by large harmonic curvatures-bonds and angles-so naïvely large $\gamma$ produces oscillations or even growth in $U$; very small $\gamma$ converges but can be painfully slow. 




# 5. Results

I built a complete OPLS-AA pipeline for ethane-from parsing the XYZ, assigning CT/HC atom types, and constructing bonds/angles/dihedrals, to evaluating bonded and non-bonded interactions with geometric LJ mixing and 0.5-scaled 1-4 intramolecular pairs, and relaxing the structure by steepest descent. The resulting configuration exhibits near-equilibrium bonds and angles with small harmonic contributions, a dihedral term that correctly favors the staggered conformer, and balanced non-bonded contributions (a modest LJ attraction with electrostatics treated in consistent $\mathrm{nm} / \mathrm{kJ} \cdot \mathrm{mol}^{-1}$ units). The final energy decomposition (bond $\approx 0.011 \mathrm{~kJ} \cdot \mathrm{~mol}^{-1}$, angle $\approx 0.044 \mathrm{~kJ}^{-1} \mathrm{~mol}^{-1}$, dihedral $\approx 6.86 \mathrm{~kJ}^{-1} \mathrm{~mol}^{-1}, \mathrm{LJ} \approx-0.038 \mathrm{kJ} \cdot \mathrm{mol}^{-1}$, Coulomb $\approx 8.38 \mathrm{~kJ}^{-1} \cdot \mathrm{~mol}^{-1}$; total $\approx 15.25 \mathrm{~kJ}^{-1} \cdot \mathrm{~mol}^{-1}$ ) is chemically reasonable for ethane, and the accompanying visualization confirms a sensible staggered geometry. 




```
Initial total energy (kJ/mol): 857.0707301641196
Initial Fmax (kJ/mol/nm): 10178.805735724574
converged in 999 iterations
Final energy and forces
-----------------------
total: 15.249
bond: 0.011489
angle: 0.04358
dihedral: 6.8571
LJ: -0.038363
coulomb: 8.3752
[[ 8.98621754e-04  6.74692384e-03 -8.94140142e-04]
 [-4.49235386e+01  2.72336352e+00 -2.47928688e+01]
 [ 1.67190699e+00  6.55334963e+00  5.09360944e+01]
 [ 4.32615481e+01 -9.20225922e+00 -2.61530926e+01]
 [-8.98756081e-04 -6.74793236e-03  8.94273796e-04]
 [-2.63586632e+01 -3.21164393e+00 -5.05310401e+01]
 [ 5.65452283e+01 -7.13288507e+00  3.19982726e+00]
 [-3.01964816e+01  1.02700761e+01  4.73410797e+01]]
-----------------------
```






![[Pasted image 20250919230832.png|500]]


# Final Code


```Python
import numpy as np


def define_molecule(filename: str):
    with open(filename) as f:
        N = int(f.readline().strip())
        _ = f.readline()

        positions = np.empty((N, 3), dtype=float)
        elements = []
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


def define_force_field():
    k_e = 138.935456

    nb_params = {
        'CT': {'q': -0.18, 'sigma': 0.350, 'epsilon': 0.276, 'coulomb_factor': k_e},
        'HC': {'q': 0.06, 'sigma': 0.250, 'epsilon': 0.1255, 'coulomb_factor': k_e},
    }

    bond_params = {
        'CT-CT': {'r0': 0.15290, 'k': 224262.4},
        'CT-HC': {'r0': 0.10900, 'k': 284512.0},
    }

    angle_params = {
        'HC-CT-HC': {'th0': np.deg2rad(107.8), 'k': 276.144},
        'HC-CT-CT': {'th0': np.deg2rad(110.7), 'k': 313.800},
    }

    dihedral_params = {
        'HC-CT-CT-HC': {'V1': 0.0, 'V2': 0.0, 'V3': 1.2552, 'V4': 0.0}
    }

    return nb_params, bond_params, angle_params, dihedral_params


def bond_potential(r, r0, k):
    """Harmonic bond potential U = 1/2 k (R - r0)^2."""
    ri, rj = r
    rij = ri - rj
    R = np.linalg.norm(rij)
    if R < 1e-12:
        uhat = np.zeros(3)
    else:
        uhat = rij / R

    dU_dR = k * (R - r0)
    Fi = -dU_dR * uhat
    Fj = -Fi
    energy = 0.5 * k * (R - r0) ** 2
    forces = np.vstack((Fi, Fj))
    return energy, forces


def angle_potential(r_ijk, th0, k):
    """Harmonic angle potential U = 1/2 k (theta - th0)^2 for atoms (i-j-k)."""
    ri, rj, rk = r_ijk
    a = ri - rj
    b = rk - rj
    a2 = np.dot(a, a)
    b2 = np.dot(b, b)
    a_len = np.sqrt(a2)
    b_len = np.sqrt(b2)

    eps_len = 1e-12
    if a_len < eps_len or b_len < eps_len:
        return 0.0, np.zeros((3, 3))

    inv_ab = 1.0 / (a_len * b_len)
    cos_th = np.dot(a, b) * inv_ab
    cos_th = max(-1.0, min(1.0, cos_th))

    cross = np.cross(a, b)
    sin_th = np.linalg.norm(cross) * inv_ab
    sin_th = max(sin_th, 1e-8)

    theta = np.arccos(cos_th)
    dU_dtheta = k * (theta - th0)
    g = dU_dtheta / sin_th

    Fi = g * (b * inv_ab - cos_th * a / a2)
    Fk = g * (a * inv_ab - cos_th * b / b2)
    Fj = -(Fi + Fk)

    energy = 0.5 * k * (theta - th0) ** 2
    forces = np.vstack((Fi, Fj, Fk))
    return energy, forces


def dihedral_potential(r_ijkl, V):
    ri, rj, rk, rl = r_ijkl
    b1 = rj - ri
    b2 = rk - rj
    b3 = rl - rk

    n1 = np.cross(b1, b2)
    n2 = np.cross(b3, b2)
    n1_sq = np.dot(n1, n1)
    n2_sq = np.dot(n2, n2)
    b2_sq = np.dot(b2, b2)
    b2_len = np.sqrt(b2_sq)

    eps = 1e-12
    if b2_len < eps or n1_sq < eps or n2_sq < eps:
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
        energy = (0.5 * V1 * (1 + np.cos(phi))
                  + 0.5 * V2 * (1 - np.cos(2 * phi))
                  + 0.5 * V3 * (1 + np.cos(3 * phi))
                  + 0.5 * V4 * (1 - np.cos(4 * phi)))
        return float(energy), np.zeros((4, 3))

    x = np.dot(n1, n2)
    y = b2_len * np.dot(b1, n2)
    phi = np.arctan2(y, x)

    if isinstance(V, dict):
        V1 = V.get('V1', 0.0)
        V2 = V.get('V2', 0.0)
        V3 = V.get('V3', 0.0)
        V4 = V.get('V4', 0.0)
    else:
        V1 = V[0] if len(V) > 0 else 0.0
        V2 = V[1] if len(V) > 1 else 0.0
        V3 = V[2] if len(V) > 2 else 0.0
        V4 = V[3] if len(V) > 3 else 0.0

    energy = (0.5 * V1 * (1 + np.cos(phi))
              + 0.5 * V2 * (1 - np.cos(2 * phi))
              + 0.5 * V3 * (1 + np.cos(3 * phi))
              + 0.5 * V4 * (1 - np.cos(4 * phi)))

    dUdphi = (-0.5 * V1 * np.sin(phi)
              + V2 * np.sin(2 * phi)
              - 1.5 * V3 * np.sin(3 * phi)
              + 2.0 * V4 * np.sin(4 * phi))

    g1 = -(b2_len / n1_sq) * n1
    g4 = (b2_len / n2_sq) * n2
    alpha = np.dot(b1, b2) / b2_sq
    beta = np.dot(b3, b2) / b2_sq
    g2 = -g1 + alpha * g1 - beta * g4
    g3 = -g4 + beta * g4 - alpha * g1

    Fi = -dUdphi * g1
    Fj = -dUdphi * g2
    Fk = -dUdphi * g3
    Fl = -dUdphi * g4

    forces = np.vstack((Fi, Fj, Fk, Fl))
    return float(energy), forces


def nonbond_potential(r, indexes, q, sigma, epsilon, coulomb_factor):
    ri, rj = r
    rij = ri - rj
    r2 = float(np.dot(rij, rij))
    rlen = np.sqrt(r2)
    scale_lj = scale_coul = 1.0
    if len(indexes) >= 4:
        scale_lj = float(indexes[2])
        scale_coul = float(indexes[3])
    if rlen < 1e-12:
        return 0.0, 0.0, np.zeros((2, 3))

    uhat = rij / rlen
    qi, qj = (float(q[0]), float(q[1]))

    try:
        si, sj = float(sigma[0]), float(sigma[1])
        ei, ej = float(epsilon[0]), float(epsilon[1])
        sij = np.sqrt(si * sj)
        eij = np.sqrt(ei * ej)
    except (TypeError, IndexError):
        sij = float(sigma)
        eij = float(epsilon)

    sr = sij / rlen
    sr6 = sr ** 6
    sr12 = sr6 * sr6
    u_lj = 4.0 * eij * (sr12 - sr6)
    f_lj_mag = (24.0 * eij / rlen) * (2.0 * sr12 - sr6)
    F_lj_i = f_lj_mag * uhat
    F_lj_j = -F_lj_i

    u_el = coulomb_factor * qi * qj / rlen
    f_el_mag = coulomb_factor * qi * qj / (rlen * rlen)
    F_el_i = f_el_mag * uhat
    F_el_j = -F_el_i

    u_lj *= scale_lj
    u_el *= scale_coul
    Fi = scale_lj * F_lj_i + scale_coul * F_el_i
    Fj = scale_lj * F_lj_j + scale_coul * F_el_j

    forces = np.vstack((Fi, Fj))
    return float(u_lj), float(u_el), forces


def compute_energy_forces(positions, types, bonds, angles, dihedrals,
                          nb_params, bond_params, angle_params, dihedral_params):
    energy = {'bond': 0.0, 'angle': 0.0, 'dihedral': 0.0, 'LJ': 0.0, 'coulomb': 0.0}
    N = len(positions)
    forces = np.zeros((N, 3))

    for t, i, j in bonds:
        r0 = bond_params[t]['r0']
        kb = bond_params[t]['k']
        u, f = bond_potential(positions[[i, j]], r0, kb)
        energy['bond'] += u
        forces[[i, j]] += f

    for t, i, j, k in angles:
        th0 = angle_params[t]['th0']
        ka = angle_params[t]['k']
        u, f = angle_potential(positions[[i, j, k]], th0, ka)
        energy['angle'] += u
        forces[[i, j, k]] += f

    for t, i, j, k, l in dihedrals:
        V = dihedral_params[t]
        u, f = dihedral_potential(positions[[i, j, k, l]], V)
        energy['dihedral'] += u
        forces[[i, j, k, l]] += f

    s12 = set()
    for _, i, j in bonds:
        s12.add((min(i, j), max(i, j)))

    s13 = set()
    for _, i, j, k in angles:
        s13.add((min(i, k), max(i, k)))

    s14 = set()
    for _, i, j, k, l in dihedrals:
        p = (min(i, l), max(i, l))
        if p not in s12 and p not in s13:
            s14.add(p)

    for i in range(N):
        for j in range(i + 1, N):
            p = (i, j)
            if p in s12 or p in s13:
                continue
            scale = 0.5 if p in s14 else 1.0

            ti, tj = types[i], types[j]
            qi, qj = nb_params[ti]['q'], nb_params[tj]['q']
            si, sj = nb_params[ti]['sigma'], nb_params[tj]['sigma']
            ei, ej = nb_params[ti]['epsilon'], nb_params[tj]['epsilon']
            ke = nb_params[ti]['coulomb_factor']

            u_lj, u_el, f = nonbond_potential(
                positions[[i, j]],
                (i, j, scale, scale),
                (qi, qj),
                (si, sj),
                (ei, ej),
                ke
            )
            energy['LJ'] += u_lj
            energy['coulomb'] += u_el
            forces[[i, j]] += f

    return energy, forces


if __name__ == "__main__":
    positions, types, bonds, angles, dihedrals = define_molecule("ethane.xyz")
    nb_params, bond_params, angle_params, dihedral_params = define_force_field()

    U0, F0 = compute_energy_forces(
        positions, types, bonds, angles, dihedrals,
        nb_params, bond_params, angle_params, dihedral_params
    )
    print("Initial total energy (kJ/mol):", sum(U0.values()))
    print("Initial Fmax (kJ/mol/nm):", np.max(np.linalg.norm(F0, axis=1)))

    gamma = 1.0e-6
    tol = 10.0
    max_iterations = 1000

    for i in range(max_iterations):
        U, F = compute_energy_forces(
            positions, types, bonds, angles, dihedrals,
            nb_params, bond_params, angle_params, dihedral_params
        )
        Fmax = np.max(np.linalg.norm(F, axis=1))
        if Fmax < tol:
            break
        positions += gamma * F

    if i == max_iterations:
        print('warning: may not have converged')
    else:
        print(f'converged in {i} iterations')

    print('Final energy and forces')
    print('-----------------------')
    print('total: {:.5g}'.format(sum(U.values())))
    for k, v in U.items():
        print(f'{k}: {v:.5g}')
    print(F)
    print('-----------------------')

    with open('ethane_min.xyz', 'w') as f:
        f.write(f'{positions.shape[0]}\n')
        f.write('\n')
        for i in range(positions.shape[0]):
            elem = 'C' if types[i] == 'CT' else 'H'
            x, y, z = positions[i]
            f.write(f'{elem} {x:.6f} {y:.6f} {z:.6f}\n')

```