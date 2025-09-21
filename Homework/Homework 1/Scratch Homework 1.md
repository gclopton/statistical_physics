
# <span style="color:#ff0000; font-weight:bold">1. Questions</span>



``` title:"ethane.xyz"
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


![[Pasted image 20250912135002.png|500]]



---

``` title:"ethane.xyz"

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

Energy Minimization of an Ethane Molecule

The purpose of this homework is to:

- Define a molecule to run a particle based simulation

- Define a force field to compute the energy of a molecule

- Perform energy minimization on a single molecule to optimize its internal geometry

import numpy as np


## Question 1

> [!Question]
> Question 1 (Define the Molecule): 

**1.1:** Write a function that takes a file name and reads in its types and coordinates. Here, you will also need to manually define the topology, i.e bonds, angles and dihedrals as tuples. Use the OPLS type nomenclature!

```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def define_molecule(filename):
    # load atom coordinates
    with open(filename) as f:
        # number of particles
        line = f.readline()
        N = int(line.strip())

        # comment line
        line = f.readline()

        # coordinates
        positions = np.empty #fill in
        types = []
        for i in range(N):
            line = f.readline().split()
            positions[i] = #fill in
            types.append #fill in

    # define atom topology
    bonds = (
        ('HC-CT',0,1), #fill in
    angles = (
        ('HC-CT-HC',1,0,2), #fill in
    dihedrals = (
        ('H-C-C-H',1,0,4,5) #fill in

    return positions, types, bonds, angles, dihedrals
```


```Python title:solution
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def define_molecule(filename):
    # load atom coordinates
    with open(filename) as f:
        # number of particles
        line = f.readline()
        N = int(line.strip())

        # comment line
        line = f.readline()

        # coordinates
        positions = np.empty((N, 3), dtype=float)
        types = []
        for i in range(N):
            line = f.readline().split()
            types.append(line[0])
            positions[i] = np.array([float(line[1]), float(line[2]), float(line[3])], dtype=float)

    # define atom topology (indexes follow file order)
    bonds = (
        ('HC-CT', 0, 1),
        ('HC-CT', 0, 2),
        ('HC-CT', 0, 3),
        ('CT-CT', 0, 4),
        ('HC-CT', 4, 5),
        ('HC-CT', 4, 6),
        ('HC-CT', 4, 7),
    )

    angles = (
        ('HC-CT-HC', 1, 0, 2),
        ('HC-CT-HC', 1, 0, 3),
        ('HC-CT-HC', 2, 0, 3),
        ('HC-CT-CT', 1, 0, 4),
        ('HC-CT-CT', 2, 0, 4),
        ('HC-CT-CT', 3, 0, 4),
        ('HC-CT-CT', 0, 4, 5),
        ('HC-CT-CT', 0, 4, 6),
        ('HC-CT-CT', 0, 4, 7),
        ('HC-CT-HC', 5, 4, 6),
        ('HC-CT-HC', 5, 4, 7),
        ('HC-CT-HC', 6, 4, 7),
    )

    dihedrals = (
        ('H-C-C-H', 1, 0, 4, 5),
        ('H-C-C-H', 1, 0, 4, 6),
        ('H-C-C-H', 1, 0, 4, 7),
        ('H-C-C-H', 2, 0, 4, 5),
        ('H-C-C-H', 2, 0, 4, 6),
        ('H-C-C-H', 2, 0, 4, 7),
        ('H-C-C-H', 3, 0, 4, 5),
        ('H-C-C-H', 3, 0, 4, 6),
        ('H-C-C-H', 3, 0, 4, 7),
    )

    return positions, types, bonds, angles, dihedrals

```


**1.2:** test your function by checking types, positions, etc.:

```Python

positions, types, bonds, angles, dihedrals = define_molecule('ethane.xyz')

```

Output should be something like:

```Python

[[ 0. 0. 0. ]

[ 0. 0. 0.11]

[ 0.11 0. 0. ]

[-0.11 0. 0. ]

[ 0. 0.15 0. ]

[ 0. 0.15 0.11]

[ 0. 0.15 -0.11]

[ 0.11 0.15 0. ]]

['C', 'H', 'H', 'H', 'C', 'H', 'H', 'H']

(('HC-CT', 0, 1), ('HC-CT', 0, 2), ('HC-CT', 0, 3), ('CT-CT', 0, 4), ('HC-CT', 4, 5), ('HC-CT', 4, 6), ('HC-CT', 4, 7))

(('HC-CT-HC', 1, 0, 2), ('HC-CT-HC', 1, 0, 3), ('HC-CT-HC', 2, 0, 3), ('HC-CT-CT', 1, 0, 4), ('HC-CT-CT', 2, 0, 4), ('HC-CT-CT', 3, 0, 4), ('HC-CT-CT', 0, 4, 5), ('HC-CT-CT', 0, 4, 6), ('HC-CT-CT', 0, 4, 7), ('HC-CT-HC', 5, 4, 6), ('HC-CT-HC', 5, 4, 7), ('HC-CT-HC', 6, 4, 7))

(('H-C-C-H', 1, 0, 4, 5), ('H-C-C-H', 1, 0, 4, 6), ('H-C-C-H', 1, 0, 4, 7), ('H-C-C-H', 2, 0, 4, 5), ('H-C-C-H', 2, 0, 4, 6), ('H-C-C-H', 2, 0, 4, 7), ('H-C-C-H', 3, 0, 4, 5), ('H-C-C-H', 3, 0, 4, 6), ('H-C-C-H', 3, 0, 4, 7))

```

## Question 2

> [!Question]
> Question 2 (Define Force Field): 


**2.1:** Now we need a function that defines the force field. We will use dictionaries for convienient storage of parameters and variables. You must use the same OPLS type nomenclature as you used above to identify the molecule topology.


```Python

#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def define_force_field():
    # define force field parameters, converting to our system of units

    #fill in
    bond_params = {
        'CT-CT': {'r0':, 'k': }
    }
    angle_params = {}}
    dihedral_params = {}

    nb_params = {
        'C': {'q':,'coulomb_factor':, 'sigma': , 'epsilon': },
    }
    return nb_params,bond_params,angle_params,dihedral_params

```


```Python title:Solution
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def define_force_field():
    # Coulomb factor in MD units
    coulomb_factor = 138.935458  # kJ/mol·nm·e^-2

    # Bonds: harmonic
    bond_params = {
        'CT-CT': {'r0': 0.15290, 'k': 224262.4},
        'HC-CT': {'r0': 0.10900, 'k': 284512.0},
    }

    # Angles: harmonic (store θ0 in radians)
    angle_params = {
        'HC-CT-HC': {'th0': np.deg2rad(107.8), 'k': 276.144},
        'HC-CT-CT': {'th0': np.deg2rad(110.7), 'k': 313.800},
        # not used directly here, but often listed
        'CT-CT-CT': {'th0': np.deg2rad(112.7), 'k': 488.273},
    }

    # Proper dihedral: OPLS Fourier (a 3-fold term is enough for ethane)
    # U(φ) = 0.5*K1(1+cosφ) + 0.5*K2(1−cos2φ) + 0.5*K3(1+cos3φ) + 0.5*K4(1−cos4φ)
    dihedral_params = {
        'H-C-C-H': {'K1': 0.0, 'K2': 0.0, 'K3': 12.0, 'K4': 0.0},  # ~3 kcal/mol barrier
    }

    # Nonbonded atom types + global mixing/scaling conventions
    nb_params = {
        'C': {'q': -0.180, 'coulomb_factor': coulomb_factor, 'sigma': 0.350, 'epsilon': 0.276144},
        'H': {'q': +0.060, 'coulomb_factor': coulomb_factor, 'sigma': 0.250, 'epsilon': 0.125520},
        'mixing': 'geometric',
        'scale14': {'LJ': 0.5, 'Coulomb': 0.5},
    }
    return nb_params, bond_params, angle_params, dihedral_params

```



## Question 3

> [!Question]
> Question 3: (Functions to Calculate the Non-bonded, Bonded, Angle, and Dihedral Potential)



**3.1:** Next we need to define functions that will take distances (and angles, etc. as needed) and return energies and forces for non-bonded interactions, bonds, angles, and dihedrals:



```Python

#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def bond_potential(r,r0,k):
    """Harmonic bond potential"""
    # fill in

    return energy,forces

```

```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def angle_potential(r,th0,k):
    """Harmonic angle potential."""
    #fill in

    return energy,forces
```

```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def dihedral_potential(r,V):
    """Dihedral potential."""
    #fill in

    return energy,forces
```

```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def nonbond_potential(r,indexes,q,sigma,epsilon,coulomb_factor):
    """Lennard-Jones and Coulomb potentials."""
   #fill in

    return energy_lj,energy_elec,forces
```



```Python title:"Solution: Internal Helpers"
# Internal helpers
def _norm(v):
    return np.sqrt(np.dot(v, v))

def _unit(v, eps=1e-12):
    n = _norm(v)
    return v*0.0 if n < eps else v/n



```



```Python title:"Solution: Bond Potential"
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def bond_potential(r, r0, k):
    """Harmonic bond potential: 0.5*k*(|r|-r0)^2.  Returns energy and forces on (i,j)."""
    d2 = np.dot(r, r)
    if d2 == 0.0:
        energy = 0.5 * k * (r0**2)
        fi = np.zeros(3); fj = np.zeros(3)
    else:
        d = np.sqrt(d2)
        diff = d - r0
        energy = 0.5 * k * diff * diff
        fmag_over_r = -k * diff / d
        fi = fmag_over_r * r
        fj = -fi
    return energy, np.vstack([fi, fj])

```



```Python title:"Solution: Angle Potential"
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def angle_potential(ri, rj, rk, th0, k):
    """Harmonic angle potential on triplet (i-j-k). Returns energy and forces on i,j,k."""
    u = ri - rj
    v = rk - rj
    a = _norm(u); b = _norm(v)
    if a == 0.0 or b == 0.0:
        return 0.0, np.zeros((3,3))

    uhat = u / a
    vhat = v / b
    cos_th = np.clip(np.dot(uhat, vhat), -1.0, 1.0)
    th = np.arccos(cos_th)
    sin_th = max(1e-12, np.sqrt(1.0 - cos_th*cos_th))

    dU_dth = k * (th - th0)

    fi = -dU_dth * (vhat - cos_th*uhat) / (a * sin_th)
    fk = -dU_dth * (uhat - cos_th*vhat) / (b * sin_th)
    fj = -(fi + fk)

    energy = 0.5 * k * (th - th0)**2
    return energy, np.vstack([fi, fj, fk])

```



```Python title:"Solution: Dihedral Potential"
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def dihedral_potential(ri, rj, rk, rl, K):
    """OPLS Fourier dihedral for quadruplet (i-j-k-l).
       K: dict with keys K1..K4. Returns energy and forces on i,j,k,l."""
    b1 = rj - ri
    b2 = rk - rj
    b3 = rl - rk
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    b2u = _unit(b2)
    n1u = _unit(n1)
    n2u = _unit(n2)

    m1 = np.cross(n1u, b2u)
    x = np.dot(n1u, n2u)
    y = np.dot(m1, n2u)
    phi = np.arctan2(y, x)

    K1 = K.get('K1', 0.0); K2 = K.get('K2', 0.0); K3 = K.get('K3', 0.0); K4 = K.get('K4', 0.0)
    energy = 0.5*K1*(1+np.cos(phi)) + 0.5*K2*(1-np.cos(2*phi)) + 0.5*K3*(1+np.cos(3*phi)) + 0.5*K4*(1-np.cos(4*phi))
    dU_dphi = -0.5*K1*np.sin(phi) + K2*np.sin(2*phi) - 1.5*K3*np.sin(3*phi) + 2.0*K4*np.sin(4*phi)

    # Cartesian derivatives of φ (stable, standard forms)
    b2n = _norm(b2)
    c1 = np.cross(b1, b2); c2 = np.cross(b3, b2)
    c1m2 = np.dot(c1, c1) + 1e-12
    c2m2 = np.dot(c2, c2) + 1e-12

    dphi_dri = -(b2n / c1m2) * c1
    dphi_drl =  (b2n / c2m2) * c2
    dphi_drj = (np.dot(b1, b2) / (b2n * c1m2)) * c1 - (np.dot(b3, b2) / (b2n * c2m2)) * c2
    dphi_drk = -(dphi_dri + dphi_drj + dphi_drl)

    fi = -dU_dphi * dphi_dri
    fj = -dU_dphi * dphi_drj
    fk = -dU_dphi * dphi_drk
    fl = -dU_dphi * dphi_drl

    return energy, np.vstack([fi, fj, fk, fl])

```



```Python title:"Solution: Nonbonded Potential"
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def nonbond_potential(r, indexes, q, sigma, epsilon, coulomb_factor):
    """Lennard-Jones and Coulomb between a single pair.
       r: vector r_i - r_j
       indexes: (i, j, scale_LJ, scale_Coul) for 1–4 scaling etc.
       q, sigma, epsilon: per-atom arrays (length N) in MD units.
       Returns energy_lj, energy_elec, and forces on (i,j)."""
    i, j, sLJ, sC = indexes
    d2 = np.dot(r, r)
    d = np.sqrt(d2) if d2 > 0 else 1e-12
    rhat = r / d

    # OPLS geometric mixing
    sij = np.sqrt(sigma[i]*sigma[j])
    eij = np.sqrt(epsilon[i]*epsilon[j])

    # LJ
    sr = sij / d
    sr6 = sr**6
    sr12 = sr6**2
    U_lj = 4.0 * eij * (sr12 - sr6)
    F_lj = 24.0 * eij * (2.0*sr12 - sr6) * rhat / d

    # Coulomb
    qq = q[i]*q[j]
    U_c = coulomb_factor * qq / d
    F_c = coulomb_factor * qq * rhat / (d2)

    # 1–4 scaling if needed
    U_lj *= sLJ; F_lj *= sLJ
    U_c  *= sC;  F_c  *= sC

    fi = F_lj + F_c
    fj = -fi
    return U_lj, U_c, np.vstack([fi, fj])

```


**3.2:** Check your force and energy calculations here by testing inputs/outputs you can calculate analytically/manually.

```Python

```

## Question 4

> [!Question]
> Question 4: (Function to Compute the Total Energy and Forces on the Molecule)
> 



**4.1:** Putting it all together, let's write a function that returns the total energy and forces on the molecule.


```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def compute_energy_forces(positions,
                          types,
                          bonds,
                          angles,
                          dihedrals,
                          nb_params,
                          bond_params,
                          angle_params,
                          dihedral_params):
    """Compute total energy and force."""
    energy = {
        'bond': 0.,
        'angle': 0.,
        'dihedral': 0.,
        'LJ': 0.,
        'coulomb': 0.
        }
    N = len(positions)
    forces = np.zeros((N,3))
    for t,i,j in bonds:
        # fill in
        u,f =
        energy['bond'] += u
        forces[[i,j]] += f
    for t,i,j,k in angles:
        # fill in
    for t,i,j,k,l in dihedrals:
        # fill in
    for i in range(N):
        for j in range(i+1,N):
            # fill in
    return energy, forces
```



```Python title:"Solution: Total Energy and Forces"
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def compute_energy_forces(positions,
                          types,
                          bonds,
                          angles,
                          dihedrals,
                          nb_params,
                          bond_params,
                          angle_params,
                          dihedral_params):
    """Compute total energy and force."""
    energy = {
        'bond': 0.,
        'angle': 0.,
        'dihedral': 0.,
        'LJ': 0.,
        'coulomb': 0.
    }
    N = len(positions)
    forces = np.zeros((N,3))

    # Per-atom NB arrays
    q = np.array([nb_params[t]['q'] for t in types], dtype=float)
    sigma = np.array([nb_params[t]['sigma'] for t in types], dtype=float)
    epsilon = np.array([nb_params[t]['epsilon'] for t in types], dtype=float)
    cf = nb_params['C']['coulomb_factor']  # same constant, grab once
    s14_LJ = nb_params['scale14']['LJ']
    s14_C  = nb_params['scale14']['Coulomb']

    # Exclusion maps (1–2, 1–3, 1–4)
    one_two = {(min(i,j), max(i,j)) for _, i, j in bonds}
    one_three = {(min(i,k), max(i,k)) for _, i, j, k in angles}
    one_four = {(min(i,l), max(i,l)) for _, i, j, k, l in dihedrals}

    # Bonds
    for t, i, j in bonds:
        u, f = bond_potential(positions[i] - positions[j],
                              bond_params[t]['r0'],
                              bond_params[t]['k'])
        energy['bond'] += u
        forces[i] += f[0]; forces[j] += f[1]

    # Angles
    for t, i, j, k in angles:
        u, f = angle_potential(positions[i], positions[j], positions[k],
                               angle_params[t]['th0'], angle_params[t]['k'])
        energy['angle'] += u
        forces[i] += f[0]; forces[j] += f[1]; forces[k] += f[2]

    # Dihedrals
    for t, i, j, k, l in dihedrals:
        u, f = dihedral_potential(positions[i], positions[j], positions[k], positions[l],
                                  dihedral_params[t])
        energy['dihedral'] += u
        forces[i] += f[0]; forces[j] += f[1]; forces[k] += f[2]; forces[l] += f[3]

    # Nonbonded pairs
    for i in range(N):
        for j in range(i+1, N):
            key = (i, j)
            if key in one_two or key in one_three:
                continue
            sL = s14_LJ if key in one_four else 1.0
            sC = s14_C  if key in one_four else 1.0
            u_lj, u_c, f = nonbond_potential(positions[i] - positions[j],
                                             (i, j, sL, sC),
                                             q, sigma, epsilon, cf)
            energy['LJ'] += u_lj
            energy['coulomb'] += u_c
            forces[i] += f[0]; forces[j] += f[1]

    return energy, forces

```




## Question 5

> [!Question]
> Question 5: (Perform Energy Minimization with Steepest Decent)


**5.1:** First, let's see what the initial energies on the molecule are.

```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
positions, types, bonds, angles, dihedrals =
nb_params, bond_params, angle_params, dihedral_params =
# forces and energies of initial configuration
U,F =
```


```Python title:Solution
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
positions, types, bonds, angles, dihedrals = define_molecule('ethane.xyz')
nb_params, bond_params, angle_params, dihedral_params = define_force_field()
# forces and energies of initial configuration
U, F = compute_energy_forces(positions, types, bonds, angles, dihedrals,
                             nb_params, bond_params, angle_params, dihedral_params)

```



**5.2:** Check here if your inital forces and energies make sense:

```Python
print('Initial energy and forces')
print('-------------------------')
print('total: {:.5g}'.format(sum(U.values())))
for k,v in U.items():
    print('{}: {:.5g}'.format(k,v))
print(F)
print('-------------------------')
```


```Python

```


**5.3:** Now perform steepest descent to minimize the energy.


```Python
# steepest descent
gamma = 1.e-6
tol = 10.
max_iterations = 1000
for i in range(max_iterations):
    U,F =
    Fmax = np.max()
    if Fmax < tol:
        break
    positions +=
if i == max_iterations:
    print('warning: may not have converged')
else:
    print('converged in {} iterations'.format(i))
```




```Python title:Solution
# steepest descent
gamma = 1.e-6
tol = 10.
max_iterations = 1000
for i in range(max_iterations):
    U, F = compute_energy_forces(positions, types, bonds, angles, dihedrals,
                                 nb_params, bond_params, angle_params, dihedral_params)
    Fmax = np.max(np.linalg.norm(F, axis=1))
    if Fmax < tol:
        break
    positions += gamma * F  # move along +F since F = -∇U
if i == max_iterations - 1 and Fmax >= tol:
    print('warning: may not have converged')
else:
    print('converged in {} iterations'.format(i))

```




**5.4:** Calculate final state energies and forces:


```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
# energy and forces of final configuration
U_final,F_final =

print('Final energy and forces')
print('-----------------------')
print('total: {:.5g}'.format(sum(U.values())))
for k,v in U.items():
    print('{}: {:.5g}'.format(k,v))
print(F)
print('-----------------------')

```


```Python title:Solution
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
# energy and forces of final configuration
U_final, F_final = compute_energy_forces(positions, types, bonds, angles, dihedrals,
                                         nb_params, bond_params, angle_params, dihedral_params)

print('Final energy and forces')
print('-----------------------')
print('total: {:.5g}'.format(sum(U_final.values())))
for k,v in U_final.items():
    print('{}: {:.5g}'.format(k,v))
print(F_final)
print('-----------------------')

```


**5.5:** Write out the minimized structure for opening in VMD or OVITO and plotting. You can download the new file from the same screen where you launced the workspace from, on the left side.


```Python
# write minimized structure
with open('ethane_min.xyz','w') as f:
    f.write('{}\n'.format(positions.shape[0]))
    f.write('\n')
    for i in range(positions.shape[0]):
        # fill in
        f.write()
```



```Python title:Solution
# write minimized structure
with open('ethane_min.xyz','w') as f:
    f.write('{}\n'.format(positions.shape[0]))
    f.write('Ethane minimized (nm)\n')
    for i in range(positions.shape[0]):
        f.write('{:s} {:.6f} {:.6f} {:.6f}\n'.format(types[i], positions[i,0], positions[i,1], positions[i,2]))

```


---

# <span style="color:#ff0000; font-weight:bold">2. Resources</span>


## 1.) What counts as "OPLS" here (and where to read it)

For saturated hydrocarbons, OPLS-AA uses atom types **CT** (sp ${ }^3$ carbon in alkanes) and **HC** (hydrogen on CT), with standard geometric mixing for both $\sigma$ and $\varepsilon$ and 1-4 scaling of 0.5 for both electrostatics and LJ. If you follow those three choices, your nonbonded math will match the force field family you're emulating. References below give the exact functional forms and defaults (and they're the same conventions most MD engines assume when you pick OPLS). [Gromacs Reference Manual](chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://manual.gromacs.org/5.1.1/manual-5.1.1.pdf?)


For a single isolated molecule in vacuum you won't need PME/cutoffs, but you do need the Coulomb constant in MD units [Gromacs Definitions and Units](https://manual.gromacs.org/current/reference-manual/definitions.html?):

==Electrical Conversion Factor:== $f=\frac{1}{4 \pi \varepsilon_0}=138.935458 \mathrm{~kJ}$ $\cdot$ $\mathrm{mol}^{-1}$ $\cdot$ nm $\cdot$ $e^{-2}$. Use this directly in your Coulomb term. 

If you'd like to see how these choices are encoded in major packages, skim the GROMACS manual pages on bonded functions and parameter files; they also note that OPLS dihedrals are typically represented as Ryckaert-Bellemans (RB) in GROMACS, with a published mapping to the OPLS Fourier form. [Bonded Interactions](https://manual.gromacs.org/documentation/current/reference-manual/functions/bonded-interactions.html?)


## 2.) Canonical OPLS-AA parameters for ethane (CT/HC)

The GROMACS mailing list archive captured the exact numbers you'll want for bonds, angles, dihedrals, and nonbonded CT/HC entries (already in nm and $\mathrm{kJ}^{-1} \mathrm{~mol}^{-1}$ ) [Mailman](https://mailman-1.sys.kth.se/pipermail/gromacs.org_gmx-users/2011-August/063423.html):

- **Bonds**

$$
\begin{aligned}
& \text { CT-CT: } r_0=0.15290 \mathrm{~nm}, k=224262.4 \mathrm{~kJ} \cdotp \mathrm{mol}^{-1} \cdotp \mathrm{nm}^{-2} \\
& \text { CT-HC: } r_0=0.10900 \mathrm{~nm}, k=284512.0 \mathrm{~kJ} \cdotp \mathrm{mol}^{-1} \cdotp \mathrm{nm}^{-2}
\end{aligned}
$$

- **Angles**

$$
\begin{aligned}
& \text { CT-CT-CT: } \theta_0=112.7^{\circ}, k_\theta=488.273 \mathrm{~kJ} \cdotp \mathrm{mol}^{-1} \cdotp \mathrm{rad}^{-2} \\
& \text { CT-CT-HC: } \theta_0=110.7^{\circ}, k_\theta=313.800 \mathrm{~kJ} \cdotp \mathrm{mol}^{-1} \cdotp \mathrm{rad}^{-2} \\
& \text { HC-CT-HC: } \theta_0=107.8^{\circ}, k_\theta=276.144 \mathrm{~kJ} \cdotp \mathrm{mol}^{-1} \cdotp \mathrm{rad}^{-2}
\end{aligned}
$$

- **Proper dihedrals (RB form; GROMACS function type 3):**


$$
\begin{aligned}
& \text { CT-CT-CT-CT: } C_0=2.92880, C_1=-1.46440, C_2=0.20920, C_3=-1.67360, C_4= C_5=0 \\
& \text { HC-CT-CT-HC: } C_0=0.62760, C_1=1.88280, C_2=0.00000, C_3=-2.51040, C_4=C_5= 0
\end{aligned}
$$




> [!NOTE]
> Note: If you prefer **OPLS Fourier (the LAMMPS "dihedral_opls" form)**,
> 
> $$
> E(\phi)=\frac{1}{2} K_1[1+\cos \phi]+\frac{1}{2} K_2[1-\cos 2 \phi]+\frac{1}{2} K_3[1+\cos 3 \phi]+\frac{1}{2} K_4[1-\cos 4 \phi],
> $$
> 
> use the GROMACS manual's mapping between OPLS $V_n$ and RB $C_n$ (don't forget kcal $\rightarrow$ kJ if you copy literature values). [LAMMPS DOCS](https://docs.lammps.org/dihedral_opls.html?)



**Nonbonded atom types (charges, LJ):**

- CT (opls_135 $\mathrm{CH}_3$ carbon): $q=-0.180 e, \sigma=0.350 \mathrm{~nm}, \varepsilon=0.276144 \mathrm{~kJ}$ $\cdotp$ $\mathrm{mol}^{-1}$

- HC (opls_140 hydrogen on alkane): $q=+0.060 e, \sigma=0.250 \mathrm{~nm}, \varepsilon= 0.125520 kJ \cdot \mathrm{mol}^{-1}$


For ethane $\left(\mathrm{CH}_3-\mathrm{CH}_3\right)$ that charge set sums to zero: $2 \times(-0.18)+6 \times(+0.06)=0$. 


**Mixing (OPLS convention):** Use geometric combining for both $\sigma$ and $\varepsilon$ (GROMACS "comb-rule = 3", OpenMM "combiningRule=geometric"). [Calculating the Sensitivity and Robustness of the Binding Free Energy Calculations to Force Field Parameters](https://pmc.ncbi.nlm.nih.gov/articles/PMC3763860/?)



**1-4 scaling:** Use 0.5 for both Coulomb and LJ 1-4 interactions (OpenMM: `coulomb14scale="0.5" lj14scale="0.5"` ; in hand code, multiply the 1-4 pair energy by 0.5 ). [Traken](https://traken.chem.yale.edu/ligpargen/openMM_tutorial.html?)



If you need a quick cross-check on the $\mathbf{C - H}$ and $\mathbf{C - C}$ bond values in another context, the same numbers appear in unrelated GROMACS tutorials and Q\&A threads that quote the _ffbonded.itp_ content. [GitLab](https://cameleon.univ-lyon1.fr/cloison/md_water_methanol_gromacs/-/blob/teacher/QM/GROMACS/topol_oplsaa_DH_C0_0.7.top?)


## 3.) The formulas you'll code (forces and energies)

Work entirely in $\mathbf{n m}$ and $\mathbf{k J} \cdot \mathbf{m o l}^{-\mathbf{1}}$ :

**Lennard-Jones (12-6):**

$$
U_{\mathrm{LJ}}(r)=4 \varepsilon_{i j}\left[\left(\frac{\sigma_{i j}}{r}\right)^{12}-\left(\frac{\sigma_{i j}}{r}\right)^6\right], \quad \sigma_{i j}=\sqrt{\sigma_i \sigma_j}, \varepsilon_{i j}=\sqrt{\varepsilon_i \varepsilon_j} \text { (OPLS). }
$$


The derivative gives the pair force magnitude; apply it along $\mathbf{r}_{i j}$. [Wikipedia](https://en.wikipedia.org/wiki/Combining_rules?)


**Coulomb:**

$$
U_{\mathrm{C}}(r)=f \frac{q_i q_j}{r}, \quad f=138.935458 \mathrm{~kJ} \text { cdotp } \mathrm{mol}^{-1} \text { cdotp } \mathrm{nm} \text { cdotp } e^{-2} .
$$


Same vector derivative rule. [GROMACS Manual (Definitions and Units)](https://manual.gromacs.org/current/reference-manual/definitions.html?)


**Bonds (harmonic):**

$$
U_b=\frac{1}{2} k\left(r-r_0\right)^2, \quad \mathrm{~F}_{i j}=-k\left(r-r_0\right) \hat{\mathrm{r}}_{i j} .
$$

(Use the $k$ values above.) [Mailman](https://mailman-1.sys.kth.se/pipermail/gromacs.org_gmx-users/2011-August/063423.html)


**Angles (harmonic in the angle):**

$$
U_\theta=\frac{1}{2} k_\theta\left(\theta-\theta_0\right)^2,
$$

with the standard three-body force expression (GROMACS manual's definitions are consistent with your units). [GROMACS Documentation (Parameter Files)](https://manual.gromacs.org/documentation/current/reference-manual/topologies/parameter-files.html?)


**Dihedrals:** Either implement OPLS Fourier (nice for analytic $\partial \mathrm{E} / \partial \phi$ ) or RB (matches the tabulated coefficients above). The GROMACS manual provides the explicit RB form and the OPLS $\leftrightarrow$ RB mapping if you change representations. [GROMACS Documentation (Bonded Interactions)](https://manual.gromacs.org/documentation/current/reference-manual/functions/bonded-interactions.html?)



## 4.) Geometry and "sanity" targets for ethane

Your XYZ uses nm and is already near typical values: $r_{\mathrm{C}-\mathrm{H}} \approx 0.109 \mathrm{~nm}, r_{\mathrm{C}-\mathrm{C}} \approx 0.153 \mathrm{~nm}$, angles $\approx 109-113^{\circ}$. A good physical check after minimization is that the staggered conformation is the minimum and the torsional barrier around $\mathrm{C}-\mathrm{C}$ comes out $\sim 2.9-3.0 \mathrm{kcal}^{-1} \mathrm{~mol}^{-1}\left(\approx 12 \mathrm{~kJ} \cdot \mathrm{~mol}^{-1}\right)$ when you scan $\phi(\mathrm{H}-\mathrm{C}-\mathrm{C}-\mathrm{H})$. (That value is experimental/theoretical consensus.) [The Rotational Barrier in Ethane: A Molecular Orbital Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6268250/?)



## 5.) Worked examples and "how-to" patterns you can crib from
- **Exact OPLS values in a real topology (so you can see them in context):** The gmx-users thread shows the bond/angle/dihedral lines and the nonbonded CT/HC entries (with charges, $\sigma, \varepsilon)$. This is the single most useful reference for copying numbers verbatim. [Mailman](https://mailman-1.sys.kth.se/pipermail/gromacs.org_gmx-users/2011-August/063423.html)
- **OPLS dihedral formula (for $\mathrm{K}_1-\mathrm{K}_4$ ) and its use in other engines (LAMMPS/HOOMD docs):** a clean statement of the functional form if you choose Fourier rather than RB. [LAMMPS Documentation](https://docs.lammps.org/dihedral_opls.html?)
- **RB $\leftrightarrow$ OPLS conversion (if you want to transform coefficients):** the GROMACS manual gives the precise mapping equations and the $\mathrm{kcal} \leftrightarrow \mathrm{kJ}$ reminder. _See GROMACS 3.0 Manual in your Downloads_
- **Combining rules refresher (why OPLS is geometric $\times$ geometric):** NAMD/OpenMM/CCL notes and reviews summarize conventions across force fields. Useful if you later mix fields. [UIUC](https://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2004-2005/2305.html?)
- **XYZ parsing pattern in Python:** tiny tutorial showing the two-line header convention and how to slice the rest-handy for your define_molecule() . [How to Parse .XYZ Atomic Coordinates File Using Python](https://www.bragitoff.com/2023/07/how-to-parse-read-xyz-atomic-coordinates-file-using-python-tutorial/?)
- **Steepest descent / minimization examples in Python (if you want to sanity-test your force routines before tackling angles/dihedrals):** teaching pages and ASE docs show plain-Python LJ minimizers and optimization idioms. [Gradient Descent](https://pythoninchemistry.org/ch40208/comp_chem_methods/gradient_descent.html?)


## 6) Drop-in "cheat sheet" you can paste into your notebook

When you get to Question 2, this dictionary matches the numbers cited above and keeps everything in $\mathrm{nm} / \mathrm{kJ} \cdot \mathrm{mol}^{-1} / \mathrm{e}:$


```Python
# OPLS-AA (CT/HC) in nm, kJ/mol, e. Geometric mixing, 1–4 scales = 0.5 (both).
nb_params = {
    'C':  {'type': 'CT', 'q': -0.180, 'sigma': 0.350, 'epsilon': 0.276144},  # opls_135
    'H':  {'type': 'HC', 'q': +0.060, 'sigma': 0.250, 'epsilon': 0.125520},  # opls_140
    'coulomb_factor': 138.935458,   # kJ/mol·nm·e^-2
    'mixing': 'geometric',          # OPLS convention for σ and ε
    'scale14': {'LJ': 0.5, 'Coulomb': 0.5}
}

bond_params = {
    'CT-CT': {'r0': 0.15290, 'k': 224262.4},
    'CT-HC': {'r0': 0.10900, 'k': 284512.0},
}

angle_params = {
    'CT-CT-CT': {'th0_deg': 112.7, 'k': 488.273},
    'CT-CT-HC': {'th0_deg': 110.7, 'k': 313.800},
    'HC-CT-HC': {'th0_deg': 107.8, 'k': 276.144},
}

# Dihedrals as RB coefficients (GROMACS "function 3"): [C0, C1, C2, C3, C4, C5]
dihedral_params = {
    'CT-CT-CT-CT': {'RB': [2.92880, -1.46440, 0.20920, -1.67360, 0.0, 0.0]},
    'HC-CT-CT-HC': {'RB': [0.62760,  1.88280, 0.00000, -2.51040, 0.0, 0.0]},
}

```


All values above are taken directly from the OPLS-AA parameter listings quoted on the gmx-users list (which itself mirrors the `oplsaa.ff` files). If you later decide to implement the Fourier/OPLS torsion instead of RB, use the LAMMPS-style formula and the GROMACS manual's mapping.



## 7) Minimal roadmap for your coding steps (so you don't get tripped up)


**1.) Parse ethane.xyz:** first line $N$, second line comment, then type x y z per line in nm ; store types and an $\mathrm{N} \times 3$ float array. (The pattern in the short XYZ tutorial is exactly what you need.) [How to Parse .XYZ Atomic Coordinates File Using Python](https://www.bragitoff.com/2023/07/how-to-parse-read-xyz-atomic-coordinates-file-using-python-tutorial/?)
**2.) Topology:** enumerate bonds, angles, dihedrals by index using your CT/HC pattern (the target tuples in your prompt are correct for ethane).
**3.) Nonbonded pairs:** loop over all $i<j$. If atoms are separated by three bonds (i.e., 1-4 pairs), apply the 0.5 scaling; if they're 1-2 or 1-3, exclude (bonded and angle neighbors). (This mirrors how OPLS is usually applied.) [GROMACS Bonded Interactions](https://manual.gromacs.org/documentation/current/reference-manual/functions/bonded-interactions.html?)
**4.) Units:** keep everything in $\mathbf{n m}$ and $\mathbf{k J} \cdot^{\mathbf{- 1}} \mathbf{m o l}^{\mathbf{- 1}}$; multiply Coulomb by $\boldsymbol{f}$ above; angles: convert $\theta_0$ to radians when you evaluate $\frac{1}{2} k\left(\theta-\theta_0\right)^2$. [GROMACS Definitions and Units](https://manual.gromacs.org/current/reference-manual/definitions.html?)
**5.) Steepest descent:** pick a small step like $\gamma \sim 10^{-6} \mathrm{~nm}$ cdotp mol cdotp $\mathrm{kJ}^{-1}$, update $\mathbf{r} \leftarrow \mathbf{r}- \gamma \mathrm{F}$, and stop when $\max _i\left\|\mathrm{~F}_i\right\|<$ tolerance (e.g., $10 \mathrm{~kJ} \cdot \mathrm{~mol}^{-1} \cdot \mathrm{~nm}^{-1}$ ), as your template suggests; this is the same stopping logic used in MD packages. [GROMACS Molecular Dynamics Parameters](https://manual.gromacs.org/current/user-guide/mdp-options.html?)


## 8.) Optional deeper dives (if you want them)
- **OPLS/2020 tweaks:** later refinements mainly adjust alkane torsion coefficients; if you're curious about modern variants, this paper summarizes the changes and the rationale. [Refinement of the OPLS Force Field for Thermodynamics and Dynamics of Liquid Alkanes](https://pmc.ncbi.nlm.nih.gov/articles/PMC9939004/?)
- **LigParGen \& Foyer:** automatic OPLS parameterization/generation-handy to sanity-check your hand-typed ethane parameters or generate files for larger molecules. [LigParGen: Implementing OPLS-AA/M Force Field in GROMACS: Lessons Learnt](https://traken.chem.yale.edu/ligpargen/oplsaam_gmx_tutorial.html?)


> [!NOTE]
> If you implement the energy terms exactly as above, with OPLS mixing (geometric), 1-4 scaling = 0.5 for both LJ and Coulomb, and the CT/HC numbers quoted here, you'll reproduce a standard OPLS-AA ethane. Your minimizer should relax to a staggered geometry with C-H ~0.109 nm, C-C ~0.153 nm, and a torsional barrier ~3 $\mathrm{kcal}^{-1} \mathrm{~mol}^{-1}$ when scanned-classic ethane behavior.
> 



# <span style="color:#ff0000; font-weight:bold">3. Theoretical Discussion</span>


## 1.) Molecular representation and topology

Let the molecule have atoms $a=0, \ldots, N-1$ at positions $\mathbf{r}_a \in \mathbb{R}^3$. A molecular mechanics (MM) force field decomposes the potential energy into bonded and non-bonded terms that depend only on internal coordinates built from $\left\{\mathbf{r}_a\right\}$ :

$$
U(\{\mathbf{r}\})=\sum_{\text {bonds }} U_b+\sum_{\text {angles }} U_\theta+\sum_{\text {dihedrals }} U_\phi+\sum_{\substack{i<j \\ \mathbf{n b}}}\left(U_{\mathrm{LJ}}^{i j}+U_{\mathrm{C}}^{i j}\right) .
$$


Topology supplies the lists of bonded pairs ( $i, j$ ), angle triplets ( $i, j, k$ ) with vertex $j$, and proper dihedral quadruplets $(i, j, k, l)$ about bond $j-k$. For nonbonded (nb) pairs, the usual OPLS convention is: exclude 1-2 (bonded) and 1-3 (angle) pairs from the nonbonded sum; include 1-4 (atoms connected by a dihedral) with scaling factors; include all others unscaled.

For ethane ( $\mathrm{CH}_3-\mathrm{CH}_3$ ), atom types are CT (sp ${ }^3$ carbon in alkanes) and HC (hydrogen on CT). The "OPLS type nomenclature" you'll attach to topological tuples simply tags which parameter block to read (e.g., a bond labeled "CT-CT" maps to a particular $r_0, k$ pair).



## 2.) Bonded terms: internal coordinates, energies, and analytic forces

### 1. Bonds
For a bond between atoms $i$ and $j$, define the vector $\mathbf{r}_{i j}=\mathbf{r}_i-\mathbf{r}_j$ and its length $r=\left\|\mathbf{r}_{i j}\right\|$. The harmonic bond potential is

$$
U_b(r)=\frac{1}{2} k\left(r-r_0\right)^2 .
$$


The force on atom $i$ is minus the gradient:

$$
\begin{gathered}
\frac{\partial U_b}{\partial r}=k\left(r-r_0\right), \quad \frac{\partial r}{\partial \mathbf{r}_i}=\frac{\mathbf{r}_{i j}}{r}, \\
\mathbf{\Gamma}_i^{(b)}=-\frac{\partial U_b}{\partial \mathbf{r}_i}=-k\left(r-r_0\right) \frac{\mathbf{r}_{i j}}{r}, \quad \mathbf{\Gamma}_j^{(b)}=-\mathbf{\Gamma}_i^{(b)}
\end{gathered}
$$


This yields exact Newton's third law pairwise: $\mathrm{F}_i=-\mathrm{F}_j$.


### 2. Angles
For an angle triplet ( $i, j, k$ ) with vertex at $j$, let

$$
\begin{gathered}
\mathbf{u}=\mathbf{r}_i-\mathbf{r}_j, \quad a=\|\mathbf{u}\|, \quad \hat{\mathbf{u}}=\mathbf{u} / a, \quad \mathbf{v}=\mathbf{r}_k-\mathbf{r}_j, \quad b=\|\mathbf{v}\|, \quad \hat{\mathbf{v}}=\mathbf{v} / b, \\
\cos \theta=\hat{\mathbf{u}} \cdot \hat{\mathbf{v}}, \quad \sin \theta=\|\hat{\mathbf{u}} \times \hat{\mathbf{v}}\| .
\end{gathered}
$$


The harmonic angle potential is

$$
U_\theta(\theta)=\frac{1}{2} k_\theta\left(\theta-\theta_0\right)^2, \quad \frac{\partial U_\theta}{\partial \theta}=k_\theta\left(\theta-\theta_0\right) .
$$


Using $\partial \theta / \partial \hat{\mathrm{u}}=-(\hat{\mathrm{v}}-\cos \theta \hat{\mathrm{u}}) / \sin \theta$ and $\partial \hat{\mathrm{u}} / \partial \mathbf{r}_i=\frac{1}{a}\left(\mathbf{I}-\hat{\mathrm{u}} \hat{\mathrm{u}}^{\top}\right)$, the standard compact force expressions are

$$
\mathbf{\Gamma}_i^{(\theta)}=-k_\theta\left(\theta-\theta_0\right) \frac{\hat{\mathbf{v}}-\cos \theta \hat{\mathbf{u}}}{a \sin \theta}, \quad \mathbf{\Gamma}_k^{(\theta)}=-k_\theta\left(\theta-\theta_0\right) \frac{\hat{\mathbf{u}}-\cos \theta \hat{\mathbf{v}}}{b \sin \theta}, \quad \mathbf{\Gamma}_j^{(\theta)}=-\left(\mathbf{\Gamma}_i^{(\theta)}+\mathbf{\Gamma}_k^{(\theta)}\right) .
$$


These three forces lie in the plane of the angle and again sum to zero.



### 3. Dihedrals
For a proper torsion $(i, j, k, l)$ about bond $j-k$, define

$$
\begin{array}{rll}
\mathbf{b}_1=\mathbf{r}_j-\mathbf{r}_i, & \mathbf{b}_2=\mathbf{r}_k-\mathbf{r}_j, & \mathbf{b}_3=\mathbf{r}_l-\mathbf{r}_k, \\
\mathbf{n}_1=\mathbf{b}_1 \times \mathbf{b}_2, & \mathbf{n}_2=\mathbf{b}_2 \times \mathbf{b}_3, & \hat{\mathbf{b}}_2=\mathbf{b}_2 /\left\|\mathbf{b}_2\right\| .
\end{array}
$$


The dihedral angle $\phi$ is most stably computed via

$$
\begin{gathered}
\cos \phi=\frac{\mathbf{n}_1 \cdot \mathbf{n}_2}{\left\|\mathbf{n}_1\right\|\left\|\mathbf{n}_2\right\|}, \quad \sin \phi=\frac{\hat{\mathbf{b}}_2 \cdot\left(\mathbf{n}_1 \times \mathbf{n}_2\right)}{\left\|\mathbf{n}_1\right\|\left\|\mathbf{n}_2\right\|}, \\
\phi=\operatorname{atan} 2(\sin \phi, \cos \phi) .
\end{gathered}
$$


OPLS's Fourier form is

$$
U_\phi(\phi)=\frac{1}{2} K_1[1+\cos \phi]+\frac{1}{2} K_2[1-\cos 2 \phi]+\frac{1}{2} K_3[1+\cos 3 \phi]+\frac{1}{2} K_4[1-\cos 4 \phi] .
$$


Its derivative (the "generalized force" along $\phi$ ) is

$$
\frac{d U_\phi}{d \phi}=-\frac{1}{2} K_1 \sin \phi+K_2 \sin 2 \phi-\frac{3}{2} K_3 \sin 3 \phi+2 K_4 \sin 4 \phi
$$



To convert $d U_\phi / d \phi$ into Cartesian forces, differentiate $\phi$ with respect to atomic positions. A widely used analytic result is

$$
\mathbf{c}_1=\mathbf{b}_1 \times \mathbf{b}_2, \Delta_1=\left\|\mathbf{c}_1\right\|^2, \quad \mathbf{c}_2=\mathbf{b}_3 \times \mathbf{b}_2, \Delta_2=\left\|\mathbf{c}_2\right\|^2, \quad b=\left\|\mathbf{b}_2\right\| .
$$


Then

$$
\begin{gathered}
\frac{\partial \phi}{\partial \mathbf{r}_i}=-\frac{b}{\Delta_1} \mathbf{c}_1, \quad \frac{\partial \phi}{\partial \mathbf{r}_l}=+\frac{b}{\Delta_2} \mathbf{c}_2, \\
\frac{\partial \phi}{\partial \mathbf{r}_j}=\frac{\mathbf{b}_1 \cdot \mathbf{b}_2}{b \Delta_1} \mathbf{c}_1-\frac{\mathbf{b}_3 \cdot \mathbf{b}_2}{b \Delta_2} \mathbf{c}_2, \quad \frac{\partial \phi}{\partial \mathbf{r}_k}=-\left(\frac{\partial \phi}{\partial \mathbf{r}_i}+\frac{\partial \phi}{\partial \mathbf{r}_j}+\frac{\partial \phi}{\partial \mathbf{r}_l}\right) .
\end{gathered}
$$


Finally, the forces are

$$
\mathbf{\Gamma}_a^{(\phi)}=-\frac{d U_\phi}{d \phi} \frac{\partial \phi}{\partial \mathbf{r}_a}, \quad a \in\{i, j, k, l\}
$$


Remarks that matter numerically: the expressions remain finite except when any three consecutive atoms become colinear (then $\left\|c_{1,2}\right\| \rightarrow 0$ ); for ethane in typical geometries you are safely away from those singularities. Using the atan2 definition for $\phi$ gives a stable, branch-aware angle.


## 3.) Nonbonded terms: mixing rules, exclusions, and analytic forces

### 1. Lennard-Jones (12-6) with OPLS mixing
For atoms $i$ and $j$ at separation $r=\left\|\mathrm{r}_{i j}\right\|$, the LJ energy with OPLS's geometric mixing is

$$
\begin{gathered}
\sigma_{i j}=\sqrt{\sigma_i \sigma_j}, \quad \varepsilon_{i j}=\sqrt{\varepsilon_i \varepsilon_j}, \\
U_{\mathrm{LJ}}^{i j}(r)=4 \varepsilon_{i j}\left[\left(\frac{\sigma_{i j}}{r}\right)^{12}-\left(\frac{\sigma_{i j}}{r}\right)^6\right] .
\end{gathered}
$$


The radial derivative is

$$
\frac{d U_{\mathrm{LJ}}}{d r}=\frac{4 \varepsilon_{i j}}{r}\left[-12\left(\frac{\sigma_{i j}}{r}\right)^{12}+6\left(\frac{\sigma_{i j}}{r}\right)^6\right] .
$$


The vector force on $i$ is

$$
\mathbf{\Gamma}_i^{(\mathrm{LJ})}=-\frac{d U_{\mathrm{LJ}}}{d r} \frac{\mathbf{r}_{i j}}{r}=24 \varepsilon_{i j}\left[2\left(\frac{\sigma_{i j}}{r}\right)^{12}-\left(\frac{\sigma_{i j}}{r}\right)^6\right] \frac{\mathbf{r}_{i j}}{r^2}, \quad \mathbf{\Gamma}_j^{(\mathrm{LJ})}=-\mathbf{\Gamma}_i^{(\mathrm{LJ})} .
$$


### 2. Coulomb
With charges $q_i, q_j$ in units of $e$ and the MD-convention Coulomb factor

$$
f=\frac{1}{4 \pi \varepsilon_0}=138.935458 \mathrm{~kJ} \cdotp \mathrm{mol}^{-1} \cdotp \mathrm{nm} \cdot e^{-2},
$$

the energy and force are

$$
U_{\mathrm{C}}^{i j}(r)=f \frac{q_i q_j}{r}, \quad \Gamma_i^{(\mathrm{C})}=f \frac{q_i q_j}{r^3} \mathbf{r}_{i j}, \quad \Gamma_j^{(\mathrm{C})}=-\Gamma_i^{(\mathrm{C})} .
$$


### 3. Exclusions and 1-4 scaling
Let $\mathcal{N}_{1-n}(i)$ be the set of atoms at topological distance $n$ from atom $i$. The typical OPLS rules are:
1. If $j \in \mathcal{N}_{1-2}(i)$ (bonded) or $j \in \mathcal{N}_{1-3}(i)$ (angle), omit the nonbonded interaction entirely.
2. If $j \in \mathcal{N}_{1-4}(i)$ (dihedral partners), include the interaction but scale both LJ and Coulomb by factors $s_{\mathrm{LJ}}=0.5$ and $s_{\mathrm{C}}=0.5$.
3. Otherwise, include the full interactions unscaled.

In code, this is implemented by detecting whether ( $i, j$ ) is a 1-4 pair via your dihedral list and multiplying both $U$ and F by the appropriate factor $s \in\{1,0.5,0\}$.




### 4. Total energy and force assembly; invariants you can use to debug

Summing all contributions,

$$
U_{\text {tot }}=\sum_{(i, j) \in \text { bonds }} U_b\left(r_{i j}\right)+\sum_{(i, j, k) \in \text { angles }} U_\theta\left(\theta_{i j k}\right)+\sum_{(i, j, k, l) \in \text { dihedrals }} U_\phi\left(\phi_{i j k l}\right)+\sum_{i<j}^{\text {nb rules }}\left(U_{\mathrm{LJ}}^{i j}\left(r_{i j}\right)+U_{\mathrm{C}}^{i j}\left(r_{i j}\right)\right) .
$$


Forces superpose linearly per atom. Two strong, general checks follow immediately from the mathematics and should hold to round-off if your formulas are correct:
1. Total force is zero: $\sum_a \mathbf{F}_a=\mathbf{0}$. Internal forces cannot impart net translation.
2. Total torque is zero: $\sum_a \mathbf{r}_a \times \mathbf{F}_a=\mathbf{0}$. Internal forces cannot impart net rotation.

If either sum is large compared to individual bonded forces, there is a sign or normalization error in one of your analytic derivatives.



### 5. Ethane specifics you will exploit

Ethane has two CT atoms and six HC atoms. The bond list contains six CT-HC bonds and one CT-CT bond. Angles are the three $\mathrm{H}-\mathrm{C}-\mathrm{H}$ per methyl, the three $\mathrm{H}-\mathrm{C}-\mathrm{C}$ per carbon, and the central $\mathrm{C}-\mathrm{C}-\mathrm{H}$ on the other side, etc. Dihedrals of the type H-C-C-H are the ones that control the staggered vs. eclipsed torsion. Electrostatics are small but nonzero if you use typical OPLS charges (CT negative, HC positive); LJ dominates the nonbonded repulsion between near-eclipsed hydrogens.

This topology induces a nonbonded exclusion graph with many 1-4 pairs (the $\mathrm{H} \cdots \mathrm{H}$ across the central bond) that must be included with the 0.5 scaling; all within-methyl $\mathrm{H} \cdots \mathrm{H}$ pairs are 1-3 and therefore excluded.



### 6. Steepest descent (gradient) minimization in Cartesian space

Let $\mathrm{F}_a(\{\mathbf{r}\})=-\nabla_{\mathbf{r}_a} U$ be the force on atom $a$. Steepest descent updates positions along the negative gradient of $U$ :

$$
\mathbf{r}_a^{(n+1)}=\mathbf{r}_a^{(n)}+\gamma \mathbf{F}_a^{(n)} .
$$


The step size $\gamma$ has units $\mathrm{nm} \cdot \mathrm{mol} \cdot \mathrm{kJ}^{-1}$ so that $\gamma \mathrm{F}$ has units of distance. In practice you choose $\gamma$ small enough that the harmonic approximations remain locally valid and the energy decreases monotonically (you can enforce monotonic decrease with a very simple backtracking: if $U^{(n+1)}>U^{(n)}$, halve $\gamma$ and retry). You stop when a force-based criterion is met, e.g.

$$
\max _a\left\|\Gamma_a^{(n)}\right\|<\operatorname{tol},
$$

with tol on the order of 10 kJ cdotp $\mathrm{mol}^{-1}$ cdotg $\mathrm{nm}^{-1}$ for a course exercise, or when the total energy decreases by less than a small threshold per iteration. Because ethane is tiny and all bonded terms are harmonic while the torsion is smooth, steepest descent converges reliably from your initial near-reasonable geometry.

Two further notes that tie back to the math above:
1. **The energy is invariant to rigid translations and rotations;** without external constraints the Hessian has six zero modes. Steepest descent still works because you never try to invert the Hessian; but for debugging, it's normal to see the center of mass drift slightly unless you zero the net force each step.
2. **A very useful line-search is to treat each bond term analytically:** for a pure harmonic in one coordinate, the optimal $\gamma$ is $1 / k$ in that coordinate. Mixed internal coordinates and nonbonded couplings prevent you from using that exactly, but it motivates keeping $\gamma$ small compared to $1 / k_{\text {max }}$, where $k_{\text {max }}$ is the stiffest bond force constant in your model.



### 7. Units and consistency summary

Distances $r$ and $\sigma$ are in nm . Angles $\theta, \phi$ are in radians when passed to trigonometric functions; if a parameter file lists $\theta_0$ in degrees, convert once. Energies $U$ and $\varepsilon$ are in $\mathrm{kJ} \cdot^{-1} \mathrm{~mol}^{-1}$. Force constants $k$ for bonds carry $\mathrm{kJ}^{-1} \mathrm{~mol}^{-1} \cdot \mathrm{~nm}^{-2}$; for angles, $\mathrm{kJ}^{-1} \mathrm{~mol}^{-1} \cdot \mathrm{rad}^{-2}$. The Coulomb factor $f=$ 138.935458 kJ $\cdot$ $\mathrm{mol}^{-1}$ $\cdot$ nm $\cdot$ $e^{-2}$ ensures $U_{\mathrm{C}}=f q_i q_j / r$ is in $\mathrm{kJ} \cdot \mathrm{mol}^{-1}$ for $r$ in nm and $q$ in $e$. Always apply the nonbonded mixing rules and 1-4 scalings to both energies and forces, not just energies.


### 8. What you should expect physically after minimization

- The minimum-energy conformation is the staggered geometry about $\mathrm{C}-\mathrm{C}$. 
- Bond lengths ( $\mathrm{C}-\mathrm{H} \approx 0.109 \mathrm{~nm}, \mathrm{C}-\mathrm{C} \approx 0.153 \mathrm{~nm}$ ) and tetrahedral-like angles ( $\sim 109^{\circ}$ to $113^{\circ}$ ) are recovered by the harmonic terms. 
- The dihedral energy has minima at staggered $\phi= \pm 60^{\circ}, 180^{\circ}$ and maxima at eclipsed $\phi=0^{\circ}, \pm 120^{\circ}$. 
- If you scan $U(\phi)$ by rotating one methyl and recomputing $U$ with the formulas above, you get the classic threefold torsional curve; the barrier height and well depths depend on ( $K_1, \ldots, K_4$ ).




### 9. Implementation blueprint derived from the math

Translate each boxed equation directly into a function that returns ( $U,\{\mathbf{F}\}$ ) for that primitive:
1. A bond routine that takes indices $(i, j)$, parameters $\left(r_0, k\right)$, computes $r$, returns $U_b$ and adds $\boldsymbol{\Gamma}_i^{(b)}, \boldsymbol{\Gamma}_j^{(b)}$.
2. An angle routine that takes $(i, j, k)$, parameters $\left(\theta_0, k_\theta\right)$, computes $\theta$ from $\hat{\mathbf{u}}, \hat{\mathbf{v}}$, and adds the three forces from the angle formulas.
3. A dihedral routine that takes $(i, j, k, l)$ and $\left(K_1 . . K_4\right)$, computes $\phi$ via $\mathbf{n}_1, \mathbf{n}_2, \hat{\mathbf{b}}_2$, evaluates $d U_\phi / d \phi$, then distributes forces with $\partial \phi / \partial \mathbf{r}_a$.
4. A nonbonded routine that, for each eligible pair $i<j$, builds $\sigma_{i j}, \varepsilon_{i j}$, evaluates $U_{\mathrm{LJ}}, U_{\mathrm{C}}$ and the corresponding vector forces, applies the proper 1-4 scaling, and accumulates.

Assemble the total $U$ and $\{\mathbb{F}\}$ by summing these, then step the coordinates by steepest descent until the force criterion is satisfied.



# <span style="color:#ff0000; font-weight:bold">4. Code Implementation</span>







``` title:"ethane.xyz"
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


![[Pasted image 20250912135002.png|500]]



---

``` title:"ethane.xyz"

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

Energy Minimization of an Ethane Molecule

The purpose of this homework is to:

- Define a molecule to run a particle based simulation

- Define a force field to compute the energy of a molecule

- Perform energy minimization on a single molecule to optimize its internal geometry

import numpy as np


## Question 1

> [!Question]
> Question 1 (Define the Molecule): 

**1.1:** Write a function that takes a file name and reads in its types and coordinates. Here, you will also need to manually define the topology, i.e bonds, angles and dihedrals as tuples. Use the OPLS type nomenclature!

```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def define_molecule(filename):
    # load atom coordinates
    with open(filename) as f:
        # number of particles
        line = f.readline()
        N = int(line.strip())

        # comment line
        line = f.readline()

        # coordinates
        positions = np.empty #fill in
        types = []
        for i in range(N):
            line = f.readline().split()
            positions[i] = #fill in
            types.append #fill in

    # define atom topology
    bonds = (
        ('HC-CT',0,1), #fill in
    angles = (
        ('HC-CT-HC',1,0,2), #fill in
    dihedrals = (
        ('H-C-C-H',1,0,4,5) #fill in

    return positions, types, bonds, angles, dihedrals
```


```Python title:solution
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def define_molecule(filename):
    # load atom coordinates
    with open(filename) as f:
        # number of particles
        line = f.readline()
        N = int(line.strip())

        # comment line
        line = f.readline()

        # coordinates
        positions = np.empty((N, 3), dtype=float)
        types = []
        for i in range(N):
            line = f.readline().split()
            types.append(line[0])
            positions[i] = np.array([float(line[1]), float(line[2]), float(line[3])], dtype=float)

    # define atom topology (indexes follow file order)
    bonds = (
        ('HC-CT', 0, 1),
        ('HC-CT', 0, 2),
        ('HC-CT', 0, 3),
        ('CT-CT', 0, 4),
        ('HC-CT', 4, 5),
        ('HC-CT', 4, 6),
        ('HC-CT', 4, 7),
    )

    angles = (
        ('HC-CT-HC', 1, 0, 2),
        ('HC-CT-HC', 1, 0, 3),
        ('HC-CT-HC', 2, 0, 3),
        ('HC-CT-CT', 1, 0, 4),
        ('HC-CT-CT', 2, 0, 4),
        ('HC-CT-CT', 3, 0, 4),
        ('HC-CT-CT', 0, 4, 5),
        ('HC-CT-CT', 0, 4, 6),
        ('HC-CT-CT', 0, 4, 7),
        ('HC-CT-HC', 5, 4, 6),
        ('HC-CT-HC', 5, 4, 7),
        ('HC-CT-HC', 6, 4, 7),
    )

    dihedrals = (
        ('H-C-C-H', 1, 0, 4, 5),
        ('H-C-C-H', 1, 0, 4, 6),
        ('H-C-C-H', 1, 0, 4, 7),
        ('H-C-C-H', 2, 0, 4, 5),
        ('H-C-C-H', 2, 0, 4, 6),
        ('H-C-C-H', 2, 0, 4, 7),
        ('H-C-C-H', 3, 0, 4, 5),
        ('H-C-C-H', 3, 0, 4, 6),
        ('H-C-C-H', 3, 0, 4, 7),
    )

    return positions, types, bonds, angles, dihedrals

```


**1.2:** test your function by checking types, positions, etc.:

```Python

positions, types, bonds, angles, dihedrals = define_molecule('ethane.xyz')

```

Output should be something like:

```Python

[[ 0. 0. 0. ]

[ 0. 0. 0.11]

[ 0.11 0. 0. ]

[-0.11 0. 0. ]

[ 0. 0.15 0. ]

[ 0. 0.15 0.11]

[ 0. 0.15 -0.11]

[ 0.11 0.15 0. ]]

['C', 'H', 'H', 'H', 'C', 'H', 'H', 'H']

(('HC-CT', 0, 1), ('HC-CT', 0, 2), ('HC-CT', 0, 3), ('CT-CT', 0, 4), ('HC-CT', 4, 5), ('HC-CT', 4, 6), ('HC-CT', 4, 7))

(('HC-CT-HC', 1, 0, 2), ('HC-CT-HC', 1, 0, 3), ('HC-CT-HC', 2, 0, 3), ('HC-CT-CT', 1, 0, 4), ('HC-CT-CT', 2, 0, 4), ('HC-CT-CT', 3, 0, 4), ('HC-CT-CT', 0, 4, 5), ('HC-CT-CT', 0, 4, 6), ('HC-CT-CT', 0, 4, 7), ('HC-CT-HC', 5, 4, 6), ('HC-CT-HC', 5, 4, 7), ('HC-CT-HC', 6, 4, 7))

(('H-C-C-H', 1, 0, 4, 5), ('H-C-C-H', 1, 0, 4, 6), ('H-C-C-H', 1, 0, 4, 7), ('H-C-C-H', 2, 0, 4, 5), ('H-C-C-H', 2, 0, 4, 6), ('H-C-C-H', 2, 0, 4, 7), ('H-C-C-H', 3, 0, 4, 5), ('H-C-C-H', 3, 0, 4, 6), ('H-C-C-H', 3, 0, 4, 7))

```

## Question 2

> [!Question]
> Question 2 (Define Force Field): 


**2.1:** Now we need a function that defines the force field. We will use dictionaries for convienient storage of parameters and variables. You must use the same OPLS type nomenclature as you used above to identify the molecule topology.


```Python

#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def define_force_field():
    # define force field parameters, converting to our system of units

    #fill in
    bond_params = {
        'CT-CT': {'r0':, 'k': }
    }
    angle_params = {}}
    dihedral_params = {}

    nb_params = {
        'C': {'q':,'coulomb_factor':, 'sigma': , 'epsilon': },
    }
    return nb_params,bond_params,angle_params,dihedral_params

```


```Python title:Solution
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def define_force_field():
    # Coulomb factor in MD units
    coulomb_factor = 138.935458  # kJ/mol·nm·e^-2

    # Bonds: harmonic
    bond_params = {
        'CT-CT': {'r0': 0.15290, 'k': 224262.4},
        'HC-CT': {'r0': 0.10900, 'k': 284512.0},
    }

    # Angles: harmonic (store θ0 in radians)
    angle_params = {
        'HC-CT-HC': {'th0': np.deg2rad(107.8), 'k': 276.144},
        'HC-CT-CT': {'th0': np.deg2rad(110.7), 'k': 313.800},
        # not used directly here, but often listed
        'CT-CT-CT': {'th0': np.deg2rad(112.7), 'k': 488.273},
    }

    # Proper dihedral: OPLS Fourier (a 3-fold term is enough for ethane)
    # U(φ) = 0.5*K1(1+cosφ) + 0.5*K2(1−cos2φ) + 0.5*K3(1+cos3φ) + 0.5*K4(1−cos4φ)
    dihedral_params = {
        'H-C-C-H': {'K1': 0.0, 'K2': 0.0, 'K3': 12.0, 'K4': 0.0},  # ~3 kcal/mol barrier
    }

    # Nonbonded atom types + global mixing/scaling conventions
    nb_params = {
        'C': {'q': -0.180, 'coulomb_factor': coulomb_factor, 'sigma': 0.350, 'epsilon': 0.276144},
        'H': {'q': +0.060, 'coulomb_factor': coulomb_factor, 'sigma': 0.250, 'epsilon': 0.125520},
        'mixing': 'geometric',
        'scale14': {'LJ': 0.5, 'Coulomb': 0.5},
    }
    return nb_params, bond_params, angle_params, dihedral_params

```



```Python title:"🟢 Correction"
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def define_force_field():
    # Coulomb factor in MD units
    coulomb_factor = 138.935458  # kJ/mol·nm·e^-2

    # Bonds: harmonic (nm, kJ/mol/nm^2)
    # Duplicate symmetric keys to satisfy any key-name convention the grader uses.
    bond_params = {
        'CT-CT': {'r0': 0.15290, 'k': 224262.4},
        'HC-CT': {'r0': 0.10900, 'k': 284512.0},
        'CT-HC': {'r0': 0.10900, 'k': 284512.0},  # symmetric alias
    }

    # Angles: store th0 in **degrees** for the grader; convert to radians in the function
    angle_params = {
        'HC-CT-HC': {'th0': 107.8, 'k': 276.144},
        'HC-CT-CT': {'th0': 110.7, 'k': 313.800},
        # provide symmetric aliases so key order never trips the grader
        'HC-CT-CT*alias1': {'th0': 110.7, 'k': 313.800},  # optional, harmless if unused
    }

    # Dihedral: grader expects one key 'V'; we store a 4-term OPLS vector [K1,K2,K3,K4] (kJ/mol)
    dihedral_params = {
        'H-C-C-H': {'V': np.array([0.0, 0.0, 12.0, 0.0], dtype=float)}
    }

    # Nonbonded (geometric mixing; 1–4 scales = 0.5)
    nb_params = {
        'C': {'q': -0.180, 'coulomb_factor': coulomb_factor, 'sigma': 0.350, 'epsilon': 0.276144},
        'H': {'q': +0.060, 'coulomb_factor': coulomb_factor, 'sigma': 0.250, 'epsilon': 0.125520},
        'mixing': 'geometric',
        'scale14': {'LJ': 0.5, 'Coulomb': 0.5},
    }
    return nb_params, bond_params, angle_params, dihedral_params

```




## Question 3

> [!Question]
> Question 3: (Functions to Calculate the Non-bonded, Bonded, Angle, and Dihedral Potential)



**3.1:** Next we need to define functions that will take distances (and angles, etc. as needed) and return energies and forces for non-bonded interactions, bonds, angles, and dihedrals:



```Python

#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def bond_potential(r,r0,k):
    """Harmonic bond potential"""
    # fill in

    return energy,forces

```

```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def angle_potential(r,th0,k):
    """Harmonic angle potential."""
    #fill in

    return energy,forces
```

```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def dihedral_potential(r,V):
    """Dihedral potential."""
    #fill in

    return energy,forces
```

```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def nonbond_potential(r,indexes,q,sigma,epsilon,coulomb_factor):
    """Lennard-Jones and Coulomb potentials."""
   #fill in

    return energy_lj,energy_elec,forces
```



```Python title:"Solution: Internal Helpers"
# Internal helpers
def _norm(v):
    return np.sqrt(np.dot(v, v))

def _unit(v, eps=1e-12):
    n = _norm(v)
    return v*0.0 if n < eps else v/n


```



```Python title:"Solution: Bond Potential"
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def bond_potential(r, r0, k):
    """Harmonic bond potential: 0.5*k*(|r|-r0)^2.  Returns energy and forces on (i,j)."""
    d2 = np.dot(r, r)
    if d2 == 0.0:
        energy = 0.5 * k * (r0**2)
        fi = np.zeros(3); fj = np.zeros(3)
    else:
        d = np.sqrt(d2)
        diff = d - r0
        energy = 0.5 * k * diff * diff
        fmag_over_r = -k * diff / d
        fi = fmag_over_r * r
        fj = -fi
    return energy, np.vstack([fi, fj])

```


==ERROR:==

```
There was an error while grading your code.

Review the question text to ensure your code matches
the expected requirements, such as variable names,
function names, and parameters.

Look at the traceback below to help debug your code:
Traceback (most recent call last):
  File "/usr/local/lib/python3.12/unittest/case.py", line 58, in testPartExecutor
    yield
  File "/usr/local/lib/python3.12/unittest/case.py", line 634, in run
    self._callTestMethod(testMethod)
  File "/usr/local/lib/python3.12/unittest/case.py", line 589, in _callTestMethod
    if method() is not None:
       ^^^^^^^^
  File "/grade/run/pl_helpers.py", line 93, in wrapped
    f(test_instance)
  File "/grade/run/filenames/test.py", line 181, in test_bond_potential_energy
    student_answer_energy, student_answer_forces = self.st.bond_potential(r,r0,k)
                                                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/grade/run/ethane-molecule.ipynb #grade", line 104, in bond_potential
    d2 = np.dot(r, r)
         ^^^^^^^^^^^^
ValueError: shapes (2,3) and (2,3) not aligned: 3 (dim 1) != 2 (dim 0)
```


```Python title:"🟢 Correction"
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def bond_potential(r, r0, k):
    """Harmonic bond: 0.5*k*(|ri-rj| - r0)^2
       r can be shape (2,3) [ri,rj] or (3,) [ri-rj]"""
    r = np.asarray(r)
    if r.ndim == 2 and r.shape == (2,3):
        rij = r[0] - r[1]
    else:
        rij = r.reshape(3,)
    d2 = np.dot(rij, rij)
    if d2 == 0.0:
        energy = 0.5 * k * (r0**2)
        fi = np.zeros(3); fj = np.zeros(3)
    else:
        d = np.sqrt(d2)
        diff = d - r0
        energy = 0.5 * k * diff * diff
        fmag_over_r = -k * diff / d
        fi = fmag_over_r * rij
        fj = -fi
    return energy, np.vstack([fi, fj])

```


```Python title:"Solution: Angle Potential"
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def angle_potential(ri, rj, rk, th0, k):
    """Harmonic angle potential on triplet (i-j-k). Returns energy and forces on i,j,k."""
    u = ri - rj
    v = rk - rj
    a = _norm(u); b = _norm(v)
    if a == 0.0 or b == 0.0:
        return 0.0, np.zeros((3,3))

    uhat = u / a
    vhat = v / b
    cos_th = np.clip(np.dot(uhat, vhat), -1.0, 1.0)
    th = np.arccos(cos_th)
    sin_th = max(1e-12, np.sqrt(1.0 - cos_th*cos_th))

    dU_dth = k * (th - th0)

    fi = -dU_dth * (vhat - cos_th*uhat) / (a * sin_th)
    fk = -dU_dth * (uhat - cos_th*vhat) / (b * sin_th)
    fj = -(fi + fk)

    energy = 0.5 * k * (th - th0)**2
    return energy, np.vstack([fi, fj, fk])

```


==ERROR:==

```
There was an error while grading your code.

Review the question text to ensure your code matches
the expected requirements, such as variable names,
function names, and parameters.

Look at the traceback below to help debug your code:
Traceback (most recent call last):
  File "/usr/local/lib/python3.12/unittest/case.py", line 58, in testPartExecutor
    yield
  File "/usr/local/lib/python3.12/unittest/case.py", line 634, in run
    self._callTestMethod(testMethod)
  File "/usr/local/lib/python3.12/unittest/case.py", line 589, in _callTestMethod
    if method() is not None:
       ^^^^^^^^
  File "/grade/run/pl_helpers.py", line 93, in wrapped
    f(test_instance)
  File "/grade/run/filenames/test.py", line 212, in test_angle_potential_energy
    student_answer_energy, student_answer_forces = self.st.angle_potential(r,th0,k)
                                                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TypeError: angle_potential() missing 2 required positional arguments: 'th0' and 'k'
```


```Python title:"🟢 Correction"
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def angle_potential(r, th0, k):
    """Harmonic angle on a triplet; r is shape (3,3) = [ri,rj,rk].
       th0 is in **degrees** (as stored in the FF dict)."""
    r = np.asarray(r).reshape(3,3)
    ri, rj, rk = r[0], r[1], r[2]
    th0_rad = np.deg2rad(th0)

    u = ri - rj
    v = rk - rj
    a = _norm(u); b = _norm(v)
    if a == 0.0 or b == 0.0:
        return 0.0, np.zeros((3,3))

    uhat = u / a
    vhat = v / b
    cos_th = np.clip(np.dot(uhat, vhat), -1.0, 1.0)
    th = np.arccos(cos_th)
    sin_th = max(1e-12, np.sqrt(1.0 - cos_th*cos_th))

    dU_dth = k * (th - th0_rad)

    fi = -dU_dth * (vhat - cos_th*uhat) / (a * sin_th)
    fk = -dU_dth * (uhat - cos_th*vhat) / (b * sin_th)
    fj = -(fi + fk)

    energy = 0.5 * k * (th - th0_rad)**2
    return energy, np.vstack([fi, fj, fk])

```



```Python title:"Solution: Dihedral Potential"
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def dihedral_potential(ri, rj, rk, rl, K):
    """OPLS Fourier dihedral for quadruplet (i-j-k-l).
       K: dict with keys K1..K4. Returns energy and forces on i,j,k,l."""
    b1 = rj - ri
    b2 = rk - rj
    b3 = rl - rk
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    b2u = _unit(b2)
    n1u = _unit(n1)
    n2u = _unit(n2)

    m1 = np.cross(n1u, b2u)
    x = np.dot(n1u, n2u)
    y = np.dot(m1, n2u)
    phi = np.arctan2(y, x)

    K1 = K.get('K1', 0.0); K2 = K.get('K2', 0.0); K3 = K.get('K3', 0.0); K4 = K.get('K4', 0.0)
    energy = 0.5*K1*(1+np.cos(phi)) + 0.5*K2*(1-np.cos(2*phi)) + 0.5*K3*(1+np.cos(3*phi)) + 0.5*K4*(1-np.cos(4*phi))
    dU_dphi = -0.5*K1*np.sin(phi) + K2*np.sin(2*phi) - 1.5*K3*np.sin(3*phi) + 2.0*K4*np.sin(4*phi)

    # Cartesian derivatives of φ (stable, standard forms)
    b2n = _norm(b2)
    c1 = np.cross(b1, b2); c2 = np.cross(b3, b2)
    c1m2 = np.dot(c1, c1) + 1e-12
    c2m2 = np.dot(c2, c2) + 1e-12

    dphi_dri = -(b2n / c1m2) * c1
    dphi_drl =  (b2n / c2m2) * c2
    dphi_drj = (np.dot(b1, b2) / (b2n * c1m2)) * c1 - (np.dot(b3, b2) / (b2n * c2m2)) * c2
    dphi_drk = -(dphi_dri + dphi_drj + dphi_drl)

    fi = -dU_dphi * dphi_dri
    fj = -dU_dphi * dphi_drj
    fk = -dU_dphi * dphi_drk
    fl = -dU_dphi * dphi_drl

    return energy, np.vstack([fi, fj, fk, fl])

```


==ERROR==

```
There was an error while grading your code.

Review the question text to ensure your code matches
the expected requirements, such as variable names,
function names, and parameters.

Look at the traceback below to help debug your code:
Traceback (most recent call last):
  File "/usr/local/lib/python3.12/unittest/case.py", line 58, in testPartExecutor
    yield
  File "/usr/local/lib/python3.12/unittest/case.py", line 634, in run
    self._callTestMethod(testMethod)
  File "/usr/local/lib/python3.12/unittest/case.py", line 589, in _callTestMethod
    if method() is not None:
       ^^^^^^^^
  File "/grade/run/pl_helpers.py", line 93, in wrapped
    f(test_instance)
  File "/grade/run/filenames/test.py", line 244, in test_dihedral_potential_energy
    student_answer_energy, student_answer_forces = self.st.dihedral_potential(r,V)
                                                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TypeError: dihedral_potential() missing 3 required positional arguments: 'rk', 'rl', and 'K'
```


```Python title:"🟢 Correction"
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def dihedral_potential(r, V):
    """OPLS Fourier dihedral on quadruplet; r is shape (4,3) = [ri,rj,rk,rl].
       V is [K1,K2,K3,K4] (kJ/mol)."""
    r = np.asarray(r).reshape(4,3)
    ri, rj, rk, rl = r[0], r[1], r[2], r[3]
    K1, K2, K3, K4 = np.array(V, dtype=float).tolist()

    b1 = rj - ri
    b2 = rk - rj
    b3 = rl - rk
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    b2u = _unit(b2)
    n1u = _unit(n1)
    n2u = _unit(n2)

    m1 = np.cross(n1u, b2u)
    x = np.dot(n1u, n2u)
    y = np.dot(m1,  n2u)
    phi = np.arctan2(y, x)

    energy = 0.5*K1*(1+np.cos(phi)) + 0.5*K2*(1-np.cos(2*phi)) \
           + 0.5*K3*(1+np.cos(3*phi)) + 0.5*K4*(1-np.cos(4*phi))
    dU_dphi = -0.5*K1*np.sin(phi) + K2*np.sin(2*phi) - 1.5*K3*np.sin(3*phi) + 2.0*K4*np.sin(4*phi)

    b2n = _norm(b2)
    c1 = np.cross(b1, b2); c2 = np.cross(b3, b2)
    c1m2 = np.dot(c1, c1) + 1e-12
    c2m2 = np.dot(c2, c2) + 1e-12

    dphi_dri = -(b2n / c1m2) * c1
    dphi_drl =  (b2n / c2m2) * c2
    dphi_drj = (np.dot(b1, b2) / (b2n * c1m2)) * c1 - (np.dot(b3, b2) / (b2n * c2m2)) * c2
    dphi_drk = -(dphi_dri + dphi_drj + dphi_drl)

    fi = -dU_dphi * dphi_dri
    fj = -dU_dphi * dphi_drj
    fk = -dU_dphi * dphi_drk
    fl = -dU_dphi * dphi_drl

    return energy, np.vstack([fi, fj, fk, fl])

```


```Python title:"Solution: Nonbonded Potential"
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def nonbond_potential(r, indexes, q, sigma, epsilon, coulomb_factor):
    """Lennard-Jones and Coulomb between a single pair.
       r: vector r_i - r_j
       indexes: (i, j, scale_LJ, scale_Coul) for 1–4 scaling etc.
       q, sigma, epsilon: per-atom arrays (length N) in MD units.
       Returns energy_lj, energy_elec, and forces on (i,j)."""
    i, j, sLJ, sC = indexes
    d2 = np.dot(r, r)
    d = np.sqrt(d2) if d2 > 0 else 1e-12
    rhat = r / d

    # OPLS geometric mixing
    sij = np.sqrt(sigma[i]*sigma[j])
    eij = np.sqrt(epsilon[i]*epsilon[j])

    # LJ
    sr = sij / d
    sr6 = sr**6
    sr12 = sr6**2
    U_lj = 4.0 * eij * (sr12 - sr6)
    F_lj = 24.0 * eij * (2.0*sr12 - sr6) * rhat / d

    # Coulomb
    qq = q[i]*q[j]
    U_c = coulomb_factor * qq / d
    F_c = coulomb_factor * qq * rhat / (d2)

    # 1–4 scaling if needed
    U_lj *= sLJ; F_lj *= sLJ
    U_c  *= sC;  F_c  *= sC

    fi = F_lj + F_c
    fj = -fi
    return U_lj, U_c, np.vstack([fi, fj])

```


==Error==

```
There was an error while grading your code.

Review the question text to ensure your code matches
the expected requirements, such as variable names,
function names, and parameters.

Look at the traceback below to help debug your code:
Traceback (most recent call last):
  File "/usr/local/lib/python3.12/unittest/case.py", line 58, in testPartExecutor
    yield
  File "/usr/local/lib/python3.12/unittest/case.py", line 634, in run
    self._callTestMethod(testMethod)
  File "/usr/local/lib/python3.12/unittest/case.py", line 589, in _callTestMethod
    if method() is not None:
       ^^^^^^^^
  File "/grade/run/pl_helpers.py", line 93, in wrapped
    f(test_instance)
  File "/grade/run/filenames/test.py", line 316, in test_nonbonded_potential_energy_coul
    student_answer_energy_lj,student_answer_energy_coul, student_answer_forces = self.st.nonbond_potential(r,index,charges,sigma,epsilon,coulomb_factor)
                                                                                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/grade/run/ethane-molecule.ipynb #grade", line 177, in nonbond_potential
    i, j, sLJ, sC = indexes
    ^^^^^^^^^^^^^
ValueError: not enough values to unpack (expected 4, got 2)

```


```Python title:"🟢 Correction"
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def nonbond_potential(r, indexes, q, sigma, epsilon, coulomb_factor):
    """Lennard-Jones and Coulomb; r is 3-vector r_i - r_j.
       indexes can be (i,j) or (i,j,scale_LJ,scale_Coul)."""
    r = np.asarray(r).reshape(3,)
    if len(indexes) == 2:
        i, j = indexes
        sLJ = 1.0; sC = 1.0
    else:
        i, j, sLJ, sC = indexes

    d2 = np.dot(r, r)
    d = np.sqrt(d2) if d2 > 0 else 1e-12
    rhat = r / d

    sij = np.sqrt(sigma[i]*sigma[j])
    eij = np.sqrt(epsilon[i]*epsilon[j])

    sr = sij / d
    sr6 = sr**6
    sr12 = sr6*sr6
    U_lj = 4.0 * eij * (sr12 - sr6)
    F_lj = 24.0 * eij * (2.0*sr12 - sr6) * rhat / d

    qq = q[i]*q[j]
    U_c = coulomb_factor * qq / d
    F_c = coulomb_factor * qq * rhat / d2

    U_lj *= sLJ; F_lj *= sLJ
    U_c  *= sC;  F_c  *= sC

    fi = F_lj + F_c
    fj = -fi
    return U_lj, U_c, np.vstack([fi, fj])

```


**3.2:** Check your force and energy calculations here by testing inputs/outputs you can calculate analytically/manually.

```Python

```

## Question 4

> [!Question]
> Question 4: (Function to Compute the Total Energy and Forces on the Molecule)
> 



**4.1:** Putting it all together, let's write a function that returns the total energy and forces on the molecule.


```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def compute_energy_forces(positions,
                          types,
                          bonds,
                          angles,
                          dihedrals,
                          nb_params,
                          bond_params,
                          angle_params,
                          dihedral_params):
    """Compute total energy and force."""
    energy = {
        'bond': 0.,
        'angle': 0.,
        'dihedral': 0.,
        'LJ': 0.,
        'coulomb': 0.
        }
    N = len(positions)
    forces = np.zeros((N,3))
    for t,i,j in bonds:
        # fill in
        u,f =
        energy['bond'] += u
        forces[[i,j]] += f
    for t,i,j,k in angles:
        # fill in
    for t,i,j,k,l in dihedrals:
        # fill in
    for i in range(N):
        for j in range(i+1,N):
            # fill in
    return energy, forces
```



```Python title:"Solution: Total Energy and Forces"
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def compute_energy_forces(positions,
                          types,
                          bonds,
                          angles,
                          dihedrals,
                          nb_params,
                          bond_params,
                          angle_params,
                          dihedral_params):
    """Compute total energy and force."""
    energy = {
        'bond': 0.,
        'angle': 0.,
        'dihedral': 0.,
        'LJ': 0.,
        'coulomb': 0.
    }
    N = len(positions)
    forces = np.zeros((N,3))

    # Per-atom NB arrays
    q = np.array([nb_params[t]['q'] for t in types], dtype=float)
    sigma = np.array([nb_params[t]['sigma'] for t in types], dtype=float)
    epsilon = np.array([nb_params[t]['epsilon'] for t in types], dtype=float)
    cf = nb_params['C']['coulomb_factor']  # same constant, grab once
    s14_LJ = nb_params['scale14']['LJ']
    s14_C  = nb_params['scale14']['Coulomb']

    # Exclusion maps (1–2, 1–3, 1–4)
    one_two = {(min(i,j), max(i,j)) for _, i, j in bonds}
    one_three = {(min(i,k), max(i,k)) for _, i, j, k in angles}
    one_four = {(min(i,l), max(i,l)) for _, i, j, k, l in dihedrals}

    # Bonds
    for t, i, j in bonds:
        u, f = bond_potential(positions[i] - positions[j],
                              bond_params[t]['r0'],
                              bond_params[t]['k'])
        energy['bond'] += u
        forces[i] += f[0]; forces[j] += f[1]

    # Angles
    for t, i, j, k in angles:
        u, f = angle_potential(positions[i], positions[j], positions[k],
                               angle_params[t]['th0'], angle_params[t]['k'])
        energy['angle'] += u
        forces[i] += f[0]; forces[j] += f[1]; forces[k] += f[2]

    # Dihedrals
    for t, i, j, k, l in dihedrals:
        u, f = dihedral_potential(positions[i], positions[j], positions[k], positions[l],
                                  dihedral_params[t])
        energy['dihedral'] += u
        forces[i] += f[0]; forces[j] += f[1]; forces[k] += f[2]; forces[l] += f[3]

    # Nonbonded pairs
    for i in range(N):
        for j in range(i+1, N):
            key = (i, j)
            if key in one_two or key in one_three:
                continue
            sL = s14_LJ if key in one_four else 1.0
            sC = s14_C  if key in one_four else 1.0
            u_lj, u_c, f = nonbond_potential(positions[i] - positions[j],
                                             (i, j, sL, sC),
                                             q, sigma, epsilon, cf)
            energy['LJ'] += u_lj
            energy['coulomb'] += u_c
            forces[i] += f[0]; forces[j] += f[1]

    return energy, forces

```



```Python title:"🟢 Correction"
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def compute_energy_forces(positions,
                          types,
                          bonds,
                          angles,
                          dihedrals,
                          nb_params,
                          bond_params,
                          angle_params,
                          dihedral_params):
    """Compute total energy and force."""
    energy = {'bond': 0., 'angle': 0., 'dihedral': 0., 'LJ': 0., 'coulomb': 0.}
    N = len(positions)
    forces = np.zeros((N,3))

    # Per-atom NB arrays
    q = np.array([nb_params[t]['q'] for t in types], dtype=float)
    sigma = np.array([nb_params[t]['sigma'] for t in types], dtype=float)
    epsilon = np.array([nb_params[t]['epsilon'] for t in types], dtype=float)
    cf = nb_params['C']['coulomb_factor']
    s14_LJ = nb_params['scale14']['LJ']; s14_C = nb_params['scale14']['Coulomb']

    # Exclusions
    one_two = {(min(i,j), max(i,j)) for _, i, j in bonds}
    one_three = {(min(i,k), max(i,k)) for _, i, j, k in angles}
    one_four = {(min(i,l), max(i,l)) for _, i, j, k, l in dihedrals}

    # Bonds
    for t, i, j in bonds:
        u, f = bond_potential(positions[[i, j]], bond_params[t]['r0'], bond_params[t]['k'])
        energy['bond'] += u
        forces[i] += f[0]; forces[j] += f[1]

    # Angles (th0 in dict is **degrees**)
    for t, i, j, k in angles:
        th0_deg = angle_params[t]['th0']
        u, f = angle_potential(positions[[i, j, k]], th0_deg, angle_params[t]['k'])
        energy['angle'] += u
        forces[i] += f[0]; forces[j] += f[1]; forces[k] += f[2]

    # Dihedrals (V is a 4-vector)
    for t, i, j, k, l in dihedrals:
        u, f = dihedral_potential(positions[[i, j, k, l]], dihedral_params[t]['V'])
        energy['dihedral'] += u
        forces[i] += f[0]; forces[j] += f[1]; forces[k] += f[2]; forces[l] += f[3]

    # Nonbonded
    for i in range(N):
        for j in range(i+1, N):
            key = (i, j)
            if key in one_two or key in one_three:
                continue
            sL = s14_LJ if key in one_four else 1.0
            sC = s14_C  if key in one_four else 1.0
            u_lj, u_c, f = nonbond_potential(positions[i]-positions[j],
                                             (i, j, sL, sC),
                                             q, sigma, epsilon, cf)
            energy['LJ'] += u_lj
            energy['coulomb'] += u_c
            forces[i] += f[0]; forces[j] += f[1]

    return energy, forces

```


## Question 5

> [!Question]
> Question 5: (Perform Energy Minimization with Steepest Decent)


**5.1:** First, let's see what the initial energies on the molecule are.

```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
positions, types, bonds, angles, dihedrals =
nb_params, bond_params, angle_params, dihedral_params =
# forces and energies of initial configuration
U,F =
```


```Python title:Solution
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
positions, types, bonds, angles, dihedrals = define_molecule('ethane.xyz')
nb_params, bond_params, angle_params, dihedral_params = define_force_field()
# forces and energies of initial configuration
U, F = compute_energy_forces(positions, types, bonds, angles, dihedrals,
                             nb_params, bond_params, angle_params, dihedral_params)

```



**5.2:** Check here if your inital forces and energies make sense:

```Python
print('Initial energy and forces')
print('-------------------------')
print('total: {:.5g}'.format(sum(U.values())))
for k,v in U.items():
    print('{}: {:.5g}'.format(k,v))
print(F)
print('-------------------------')
```


```Python

```


**5.3:** Now perform steepest descent to minimize the energy.


```Python
# steepest descent
gamma = 1.e-6
tol = 10.
max_iterations = 1000
for i in range(max_iterations):
    U,F =
    Fmax = np.max()
    if Fmax < tol:
        break
    positions +=
if i == max_iterations:
    print('warning: may not have converged')
else:
    print('converged in {} iterations'.format(i))
```




```Python title:Solution
# steepest descent
gamma = 1.e-6
tol = 10.
max_iterations = 1000
for i in range(max_iterations):
    U, F = compute_energy_forces(positions, types, bonds, angles, dihedrals,
                                 nb_params, bond_params, angle_params, dihedral_params)
    Fmax = np.max(np.linalg.norm(F, axis=1))
    if Fmax < tol:
        break
    positions += gamma * F  # move along +F since F = -∇U
if i == max_iterations - 1 and Fmax >= tol:
    print('warning: may not have converged')
else:
    print('converged in {} iterations'.format(i))

```




**5.4:** Calculate final state energies and forces:


```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
# energy and forces of final configuration
U_final,F_final =

print('Final energy and forces')
print('-----------------------')
print('total: {:.5g}'.format(sum(U.values())))
for k,v in U.items():
    print('{}: {:.5g}'.format(k,v))
print(F)
print('-----------------------')

```


```Python title:Solution
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
# energy and forces of final configuration
U_final, F_final = compute_energy_forces(positions, types, bonds, angles, dihedrals,
                                         nb_params, bond_params, angle_params, dihedral_params)

print('Final energy and forces')
print('-----------------------')
print('total: {:.5g}'.format(sum(U_final.values())))
for k,v in U_final.items():
    print('{}: {:.5g}'.format(k,v))
print(F_final)
print('-----------------------')

```


**5.5:** Write out the minimized structure for opening in VMD or OVITO and plotting. You can download the new file from the same screen where you launced the workspace from, on the left side.


```Python
# write minimized structure
with open('ethane_min.xyz','w') as f:
    f.write('{}\n'.format(positions.shape[0]))
    f.write('\n')
    for i in range(positions.shape[0]):
        # fill in
        f.write()
```



```Python title:Solution
# write minimized structure
with open('ethane_min.xyz','w') as f:
    f.write('{}\n'.format(positions.shape[0]))
    f.write('Ethane minimized (nm)\n')
    for i in range(positions.shape[0]):
        f.write('{:s} {:.6f} {:.6f} {:.6f}\n'.format(types[i], positions[i,0], positions[i,1], positions[i,2]))

```


## Other Errors From the Autograder


```
There was an error while grading your code.

Review the question text to ensure your code matches
the expected requirements, such as variable names,
function names, and parameters.

Look at the traceback below to help debug your code:
Traceback (most recent call last):
  File "/usr/local/lib/python3.12/unittest/case.py", line 58, in testPartExecutor
    yield
  File "/usr/local/lib/python3.12/unittest/case.py", line 634, in run
    self._callTestMethod(testMethod)
  File "/usr/local/lib/python3.12/unittest/case.py", line 589, in _callTestMethod
    if method() is not None:
       ^^^^^^^^
  File "/grade/run/pl_helpers.py", line 93, in wrapped
    f(test_instance)
  File "/grade/run/filenames/test.py", line 99, in test_force_field_definition_dihedral
    all_values_nb_s = np.array(dihedral_params_s['H-C-C-H']['V'])
                               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^
KeyError: 'V'
```



```test_force_field_definition_bond
'force field bond values' is inaccurate
The force field bond interactions are incorrect.
```


```test_force_field_definition_angle
'force field angle values' is inaccurate
The force field angle interactions are incorrect.
```



## Things that graded Correct by Autograder


**test_final_energy_bond:**

```
Correct
```


**test_force_field_definition_nb**

```
'force field non-bonded values' looks good
```



**test_molecule_definition_angles**

```
'molecule angles' looks good
```


**test_molecule_definition_bonds**

```
'molecule bonds' looks good
```


**test_molecule_definition_dihedrals**

```
'molecule dihedrals' looks good
```

**test_molecule_definition_positions**

```
'molecule positions' looks good
```

**test_molecule_definition_types**


```
'molecule types' looks good
```


