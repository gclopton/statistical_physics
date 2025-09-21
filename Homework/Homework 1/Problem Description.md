

I am solving a problem involving the energy Minimization of an ethane molecule. Here is a description of the problem:

I am given an .xzy file. 


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


I need to do the following:

## 1.) Using OPLS nomenclature, define the topology of the molecule, including the bonds, angles, and dihedrals as tuples.

```Python
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


## 2.) I need to define the force field using dictionaries to store the following information:

```Python
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





## 3.) Write functions to calculate the non-bonded, bonded, angle, and dihedral potentials:**

```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def bond_potential(r,r0,k):
    """Harmonic bond potential"""
    # fill in

    return energy,forces
    
    
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def angle_potential(r,th0,k):
    """Harmonic angle potential."""
    #fill in

    return energy,forces
    
    
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def dihedral_potential(r,V):
    """Dihedral potential."""
    #fill in

    return energy,forces
    
    
def nonbond_potential(r,indexes,q,sigma,epsilon,coulomb_factor):
    """Lennard-Jones and Coulomb potentials."""
   #fill in

    return energy_lj,energy_elec,forces
```


## 4.) We need to then write a function to compute the total energy and forces on the molecule:

```Python
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


## 5.) Finally we need to perform energy minimization with steepest descent. 

We will calculate the initial energies of the molecule:


```Python
positions, types, bonds, angles, dihedrals =
nb_params, bond_params, angle_params, dihedral_params =
# forces and energies of initial configuration
U,F =
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


