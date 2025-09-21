

# 4.1 Particle Coordinates


You want to confine a particle to a box of side length $L$ centered about the origin (i.e. the simulation box is defined from $-L / 2 \leq x, y, z<L / 2$.

If the particle's literal position is outside the box, each coordinate should be wrapped at its respective boundary so that the position is mapped back into the box. The way to do this for a given coordinate is to subtract/add some integer multiple of $L$ tha puts the coordinate between $-L / 2$ and $L / 2$.

Write a Python function my_pos_in_box that takes a list of atomic positions in open boundary conditions and return them unde PBC (i.e. wrapped back into the simulation box).

It is highly recommended to use modulo operations instead of while loops in this function. Given two numbers $x, y$ (int or float, $\mathrm{y}>0$ ), it is possible to find $n$ and $r$ such that

$$
x=n * y+r, \quad n \in \mathbb{Z}, \quad r \in[0, y) .
$$


In python, one can get $n$ and $r$ using the integer divide " $/ /$ " and modulo " $\%$ " operators, respectively.
hint: think of the input (pos) as a float rather than an array. numpy will automatically apply your operation to all elements of the array.


```Python
import numpy as np
def my_pos_in_box(pos, lbox):
  """ wrap positions inside simulation box

  Args:
    pos (np.array): positions, shape (natom, ndim)
    lbox (float): box side length
  Returns:
    np.array: pos in box
  """
  pass
```




# 4.2 Discpacement Vector


In this problem, you will implement the minimum image convention (MIC). Consider two particles at $R_1=[-0.4,0.0,0.0]$ and $R_2=[0.4,0.0,0.0]$ in a box containing $[-0.5,0.5)$ (i.e. $L=1.0$ ) in all 3 spatial dimensions. With open boundary condition, $R_1-R_2=[-0.8,0.0,0.0]$. However, with the MIC, $R_1-R_2=[0.2,0.0,0.0]$.

In periodic boundary conditions, there are copies of $R_2$ at $R_2+L *\left[n_1, n_2, n_3\right]$ for all integers $n_1, n_2, n_3$. In the MIC, when asked for pair differences $R_1-R_2$, we only consider the difference between $R_1$ and the nearest image of $R_2$. For a cubic box, this means that pair distances will never be larger than $\sqrt{3} L / 2$.

Write a Python function my_disp_in_box that imposes the minimum image convention on a displacement vector drij = R_i-R_j in open boundary condition.

Your function should be consistent about the boundary case such that $f\left(R_i-R_j\right)=-f\left(R_j-R_i\right)$. "np.round" may be of help. Depending on how you wrote the previous exercise, the same exact routine might work here!
hint: think of the input (drij) as a float rather than an array. numpy will automatically apply your operation to all elements of the array.

IMPORTANT: Your function will receive full credit even if it does NOT pass all tests. Please read through the grading output. Pay special attention to test cases 5 and 6 . If they are not handled consistently, then you must manually symmetrize the displacement table in your MD code.


```Python
import numpy as np
def my_disp_in_box(drij, lbox):
    """ 
    Impose minimum image condition on displacement vector drij=ri-rj

    Args:
      drij (np.array): length-3 displacement vector ri-rj
      lbox (float): length of cubic cell
    Returns:
      np.array: drij under MIC
    """
    pass
```




# 4.3 Kinetic Energy


Given a list of velocities of particles, compute the total kinetic energy

$$
K=\sum_i \frac{1}{2} m v_i^2
$$

where the mass $m$ is the same for every particle, and passed in as a parameter.
For the test input, vel is a (nparticle, ndim) numpy array. For example, vel [0], the velocity vector of particle 0 , is a length- 3 numpy array. In turn, $\mathrm{V}[0][0], \mathrm{V}[0][1], \mathrm{V}[0][2]$ are the $\mathrm{x}, \mathrm{y}, \mathrm{z}$ components of that particle's velocity.


```Python
import numpy as np
def my_kinetic_energy(vel, mass):
  """ Calculate total kinetic energy.

  Args:
    vel (np.array): particle velocities, shape (nparticle, ndim)
    mass (float): particle mass
  Return:
    float: total kinetic energy
  """
  return 0.0
```



# 4.4 Potential Energy


Given a distance table, compute the total Lennard-Jones potential energy

$$
P E_{\text {total }}=\sum_{i<j} V\left(r_{i j}\right),
$$

where $r_{i j}$ is the distance particles ( $\mathrm{i}, \mathrm{j}$ ). $V(r)$ is the Lennard-Jones pair potential

$$
V(r)=4 \epsilon\left(\frac{\sigma}{r}\right)^6\left[\left(\frac{\sigma}{r}\right)^6-1\right]
$$


Use Lennard-Jones units $\epsilon=\sigma=1$.
Do NOT apply any cutoff or shift to the potential.
For the test input, rij is a distance table. That is, $\mathrm{rij}[0,1]$ gives the distance between particles 0 and 1 in the minimum image convention.



```Python
import numpy as np
def my_potential_energy(rij):
  """ Calculate total potential energy.

  Args:
    rij (np.array): distance table, shape (nparticles, nparticles)
  Return:
    float: total potential energy
  """
  return 0.0
```



# 4.5 Potential Energy Cut and Shifted


Given a distance table and cutoff, compute the total Lennard-Jones potential energy using a cut-and-shifted potential

$$
v_{c s}(r)= \begin{cases}v(r)-v\left(r_c\right) & r<r_c \\ 0 & r \geq r_c\end{cases}
$$


Use Lennard-Jones units $\epsilon=\sigma=1$.
For the test input, rij is a distance table. That is, $\mathrm{rij}[0,1]$ gives the distance between particles 0 and 1 in the minimum image convention.



```Python
import numpy as np
def my_potential_energy(rij, rc):
  """ Calculate total potential energy using
   cut-and-shifted Lennard-Jones potential.

  Args:
    rij (np.array): distance table, shape (nparticle, nparticle)
    rc (float): cutoff radius
  Return:
    float: total potential energy
  """
  return 0.0
```




# 4.6 Lennard-Jones Force

Given two particles at $\vec{r}_i$ and $\vec{r}_j$, define $\vec{r}_{i j} \equiv \vec{r}_i-\vec{r}_j$, and $r_{i j} \equiv\left|\vec{r}_{i j}\right|$. The Lennard-Jones force on particle $i$ is

$$
\vec{F}_i=-\sum_{j \neq i} \vec{\nabla}_i V\left(r_{i j}\right)=24 \epsilon \sum_{j \neq i} \frac{\sigma}{r_{i j}^2}\left(\frac{\sigma}{r_{i j}}\right)^6\left[2\left(\frac{\sigma}{r_{i j}}\right)^6-1\right] \vec{r}_{i j} .
$$


Recall that the Lennard-Jones potential is

$$
v(r)=4 \epsilon\left(\frac{\sigma}{r}\right)^6\left[\left(\frac{\sigma}{r}\right)^6-1\right] .
$$


Use Lennard-Jones units $\epsilon=\sigma=1$.

Write a function my_force_on that takes an index i, a list of positions pos, the box side length lbox and returns the total Lennard-Jones force on particle i.

For the test input, pos is a $n-$ by -3 numpy array, where $n$ is the number of particles in the list. For example, pos [0], the position vector of particle 0 , is a length-3 numpy array. So pos [0] [0], pos [0] [1], pos [0] [2] are the $x, y, z$ coordinates of that particle. pos [i] is the position vector of particle $i$, on which we are going to calculate the total Lennard-Jones force exerted by all the other particles.

Remember to apply the minimum image convention to displacement vectors.

<span style="color:#ff0000; font-weight:bold">Do NOT cut off the force in this problem</span>


```Python
import numpy as np
def my_force_on(i, pos, lbox):
    """
    Compute force on particle i

    Args:
      i (int): particle index
      pos (np.array) : particle positions, shape (nparticles, ndim)
      lbox (float): side length of cubic box
    Returns:
      np.array: force on particle i, a length-3 vector
    """
    force = np.zeros(3)
    return force
```


# 4.7 Lennard-Jones Force Cut


Add the cutoff radius `rc` to `my_force_on` implemented in the previous question. Do not apply any smoothing scheme.
Remember to apply the minimum image convention to displacement vectors.



```Python
import numpy as np
def my_force_on(i, pos, lbox, rc):
    """
    Compute force on atom i, cut off pair interaction beyond rc

    Args:
      i (int): particle index
      pos (np.array) : particle positions, shape (natom, ndim)
      lbox (float): side length of cubic box
      rc (float): cutoff radius
    Returns:
      np.array: force on atom i, a length-3 vector
    """
    force = np.zeros(3)
    return force
```



