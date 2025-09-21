
# 3.1 Integration of Harmonic Oscillator


![[Screenshot 2025-09-21 at 12.42.35 AM.png]]



Here, the goal is to perform a NVE integration of the position $x$ and momentum $p$ of an harmonic oscillator and determine its phase space location ( $x, p$ ). Then we will compare the nummerical solution to the analytical one.


## Set-Up Initial Conditions


```Python 
import numpy as np
import matplotlib.pyplot as plt
```



```Python 
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
# initial condition
x =
p =
m =
```



```Python 
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
k =
dt =
```


## Force on Harmonic Oscillator


```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def harmonic_force(x,k):
    # fill in
    return force
```



## Calculate Number of Steps to Run


```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
#natural period of oscillator
T =

# Number of steps for 10 periods
num_steps =

```


## Define Velocity Verlet Integration


```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def velocity_verlet(dt,k,x,p):
    # fill in
    return x,p

```



## Run Velocity Verlet to Create Trajectory


```Python
outfile = 'nve_%s.dat'%dt
traj = open(outfile,'w')

for i in range(num_steps):
    traj.write() # fill in
    x,p = # fill in

traj.close()
```



## Define Analytical Solution


```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
# analytical solution
t =
x_exact =
p_exact =
```



# 3.2 Integration of Harmonic Oscillator Langevin Dynamics

Here, the goal is to perform a NVT integration of the position $x$ and momentum $p$ of an harmonic oscillator and determine its phase space location ( $x, p$ ). We will use Langevin dynamics.

> [!Question]
> ![[Screenshot 2025-09-21 at 12.47.30 AM.png]] 




## Harmonic Oscilator


```Python
import numpy as np
import matplotlib.pyplot as plt
```

## Initial Condition

```Python
# initial condition
x =
p =
m =
k =
```


## Force on Harmonic Oscillator

```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def harmonic_force(x,k):
    # fill in
```


## Integration Parameters

```Python
# integration parameters
dt =
total_time =
sample_time =
num_steps =
num_sample_steps =
gamma =
kT =
```



## Langevin Integrator


```Python
rng = np.random.default_rng(seed=42)
```


```Python
#grade (enter your code in this cell - DO NOT DELETE THIS LINE)
def langevin(dt,k,x,p,gamma,kT,m,rng):
    # fill in
    return x,p
```


```Python
# open file for writing
traj = open('ld.dat','w')

for i in range(num_steps):
    if i % num_sample_steps == 0:
        traj.write()
    x,p = # fill in

traj.close()
```



