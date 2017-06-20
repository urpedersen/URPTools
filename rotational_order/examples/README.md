This folder contains some files with example configuration that will help the user to get familiar with the program.
Below we give some usage examples on these configurations. For a summary on using the program in general try
```
../bin/Q6 -h
```

# Global bond order
The following command will run the program on a configuration with a fcc crystal
```
../bin/Q6 -i fcc.xyz
```
and gived the output
```
Total number of particles:    500
Number of types:              1
Number of particles of types: 500
Lengths of box vectors:       8 8 8
Box volume:                   512
Number density:               0.976562
degree:                    l= 6
rcut:                     rc= 1.4
Average number of neighbors:  12
Average qi:               Ql= 0.574524
Average qiAvg:         QlAvg= 0.574524
Q minimum value:        Qmin= 0
Q maximum value:        Qmax= 1
Number of selected par.:      500

```
On line 10 we read that the Steainhard Q6 orderparameter is 0.574, and on line 11 that the Lechner-Dellago average version is  also 0.575. 
## Size of periodic cell
The size of the periodic simulation cell is read from the header of the file (try: `head fcc.xyz`) in the format `sim_box=RectangularSimulationBox,8,8,8'` (using the convention of the RUMD program). The size can also be specified from the commandline with the `-L` flag:
```
../bin/Q6 -i fcc.xyz -L8,8,8
```
## Changing cut-off
Here we used the default cut-off of 1.4, this can be change with the `-r` flag like this:
```
../bin/Q6 -i fcc.xyz -r1.3
```

## LAMMPS files
The program can also read output in the LAMMPS format: 
```
../bin/Q6 -i liquid.atom
```
## Computing Q4 (or other degrees)
The degree of the speherical harmonics used for the bond-order can be changed with the `-l` flag like this
```
../bin/Q6 -i fcc.xyz -l4
```
## Scan a trajectory
The program can analyse concatenated configurations if the `-s` flag is used:
```
../bin/Q6 -s liquid.atom
```
# Local bond order
The single-particle bond-order parameters can be written to a xyz-configuration like this
```
../bin/Q6 -i liq_fcc.xyz -o Q6.xyz
```
The first part of the produced Q6.xyz file look like this (try: `head Q6.xyz`):
```
8000
 numTypes=2 sim_box=RectangularSimulationBox,15.2274,15.2274,31.977 columns=type,x,y,z,Ni,ql,qlAvg,Ns Navg=12.4682 Ql=0.41153 QlAvg=0.295578
1 0.0846575 15.1051 0.205621 12 0.497299 0.379776
1 0.914088 0.556652 31.9627 13 0.331681 0.394229
1 0.888056 15.1789 0.785834 12 0.522974 0.37958
1 0.187222 0.867201 0.832071 14 0.272439 0.408504
1 1.70258 15.085 31.8381 14 0.328237 0.363577
1 2.39801 0.720773 31.9501 13 0.480408 0.371095
1 2.29595 0.0184157 0.824872 13 0.419481 0.370507
1 1.54915 0.65325 0.771499 13 0.475937 0.383607
```
This is an extended xyz-file where each line refere to a particle. The columns are the particle `id`, the `x`, `y` and `z` coordinate, the number of neighbours within the cut-off `Ni`, the local Steinhard order-parameter `ql`, and the average Lechner-Dellago version `qlAvg`.

We can plot the local bond order parameter along the z-coordinate with xmgrace:
```
tail -n 8000 Q6.xyz > Q6.dat
xmgrace -block Q6.dat -bxy 4:6 -bxy 4:7
```
The plot suggest that we have a liqud and a crystaline fcc slap (henche the name liq_fcc.xyz), and that crystaline particles are the ones with a Q6 larger that 0.30. We will use this information in the next section.
## Selecting particles
We can select particles that are within a range with the following commands:
```
../bin/Q6 -i liq_fcc.xyz -Q0.30,1.00 -o slab_cry.xyz
../bin/Q6 -i liq_fcc.xyz -Q0.00,0.25 -o slab_liq.xyz
../bin/Q6 -i liq_fcc.xyz -Q0.25,0.30 -o slab_between.xyz
```

# Cluster analysis
The dotprocuct Sij=qi.qj between the complex qi vectors of neighboring particles can be used to assing connections. The program can make a cluster analysis using the `-S` flag. First, we include all particles by using the value 0 for the `-S` flag (and use gnuplot)
```
../bin/Q6 -i liq_fcc.xyz -S0 -o Q6.xyz
awk '{if(NR>2){print $9,$6}}' node_connections.dat > S.dat
echo plot \'S.dat\'|gnuplot -p
```
The output file named `node_connections.dat` is a space matrix with particle connections. 
The first part of the file looks like this (try: `head node_connections.dat`)
```
8000
List of connected particles using Sij > 0
0 1993 # S = 0.166418  r = 1.03138
0 6673 # S = 0.14835  r = 1.15086
0 1423 # S = 0.173737  r = 1.15723
0 6008 # S = 0.154666  r = 1.19377
0 7427 # S = 0.162616  r = 1.00511
0 6007 # S = 0.132918  r = 1.18183
0 667 # S = 0.142214  r = 1.20957
0 2 # S = 0.146211  r = 0.993752
```
This file can be anaysed further using the program `cluster_analysis`. The plot where we plot S agains r suggest that we should use rc=1.3 for the distance cut-off and Sc=0.1 for the Sij cut-off.
```
../bin/Q6 -i liq_fcc.xyz -r1.3 -S0.1 -o Q6.xyz
```
The largest cluster of connected particles are witten to `largest_cluster.xyz`.
