# rotational_order

## Q6
Compute the Q6 rotational bond order-parameters by Steinhard (ql),
and the Lechner-Dellago averaged versions (qlAvl).
Local bond-orders parameters are printed to the last columns of an xyz-file.

Written by Ulf R. Pedersen (2016), www.urp.dk.

### Build
This program uses the boost c++ library. Get it on Debian or Ubuntu with:
```
sudo apt-get install libboost-all-dev
```
Then type
```
make
```
to build. The executable program is
```
./bin/Q6
```
Type
```
./bin/Q6 -h
```
for more infomation on using the program.
### Usage:
```
Usage examples:
./bin/Q6 --input=input.xyz --Lenghts=15.0,15.0.30.0
./bin/Q6 -i bcc.xyz -L 15.0 -r 1.5 -l 4 -o Q4.xyz

Optional flags:     [default]
 -h, --help                     Print information on using this program.
 -q, --quiet                    Hide program output.
 -l, --degree=INT   [6]         The degree of the bond rotational-order.
 -r, --rcut=NUM     [1.4]       Neighbour cutoff distance.
 -i, --input=FILE   [traj.xyz]  Input file (*.xyz, *.xyz.gz or *.atom).
 -f, --frame=INT    [0]         Frame of input file (0=first frame).
 -s, --scan                     Scan frames in file (only std output).
 -L, --Lenghts=NUM  [10]        Size of periodic box,
     --Lenghts=NUM,NUM,NUM        unless it is provided in the input file.
 -Q  --QminQmax=NUM,NUM         Minimum and maximum q6 limits
                    [0.0,1.0]     Default values includes all particles.
 -S, --Sij=NUM      [-1.0]      Min threshold values for the Sij connection matrix.
                                  The default value of -1.0 skips computation.
                                  A sparse matrix is written to node_connections.dat and
                                  the largest cluster is written to largest_cluster.xyz.
                                  The Lechner-Dellago vectors are used to compute Sij.
 -o, --output=FILE  [none]      Output file (*.xyz or *.xyz.gz).
                                  Default value does not produce an output file.
                                  The 5th column gives number of neighbours. 
                                  The 6th column gives Steinhard bond-order. 
                                  The 7th column gives Lechner-Dellago bond-order. 
                                  The 8th column gives the number of Sij connections, 
                                  if the -S flag is used. 
```
### References
* P. Steinhardt, D. Nelson, and M. Ronchetti, Phys. Rev. B 28, 784 (1983), https://doi.org/10.1103/PhysRevB.28.784
* W. Lechner and C. Dellago, J. Chem. Phys. 129, 114707 (2008), http://dx.doi.org/10.1063/1.2977970

### Table of local orderparameters
For reference, the table below list local bond-order parameters using a cut-off including N neighbours.
The last two coulmns contains the average qi values of the Lennard-Jones liquid at the coexistence state-point at T = 2 and rho = 1.062 (see J. Chem. Phys. 139, 104102 (2013); http://dx.doi.org/10.1063/1.4818747 for details on this state-point; the cut-off is 1.38). (ste) referer to the Steinhard definition, while (lec) is the Lechner-Dellago averaged version. 

|         | bcc   | bcc   | fcc   | hcp   | Icosahedral | sc    | liq(ste) | liq(lec) |
| ------- |:-----:|:-----:|:-----:|:-----:|:-----------:|:-----:|:--------:|---------:|
| N       | 8     | 14    | 12    | 12    | 12          | 6     | 12       | 12       |
| l=2     | 0     | 0     | 0     | 0     | 0           | 0     | 0.094    | 0.029    |
| l=3     | 0     | 0     | 0     | 0.076 | 0           | 0     | 0.110    | 0.027    |
| l=4     | 0.509 | 0.036 | 0.190 | 0.097 | 0           | 0.764 | 0.154    | 0.051    |
| l=5     | 0     | 0     | 0     | 0.252 | 0           | 0     | 0.254    | 0.064    |
| l=6     | 0.629 | 0.511 | 0.575 | 0.484 | 0.663       | 0.354 | 0.376    | 0.155    |
| l=7     | 0     | 0     | 0     | 0.311 | 0           | 0     | 0.361    | 0.087    |
| l=8     | 0.213 | 0.429 | 0.404 | 0.317 | 0           | 0.718 | 0.289    | 0.111    |
| l=9     | 0     | 0     | 0     | 0.138 | 0           | 0     | 0.247    | 0.060    |
| l=10    | 0.650 | 0.195 | 0.013 | 0.010 | 0.363       | 0.411 | 0.241    | 0.074    |
| l=11    | 0     | 0     | 0     | 0.123 | 0           | 0     | 0.264    | 0.062    |
| l=12    | 0.415 | 0.405 | 0.600 | 0.565 | 0.585       | 0.696 | 0.295    | 0.103    |

A part this table is from J. Chem. Phys. 138, 044501 (2013); http://dx.doi.org/10.1063/1.4774084

