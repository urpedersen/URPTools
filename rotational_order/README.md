# rotational_order

## Q6
    Compute the Q6 rotational bond order-parameters by Steinhard (ql),
            and the Lechner-Dellago averaged versions (qlAvl).
     Local bond orders are printed to the last columns of an xyz-file.

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
to build. The executable program named Q6 is located in
```
./bin/
```

### Usage examples:
```
./bin/Q6 --input=input.xyz --Lenghts=15.0,15.0.30.0
./bin/Q6 -i bcc.xyz -L 15.0 -r 1.5 -l 4 -o Q4.xyz

Optional flags:     [default]
 -h, --help                     Print information on using this program.
 -q, --quiet                    Hide program output.
 -l, --l=INT        [6]         The degree of the bond rotational-order.
 -r, --rcut=NUM     [1.4]       Neighbour cutoff distance.
 -i, --input=FILE   [input.xyz] Input file (*.xyz or *.xyz.gz).
 -L, --Lenghts=NUM  [10]        Size of periodic box.
     --Lenghts=NUM,NUM,NUM 
 -Q  --QminQmax=NUM,NUM         Minimum and maximum Q6 limits
                    [0.0,1.0]   Default values includes all particles.
 -S, --Sij=NUM      [-1.0]      Min threshold values for the Sij connection matrix.
                                Default value -1.0 to skip computation. A sparse matrix is written to node_connections.dat
 -o, --output=FILE  [Q6.xyz]    Output file (*.xyz or *.xyz.gz).
```
## References
Accurate determination of crystal structures based on averaged local bond order parameters
Wolfgang Lechner and Christoph Dellago
J. Chem. Phys. 129, 114707 (2008)
http://dx.doi.org/10.1063/1.2977970
