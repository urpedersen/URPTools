# rotational_order

## Q6
Compute the Q6 rotational bond order-parameters by Steinhard (ql),
and the Lechner-Dellago averaged versions (qlAvl).
Bond are written to the last columns of an xyz-file.

This program is written by Ulf R. Pedersen (2016), www.urp.dk.

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
./bin/Q6 --input=in.xyz --Lenghts=15.0,15.0.30.0
./bin/Q6 -iin.xyz -L15.0,15.0.30.0

Optional flags:
 -h, --help                    Print information on using this program.
 -q, --quiet                   Hide program output.
 -i, --input=FILE              *.xyz or *.xyz.gz input file.
 -L, --Lenghts=NUM             Size of periodic box
     --Lenghts=NUM,NUM,NUM 
 -o, --output=FILE             NOT IMPLIMENTET
                               Current version write q6.xyz and 14.xyz.
```
## References
Accurate determination of crystal structures based on averaged local bond order parameters
Wolfgang Lechner and Christoph Dellago
J. Chem. Phys. 129, 114707 (2008)
http://dx.doi.org/10.1063/1.2977970
