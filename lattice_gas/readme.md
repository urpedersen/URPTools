# lattice_gas
A program to perform a 2D computations of a lattice gas by [Ulf R. Pedersen](http://urp.dk) (2012).

## Build 
This program uses the [boost c++ library](http://www.boost.org). Get it on Debian with:
```
    apt-get install libboost-all-dev
```

Type
``` 
    make 
```
to generate an executable into the bin library.

## Usage
```
make
cd Examples
nano input.cfg
../bin/lattice_gas | tee lattice_gas.log 
../Tools/make_movie.sh
../Tools/cleanup.sh
```
Use the -h flag for additional help.
