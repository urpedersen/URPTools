# lattice_gas
A program to perform a 2D computations of a lattice gas by [Ulf R. Pedersen](http://urp.dk) (2012).

## Build 
This program uses the [boost c++ library](http://www.boost.org). Get it on Debian or Ubuntu with:
```
    sudo apt-get install libboost-all-dev
```

Type
``` 
    make 
```
to generate an executable into the bin library.

## Usage example
```
make
cd Examples
../bin/lattice_gas | tee lattice_gas.log 
sh ../Tools/make_movie.sh
sh ../Tools/cleanup.sh
```
The file input.cfg contains the setting for this computations. Use the -h flag for additional help.
