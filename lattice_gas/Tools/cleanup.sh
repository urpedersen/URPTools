#!/bin/bash


# Clean up stuff and make a gif movie of ppm image files
mkdir ppm
mv conf*.ppm ppm
convert -delay 10 ppm/conf*ppm movie.gif
tar -czf ppm.tar.gz ppm
rm -r ppm

# Extract energy
awk '/ener:/{print $2,$3}' mc_lattice.out > ener.dat 

# Extract number of particles
awk '/ener:/{print $2,$5}' mc_lattice.out > N.dat 


