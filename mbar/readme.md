# mbar
Implementation of the MBAR algorithm for reweighing of probability distributions.
Please read [Shirts and Chodera, J. Chem. Phys. 129, 124105 (2008); doi:10.1063/1.2978177] for details on the MBAR algorithm.

## Build
```
make
```
(this program has no dependencies)

## Usage
The program reads files from a set of umbrella simulations with harmonic bias potentials U_umbrella=0.5*kappa*(Q-a). The meta information for the umbrella windows is in a file named mbar_window_info.txt that should look something like this:
```
NOumbrella.dat 2.5 0.0 0.0 
umbrella01.dat 2.5 1.0 0.2
umbrella02.dat 2.5 2.0 0.2
```
The 1st column is the filename of the simulation data, the 2nd column is the inverse temperature, the 3rd coulmn is the anchor point of the umbrella, and the 4th coulmn is the spring constant. The files with simulation data should look lile this
```
4 2
6 3
7 5
9 6
3 7
3 8
```
The first column is the Q value of the parameter that the umbrella has been applies to, 
and the 2nd column is value of the observable that need tobe reveighed (note that these can be the same).

The log of the reweighted distribution function is in log_probability.dat.

Usage example:
```
make
cd Examples
../bin/mbar --verbose=10 | mbar.log
xmgrace log_probability.dat.
```
The output has more information on usage of the program (set --verbose=5 or --verbose=0 to hide some or all of this information).
