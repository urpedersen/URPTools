# URPTools

This is a bunch of stand-alone (linux) tools that I during the years have made for numerical computations in the field of condensed matter physics/chemistry. Feel free to use the programs, write me about bugs, suggestions or questions (e-mail: ulf AT urp.dk). I would also be happy to hear if you make use of these programs.

[Ulf R. Pedersen](http://urp.dk)

## List of programs

| Program Name      | Brief Description                                                            |
|-------------------|------------------------------------------------------------------------------|
| cluster_analysis  | Perform a cluster analysis of connected nodes                                |
| dat2wav           | Convert an input of float values to an WAV audio file                        |
| hard_sphere       | Simulation of hard spheres                                                   |
| lattice_gas       | Monte Carlo simulations of a 2D lattice gas                                  |
| make_lattice      | Setup crystal lattices like FCC, BCC or NaCl                                 |
| mbar              | Implimentation of the MBAR algorithm for reweighting an umbrella sampling    |
| rotational_order  | Computes the Steinhard Q6 order-parameter, and the Lechner-Dellago version   |
| traj              | Programs for analysing molecular dynamics trajectories                       |


## Download and build (linux)
1. Clone and navigate to the git repository
```sh
git clone https://github.com/urpedersen/URPTools.git
cd URPTools
```
Later, the newest branch can be downloaded with
```
git pull
```
2. Move into a program directory and read about how to build and usage
```sh
cd PROGRAM
cat readme.md
```

