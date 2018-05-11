Program for simulation of Hard Spheres
======================================

This program is currently in development ...

  - Ulf, May, 2018.

## Input parameters
The below information is available with the `-h` switch:

```
      Simulation of hard-spheres, by Ulf R. Pedersen (http://urp.dk).

 -h, --help                       Prints this help.
 -q, --quiet                      Hide program output.
 -r, --rcut=NUM        [1.7]      Neighbour cut-off distance.
 -i, --input=FILE      [none]     Input file (*.xyz, *.xyz.gz or *.atom).
                                    Default [none]: an ideal gas configuration.
 -t, --time_steps=INT  [20]       Time steps per frame.
 -f, --frames=INT      [250]      Number of frames.
 .s, --step_length=NUM [0.05]      The maximum particle step length.
 -L, --Lenghts=NUM     [10]          Size of periodic box,
     --Lenghts=NUM,NUM,NUM           used when not provided in the input file.
 -p, --pressure=NUM    [1.0]       Pressure for barostat (if applied).
 -v, --volume_step=NUM [0.0]       Volume step for barostat (if applied).
                                     The default is 0.0 resulting in a NVT simulation
 -o, --output=FILE     [traj.xyz]  Output file for trajectory (*.xyz or *.xyz.gz).

  Documentation and source code on github: https://github.com/urpedersen/URPTools 
```

## Error handling

If you get `error: Overlap (befor NL update)` then the input configuration have may have overlapping particles.

If you get `error: Overlap (after NL update)` during a run, then increase the neighbor-list cut-off (`-r`).

