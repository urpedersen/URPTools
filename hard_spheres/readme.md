Program for simulation of Hard Spheres
======================================

This program is currently in development ...

  - Ulf, May, 2018.

## Input parameters
The below information is available with the `-h` switch:

```

      Simulation of hard-spheres, by Ulf R. Pedersen (http://urp.dk, 2018).

 -h, --help                        Prints this help.
 -q, --quiet                       Hide program output.
 -r, --rcut=NUM        [1.7]       Neighbour list cut-off distance.
 -i, --input=FILE      [none]      Input file (*.xyz, *.xyz.gz or *.atom).
                                     Default [none]: an ideal gas configuration.
 -t, --time_steps=INT  [20]        Monte-Carlo time steps per frame.
 -f, --frames=INT      [250]       Number of frames.
 -s, --step_length=NUM [0.05]      The maximum particle Monte-Carlo step length.
 -L, --Lenghts=NUM     [10]        Size of periodic box,
     --Lenghts=NUM,NUM,NUM           used when not provided in the input file.
 -p, --pressure=NUM    [1.0]       Pressure for Monte-Carlo barostat (if applied),
                                     by default the program runs a NVT simulation.
 -v, --volume_step     [0.01]      Max step size in log(V) for Monte-Carlo barostat.
 -o, --output=FILE     [traj.xyz]  Output file for trajectory (*.xyz or *.xyz.gz).

  Documentation and source code at: http://urp.dk/tools 


```

## Error handling

If you get `error: Overlap (befor NL update)` then the input configuration have may have overlapping particles.

If you get `error: Overlap (after NL update)` during a run, then increase the neighbor-list cut-off (`-r`).

