# Examples of using the hard-sphere program

Perform a equnriliation simulation
```
  ../../make_lattice/bin/make_lattice -lfcc -c6,6,6 -r 0.5
  ../bin/hard_spheres -i start.xyz.gz | tee hard_spheres.log
```

Clean-up, and run NVT simulation
```
  mv final.xyz first.xyz
  rm traj.xyz
  ../bin/hard_spheres -i first.xyz | tee hard_spheres.log
```

Compute the radial distribution function
```
  ../../traj/bin/traj_gr --load_traj_type=3 --bboxX=12
```

Clean-up, and run NpT simulation
```
  mv final.xyz first.xyz
  rm traj.xyz
  ../bin/hard_spheres -i first.xyz -p 1.0 -v 0.01 -t 50 | tee hard_spheres.log
```

Plot time-trajectory of packing fraction
```
  awk '/ener:/&&!/frame/{print $2,$4}' hard_spheres.log > thermo.dat
  xmgrace thermo.dat
```

Compute the radial distribution function.
```
  ../../traj/bin/traj_gr --load_traj_type=3
```


Perform equlibrium simulation
```
  ../bin/hard_spheres -i first.xyz
```

Compute the radial distribution function
```
  ../../traj/bin/traj_gr --load_traj_type=3 -o Q6.xyz
```

Conpute rotational order
```
../../rotational_order/bin/Q6 -r 1.3 -o Q6.xyz -Q 0.2,1.0
awk '{print $7}' Q6.xyz > tmp.dat
```

