# Examples of using the hard-sphere program

Perform a simulation from an ideal gas configuration
```
  ../bin/hard_spheres -i test.xyz -p 1.0 -v 0.1 
```

Compute the radial distribution function
```
  ../../traj/bin/traj_gr --load_traj_type=3 --bboxX=16
```

Clean-up
```
  mv final.xyz first.xyz
  rm traj.xyz
```

Perform equlibrium simulation
```
  ../bin/hard_spheres -i first.xyz
```

Re-compute the radial distribution function
```
  ../../traj/bin/traj_gr --load_traj_type=3 --bboxX=14.5097 
```

