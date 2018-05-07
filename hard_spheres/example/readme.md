# Examples of using the hard-sphere program

Perform a simulation from an ideal gas configuration
```
  ../bin/hard_spheres -L 16
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
  ../bin/hard_spheres -L 16
```

Re-compute the radial distribution function
```
  ../../traj/bin/traj_gr --load_traj_type=3 --bboxX=16
```

