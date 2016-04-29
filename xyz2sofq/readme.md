# xyz2sofq
Compute the intermediate scattering S(q) of a xyz-formated configuration

by Ulf R. Pedersen (2015), http://urp.dk 

## Build 
    make

## Usage
```
Usage 2D-mode:              ./bin/xyz2sofq [XYZ file] [X] [Y] [X_pixels] [Y_pixels]
Usage 3D-mode (cubic box):  ./bin/xyz2sofq [XYZ file] [X] [X_pixels] 
Usage 3D-mode:              ./bin/xyz2sofq [XYZ file] [X] [Y] [Z] [X_pixels] [Y_pixels] [Z_pixels] 
          * [XYZ file] is a file with the input trajector in the xyz-format.
          * X, Y and Z gives the size of the periodic box.
          * X_pixels, Y_pixels and Z_pixels determine largest q-vector computed.
          * In 3D-mode, all possible q-vectors are computed,
              - a scattering image is not written,
              - but a Sq_binned.dat file is written if the box is cubic.
Notes:
         The type in xyz file must be a integer
```

## Usage examples
```
           ./bin/xyz2sofq traj.xyz 20.0 30.0 201 301
           ./bin/xyz2sofq traj.xyz 20.0 201
           ./bin/xyz2sofq traj.xyz 20.0 19.0 18.0 201 191 181
```



