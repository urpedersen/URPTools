This folder contains the tools
  * xyz2gro
  * gro2xyz
  * xyz2pov

xyz2gro and gro2xyz
===========================
Two script that converts between xyz-files and the [Gromacs format](http://manual.gromacs.org/current/online/gro.html).

xyz2pov
===========================

A bash script that make an image of a xyz-file using [POV-ray](www.povray.org).
The below information is given with the `-h` switch. 

```

  Convert an xyz configuration to a image using POV-ray [www.povray.org].
    by Ulf R. Pedersen (http://urp.dk)

Usage: ./Tools/xyz2pov [-i start.xyz.gz] [-o start]

  -h,         Print this help message.
  -i FILE     Input configuration.
                The default is start.xyz.gz. Also accept *.xyz files.
  -B NUM,NUM,NUM
              Set the size of box: Lx,Ly,Lz.
              By default it is assumed that that the header of the xyz-file
                has the info like this: sim_box=RectangularSimulationBox,10,10,10
  -f INT      Frame index.
                The default is 0 (the first frame).
  -p FILE     Povray input header file.
                The default is myPovray.ini.
                A default file is generated if the file does not exist.
  -o STR      Output file name.
                Will generate STR.pov and STR.png.
                STR=INPUT_NAME gives name of the input configuration.
  -a NUM      Camera angle position (radians).
                The default is -pi/8
  -y NUM      Camera height
                The default is 0.66*Ly
  -r NUM      Camera radial distance
                The default is 2.0*Lx
  -L NUM,NUM,NUM
              Camera look_at position
                The default is 0,-0.15*Ly,0
  -W INT,     Width of output image. Height is set to 3/4 of width.
                The default is 800 (height of 600).
  -d,         Disable running povray (but do generate *.pov file).

Usage example:
  ./Tools/xyz2pov -d -i start.xyz.gz -o start
  povray +P +W1600 +H1200 +A +HImyPovray.ini +Istart.pov

Caveat: Default myPovray.ini file can only handle up to 16 particle types. 

Dependency: POV-ray [www.povray.org].
```


### How to generate movie of a configuration

A movie of a configuration can be generated using [ffmpeg](http://www.ffmpeg.org/)
together with the `xyz2pov` script and [POV-ray](http://www.povray.org/). 
We will generate a POV-ray scene where the position of the camera is given by 
the value `clock` (that will vary from 0-1 during the movie).

1. Generate an example configuration (`start.xyz.gz`)

```
  make_lattice
```

2. Generate a POV-ray scene where the angular position of the camera
   is set using the `clock` value that will range from 0-1.

```
  xyz2pov -d -a 2*pi*clock -o scene  `
```

3. Generate 100 images for the movie using [POV-ray](http://www.povray.org/). 
   The option `+KFF100` dictates that 100 images should be made,
   `+W400 +H300` set the pixel size, and `+KC` that the movie can be looped.

```
    povray -D +W400 +H300 +HImyPovray.ini +Iscene.pov +KFF100 +KC +Oframe.png
```

4. Use [ffmpeg](https://www.ffmpeg.org/) to generate a movie from images

```
    ffmpeg -i frame%03d.png movie.mp4 
```

5. Clean up by deleting images

```
    rm frame???.png
```

6. View movie, quit by pressing `q` (the `-loop 0` option makes the movie loop indefinitely)

```
    ffplay -loop 0 movie.mp4
```

Hint: You can edit `myPovray.ini` and `scene.pov` between step 2 and 3 to modify the scene.

Implementation
-----------------------
The `xyz2pov` tool is a bash script that read a NAME.xyz.gz or NAME.xyz file
and writes a file for [POV-ray](http://www.povray.org/).
