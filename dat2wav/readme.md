# data2wav
Convert data stream into an audio file in the uncompressed PCM WAV format (out.wav)

## Usage:

The program expects a standard input with one float value on each line. The stream should be ended with EOF (ctrl-d). The input data is rescaled so the volume of the wav-output is adjusted to the range of the input. Optional flags are:

 -h or --help
     Print this help message.

 -o [output filename]
     Default: out.wav

 -r [sample rate in Hz]
     Default: 44100. Recommended values are 96000, 48000, 44100, 32000, 22050, 16000, 11025 or 8000.

 -b [bits per sample]
     Default: 16. Accepted values are 8, 16 or 32.

 -c [number of channels]
     Default: 1. Use 2 for stereo. Stereo input is given as L1 R1 L2 R2 L3 R3 ... etc.


## Usage examples: 
```
   echo -e "0.1\n0.2\n0.3" | ./bin/dat2wav

   cat myData.ascii | ./bin/dat2wav -o myData.wav

   paste left.dat right.dat | awk '{print $1;print $2}'  | ./bin/dat2wav -c 2 -o stereo.wav

   awk -F, '{print $2}' data.csv | ./bin/dat2wav

   awk 'BEGIN{twopi=8.*atan2(1.,1.);for(n=0;n<48000;n++){print sin(twopi*n/48000*110*2^(7/12))}}' | ./bin/dat2wav -r 48000 -o E3.wav

   awk 'BEGIN{twopi=8.*atan2(1.,1.);for(n=0;n<32000;n++){print sin(twopi*n/32000*110*2^(7/12));print sin(twopi*n/32000*220*2^(7/12))}}' | ./bin/dat2wav -r 32000 -c 2 -o E3E4_stereo.wav
```

## Copying
Copyright (C) 2013 Ulf R. Pedersen (http://urp.dk). This program comes with ABSOLUTELY NO WARRANTY. 
This is open-source free software (GNU GPL 3.0), and you are welcome to redistribute it under certain conditions (please read LICENSE details).


