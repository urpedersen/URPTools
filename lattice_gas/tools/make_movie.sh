#!/bin/bash

mkdir ppm
mv conf*.ppm ppm

convert -delay 10 ppm/conf*ppm movie.gif

tar -czf ppm.tar.gz ppm
rm -r ppm
