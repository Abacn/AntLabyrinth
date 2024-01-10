#!/bin/sh

# cp ant ../bin/ant
# cp msdlim ../bin/msdlim
make DIM=2D $@
make DIM=3D $@
make DIM=4D $@
make DIM=5D $@
make DIM=6D $@
make DIM=7D $@
