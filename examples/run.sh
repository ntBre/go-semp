#!/bin/bash

semp -atoms         "" \
     -geoms     file07 \
     -energies rel.dat \
     -params   opt.out \
     -debug=false      \
     -cpu           "" \
     -gauss        g16 \
     -lambda      1e-8 \
     -maxit        250 \
     -one           "" \
     semp & disown -h
