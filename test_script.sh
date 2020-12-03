#!/usr/bin/env bash
ifort interaction_module3.f90 4th_prog_test.f90 -mkl
rm rndom
./a.out
