#!/usr/bin/env bash
rm a.out
ifort gaussquad.f90
./a.out
rm -r running
mkdir running

ifort -c interaction_module3.f90
ifort interaction_module3.f90 main_Tully_2009.f90 -mkl

cp fort.23 running
cp 528atom.dat running
cp q_number.dat running
cp gold_neigh_pos.dat running
cp gold_positions_fort.20 running
cp gold_velocities_fort.21 running

cp a.out running
cp parallel_script_v4.py running
cp raw_x.txt running
cp raw_w.txt running
cp job.sh running
cd running


python3 parallel_script_v4.py



