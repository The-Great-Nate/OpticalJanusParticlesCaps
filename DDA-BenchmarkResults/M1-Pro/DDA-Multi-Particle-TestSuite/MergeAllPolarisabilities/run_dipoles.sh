#!/bin/bash

for i in 4224 552 160 56 32 8 
do
    python DipolesMulti2025Eigen.py bigjanus_${i}dip
done

#################### NOTES #####################
# 552 took up more than 64GB memory lol (if using 32 particles)
# 4224 Dipoles on 1 particle eats 9GB.