#!/bin/bash

tests=(
    BooleanMask
    LinearSearch
    MergeAllPolarisabilities
)

for test in "${tests[@]}"
do
    echo "Entering $test"
    cd "$test" || exit 1
    cd "pyoptics" || exit 1
    ./compile_parallel.sh
    cd ../..
done


dipoles=(8 32 56 160 552 4224)

for test in "${tests[@]}"
do
    echo "Entering $test"
    cd "$test" || exit 1

    for n in "${dipoles[@]}"
        do
            echo "  Running bigjanus_${n}dip"
            python3 DipolesMulti2025Eigen.py "bigjanus_${n}dip"
        done

    cd ..
done

