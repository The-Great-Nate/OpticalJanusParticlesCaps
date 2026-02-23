#!/bin/bash
set -e  # stops on first error

g++ -std=c++14 -O3 -fopenmp -Wall -Wextra -pedantic -c -fPIC beamslib.cpp -o beamslib.o
g++ -O3 -shared -fopenmp beamslib.o -o beamslib.so
g++ -std=c++14 -O3 -fopenmp -Wall -Wextra -pedantic -c -fPIC dipoleslib.cpp -o dipoleslib.o
g++ -O3 -shared -fopenmp dipoleslib.o beamslib.o -o dipoleslib.so

echo "Build complete."
