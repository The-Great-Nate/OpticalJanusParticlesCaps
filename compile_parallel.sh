#!/bin/bash
set -e  # stops on first error



cd pyoptics

# g++ -std=c++14 -O3 -fopenmp -Wall -Wextra -pedantic -c -fPIC beamslib.cpp -o beamslib.o
# g++ -O3 -shared -fopenmp beamslib.o -o beamslib.dylib
# g++ -std=c++14 -O3 -fopenmp -Wall -Wextra -pedantic -c -fPIC dipoleslib.cpp -o dipoleslib.o
# g++ -O3 -shared -fopenmp dipoleslib.o beamslib.o -o dipoleslib.dylib

g++-15 -std=c++17 -O3 -fopenmp -Wall -Wextra -pedantic -c -fPIC beamslib.cpp -o beamslib.o
g++-15 -O3 -shared -fopenmp beamslib.o -o beamslib.dylib
g++-15 -std=c++17 -O3 -fopenmp -Wall -Wextra -pedantic -c -fPIC dipoleslib.cpp -o dipoleslib.o
g++-15 -O3 -shared -fopenmp dipoleslib.o beamslib.o -o dipoleslib.dylib

cd ..

echo "Build complete."
