
To build the C++ libraries, you will need to have the Eigen library available.  This is handled via an include file available from:

https://eigen.tuxfamily.org/index.php?title=Main_Page

clang++ -std=c++14 -O3 -Wall -Wextra -pedantic -c -fPIC Beams.cpp -o Beams.o
clang++ -O3 -shared Beams.o -o Beams.dylib
python DipolesMulti2024Eigen.py testfile

g++-14 also works, similar speed, but need to ensure the abs are changed to fabs.

g++-14 -std=c++14 -O3 -Wall -Wextra -pedantic -c -fPIC Beams.cpp -o Beams.o
g++-14 -O3 -shared Beams.o -o Beams.dylib
g++-14 -std=c++14 -O3 -Wall -Wextra -pedantic -c -fPIC Dipoles.cpp -o Dipoles.o
g++-14 -O3 -shared Dipoles.o Beams.o -o Dipoles.dylib
python DipolesMulti2024Eigen.py testfile

Parallel version:

g++-15 -std=c++14 -O3 -fopenmp -Wall -Wextra -pedantic -c -fPIC beamslib.cpp -o beamslib.o
g++-15 -O3 -shared -fopenmp beamslib.o -o beamslib.dylib
g++-15 -std=c++14 -O3 -fopenmp -Wall -Wextra -pedantic -c -fPIC dipoleslib.cpp -o dipoleslib.o
g++-15 -O3 -shared -fopenmp dipoleslib.o beamslib.o -o dipoleslib.dylib
python DipolesMulti2025Eigen.py testfile

testfile.yml is a configuration file.  Results will be returned in testfile.xlsx etc.

On linux, name shared libraries .so rather than .dylib.  Some other renaming of files is needed - ask Simon for details.


Latest version with some boost:

g++-15 -std=c++14 -fopenmp -O3 -I /opt/homebrew/Cellar/boost/1.88.0/include -L /opt/homebrew/Cellar/boost/1.88.0/lib -Wall -Wextra -pedantic -c -fPIC beamslib.cpp -o beamslib.o
g++-15 -O3 -shared -fopenmp beamslib.o -o beamslib.dylib
g++-15 -std=c++14 -fopenmp -O3 -I /opt/homebrew/Cellar/boost/1.88.0/include -L /opt/homebrew/Cellar/boost/1.88.0/lib -Wall -Wextra -pedantic -c -fPIC dipoleslib.cpp -o dipoleslib.o
g++-15 -O3 -shared -fopenmp dipoleslib.o beamslib.o -o dipoleslib.dylib
