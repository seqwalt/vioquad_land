#!/bin/bash

# Compile code generator, named <project name>_codegen
cd build
cmake ..
make

# Generate the code. This generates a directory called autogen_code
cd ..
./*_codegen

# Move the test file and qpoases into generated code directory
cd autogen_code
cp -r ../test.cpp ../qpoases .

# Compile the test script
make

# While in directory autogen_code, the example can
# now be run by executing ./test

echo "



If no compilation errors occured, run the compiled test with:

$ cd autogen_code
$ ./test
"
