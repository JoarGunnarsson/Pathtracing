#!/bin/bash

# Compile the C file
g++ -std=c++11 main.cpp -o main -O3

# Execute the compiled program and redirect output
./main > Images/result.ppm