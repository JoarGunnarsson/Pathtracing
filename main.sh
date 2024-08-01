#!/bin/bash

# Compile the C file
echo "Compiling."
g++ -std=c++11 main.cpp -o main -O3
echo "Finished compiling."
echo "Running program."
./main | python to_png.py