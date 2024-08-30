#!/bin/bash

# Compile the C file
echo "Compiling."
g++ -std=c++11 src/main.cpp -o main -O3
echo "Finished compiling."
echo "Running program."
./main > temp/log.txt && python python_utils/to_png.py