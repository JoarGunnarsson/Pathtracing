#!/bin/bash

show_help() {
    echo "Usage: ./main.sh [OPTIONS] [ARGS]"
    echo ""
    echo "Options:"
    echo "  -c, --compile           Compiles the project before running. (optional)"
    echo "  -n, --name <name>       Specify a name that ends in '.png', default 'result.png' (optional)"
    echo "  -h, --help              Show this message (optional)"
    echo ""
    echo "Example:"
    echo "  ./main.sh --compile --name 'result.png'"
}


compile(){
    echo "Compiling."
    clang++ -std=c++11 src/*.cpp -o main -O3
    echo "Finished compiling."
}


name="result.png"
echo ""

# Parse arguments. Compiles project if compile flag is set. Also sets resulting image name if provided.
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -c|--compile)
            compile
            ;;
        -n|--name)
            name="$2"
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 0
            ;;
    esac
    shift
done

if [[ "$name" != *.png ]]; then
        echo "$name is not a valid image name. Image name must end in '.png'. "
        exit 1
fi

echo "Running program. The result can be found in Images/$name"
width=$(./main)
python python_utils/to_png.py --name $name --width $width
