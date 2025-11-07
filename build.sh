#!/bin/bash

show_help() {
    echo "usage: main.sh [-h] [--clean]"
    echo ""
    echo "options:"
    echo "  -h, --help              show this message (optional)"
    echo "  -c, --clean             cleans the build directory before building"
}

clean(){
    if [ -d build ]; then
        rm -r build/*
    fi
}

# Parse arguments. Compiles project if compile flag is set. Also sets resulting image name if provided.
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -c|--clean)
            clean
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


cmake -B build
make -C build
