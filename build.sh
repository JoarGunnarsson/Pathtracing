#!/bin/bash

SOURCE_DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
BUILD_DIR=$SOURCE_DIR/build


show_help() {
    echo "usage: build.sh [-h] [--clean]"
    echo ""
    echo "options:"
    echo "  -h, --help              show this message"
    echo "  -c, --clean             cleans the build directory before building"
}


clean(){
    if [ -d $BUILD_DIR ]; then
        rm -rf $BUILD_DIR/*
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


cmake -S $SOURCE_DIR -B $BUILD_DIR
make -C $BUILD_DIR
