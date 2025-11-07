#!/bin/bash

SOURCE_DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
cd $SOURCE_DIR

show_help() {
    echo "usage: main.sh [-h] [--name NAME] [--compile]"
    echo ""
    echo "options:"
    echo "  -h, --help              show this message (optional)"
    echo "  -n, --name <name>       the filename of the generated image, default 'result.png' (optional)"
}


name="result.png"


# Parse arguments. Compiles project if compile flag is set. Also sets resulting image name if provided.
while [[ "$#" -gt 0 ]]; do
    case "$1" in
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
python ./python_utils/to_png.py --name $name --width $width
