#!/bin/bash
set -e

SOURCE_DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
cd $SOURCE_DIR

if [ ! -d "venv" ]; then
    echo "Virtual environment could not be found. Run 'setup.sh' to set up the project environment first.";
    exit;
fi

if [ ! -f "main" ]; then
    echo "Main executable could not be found. Run 'setup.sh' to set up the project environment first.";
    exit;
fi

source venv/bin/activate

show_help() {
    echo "usage: main.sh [-h] [--name NAME] [--scene_file FILE_NAME] [--settings_file FILE_NAME]"
    echo ""
    echo "options:"
    echo "  -h, --help                            show this message"
    echo "  -n, --name <name>                     the name of the generated image, default 'result.png'"
    echo "  -c, --scene_file <file name>          the path to the scene file, default 'scenes/example/scene.json'"
    echo "  -s, --settings_file <file name>       the path to the scene file, default 'scenes/example/settings.json'"
}

IMAGE_NAME="result.png"
SCENE_FILE="scenes/example/scene.json"
SETTINGS_FILE="scenes/example/settings.json"

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -n|--name)
            IMAGE_NAME="$2"
            shift
            ;;
        -c|--scene_file)
            SCENE_FILE="$2"
            shift
            ;;
        -s|--settings_file)
            SETTINGS_FILE="$2"
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

if [[ "$IMAGE_NAME" != *.png ]]; then
        echo "$IMAGE_NAME is not a valid image name. Image name must end in '.png'. "
        exit 1
fi

echo "Rendering scene '$SCENE_FILE', using settings file '$SETTINGS_FILE'. The result can be found in images/$IMAGE_NAME"

./main $SCENE_FILE $SETTINGS_FILE
python python_utils/to_png.py --name $IMAGE_NAME --settings_file $SETTINGS_FILE
