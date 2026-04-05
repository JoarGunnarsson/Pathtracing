#! /bin/bash
set -e

PROJECT_DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
cd $PROJECT_DIR

if test ! -d "venv" ; then
    echo "Virtual environment could not be found. Run 'setup.sh' to set up the project environment first.";
    exit;
fi

if test ! -f "main" ; then
    echo "Main executable could not be found. Run 'setup.sh' to set up the project environment first.";
    exit;
fi

source venv/bin/activate

show_help() {
    echo "usage: main.sh [-h] [--name NAME] [--scene_directory DIRECTORY]"
    echo ""
    echo "options:"
    echo "  -h, --help                            show this message"
    echo "  -n, --name <name>                     the name of the generated image, default 'result.png'"
    echo "  -c, --scene_directory <directory>     path to the scene directory, default scenes/example"
}

IMAGE_NAME="result.png"
SCENE_DIR="scenes/example"

while test "$#" -gt 0; do
    case "$1" in
        -n|--name)
            IMAGE_NAME="$2"
            shift
            ;;
        -s|--scene_directory)
            SCENE_DIR="$2"
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
    shift
done

if [[ "$IMAGE_NAME" != *.png ]]; then
        echo "$IMAGE_NAME is not a valid image name. Image name must end in '.png'. "
        exit 1
fi

echo "Rendering scene '$SCENE_DIR'. The result can be found in images/$IMAGE_NAME"

./main $SCENE_DIR
python python_utils/to_png.py --name $IMAGE_NAME --scene_dir $SCENE_DIR
