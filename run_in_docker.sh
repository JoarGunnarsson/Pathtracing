#! /bin/bash
set -e

docker run --rm \
    --mount type=bind,src=./maps,dst=/pathtracing/maps,readonly \
    --mount type=bind,src=./models,dst=/pathtracing/models,readonly \
    --mount type=bind,src=./images,dst=/pathtracing/images \
    --name pathtracer pathtracing /pathtracing/main.sh $@
