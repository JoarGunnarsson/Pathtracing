FROM alpine:latest

WORKDIR /pathtracing

RUN apk --no-cache add bash build-base make cmake python3

COPY ./*.sh /pathtracing/
COPY ./CMakeLists.txt /pathtracing/
COPY ./requirements.txt /pathtracing/
RUN cd /pathtracing && ./setup_venv.sh

COPY ./include /pathtracing/include
COPY ./src /pathtracing/src
COPY ./app /pathtracing/app
COPY ./python_utils /pathtracing/python_utils

RUN ./build.sh

RUN mkdir /pathtracing/temp

COPY ./scenes /pathtracing/scenes
