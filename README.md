# C++ Pathtracing

This is a path tracer, written in C++, aiming to create realistic renders of a scene by simulating the path light takes. Primitives such as spheres, planes, rectangles and triangles have been implemented. Objects can be specified in a .obj file, and will be converted to a series of triangle primitives. Multiple different material types have also been implemented; diffuse materials, reflective materials, transparent materials, and microfacet materials. In order to enable rendering complex 3D models, a Bounding Volume Hierarcy has also been implemented.


### Example scene
Below are a few different example scenes, showcasing different objects and materials.


| ![Example 1](Images/Example1.png) | ![Example 2](Images/Example2.png) |
|:----------------------------------:|:----------------------------------:|
|             Example 1              |             Example 2              |
| ![Example 3](Images/Example3.png) | ![Example 4](Images/Example4.png) |
|             Example 3              |             Example 4              |


## Usage

#### Main program
To run the ray tracing simulation and generate an image, simply execute the shell script `main.sh` file:

```
usage: main.sh [-h] [--name NAME] [--compile]

options:
  -h, --help              show this message (optional)
  -c, --compile           compiles the project before running (optional)
  -n, --name <name>       the filename of the generated image, default 'result.png' (optional)
```

#### Utilities

##### get_map.py

`get_map.py` is a utility program used for turning images into UV-maps for different kinds of material properties, for example albedo.

```
Usage: get_map.py [-h] [-m MODE] in_file out_file

positional arguments:
  in_file               input file path, relative to maps/
  out_file              output file path, relative to maps/

options:
  -h, --help            show this help message and exit
  -m MODE, --mode MODE  which mode to use. Available modes are: 'albedo'
```


#### to_png.py

`to_png.png`is a utility program used to convert the raw image data from the Pathtracing program into a proper image format. Can be used while the Patracing program is running in order to view progress.

```
usage: to_png.py [-h] [--name NAME] [--width WIDTH]

options:
  -h, --help     show this help message and exit
  --name NAME    name of the resulting image.
  --width WIDTH  the width of the image.
```


### Notes

This project does not explicitly support objects intersecting other objects, and can result in inaccurate results in regard to transparent.

Furthermore, this project was tested and compiled on MacOS using the clang compiler, and is not guaranteed to work on other platforms or with other compilers.
