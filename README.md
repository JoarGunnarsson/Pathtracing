# C++ Pathtracing

This project implements Monte Carlo Pathtracing and is written in C++, with the aim of creating realistic renders.


### Objects

Multiple kinds of primitives such as spheres, planes, rectangles, and triangles are implemented. Additionally, complex objects consisting of multiple triangles/quads can be specified in a .obj file using the Wavefront OBJ format. These objects can consist of a large amount of primitives, enabled by the implementation of a Bounding Volume Hierarchy (BVH) datastructure.


### Materials

This project includes multiple different material models, including but not limited to diffuse materials, reflective materials, transparent materials, and microfacet materials. Some simple media are also implemented, allowing the modelling of homogenous scattering media such as clouds, as well as clear media like coloured glass. UV mapping is also supported, allowing the use of complex textures for albedo, roughness, as well as light emission colour & strength.


### Denoising

In order to reduce the noise level of the resulting image, À-trous wavelet denoising has been implemented, with configurable parameters.


### Example scenes
Below are a few different example scenes, showcasing different objects and materials in a cornell box-like scene.

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


### Limitations

This project does not explicitly support transparent objects intersecting eachother, and their inclusion can result in inaccurate results.

Furthermore, this project was tested and compiled on MacOS using the clang compiler, and is not guaranteed to work on other platforms or with other compilers.
