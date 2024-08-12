# C++ Path traced

This is a path tracer, written in C++, aiming to create realistic renders of a scene by simulating the path light takes. Primitives such as spheres, planes and rectangles have been implemented. Multiple different material types have also been implemented; diffuse materials, reflective materials, and transparent materials.


### Example scene
Below are a few different an example scenes, showcasing the different material types and light sources.


| ![Example 1](Images/Example1.png) | ![Example 2](Images/Ball.png) |
|:----------------------------------:|:----------------------------------:|
|             Example 1              |             Example 2              |
| ![Example 3](Images/NoNee.png) | ![Example 4](Images/example_4.png) |
|             Example 3              |             Example 4              |


## Usage

To run the ray tracing simulation and generate an image, simply execute the shell script `main.sh` file:

```
./main.sh
```

### Notes

This project does not explicitly support objects intersecting other objects, and can result in inaccurate results in regard to transparent and refractive materials. As a consequence of this fact, transparent planar objects should only have a refractive index equal to 1.

Furthermore, this project was tested and compiled on MacOS using the g++ compiler, and is not guaranteed to work on other platforms.