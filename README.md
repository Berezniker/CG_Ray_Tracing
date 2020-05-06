# Ray tracing
![Image](https://github.com/Berezniker/CG_Ray_Tracing/blob/master/out/out_1_scene.jpg)

### Load:
```sh
$ git clone https://github.com/Berezniker/CG_Ray_Tracing.git
$ cd CG_Ray_Tracing
```
### Compile:
```sh
$ mkdir build && cd build
$ cmake .. && make
```
### Run:
```sh
$ ./rt −out <output_path> −scene <scene_number> −threads <threads>
```

### implemented:
- Phong reflection model
- Shadows
- Mirror reflection
- Glass refraction
- Spherical texture
- Triangular meshes
- Anti-aliasing
- Tone mapping
- Distance Aided Ray Marching
- Cubemap

### [geometric primitives](https://github.com/Berezniker/CG_Ray_Tracing/blob/master/objects.h):
- Triangle
- Sphere
- Plane
- Disk
- Ring
- Cylinder
- Hyperboloid
- *Fractal*
