#ifndef SCENE_H
#define SCENE_H

#include "objects.h"
#include "materials.h"
#include "medium.h"
#include "camera.h"

struct Scene {
    Object** objects;
    int number_of_objects;
    Camera* camera;
    MaterialManager* material_manager;
    Medium* medium;
};

void load_settings();
Scene create_scene();

#endif
