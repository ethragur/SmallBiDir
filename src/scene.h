#ifndef SCENE_H
#define SCENE_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>


#include "sphere.h"
#include "triangle.h"
#include "ray.h"

class Scene
{
    private:
        std::vector<Shape*> scene_shapes;

    public:
        const std::vector<Shape*> & get_shapes() const;

        void load_scene( std::string filename, std::string matfiles );
        void add_sphere( Sphere *new_sphere );
};

#endif
