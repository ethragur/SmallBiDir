#include "scene.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

/******************************************************************
 * Load a OBJ file and build the scene from triangles out of it
 ******************************************************************/
void Scene::load_scene(std::string filename, std::string matfiles)
{
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;

    std::cout << "Loading " << filename << "..." << std::endl;

    // Load OBJ elements with tinyOBJLoader
    std::string err;
    bool ret = tinyobj::LoadObj(shapes, materials, err, filename.c_str(), matfiles.c_str());

    if (!err.empty())
    {
        std::cerr << err << std::endl;
    }

    if (!ret)
    {
        exit(1);
    }

    std::cout << "# of shapes: " << shapes.size() << std::endl;
    std::cout << "# of materials: " << materials.size() << std::endl;

    // Extract elements and build triangles for our scene from them
    std::vector<Vector> vertices;
    for (size_t i = 0; i < shapes.size(); i++)
    {
        // Collect the vertices for this object
        for (size_t v = 0; v < shapes[i].mesh.positions.size() / 3; v++)
        {
            /*
            printf(" v[%ld][%ld] = (%f, %f, %f) Emit: %f\n", i, v,
               shapes[i].mesh.positions[3*v+0],
               shapes[i].mesh.positions[3*v+1],
               shapes[i].mesh.positions[3*v+2],
               materials[shapes[i].mesh.material_ids[0]].emission[0]);
               */

            vertices.push_back(
                    Vector(shapes[i].mesh.positions[3 * v + 0] * 50,
                           shapes[i].mesh.positions[3 * v + 1] * 50,
                           shapes[i].mesh.positions[3 * v + 2] * 50)
            );
        }

        //TODO(max): Add translucency and reflection too
        // Determine color of object
        Color triColor = Color(materials[shapes[i].mesh.material_ids[0]].diffuse[0],
                               materials[shapes[i].mesh.material_ids[0]].diffuse[1],
                               materials[shapes[i].mesh.material_ids[0]].diffuse[2]);

        // Check if its supposed to be a lightsource
        float lightmodifier = 0.0f;
        if (materials[shapes[i].mesh.material_ids[0]].emission[0])
        {
            lightmodifier = materials[shapes[i].mesh.material_ids[0]].emission[0];
        }

        for (size_t b = 0; b < shapes[i].mesh.indices.size(); b += 3)
        {
            //cout << "Vertices left: " << shapes[i].mesh.indices.size() - b << endl;

            Triangle *newtri = new Triangle(vertices[shapes[i].mesh.indices[b]],
                                       vertices[shapes[i].mesh.indices[b + 1]],
                                       vertices[shapes[i].mesh.indices[b + 2]],
                                       triColor * lightmodifier, triColor, DIFF);

            scene_shapes.push_back(newtri);
        }
        // Clear stack for next triangles
        vertices.clear();
    }

    std::cout << "Done Loading!" << std::endl;
}

/*********************************************************
 * Since LoadScene can only load triangles, add a sphere
 * with this function!
 *********************************************************/
void Scene::add_sphere( Sphere *new_sphere )
{
     scene_shapes.push_back(new_sphere);
}

/*********************************************************
 * Returns the scene as an array of shapes.
 *********************************************************/
const std::vector<Shape*> & Scene::get_shapes() const
{
    return scene_shapes;
}

/*
Sphere spheres[] =
        {
                Sphere(1e5, Vector(1e5 + 1, 40.8, 81.6), Vector(), Vector(.75, .25, .25), DIFF),  Left wall
                Sphere(1e5, Vector(-1e5 + 99, 40.8, 81.6), Vector(), Vector(.25, .25, .75), DIFF),  Rght wall
                Sphere(1e5, Vector(50, 40.8, 1e5), Vector(), Vector(.75, .75, .75), DIFF),  Back wall
                Sphere(1e5, Vector(50, 40.8, -1e5 + 170), Vector(), Vector(), DIFF),  Front wall
                Sphere(1e5, Vector(50, 1e5, 81.6), Vector(), Vector(.75, .75, .75), DIFF),  Floor
                Sphere(1e5, Vector(50, -1e5 + 81.6, 81.6), Vector(), Vector(.75, .75, .75), DIFF),  Ceiling

                Sphere(16.5, Vector(27, 16.5, 47), Vector(), Vector(1, 1, 1) * .999, SPEC),  Mirror sphere
                Sphere(16.5, Vector(73, 16.5, 78), Vector(), Vector(1, 1, 1) * .999, REFR),  Glas sphere

                Sphere(1.5, Vector(50, 81.6 - 16.5, 81.6), Vector(4, 4, 4) * 100, Vector(), DIFF),  Light
        };
*/

