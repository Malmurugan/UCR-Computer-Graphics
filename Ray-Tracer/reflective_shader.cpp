#include "reflective_shader.h"
#include "ray.h"
#include "render_world.h"


vec3 Reflective_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& normal,int recursion_depth) const 
    {
    vec3 color = shader->Shade_Surface(ray, intersection_point,normal, recursion_depth);
    vec3 direc = ray.direction;
    if (recursion_depth >= world.recursion_depth_limit) {
        return (1 - reflectivity) * color;
    }
    Ray r(intersection_point, direc - 2 * dot(direc, normal) * normal);
    return (1 - reflectivity) * color + reflectivity * world.Cast_Ray(r, ++recursion_depth);
}
 