#include "light.h"
#include "phong_shader.h"
#include "ray.h"
#include "render_world.h"
#include "object.h"


using namespace std;

vec3 Phong_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& normal,int recursion_depth) const
{
    vec3 color;
    Hit hit;
    Ray rayz;
    color = world.ambient_color * world.ambient_intensity * color_ambient;
    int count = 0;
    for(unsigned i = 0; i < world.lights.size(); i++) {
        vec3 l = world.lights.at(i)->position - intersection_point;
        if(world.enable_shadows) {
            rayz.endpoint = intersection_point;
            rayz.direction = l.normalized();
            hit = world.Closest_Intersection(rayz);
            if(hit.object != NULL) {
		        if(hit.dist < l.magnitude()){
                    count ++;
            }}
        }
    if (count == 0){
        vec3 ref = (2 * dot(l, normal) * normal - l).normalized();
        vec3 l_color = world.lights.at(i)->Emitted_Light(l.normalized()) / (l.magnitude_squared());
        double diff_int = std::max(0.0, dot(l.normalized(), normal));
        vec3 inter = ray.endpoint - intersection_point;
        double spec_int = pow(std::max(0.0, dot(ref, inter.normalized())), specular_power);
	    color = l_color * (color_diffuse * diff_int + color_specular * spec_int) + color;
    }
    }
    return color;
}