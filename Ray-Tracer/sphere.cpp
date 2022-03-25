#include "sphere.h"
#include "ray.h"

// Determine if the ray intersects with the sphere
Hit Sphere::Intersection(const Ray& ray, int part) const
{
    
    Hit intersection = {0, 0, 0};
    vec3 end_point = ray.endpoint - center;
 
    double a = dot(ray.direction,ray.direction);
    double b = 2 * (dot(end_point,ray.direction));
    double c = dot(end_point,end_point);
    c -= radius*radius;
    
    double disc = (b*b) - (4*a*c);
    if(disc>0){
    	double posdis = (-b + sqrt(disc))/(2*a);
    	double negdis = (-b - sqrt(disc))/(2*a);
    	if(posdis>0 && negdis>0){
			intersection.object = this;
			intersection.dist = negdis;
			intersection.part = part;
		}
	}
    return intersection;
}

vec3 Sphere::Normal(const vec3& point, int part) const
{
    vec3 normal;
    vec3 temp = (point - center) / radius;
    normal = temp.normalized(); // compute the normal direction
    return normal;
}

Box Sphere::Bounding_Box(int part) const
{
    Box box;
    vec3 radi(radius,radius,radius);
    box.lo = center - radi;
    box.hi = center + radi;
    return box;
}
