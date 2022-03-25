#include "mesh.h"
#include "plane.h"
#include <fstream>
#include <string>
#include <limits>


// Consider a triangle to intersect a ray if the ray intersects the plane of the
// triangle with barycentric weights in [-weight_tolerance, 1+weight_tolerance]
static const double weight_tolerance = 1e-4;

// Read in a mesh from an obj file.  Populates the bounding box and registers
// one part per triangle (by setting number_parts).
void Mesh::Read_Obj(const char* file)
{
    std::ifstream fin(file);
    if(!fin)
    {
        exit(EXIT_FAILURE);
    }
    std::string line;
    ivec3 e;
    vec3 v;
    box.Make_Empty();
    while(fin)
    {
        getline(fin,line);

        if(sscanf(line.c_str(), "v %lg %lg %lg", &v[0], &v[1], &v[2]) == 3)
        {
            vertices.push_back(v);
            box.Include_Point(v);
        }

        if(sscanf(line.c_str(), "f %d %d %d", &e[0], &e[1], &e[2]) == 3)
        {
            for(int i=0;i<3;i++) e[i]--;
            triangles.push_back(e);
        }
    }
    number_parts=triangles.size();
}

Hit Mesh::Intersection(const Ray& ray, int part) const
{
    Hit hit;
    if (part >= 0) {
        if (Intersect_Triangle(ray, part, hit.dist)) {
           return hit;
        }
    } else {
        hit.dist = std::numeric_limits<double>::max();
        for (unsigned i = 0; i < triangles.size(); i++) {
            double dist;
            if (Intersect_Triangle(ray, i, dist)) {
                if (dist < hit.dist) {
                    hit.object = this;
                    hit.dist = dist;
                    hit.part = i;
                }
            }
        }        
    }
    return hit;
}

// Compute the normal direction for the triangle with index part.
vec3 Mesh::Normal(const vec3& point, int part) const
{
    assert(part>=0);
    ivec3 tri = triangles[part];
    vec3 a = vertices[tri[1]] - vertices[tri[0]];
    vec3 b = vertices[tri[2]] - vertices[tri[0]];
    return cross(a, b).normalized();
}

bool Mesh::Intersect_Triangle(const Ray& ray, int tri, double& dist) const
{
   
    Hit hit = Plane(vertices[triangles[tri][0]], Normal(vertices[triangles[tri][0]], tri)).Intersection(ray, tri);

    if (!hit.object) {
        return false;
    }

    vec3 u = ray.direction;
    vec3 v =  vertices[triangles[tri][1]] -  vertices[triangles[tri][0]];
    vec3 w =  vertices[triangles[tri][2]] -  vertices[triangles[tri][0]];
    vec3 y = ray.Point(dist) -  vertices[triangles[tri][0]];

    double denom = dot(cross(u, v), w);
    
    if (!denom) {
        return false;}

    double beta = dot(cross(w, u), y) / denom;
    double gamma = dot(cross(u, v), y) / denom;
    double alpha = 1 - (gamma + beta);

    if (gamma > -weight_tol && beta > -weight_tol && alpha > -weight_tol) 
    {
        dist = hit.dist;
        return true;
    }
    else{
    return false;}
}

// Compute the bounding box.  Return the bounding box of only the triangle whose
// index is part.
Box Mesh::Bounding_Box(int part) const
{
    Box b;
    TODO;
    return b;
}
