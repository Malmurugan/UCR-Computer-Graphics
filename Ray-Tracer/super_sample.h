#ifndef __SUPER_SAMPLE_H__
#define __SUPER_SAMPLE_H__

#include "vec.h"
#include "camera.h"
#include "render_world.h"

class Render_World;

class Super_Sample
{
public:
    unsigned samples;
    Camera& camera;
    Render_World& world;
    
    Super_Sample(Render_World& world_input, Camera& camera_input);

    vec3 Sample(const ivec2& pixel_index);
    void Pixel_Positions(const ivec2& pixel_index, vec3* pixels);
    void Jitter(vec2& pixel_position);

private:
    void Fill_Pixel_Positions(vec2 * pix_pos, const ivec2& pixel_index);
};

#endif
