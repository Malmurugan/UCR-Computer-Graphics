/*
WARNING: the code that follows will make you cry; a safety pig is provided 
below for your benefit.

                         _
 _._ _..._ .-',     _.._(`))
'-. `     '  /-._.-'    ',/
   )         \            '.
  / _    _    |             \
 |  a    a    /              |
 \   .-.                     ;  
  '-('' ).-'       ,'       ;
     '-;           |      .'
        \           \    /
        | 7  .__  _.-\   \
        | |  |  ``/  /`  /
       /,_|  |   /,_/   /
          /,_/      '`-'
*/

#include "super_sample.h"

#include "globals.h"
#include "assert.h"

#include <cstdlib>

Super_Sample::Super_Sample(Render_World& world_input, Camera& camera_input) 
    :camera(camera_input), world(world_input) {

    samples = pow(2, world.enable_anti_aliasing * 2);
}

vec3 Super_Sample::Sample(const ivec2& pixel_index) {
    vec3 total_color = {0, 0, 0};
    vec3 pixels[samples];
    
    Pixel_Positions(pixel_index, pixels);

    for (unsigned i = 0; i < samples; i++) {
        total_color += world.Cast_Ray(pixels[i]);
    }

    return total_color / samples;
}

void Super_Sample::Pixel_Positions(const ivec2& pixel_index, vec3* pixels) {
    vec2 pix_pos[samples];

    Fill_Pixel_Positions(pix_pos, pixel_index);
    
    for (unsigned i = 0; i < samples; i++) {
        Jitter(pix_pos[i]);
        pix_pos[i] = camera.min + pix_pos[i] * camera.pixel_size;
        pixels[i] = camera.film_position
                  + (camera.horizontal_vector * pix_pos[i][X])
                  + (camera.vertical_vector * pix_pos[i][Y]);
    }

}

void Super_Sample::Jitter(vec2& pixel_position) {
    bool add;
    float amount;
    for (unsigned i = 0; i < 2; i++) {
        add = rand() % 2;
        amount = (rand() % 25) / 100.0f;
        pixel_position[i] += add ? amount : -amount;
    }
}

void Super_Sample::Fill_Pixel_Positions(vec2 * pix_pos,
    const ivec2& pixel_index) {

    // This only works if we use 4 or 16 sampling because I hardcoded
    // the pixel positions
    assert(samples == 4 || samples == 16);
    
    vec2 c = vec2(pixel_index) + vec2(.5f, .5f);

    if (samples == 4) {
        pix_pos[0] = vec2(c[X] - .25f, c[Y] + .25f);
        pix_pos[1] = vec2(c[X] + .25f, c[Y] + .25f);
        pix_pos[2] = vec2(c[X] - .25f, c[Y] - .25f);
        pix_pos[3] = vec2(c[X] + .25f, c[Y] - .25f);
    
    } else if (samples == 16) {
        // Subsection: A
        pix_pos[0] = vec2(c[X] - .125f * 2, c[Y] + .125f * 2);
        pix_pos[1] = vec2(c[X] - .125f * 2, c[Y] + .125f);
        pix_pos[2] = vec2(c[X] - .125f, c[Y] + .125f * 2);
        pix_pos[3] = vec2(c[X] - .125f, c[Y] + .125f);

        // Subsection: B
        pix_pos[4] = vec2(c[X] + .125f * 2, c[Y] + .125f * 2);
        pix_pos[5] = vec2(c[X] + .125f * 2, c[Y] + .125f);
        pix_pos[6] = vec2(c[X] + .125f, c[Y] + .125f * 2);
        pix_pos[7] = vec2(c[X] + .125f, c[Y] + .125f);

        // Subsection: C
        pix_pos[8] = vec2(c[X] - .125f * 2, c[Y] - .125f * 2);
        pix_pos[9] = vec2(c[X] - .125f * 2, c[Y] - .125f);
        pix_pos[10] = vec2(c[X] - .125f, c[Y] - .125f * 2);
        pix_pos[11] = vec2(c[X] - .125f, c[Y] - .125f);

        // Subsection: D
        pix_pos[12] = vec2(c[X] + .125f * 2, c[Y] - .125f * 2);
        pix_pos[13] = vec2(c[X] + .125f * 2, c[Y] - .125f);
        pix_pos[14] = vec2(c[X] + .125f, c[Y] - .125f * 2);
        pix_pos[15] = vec2(c[X] + .125f, c[Y] - .125f);
    }
}
