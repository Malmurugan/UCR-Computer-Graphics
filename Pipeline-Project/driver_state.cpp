#include "driver_state.h"
#include <cstring>

using namespace std;

float get_tri_area(vec2 A, vec2 B, vec2 C);
void ver_new(driver_state& state, data_geometry* triangle, const data_geometry* x, const data_geometry* y, int plane, bool pos, float* newd);

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    int size = width * height;
	state.image_color = new pixel[size];
	state.image_depth = new float[size];

	for (int i = 0; i < (size); i++) {
		state.image_color[i] = make_pixel(0, 0, 0);
        state.image_depth[i] = 1;
    }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    
    switch (type) {           
        case(render_type::triangle): {

            data_vertex* v_data[state.num_vertices];
            data_geometry* v_geo[state.num_vertices];
            for (int i = 0, j = 0; i < state.num_vertices* state.floats_per_vertex; i += state.floats_per_vertex, j++) {
                v_data[j] = new data_vertex[MAX_FLOATS_PER_VERTEX];
                v_data[j]->data = &state.vertex_data[i];
            }
            
            for (int i = 0; i < state.num_vertices; i++) {
                v_geo[i] = new data_geometry();
                v_geo[i]->data = v_data[i]->data;
                state.vertex_shader(*v_data[i], *v_geo[i], state.uniform_data);
            }

            for (int i = 0; i < state.num_vertices / 3; i++) {
                const data_geometry** pair_vertices = const_cast<const data_geometry**>(v_geo + (3 * i));
                clip_triangle(state, *pair_vertices[0], *pair_vertices[1], *pair_vertices[2], 0);		
            }

            break;
        }
        case(render_type::indexed): {
		
		data_vertex* v_data[state.num_vertices];
		for (int i = 0, j = 0; i < state.num_vertices* state.floats_per_vertex; i += state.floats_per_vertex, j++) {
			v_data[j] = new data_vertex[MAX_FLOATS_PER_VERTEX];
			v_data[j]->data = &state.vertex_data[i];
		}

		data_geometry* v_geo[state.num_triangles * 3];
		for (int i = 0; i < state.num_triangles * 3; i++) {
			v_geo[i] = new data_geometry();
			v_geo[i]->data = v_data[state.index_data[i]]->data;
			state.vertex_shader(*v_data[state.index_data[i]], *v_geo[i], state.uniform_data);
		}
		for (int i = 0; i < state.num_triangles; i++) {
			const data_geometry** pair_vertices = const_cast<const data_geometry**>(v_geo + (3 * i));
			clip_triangle(state, *pair_vertices[0], *pair_vertices[1], *pair_vertices[2], 0);
		}

		break;
	}
        case(render_type::fan): {
	
		data_vertex* v_data[state.num_vertices];
		for (int i = 0, j = 0; i < state.num_vertices* state.floats_per_vertex; i += state.floats_per_vertex, j++) {
			v_data[j] = new data_vertex[MAX_FLOATS_PER_VERTEX];
			v_data[j]->data = &state.vertex_data[i];
		}

		data_geometry* v_geo[state.num_vertices];
		for (int i = 0; i < state.num_vertices; i++) {
			v_geo[i] = new data_geometry();
			v_geo[i]->data = v_data[i]->data;
			state.vertex_shader(*v_data[i], *v_geo[i], state.uniform_data);
		}

		for (int i = 1, j = 2; i < state.num_vertices - 1; i++, j++) {
			data_geometry* tri_grp[3];
			tri_grp[0] = v_geo[0];
			tri_grp[1] = v_geo[i];
			tri_grp[2] = v_geo[j];

			const data_geometry** pair_vertices = const_cast<const data_geometry**>(tri_grp);
			clip_triangle(state, *pair_vertices[0], *pair_vertices[1], *pair_vertices[2], 0);
		}


		break;
	}
        case(render_type::strip): {

		data_vertex* v_data[state.num_vertices];
		for (int i = 0, j = 0; i < state.num_vertices* state.floats_per_vertex; i += state.floats_per_vertex, j++) {
			v_data[j] = new data_vertex[MAX_FLOATS_PER_VERTEX];
			v_data[j]->data = &state.vertex_data[i];
		}

		data_geometry* v_geo[state.num_vertices];
		for (int i = 0; i < state.num_vertices; i++) {
			v_geo[i] = new data_geometry();
			v_geo[i]->data = v_data[i]->data;
			state.vertex_shader(*v_data[i], *v_geo[i], state.uniform_data);
		}

		for (int i = 0, j = 1, k = 2; i < state.num_vertices - 2; i++, j++, k++) {
			data_geometry* tri_grp[3];
			tri_grp[0] = v_geo[i];
			tri_grp[1] = v_geo[j];
			tri_grp[2] = v_geo[k];

			const data_geometry** pair_vertices = const_cast<const data_geometry**>(tri_grp);
			clip_triangle(state, *pair_vertices[0], *pair_vertices[1], *pair_vertices[2], 0);
		}

		break;
	}
        default:;
	}

}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,const data_geometry& v1, const data_geometry& v2,int face)
{   int plane = 0;
    bool a_inter = false, b_inter = false, c_inter = false, positive = true;
	
	if (face == 6)
	{
		rasterize_triangle(state,v0,v1,v2);
		return;
	}

	switch (face)
	{
	case 0:
		positive = true;
		plane = 0;
		a_inter = ( v0.gl_Position[plane] <=  v0.gl_Position[3]) ? true : false;
		b_inter = ( v1.gl_Position[plane] <=  v1.gl_Position[3]) ? true : false;
		c_inter = ( v2.gl_Position[plane] <=  v2.gl_Position[3]) ? true : false;
		break;
	case 1:
		positive = false;
		plane = 0;
		a_inter = ( v0.gl_Position[plane] >= - v0.gl_Position[3]) ? true : false;
		b_inter = ( v1.gl_Position[plane] >= - v1.gl_Position[3]) ? true : false;
		c_inter = ( v2.gl_Position[plane] >= - v2.gl_Position[3]) ? true : false;
		break;
	case 2:
		positive = true;
		plane = 1;
		a_inter = ( v0.gl_Position[plane] <=  v0.gl_Position[3]) ? true : false;
		b_inter = ( v1.gl_Position[plane] <=  v1.gl_Position[3]) ? true : false;
		c_inter = ( v2.gl_Position[plane] <=  v2.gl_Position[3]) ? true : false;
		break;
	case 3:
		positive = false;
		plane = 1;
		a_inter = ( v0.gl_Position[plane] >= - v0.gl_Position[3]) ? true : false;
		b_inter = ( v1.gl_Position[plane] >= - v1.gl_Position[3]) ? true : false;
		c_inter = ( v2.gl_Position[plane] >= - v2.gl_Position[3]) ? true : false;
		break;
	case 4:
		positive = true;
		plane = 2;
		a_inter = ( v0.gl_Position[plane] <=  v0.gl_Position[3]) ? true : false;
		b_inter = ( v1.gl_Position[plane] <=  v1.gl_Position[3]) ? true : false;
		c_inter = ( v2.gl_Position[plane] <=  v2.gl_Position[3]) ? true : false;
		break;
	case 5:
		positive = false;
		plane = 2;
		a_inter = ( v0.gl_Position[plane] >= - v0.gl_Position[3]) ? true : false;
		b_inter = ( v1.gl_Position[plane] >= - v1.gl_Position[3]) ? true : false;
		c_inter = ( v2.gl_Position[plane] >= - v2.gl_Position[3]) ? true : false;
		break;
	default:
		break;
	}

	data_geometry* tri_1[3];
	data_geometry* tri_2[3];
	for (int i = 0; i < 3; i++) {
		tri_1[i] = new data_geometry();
		tri_2[i] = new data_geometry();
	}
	float *d1 = new float[MAX_FLOATS_PER_VERTEX], *d2 = new float[MAX_FLOATS_PER_VERTEX];

	if (a_inter && b_inter && c_inter) {
		clip_triangle(state,  v0,v1,v2, face + 1);
	}
	else if (a_inter && b_inter && !c_inter) {

		tri_1[0]->gl_Position =  v0.gl_Position; 
		tri_1[0]->data = v0.data;				

		tri_1[1]->gl_Position =  v1.gl_Position;
		tri_1[1]->data = v1.data;				

		ver_new(state, tri_1[2], &v0, &v2, plane, positive, d1); 

		tri_2[0]->gl_Position =  v1.gl_Position; 
		tri_2[0]->data = v1.data;				

		ver_new(state, tri_2[1], &v1, &v2, plane, positive, d2); 

		tri_2[2]->gl_Position = tri_1[2]->gl_Position;  
		tri_2[2]->data = tri_1[2]->data;			


		clip_triangle(state, *tri_1[0], *tri_1[1], *tri_1[2], face + 1);
		clip_triangle(state, *tri_2[0], *tri_2[1], *tri_2[2], face + 1);
	}
	else if (a_inter && !b_inter && c_inter) {
		tri_1[0]->gl_Position =  v2.gl_Position;	
		tri_1[0]->data = v2.data;				

		tri_1[1]->gl_Position =  v0.gl_Position; 
		tri_1[1]->data = v0.data;				

		ver_new(state, tri_1[2], &v2, &v1, plane, positive, d1); 

		tri_2[0]->gl_Position =  v0.gl_Position; 
		tri_2[0]->data = v0.data;				

		ver_new(state, tri_2[1], &v0, &v1, plane, positive, d2); 

		tri_2[2]->gl_Position = tri_1[2]->gl_Position;  
		tri_2[2]->data = tri_1[2]->data;				


		clip_triangle(state, *tri_1[0], *tri_1[1], *tri_1[2], face + 1);
		clip_triangle(state, *tri_2[0], *tri_2[1], *tri_2[2], face + 1);
	}
	else if (!a_inter && b_inter && c_inter) {
		tri_1[0]->gl_Position =  v1.gl_Position;	
		tri_1[0]->data = v1.data;				
		tri_1[1]->gl_Position =  v2.gl_Position; 
		tri_1[1]->data = v2.data;				

		ver_new(state, tri_1[2], &v1, &v0, plane, positive, d1); 

		tri_2[0]->gl_Position =  v2.gl_Position; 
		tri_2[0]->data = v2.data;				

		ver_new(state, tri_2[1], &v2, &v0, plane, positive, d2); 

		tri_2[2]->gl_Position = tri_1[2]->gl_Position;	
		tri_2[2]->data = tri_1[2]->data;				


		clip_triangle(state, *tri_1[0],*tri_1[1],*tri_1[2], face + 1);
		clip_triangle(state, *tri_2[0],*tri_2[1],*tri_2[2], face + 1);
	}
	else if (a_inter && !b_inter && !c_inter) {
		tri_1[0]->gl_Position = v0.gl_Position; 
		tri_1[0]->data = v0.data;				
		ver_new(state, tri_1[1], &v0, &v1, plane, positive, d1); 
		ver_new(state, tri_1[2], &v0, &v2, plane, positive, d2); 
		clip_triangle(state, *tri_1[0], *tri_1[1], *tri_1[2], face + 1);
	}
	else if (!a_inter && b_inter && !c_inter) {
		tri_1[0]->gl_Position = v1.gl_Position; 
		tri_1[0]->data = v1.data;				
		ver_new(state, tri_1[1], &v1, &v2, plane, positive, d1);
		ver_new(state, tri_1[2], &v1, &v0, plane, positive, d2); 
		clip_triangle(state, *tri_1[0], *tri_1[1], *tri_1[2], face + 1);
	}
	else if (!a_inter && !b_inter && c_inter) {
		tri_1[0]->gl_Position = v2.gl_Position; 
		tri_1[0]->data = v2.data;				
		ver_new(state, tri_1[1], &v2, &v0, plane, positive, d1); 
		ver_new(state, tri_1[2], &v2, &v1, plane, positive, d2); 
		clip_triangle(state, *tri_1[0], *tri_1[1], *tri_1[2], face + 1);
	}

}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    int width = state.image_width;
	int height = state.image_height;

  
    data_geometry v_array[3] = {v0,v1,v2};
    vec2 p_ind[3];
    for (int i = 0; i < 3; i++) {
		p_ind[i][0] = ((width / 2) * v_array[i].gl_Position[0] / v_array[i].gl_Position[3]) + ((width / 2) - (0.5));
		p_ind[i][1] = ((height / 2) * v_array[i].gl_Position[1] / v_array[i].gl_Position[3]) + ((height / 2) - (0.5));
	}
    

	int x_min = min(min(p_ind[0][0], p_ind[1][0]), p_ind[2][0]);
	int y_min = min(min(p_ind[0][1], p_ind[1][1]), p_ind[2][1]);
	int x_max = max(max(p_ind[0][0], p_ind[1][0]), p_ind[2][0]);
	int y_max = max(max(p_ind[0][1], p_ind[1][1]), p_ind[2][1]);

	float abc_area = get_tri_area(p_ind[0], p_ind[1], p_ind[2]);

	for (int i = x_min; i <= x_max; i++) {
		for (int j = y_min; j <= y_max; j++) {

			float alpha = 0;
            float beta = 0;
            float gamma = 0;

			vec2 curr_p(i, j);
			
			alpha = get_tri_area(curr_p, p_ind[1], p_ind[2]) / abc_area;
			beta = get_tri_area(curr_p, p_ind[2], p_ind[0]) / abc_area;
			gamma = get_tri_area(curr_p, p_ind[0], p_ind[1]) / abc_area;

			if ((alpha >= 0) && (beta >= 0) && (gamma >= 0)) {

				float z_point = (alpha * v0.gl_Position[2] / v0.gl_Position[3]) +(beta * v1.gl_Position[2] / v1.gl_Position[3]) + (gamma * v2.gl_Position[2] / v2.gl_Position[3]);
				if (z_point < state.image_depth[(width * j) + i]) {

					state.image_depth[(width * j) + i] = z_point;

					data_output* color_f = new data_output();
					data_fragment* color_data = new data_fragment();
					float* inter_color_data = new float[MAX_FLOATS_PER_VERTEX];

					for (int n = 0; n < state.floats_per_vertex; n++) {

						switch (state.interp_rules[n]) {
						case(interp_type::noperspective): { 
							inter_color_data[n] = (alpha * v0.data[n] + beta * v1.data[n] + gamma * v2.data[n]);
							break;
						}
						case(interp_type::smooth): { 
							float new_alpha = 0, new_beta = 0, new_gamma = 0, c = 0;

							c = (alpha / v0.gl_Position[3]) + (beta / v1.gl_Position[3]) + (gamma / v2.gl_Position[3]);

							new_alpha = alpha / (v0.gl_Position[3] * c);
							new_beta = beta / (v1.gl_Position[3] * c);
							new_gamma = gamma / (v2.gl_Position[3] * c);

							inter_color_data[n] = (new_alpha * v0.data[n] + new_beta * v1.data[n] + new_gamma * v2.data[n]);
							break;
						}
						case(interp_type::flat): { 

							inter_color_data[n] = v0.data[n];
							break;
						}
						default:;
						}
					}

					color_data->data = inter_color_data;
					const data_fragment* c_color_data = const_cast<const data_fragment*>(color_data);
					state.fragment_shader(*c_color_data, *color_f, state.uniform_data);
					state.image_color[(width * j) + i] = make_pixel(color_f->output_color[0] * 255, color_f->output_color[1] * 255, color_f->output_color[2] * 255);

				}
			}

		}
	}
 
}
void ver_new(driver_state& state, data_geometry* tri, const data_geometry* x, const data_geometry* y, int plane, bool pos, float* newd) {

	float sm_a = 0;
	float np_a = 0;

	if (pos)
		{
			sm_a = (y->gl_Position[3] - y->gl_Position[plane]);
			sm_a = sm_a / (x->gl_Position[plane] - x->gl_Position[3] + y->gl_Position[3] - y->gl_Position[plane]);
		}
	else
		{
			sm_a = (-y->gl_Position[3] - y->gl_Position[plane]);
			sm_a = sm_a / (x->gl_Position[plane] + x->gl_Position[3] - y->gl_Position[3] - y->gl_Position[plane]);
		}

	tri->gl_Position = sm_a * x->gl_Position + (1 - sm_a) * y->gl_Position;

	np_a = sm_a * x->gl_Position[3] / (sm_a * x->gl_Position[3] + (1 - sm_a) * y->gl_Position[3]);

	for (int i = 0; i < state.floats_per_vertex; i++) {
		switch (state.interp_rules[i]) {
		case(interp_type::noperspective): {
			newd[i] = np_a * x->data[i] + (1 - np_a) * y->data[i];
			break;
		}
		case(interp_type::flat): {
			newd[i] = x->data[i];
			break;
		}
		case(interp_type::smooth): {
			newd[i] = sm_a * x->data[i] + (1 - sm_a) * y->data[i];

			break;
		}
		default: ;
		}

	}
	tri->data = newd;

}

float get_tri_area(vec2 A, vec2 B, vec2 C) {
    float Area = 0.5 * (((B[0] * C[1]) - (C[0] * B[1]))  + ((A[0] * B[1]) - (B[0] * A[1])) - ((A[0] * C[1]) - (C[0] * A[1])));
    return Area;
}

