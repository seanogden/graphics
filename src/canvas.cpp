/*
 * canvas.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: nbingham
 */

#include "canvas.h"
#include "core/geometry.h"
#include "light.h"
#include "material.h"

canvashdl::canvashdl(int w, int h)
{
	last_reshape_time = -1.0;
	width = w;
	height = h;
	reshape_width = w;
	reshape_height = h;

	matrices[viewport_matrix] = mat4f((float)width/2.0, 0.0, 0.0, (float)width/2.0,
									  0.0, (float)height/2.0, 0.0, (float)height/2.0,
									  0.0, 0.0, (float)depth/2.0, (float)depth/2.0,
									  0.0, 0.0, 0.0, 1.0);

	initialized = false;

	color_buffer = new unsigned char[width*height*3];
	depth_buffer = new unsigned short[width*height];

	screen_texture = 0;
	screen_geometry = 0;
	screen_shader = 0;

	active_matrix = modelview_matrix;

	for (int i = 0; i < 4; i++)
		matrices[i] = identity<float, 4, 4>();

	polygon_mode = fill;
	shade_model = smooth;
	culling = backface;
}

canvashdl::~canvashdl()
{
	if (color_buffer != NULL)
	{
		delete [] color_buffer;
		color_buffer = NULL;
	}

	if (depth_buffer != NULL)
	{
		delete [] depth_buffer;
		depth_buffer = NULL;
	}
}

void canvashdl::clear_color_buffer()
{
	memset(color_buffer, 0, width*height*3*sizeof(unsigned char));
}

void canvashdl::clear_depth_buffer()
{
	memset(depth_buffer, 255, width*height*sizeof(unsigned short));
}

void canvashdl::reallocate(int w, int h)
{
	last_reshape_time = -1.0;

	if (color_buffer != NULL)
	{
		delete [] color_buffer;
		color_buffer = NULL;
	}

	if (depth_buffer != NULL)
	{
		delete [] depth_buffer;
		depth_buffer = NULL;
	}

	width = w;
	height = h;

	color_buffer = new unsigned char[w*h*3];
	depth_buffer = new unsigned short[w*h];

	glActiveTexture(GL_TEXTURE0);
	check_error(__FILE__, __LINE__);
	glBindTexture(GL_TEXTURE_2D, screen_texture);
	check_error(__FILE__, __LINE__);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, color_buffer);
	check_error(__FILE__, __LINE__);
}

/* set_matrix
 *
 * Change which matrix is active.
 */
void canvashdl::set_matrix(matrix_id matid)
{
	active_matrix = matid;
}

/* load_identity
 *
 * Set the active matrix to the identity matrix.
 * This implements: https://www.opengl.org/sdk/docs/man2/xhtml/glLoadIdentity.xml
 */
void canvashdl::load_identity()
{
	matrices[active_matrix] = identity<float, 4, 4>();
}

/* rotate
 *
 * Multiply the active matrix by a rotation matrix.
 * This implements: https://www.opengl.org/sdk/docs/man2/xhtml/glRotate.xml
 */
void canvashdl::rotate(float angle, vec3f axis)
{
	vec3f u = norm(axis);
	mat3f S(0.0, -u[2], u[1],
			u[2], 0.0, -u[0],
			-u[1], u[0], 0.0);
	mat3f uut(axis[0]*axis[0], axis[0]*axis[1], axis[0]*axis[2],
			  axis[0]*axis[1], axis[1]*axis[1], axis[1]*axis[2],
			  axis[0]*axis[2], axis[1]*axis[2], axis[2]*axis[2]);

	mat4f R = identity<float, 4, 4>();
	R.set(0,3,0,3, uut + (float)cos(angle)*(identity<float, 3, 3>() - uut) + (float)sin(angle)*S);

	matrices[active_matrix] = matrices[active_matrix]*R;
}

/* translate
 *
 * Multiply the active matrix by a translation matrix.
 * This implements: https://www.opengl.org/sdk/docs/man2/xhtml/glTranslate.xml
 */
void canvashdl::translate(vec3f direction)
{
	matrices[active_matrix] = matrices[active_matrix]*mat4f(1.0, 0.0, 0.0, direction[0],
									 0.0, 1.0, 0.0, direction[1],
									 0.0, 0.0, 1.0, direction[2],
									 0.0, 0.0, 0.0, 1.0);
}

/* scale
 *
 * Multiply the active matrix by a scaling matrix.
 * This implements: https://www.opengl.org/sdk/docs/man2/xhtml/glScale.xml
 */
void canvashdl::scale(vec3f size)
{
	matrices[active_matrix] = matrices[active_matrix]*mat4f(size[0], 0.0, 0.0, 0.0,
									 0.0, size[1], 0.0, 0.0,
									 0.0, 0.0, size[2], 0.0,
									 0.0, 0.0, 0.0, 1.0);
}

/* perspective
 *
 * Multiply the active matrix by a perspective projection matrix.
 * This implements: https://www.opengl.org/sdk/docs/man2/xhtml/gluPerspective.xml
 */
void canvashdl::perspective(float fovy, float aspect, float n, float f)
{
	float x = tan(m_pi/2.0 - fovy/2);
	matrices[active_matrix] = mat4f(x/aspect, 0.0, 0.0, 0.0,
									0.0, x, 0.0, 0.0,
									0.0, 0.0, -(f+n)/(f-n), -2.0*f*n/(f-n),
									0.0, 0.0, -1.0, 0.0)*matrices[active_matrix];
}

/* frustum
 *
 * Multiply the active matrix by a frustum projection matrix.
 * This implements: https://www.opengl.org/sdk/docs/man2/xhtml/glFrustum.xml
 */
void canvashdl::frustum(float l, float r, float b, float t, float n, float f)
{
	matrices[active_matrix] = mat4f(2.0*n/(r - l), 0.0, (r + l)/(r - l), 0.0,
								    0.0, 2.0*n/(t - b), (t + b)/(t - b), 0.0,
								    0.0, 0.0, -(f + n)/(f - n), -2.0*f*n/(f - n),
								    0.0, 0.0, -1.0, 0.0)*matrices[active_matrix];
}

/* ortho
 *
 * Multiply the active matrix by an orthographic projection matrix.
 * This implements: https://www.opengl.org/sdk/docs/man2/xhtml/glOrtho.xml
 */
void canvashdl::ortho(float l, float r, float b, float t, float n, float f)
{
	matrices[active_matrix] = mat4f(2.0/(r - l), 0.0, 0.0, -(r + l)/(r - l),
						   0.0, 2.0/(t - b), 0.0, -(t + b)/(t - b),
						   0.0, 0.0, -2.0/(f - n), -(f + n)/(f - n),
						   0.0, 0.0, 0.0, 1.0)*matrices[active_matrix];
}

void canvashdl::viewport(int left, int bottom, int right, int top)
{
	matrices[viewport_matrix] = mat4f((float)(right - left)/2.0, 0.0, 0.0, (float)(right + left)/2.0,
									  0.0, (float)(top - bottom)/2.0, 0.0, (float)(top + bottom)/2.0,
									  0.0, 0.0, (float)depth/2.0, (float)depth/2.0,
									  0.0, 0.0, 0.0, 1.0);

	resize(right - left, top - bottom);
}

/* look_at
 *
 * Move and orient the modelview so the camera is at the 'at' position focused on the 'eye'
 * position and rotated so the 'up' vector is up
 * This implements: https://www.opengl.org/sdk/docs/man2/xhtml/gluLookAt.xml
 */
void canvashdl::look_at(vec3f eye, vec3f at, vec3f up)
{
	vec3f f = norm(at - eye);
	up = norm(up);
	vec3f s = cross(f, up);
	vec3f u = cross(norm(s), f);
	matrices[active_matrix] = mat4f(s[0], s[1], s[2], 0.0,
									u[0], u[1], u[2], 0.0,
									-f[0], -f[1], -f[2], 0.0,
									0.0, 0.0, 0.0, 1.0)*matrices[active_matrix];
	translate(-eye);
}

void canvashdl::update_normal_matrix()
{
	matrices[normal_matrix] = transpose(inverse(matrices[modelview_matrix]));
}

/* to_window
 *
 * Given a pixel coordinate (x from 0 to width and y from 0 to height),
 * convert it into window coordinates (x from -1 to 1 and y from -1 to 1).
 */
vec3f canvashdl::to_window(vec2i pixel)
{
	vec3f result(2.0*(float)pixel[0]/(float)(width-1) - 1.0,
				 2.0*(float)(height - 1 - pixel[1])/(float)(height-1) - 1.0,
				 1.04); // <-- I don't actually understand why I need this magic number... Seems strange to me... Shouldn't it be 1.0?

	return result;
}

/* unproject
 *
 * Unproject a window coordinate into world coordinates.
 * This implements: https://www.opengl.org/sdk/docs/man2/xhtml/gluUnProject.xml
 */
vec3f canvashdl::unproject(vec3f window)
{
	return inverse(matrices[modelview_matrix])*inverse(matrices[projection_matrix])*homogenize(window);
}

/* shade_vertex
 *
 * This is the vertex shader.
 * v[0] to v[2] is position
 * v[3] to v[5] is normal
 * v[7] to v[8] is texture coordinates
 * The result from this function is interpolated and passed on to the fragment shader
 * (its also used to draw the geometry during rasterization)
 * Note that the only requirements for the returned values are that the first 3 components
 * be a projected position. The rest are yours to decide what to do with.
 */
vec3f canvashdl::shade_vertex(vec8f v, vector<float> &varying)
{
    
    const materialhdl *m;
    get_uniform("material", m);


    if (m != NULL)
    {
        vec3f x = m->shade_vertex(this, v(0,3), v(3, 6), varying);
        return x;
    }

    /* There might be null material, so we should pass default */
    /* only construct whitehdl if we really need to */
    whitehdl w;
    return w.shade_vertex(this, v(0,3), v(3,6), varying);
}

/* shade_fragment
 *
 * This is the fragment shader. The pixel color is determined here.
 * the values for v are the interpolated result of whatever you returned from the vertex shader
 */
vec3f canvashdl::shade_fragment(vector<float> varying)
{
    //same as above but with shade fragment
    const materialhdl *m;
    get_uniform("material", m);


    if (m != NULL)
        return m->shade_fragment(this, varying);

    /* There might be null material, so we should pass default */
    /* only construct whitehdl if we really need to */
    whitehdl w;
    return w.shade_fragment(this, varying);

}

/* Translates from fp to int color
 * representations.
 */
inline vec3i translate_color(vec3f p)
{
    return vec3i(p)*255; 
}

/* plot
 *
 * Plot a pixel and check it against the depth buffer.
 */
void canvashdl::plot(vec3i xyz, vector<float> varying)
{
	int idx = xyz[1]*width + xyz[0];
	if (xyz[0] >= 0 && xyz[0] < width && xyz[1] >= 0 && xyz[1] < height && xyz[2] < depth_buffer[idx])
	{
		depth_buffer[idx] = xyz[2];
		vec3f color = shade_fragment(varying);
		color_buffer[idx*3 + 0] = (int)(color[0]*255.0);
		color_buffer[idx*3 + 1] = (int)(color[1]*255.0);
		color_buffer[idx*3 + 2] = (int)(color[2]*255.0);
	}
}

/* plot_point
 *
 * Plot a point given in window coordinates.
 */
void canvashdl::plot_point(vec3f v, vector<float> varying)
{
	plot((vec3i)v, varying);
}

/* plot_line
 *
 * Plot a line defined by two points in window coordinates. Use Bresenham's
 * Algorithm for this. Don't forget to interpolate the normals and texture
 * coordinates as well.
 */
void canvashdl::plot_line(vec3f v1, vector<float> v1_varying, vec3f v2, vector<float> v2_varying)
{
	vec3i s1 = (vec3i)v1;
	vec3i s2 = (vec3i)v2;

	int b = (abs(s2[1] - s1[1]) > abs(s2[0] - s1[0]));

	vec2i dv = (vec2i)s2 - (vec2i)s1;

	vec2i step((int)(dv[0] > 0) - (int)(dv[0] < 0),
			   (int)(dv[1] > 0) - (int)(dv[1] < 0));
	dv[0] *= step[0];
	dv[1] *= step[1];

	int D = 2*dv[1-b] - dv[b];

    std::vector<float> varying(v1_varying.size());
    if (shade_model != flat)
    {
        varying = v1_varying;
    }
    else
    {
        for (int i = 0; i < v1_varying.size(); ++i)
        {
            varying[i] = (v1_varying[i] + v2_varying[i])/2.0;
        }
    }
	plot(s1, varying);

	vec3i p = s1;
	float increment = (float)step[b]/(float)(s2[b] - s1[b]);
	float interpolate = 0.0f;
	float zdiff = v2[2] - v1[2];
	p[b]+=step[b];
	while (step[b]*p[b] < step[b]*s2[b])
	{
		if (D > 0)
		{
			p[1-b] += step[1-b];
			D += 2*(dv[1-b] - dv[b]);
		}
		else
			D += 2*dv[1-b];

		p[2] = (int)(zdiff*interpolate) + s1[2];

        if (shade_model != flat)
        {
            for (int i = 0; i < varying.size(); ++i)
            {
                float inc = v2_varying[i] - v1_varying[i];
                varying[i] = v1_varying[i] + inc * interpolate;
            }

            plot(p, varying);
        }
        else
        {
            plot(p, varying);
        }
     
		p[b] += step[b];
		interpolate += increment;
	}
}

/* plot_half_triangle
 *
 * Plot half of a triangle defined by three points in window coordinates (v1, v2, v3).
 * The remaining inputs are as follows (s1, s2, s3) are the pixel coordinates of (v1, v2, v3),
 * and (ave) is the average value of the normal and texture coordinates for flat shading.
 * Use Bresenham's algorithm for this. You may plot the horizontal half or the vertical half.
 */
void canvashdl::plot_half_triangle(vec3i s1, vector<float> v1_varying, vec3i s2, vector<float> v2_varying, vec3i s3, vector<float> v3_varying, vector<float> ave_varying)
{
	int b12 = (abs(s2[1] - s1[1]) > abs(s2[0] - s1[0]));
	int b13 = (abs(s1[1] - s3[1]) > abs(s1[0] - s3[0]));

	vec2i dv12 = s2 - s1;
	vec2i dv13 = s3 - s1;

	vec2i step12, step13;
	for (int i = 0; i < 2; i++)
	{
		step12[i] = (int)(dv12[i] > 0) - (int)(dv12[i] < 0);
		dv12[i] *= step12[i];
		step13[i] = (int)(dv13[i] > 0) - (int)(dv13[i] < 0);
		dv13[i] *= step13[i];
	}

	int D12 = 2*dv12[1-b12] - dv12[b12];
	int D13 = 2*dv13[1-b13] - dv13[b13];

    plot(s1, v1_varying);

	vec3i p12 = s1, p13 = s1;
	float increment12 = (float)step12[b12]/(float)(s2[b12] - s1[b12]);
	float increment13 = (float)step13[b13]/(float)(s3[b13] - s1[b13]);
	float interpolate12 = 0.0f;
	float interpolate13 = 0.0f;
	while (step12[b12] != 0 && step13[b13] != 0 && step12[b12]*p12[b12] < step12[b12]*s2[b12] && step13[b13]*p13[b13] < step13[b13]*s3[b13])
	{
		vec2i old = p12;
		do
		{
			if (D12 > 0)
			{
				p12[1-b12] += step12[1-b12];
				D12 += 2*(dv12[1-b12] - dv12[b12]);
			}
			else
				D12 += 2*dv12[1-b12];

			p12[b12] += step12[b12];
			interpolate12 += increment12;
		} while ((step12[b12] != 0 || step12[1-b12] != 0) && (vec2i)p12 == old && step12[b12]*p12[b12] < step12[b12]*s2[b12]);

		while ((step13[b13] != 0 || step13[1-b13] != 0) && p13[0] != p12[0] && step13[b13]*p13[b13] < step13[b13]*s3[b13])
		{
			p13[b13] += step13[b13];
			interpolate13 += increment13;

			if (D13 > 0)
			{
				p13[1-b13] += step13[1-b13];
				D13 += 2*(dv13[1-b13] - dv13[b13]);
			}
			else
				D13 += 2*dv13[1-b13];
		}

		p12[2] = (int)((float)(s2[2] - s1[2])*interpolate12 + (float)s1[2]);
		p13[2] = (int)((float)(s3[2] - s1[2])*interpolate13 + (float)s1[2]);

        vector<float> v2(v1_varying);
        vector<float> v3(v1_varying);

		if (shade_model != flat)
		{
			for (int i = 0; i < v1_varying.size(); ++i)
			{
				v2[i] += (v2_varying[i] - v1_varying[i])*interpolate12;
				v3[i] += (v3_varying[i] - v1_varying[i])*interpolate13;
			}
		}

		vec3i Ai = p12;
		vec3i Bi = p13;

		if (Ai[1] > Bi[1])
		{
			swap(Ai, Bi);
			swap(v2, v3);
		}

		vec3i p = Ai;
		float increment = 1.0f/(float)(Bi[1] - Ai[1]);
		float interpolateab = 0.0f;

		float zdiff = (float)(Bi[2] - Ai[2]);

		for (; p[1] < Bi[1]; p[1]++)
		{
			p[2] = (int)(zdiff*interpolateab) + Ai[2];

			if (shade_model != flat)
            {
				for (int i = 0; i < v1_varying.size(); ++i)
                {
					ave_varying[i] = v2[i] + (v3[i] - v2[i])*interpolateab;
                }
            }

			plot(p, ave_varying);
			interpolateab += increment;
		}
		plot(Bi, v3);
	}
}

/* plot_triangle
 *
 * Use the above functions to plot a whole triangle. Don't forget to
 * take into account the polygon mode. You should be able to render the
 * triangle as 3 points, 3 lines, or a filled in triangle. (v1, v2, v3)
 * are given in window coordinates.
 */
void canvashdl::plot_triangle(vec3f v1, vector<float> v1_varying, vec3f v2, vector<float> v2_varying, vec3f v3, vector<float> v3_varying)
{
	if (polygon_mode == point)
	{
		plot_point(v1, v1_varying);
		plot_point(v2, v2_varying);
		plot_point(v3, v3_varying);
	}
	else if (polygon_mode == line)
	{
		plot_line(v1, v1_varying, v2, v2_varying);
		plot_line(v2, v2_varying, v3, v3_varying);
		plot_line(v3, v3_varying, v1, v1_varying);
	}
	else if (polygon_mode == fill)
	{
		plot_line(v1, v1_varying, v2, v2_varying);
		plot_line(v2, v2_varying, v3, v3_varying);
		plot_line(v3, v3_varying, v1, v1_varying);

		if (v1[0] > v2[0])
		{
			swap(v1, v2);
			swap(v1_varying, v2_varying);
		}
		if (v1[0] > v3[0])
		{
			swap(v1, v3);
			swap(v1_varying, v3_varying);
		}
		if (v2[0] > v3[0])
		{
			swap(v2, v3);
			swap(v2_varying, v3_varying);
		}

		vector<float> ave_varying(v1_varying.size());

		for (int i = 0; i < (int)v1_varying.size(); i++)
			ave_varying[i] = (v1_varying[i] + v2_varying[i] + v3_varying[i])/3.0f;

		plot_half_triangle((vec3i)v1, v1_varying, (vec3i)v2, v2_varying, (vec3i)v3, v3_varying, ave_varying);
		plot_half_triangle((vec3i)v3, v3_varying, (vec3i)v1, v1_varying, (vec3i)v2, v2_varying, ave_varying);
	}
}

/* draw_points
 *
 * Draw a set of 3D points on the canvas. Each point in geometry is
 * formatted (vx, vy, vz, nx, ny, nz, s, t). Don't forget to test the
 * points against the clipping plains of the projection. If you don't
 * you'll get weird behavior (especially when objects behind the camera
 * are being rendered).
 */
void canvashdl::draw_points(const vector<vec8f> &geometry)
{
	update_normal_matrix();
	mat4f transform = matrices[projection_matrix]*matrices[modelview_matrix];
	vec4f planes[6];
	for (int i = 0; i < 6; i++)
		planes[i] = transform.row(3) + (float)pow(-1.0, i)*(vec4f)transform.row(i/2);

	vector<pair<vec3f, vector<float> > > processed_geometry;
	processed_geometry.reserve(geometry.size());

	for (int i = 0; i < geometry.size(); i += 3)
	{
		bool keep = true;
		for (int j = 0; j < 6 && keep; j++)
			if (dot(homogenize((vec3f)geometry[i]), planes[j]) <= 0)
				keep = false;

		if (keep)
		{
			vector<float> varying;
			vec3f position = matrices[viewport_matrix]*homogenize(shade_vertex(geometry[i], varying));
			processed_geometry.push_back(pair<vec3f, vector<float> >(position, varying));
		}
	}

	for (int i = 0; i < processed_geometry.size(); i++)
		plot_point(processed_geometry[i].first, processed_geometry[i].second);
}

/* draw_lines
 *
 * Draw a set of 3D lines on the canvas. Each point in geometry
 * is formatted (vx, vy, vz, nx, ny, nz, s, t). Don't forget to clip
 * the lines against the clipping planes of the projection. You can't
 * just not render them because you'll get some weird popping at the
 * edge of the view.
 */
void canvashdl::draw_lines(const vector<vec8f> &geometry, const vector<int> &indices)
{
	update_normal_matrix();
	mat4f transform = matrices[projection_matrix]*matrices[modelview_matrix];
	vec4f planes[6];
	for (int i = 0; i < 6; i++)
		planes[i] = transform.row(3) + (float)pow(-1.0, i)*(vec4f)transform.row(i/2);

	vector<pair<vec3f, vector<float> > > processed_geometry;
	vector<int> processed_indices;
	processed_geometry.reserve(geometry.size());
	processed_indices.reserve(indices.size());
	vector<int> index_map;
	index_map.resize(geometry.size(), -1);

	for (int i = 0; i < indices.size(); i += 2)
	{
		pair<vec8f, int> x0 = pair<vec8f, int>(geometry[indices[i]], indices[i]);
		pair<vec8f, int> x1 = pair<vec8f, int>(geometry[indices[i+1]], indices[i+1]);
		bool keep = true;

		for (int j = 0; j < 6 && keep; j++)
		{
			float d0 = dot(homogenize(x0.first), planes[j]);
			float d1 = dot(homogenize(x1.first), planes[j]);
			float del = dot(homogenize(x1.first) - homogenize(x0.first), planes[j]);
			float p = -d0/del;

			if (d0 > 0.0 && d1 <= 0.0)
			{
				x1.first = (1-p)*x0.first + p*x1.first;
				x1.second = -1;
			}
			else if (d0 <= 0.0 && d1 > 0.0)
			{
				x0.first = (1-p)*x0.first + p*x1.first;
				x0.second = -1;
			}
			else if (d0 <= 0.0 && d1 <= 0.0)
				keep = false;
		}

		if (keep)
		{
			vector<float> varying;
			vec3f position;
			if (x0.second == -1)
			{
				x0.second = processed_geometry.size();
				position = matrices[viewport_matrix]*homogenize(shade_vertex(x0.first, varying));
				processed_geometry.push_back(pair<vec3f, vector<float> >(position, varying));
			}
			else if (index_map[x0.second] == -1)
			{
				index_map[x0.second] = processed_geometry.size();
				x0.second = processed_geometry.size();
				position = matrices[viewport_matrix]*homogenize(shade_vertex(x0.first, varying));
				processed_geometry.push_back(pair<vec3f, vector<float> >(position, varying));
			}
			else
				x0.second = index_map[x0.second];

			if (x1.second == -1)
			{
				x1.second = processed_geometry.size();
				position = matrices[viewport_matrix]*homogenize(shade_vertex(x1.first, varying));
				processed_geometry.push_back(pair<vec3f, vector<float> >(position, varying));
			}
			else if (index_map[x1.second] == -1)
			{
				index_map[x1.second] = processed_geometry.size();
				x1.second = processed_geometry.size();
				position = matrices[viewport_matrix]*homogenize(shade_vertex(x1.first, varying));
				processed_geometry.push_back(pair<vec3f, vector<float> >(position, varying));
			}
			else
				x1.second = index_map[x1.second];

			processed_indices.push_back(x0.second);
			processed_indices.push_back(x1.second);
		}
	}

	for (int i = 1; i < processed_indices.size(); i+=2)
		plot_line(processed_geometry[processed_indices[i-1]].first, processed_geometry[processed_indices[i-1]].second,
				  processed_geometry[processed_indices[i]].first, processed_geometry[processed_indices[i]].second);
}

/* draw_triangles
 *
 * Draw a set of 3D triangles on the canvas. Each point in geometry is
 * formatted (vx, vy, vz, nx, ny, nz, s, t). Don't forget to clip the
 * triangles against the clipping planes of the projection. You can't
 * just not render them because you'll get some weird popping at the
 * edge of the view. Also, this is where font/back face culling is implemented.
 */
void canvashdl::draw_triangles(const vector<vec8f> &geometry, const vector<int> &indices)
{
	update_normal_matrix();
	mat4f transform = matrices[projection_matrix]*matrices[modelview_matrix];
	vec4f planes[6];
	for (int i = 0; i < 6; i++)
		planes[i] = transform.row(3) + (float)pow(-1.0, i)*(vec4f)transform.row(i/2);
	vec3f eye = matrices[modelview_matrix].col(3)(0,3);

	vector<pair<vec3f, vector<float> > > processed_geometry;
	vector<int> processed_indices;
	processed_geometry.reserve(geometry.size());
	processed_indices.reserve(indices.size());
	vector<int> index_map;
	index_map.resize(geometry.size(), -1);

	for (int i = 0; i < indices.size(); i += 3)
	{
		vector<pair<vec8f, int> > polygon;
		for (int j = 0; j < 3; j++)
			polygon.push_back(pair<vec8f, int>(geometry[indices[i+j]], indices[i+j]));

		vector<pair<vec8f, int> > clipped;
		for (int j = 0; j < 6; j++)
		{
			pair<vec8f, int> x0 = polygon[polygon.size()-1];
			float d0 = dot(homogenize(x0.first), planes[j]);
			for (int k = 0; k < polygon.size(); k++)
			{
				pair<vec8f, int> x1 = polygon[k];
				float d1 = dot(homogenize(x1.first), planes[j]);
				float del = dot(homogenize(x1.first) - homogenize(x0.first), planes[j]);
				float p = -d0/del;

				if (d0 >= 0.0 && d1 >= 0.0)
					clipped.push_back(x1);
				else if (d0 >= 0.0 && d1 < 0.0)
					clipped.push_back(pair<vec8f, int>((1-(p+0.001f))*x0.first + (p+0.001f)*x1.first, -1));
				else if (d0 < 0.0 && d1 >= 0.0)
				{
					clipped.push_back(pair<vec8f, int>((1-(p-0.001f))*x0.first + (p-0.001f)*x1.first, -1));
					clipped.push_back(x1);
				}

				x0 = x1;
				d0 = d1;
			}

			polygon = clipped;
			clipped.clear();
		}

		if (polygon.size() > 2)
		{
			for (int i = 0; i < polygon.size(); i++)
			{
				vector<float> varying;
				vec3f position;
				if (polygon[i].second == -1)
				{
					polygon[i].second = processed_geometry.size();
					position = matrices[viewport_matrix]*homogenize(shade_vertex(polygon[i].first, varying));
					polygon[i].first = position;
					processed_geometry.push_back(pair<vec3f, vector<float> >(position, varying));
				}
				else if (index_map[polygon[i].second] == -1)
				{
					index_map[polygon[i].second] = processed_geometry.size();
					polygon[i].second = processed_geometry.size();
					position = matrices[viewport_matrix]*homogenize(shade_vertex(polygon[i].first, varying));
					polygon[i].first = position;
					processed_geometry.push_back(pair<vec3f, vector<float> >(position, varying));
				}
				else
				{
					polygon[i].second = index_map[polygon[i].second];
					polygon[i].first = processed_geometry[polygon[i].second].first;
				}
			}

			for (int i = 2; i < polygon.size(); i++)
			{
				vec3f normal = cross(norm((vec3f)polygon[0].first - (vec3f)polygon[i-1].first),
								     norm((vec3f)polygon[i].first - (vec3f)polygon[i-1].first));

				if (culling == disable || (normal[2] >= 0.0 && culling == backface) || (normal[2] <= 0.0 && culling == frontface))
				{
					processed_indices.push_back(polygon[0].second);
					processed_indices.push_back(polygon[i-1].second);
					processed_indices.push_back(polygon[i].second);
				}
			}
		}
	}

	for (int i = 2; i < processed_indices.size(); i+=3)
		plot_triangle(processed_geometry[processed_indices[i-2]].first, processed_geometry[processed_indices[i-2]].second,
					  processed_geometry[processed_indices[i-1]].first, processed_geometry[processed_indices[i-1]].second,
					  processed_geometry[processed_indices[i]].first, processed_geometry[processed_indices[i]].second);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Do not edit anything below here, that code just sets up OpenGL to render a single
 * quad that covers the whole screen, applies the color_buffer as a texture to it, and
 * keeps the buffer size and texture up to date.
 */
void canvashdl::load_texture()
{
	glGenTextures(1, &screen_texture);
	check_error(__FILE__, __LINE__);
	glActiveTexture(GL_TEXTURE0);
	check_error(__FILE__, __LINE__);
	glBindTexture(GL_TEXTURE_2D, screen_texture);
	check_error(__FILE__, __LINE__);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	check_error(__FILE__, __LINE__);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	check_error(__FILE__, __LINE__);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	check_error(__FILE__, __LINE__);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	check_error(__FILE__, __LINE__);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	check_error(__FILE__, __LINE__);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, color_buffer);
	check_error(__FILE__, __LINE__);
}

void canvashdl::load_geometry()
{
	// x, y, s, t
	const GLfloat geometry[] = {
			-1.0, -1.0, 0.0, 0.0,
			 1.0, -1.0, 1.0, 0.0,
			-1.0,  1.0, 0.0, 1.0,
			-1.0,  1.0, 0.0, 1.0,
			 1.0, -1.0, 1.0, 0.0,
			 1.0,  1.0, 1.0, 1.0
	};

	glGenBuffers(1, &screen_geometry);
	glBindBuffer(GL_ARRAY_BUFFER, screen_geometry);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*4*6, NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(GLfloat)*4*6, geometry);
}

void canvashdl::load_shader()
{
	GLuint vertex = load_shader_file(working_directory + "res/canvas.vx", GL_VERTEX_SHADER);
	GLuint fragment = load_shader_file(working_directory + "res/canvas.ft", GL_FRAGMENT_SHADER);

	screen_shader = glCreateProgram();
	glAttachShader(screen_shader, vertex);
	glAttachShader(screen_shader, fragment);
	glLinkProgram(screen_shader);
}

void canvashdl::init_opengl()
{
	glEnable(GL_TEXTURE_2D);
	glViewport(0, 0, width, height);

	load_texture();
	load_geometry();
	load_shader();
	initialized = true;
}

void canvashdl::check_error(const char *file, int line)
{
	GLenum error_code;
	error_code = glGetError();
	if (error_code != GL_NO_ERROR)
		cerr << "error: " << file << ":" << line << ": " << gluErrorString(error_code) << endl;
}

double canvashdl::get_time()
{
	timeval gtime;
	gettimeofday(&gtime, NULL);
	return gtime.tv_sec + gtime.tv_usec*1.0E-6;
}

void canvashdl::resize(int w, int h)
{
	glViewport(0, 0, w, h);
	last_reshape_time = get_time();
	reshape_width = w;
	reshape_height = h;
}

void canvashdl::swap_buffers()
{
	if (!initialized)
		init_opengl();

	if (last_reshape_time > 0.0 && get_time() - last_reshape_time > 0.125)
		resize(reshape_width, reshape_height);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(screen_shader);
	check_error(__FILE__, __LINE__);

	glActiveTexture(GL_TEXTURE0);
	check_error(__FILE__, __LINE__);
	glBindTexture(GL_TEXTURE_2D, screen_texture);
	check_error(__FILE__, __LINE__);
	glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, color_buffer);
	check_error(__FILE__, __LINE__);
	glUniform1i(glGetUniformLocation(screen_shader, "tex"), 0);
	check_error(__FILE__, __LINE__);

	glBindBuffer(GL_ARRAY_BUFFER, screen_geometry);
	check_error(__FILE__, __LINE__);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	check_error(__FILE__, __LINE__);
	glEnableClientState(GL_VERTEX_ARRAY);
	check_error(__FILE__, __LINE__);

	glTexCoordPointer(2, GL_FLOAT, sizeof(GLfloat)*4, (float*)(sizeof(GLfloat)*2));
	check_error(__FILE__, __LINE__);
	glVertexPointer(2, GL_FLOAT, sizeof(GLfloat)*4, NULL);
	check_error(__FILE__, __LINE__);

	glDrawArrays(GL_TRIANGLES, 0, 6);
	check_error(__FILE__, __LINE__);

	glDisableClientState(GL_VERTEX_ARRAY);
	check_error(__FILE__, __LINE__);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	check_error(__FILE__, __LINE__);

	glutSwapBuffers();
	check_error(__FILE__, __LINE__);
}

int canvashdl::get_width()
{
	return width;
}

int canvashdl::get_height()
{
	return height;
}
