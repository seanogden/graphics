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
#include <algorithm>

#define ABS(X) ((X) < 0 ? -(X) : (X))

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

	polygon_mode = line;
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
    float c = cos(angle);
    float s = sin(angle);
    float x = axis[0];
    float y = axis[1];
    float z = axis[2];
    float x_2 = x*x;
    float y_2 = y*y;
    float z_2 = z*z;
    
    //Genereal rotation matrix is from the opengl glRotate documentation.
    mat4f R = mat4f(x_2*(1-c)+c  , x*y*(1-c)-z*s, x*z*(1-c)+y*s, 0.0,
                    y*x*(1-c)+z*s, y_2*(1-c)+c  , y*z*(1-c)-x*s, 0.0,
                    x*z*(1-c)-y*s, y*z*(1-c)+x*s, z_2*(1-c)+c  , 0.0,
                    0.0          , 0.0          , 0.0          , 1.0); 

    matrices[active_matrix] = matrices[active_matrix]*R;
}

/* translate
 *
 * Multiply the active matrix by a translation matrix.
 * This implements: https://www.opengl.org/sdk/docs/man2/xhtml/glTranslate.xml
 */
void canvashdl::translate(vec3f direction)
{
    mat4f T = mat4f(1.0, 0.0, 0.0, direction[0],
                    0.0, 1.0, 0.0, direction[1],
                    0.0, 0.0, 1.0, direction[2],
                    0.0, 0.0, 0.0, 1.0); 

    //Genereal translation matrix is from the opengl glTranslate documentation.
    matrices[active_matrix] = matrices[active_matrix]*T;
}

/* scale
 *
 * Multiply the active matrix by a scaling matrix.
 * This implements: https://www.opengl.org/sdk/docs/man2/xhtml/glScale.xml
 */
void canvashdl::scale(vec3f size)
{
    mat4f T = mat4f(size[0], 0.0    , 0.0    , 0.0,
                    0.0    , size[1], 0.0    , 0.0,
                    0.0    , 0.0    , size[2], 0.0, 
                    0.0    , 0.0    , 0.0    , 1.0); 

    //Genereal translation matrix is from the opengl glTranslate documentation.
    matrices[active_matrix] = matrices[active_matrix]*T;
}

/* perspective
 *
 * Multiply the active matrix by a perspective projection matrix.
 * This implements: https://www.opengl.org/sdk/docs/man2/xhtml/gluPerspective.xml
 */
void canvashdl::perspective(float fovy, float aspect, float n, float f)
{
    float A = (f+n)/(n-f);
    float B = 2*f*n/(n-f);
    float g = 1.0/tan(fovy/2);

    mat4f F = mat4f(g/aspect, 0.0, 0.0 , 0.0,
                    0.0     , g  , 0.0 , 0.0,
                    0.0     , 0.0, A   , B  , 
                    0.0     , 0.0, -1.0, 1.0); 

    //Genereal translation matrix is from the opengl glTranslate documentation.
    matrices[active_matrix] = matrices[active_matrix]*F;
}

/* frustum
 *
 * Multiply the active matrix by a frustum projection matrix.
 * This implements: https://www.opengl.org/sdk/docs/man2/xhtml/glFrustum.xml
 */
void canvashdl::frustum(float l, float r, float b, float t, float n, float f)
{
    float A = (r + l)/(r-l);
    float B = (t + b)/(t - b);
    float C = -(f + n)/(f - n);
    float D = -(2*f*n)/(f-n);

    mat4f F = mat4f(2*n/(r-l), 0.0      , A   , 0.0,
                    0.0      , 2*n/(t-b), B   , 0.0,
                    0.0      , 0.0      , C   , D  , 
                    0.0      , 0.0      , -1.0, 1.0); 

    //Genereal translation matrix is from the opengl glTranslate documentation.
    matrices[active_matrix] = matrices[active_matrix]*F;
}

/* ortho
 *
 * Multiply the active matrix by an orthographic projection matrix.
 * This implements: https://www.opengl.org/sdk/docs/man2/xhtml/glOrtho.xml
 */
void canvashdl::ortho(float l, float r, float b, float t, float n, float f)
{
    float tx = -(r + l)/(r - l);
    float ty = -(t + b)/(t - b);
    float tz = -(f + n)/(f - n);
    float A  = 2/(r - l);
    float B  = 2/(t - b);
    float C  = -2/(f - n);

    mat4f F = mat4f(A  , 0.0, 0.0, tx,
                    0.0, B  , 0.0, ty,
                    0.0, 0.0, C  , tz, 
                    0.0, 0.0, 0.0, 1.0); 

    matrices[active_matrix] = matrices[active_matrix]*F;
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
	// TODO Assignment 1: Emulate the functionality of gluLookAt
    vec3f f = norm(at - eye);
    vec3f UP = norm(up);
    vec3f s = cross(f, UP);
    vec3f u = cross(norm(s), f);

    std::cout << eye << std::endl;
    std::cout << at << std::endl;
    std::cout << up << std::endl;
    //std::cout << s << std::endl;
    //std::cout << u << std::endl;

    mat4f M = mat4f( s[0],  s[1],  s[2], 0.0,
                     u[0],  u[1],  u[2], 0.0,
                    -f[0], -f[1], -f[2], 0.0,
                      0.0,   0.0,   0.0, 1.0);

    matrices[active_matrix] = matrices[active_matrix] * M;
    translate(-eye);
}

void canvashdl::update_normal_matrix()
{
	// TODO Assignment 3: calculate the normal matrix
}

/* to_window
 *
 * Given a pixel coordinate (x from 0 to width and y from 0 to height),
 * convert it into window coordinates (x from -1 to 1 and y from -1 to 1).
 */
vec3f canvashdl::to_window(vec2i pixel)
{
    vec2f w = vec2f(2.0, 2.0) * vec2f(pixel) / vec2f(width, height) - vec2f(1.0, 1.0);
    return vec3f(w[0], -w[1], 1.0);
}

vec3i canvashdl::to_pixel(vec3f p)
{
    return vec3i(vec3f(width, height, 1.0) * (p + vec3f(1.0, 1.0, 0.0))/vec3f(2.0, 2.0, 1.0));
}


/* unproject
 *
 * Unproject a window coordinate into world coordinates.
 * This implements: https://www.opengl.org/sdk/docs/man2/xhtml/gluUnProject.xml
 */
vec3f canvashdl::unproject(vec3f window)
{
    return inverse(matrices[modelview_matrix])*
           (inverse(matrices[projection_matrix]))*homogenize(window);
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
    vec4f vertex = vec4f(v[0], v[1], v[2], 1.0);
    vertex = matrices[projection_matrix] * matrices[modelview_matrix] * vertex;
    vertex = vertex/vertex[3];

	/* TODO Assignment 3: Get the material from the list of uniform variables and
	 * call its vertex shader.
	 */
	return vec3f(vertex[0], vertex[1], vertex[2]);
}

/* shade_fragment
 *
 * This is the fragment shader. The pixel color is determined here.
 * the values for v are the interpolated result of whatever you returned from the vertex shader
 */
vec3f canvashdl::shade_fragment(vector<float> varying)
{
	return white; //white!

	/* TODO Assignment 3: Get the material from the list of uniform variables and
	 * call its fragment shader.
	 */
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
    vec3i color = translate_color(shade_fragment(varying));
    size_t offset = 3*(width*xyz[1] + xyz[0]);

    if (xyz[0] >= 0 && xyz[0] < width && xyz[1] >= 0 && xyz[1] < height)
    { 
        unsigned char* p = color_buffer + offset;
        p[0] = (unsigned char)color[0];
        p[1] = (unsigned char)color[1];
        p[2] = (unsigned char)color[2];

    }
	/* TODO Assignment 3: Compare the z value against the depth buffer and
	 * only render if its less. Then set the depth buffer.
	 */
}

/* plot_point
 *
 * Plot a point given in window coordinates.
 */
void canvashdl::plot_point(vec3f v, vector<float> varying)
{
    plot(to_pixel(v), varying);
}


void canvashdl::plot_line(vec3f v1, vector<float> v1_varying, vec3f v2, vector<float> v2_varying)
{
    //algorithm explained at http://tech-algorithm.com/articles/drawing-line-using-bresenham-algorithm/ 
    int x, y, w, h, dx_shortest, dy_shortest, dx_longest, dy_longest, longest, shortest, i, numerator;
    vec3i p1, p2;
    p1 = to_pixel(v1);
    p2 = to_pixel(v2);

    w = p2[0] - p1[0];
    h = p2[1] - p1[1];
    x = p1[0];
    y = p1[1];
    
    if (w < 0) 
    {
        //first point is right of second point, so x decreases
        dx_shortest = -1;
        dx_longest = -1;
    }
    else if (w > 0)
    {
        //first point is left of second point, so x increases
        dx_shortest = 1;
        dx_longest = 1;
    }
    else
    {
        //Width = 0, this is a vertical line with no change in x.
        dx_shortest = 0;
        dx_longest = 0;
    }

    if (h < 0) 
    {
        //first point is higher than second point, so y decreases.
        dy_shortest = -1;
    }
    else if (h > 0)
    {
        //first point is lower than second point, so y increases.
        dy_shortest = 1;
    }
    else
    {
        //horizontal line, no change in y
        dy_shortest = 0;
    }

    longest = ABS(w) ;
    shortest = ABS(h) ;

    if (longest < shortest) 
    {
        //line is steep (m > 0.5)
        std::swap(longest, shortest);

        //extra steps vertically for steep lines.
        if (h < 0) 
        {
            dy_longest = -1; 
        }
        else if (h > 0) 
        {
            dy_longest = 1;
        }

        //for steep lines, don't move horizontally until we've moved vertically a few times.
        dx_longest = 0; 
    }
    else
    {
        dy_longest = 0;
    }

    numerator = longest / 2;

    for (i=0; i<=longest; ++i) 
    {
        plot(vec3i(x,y,0), vector<float>());
        
        //add shortest to numerator each time we increment in the long direction.
        numerator += shortest ;

        if (numerator >= longest) {
            //reset numerator when we've moved enough steps in the long direction. 
            numerator -= longest ;
            x += dx_shortest ;
            y += dy_shortest ;
        } 
        else 
        {
            x += dx_longest ;
            y += dy_longest ;
        }
    }
}

void canvashdl::plot_horizontal(int x1, int x2, int z1, int z2, int y, vector<float> varying)
{
    int inc = 1;
    int width = ABS(x1 - x2);
    if (x2 < x1) inc = -1;

    for (int i = 0; i < width; ++i)
    {
        plot(vec3i(x1 + i, y, z1), varying);
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
	// TODO Assignment 2: Implement Bresenham's half triangle fill algorithm
    // 1.  Find highest y-coordinate point. Call this A.
    // 2.  Find middle y-coordinate point.  Call this B.
    // 3.  Find intersection of horizontal line between middle y-coordinate point and line between other 2 points.  Call this C
    // 4.  Draw lines from A to B and A to C simultaneously.
    //     -When you move y coordinates, interrupt line drawing to call plot_line between current coordinates on each line.
    //     -resume line drawing. 
	// TODO Assignment 3: Interpolate the varying values before passing them into plot.
    int x1, y1, w1, h1, numerator1;
    int x2, y2, w2, h2, numerator2;

    w1 = s2[0] - s1[0];
    h1 = s2[1] - s1[1];
    x1 = s1[0];
    y1 = s1[1];
 
    int dx_shortest1 = 0, dy_shortest1 = 0, dx_longest1 = 0, dy_longest1 = 0 ;
    if (w1<0) dx_shortest1 = -1 ; else if (w1>0) dx_shortest1 = 1 ;
    if (h1<0) dy_shortest1 = -1 ; else if (h1>0) dy_shortest1 = 1 ;
    if (w1<0) dx_longest1 = -1 ; else if (w1>0) dx_longest1 = 1 ;
    int longest1 = ABS(w1) ;
    int shortest1 = ABS(h1) ;

    if (!(longest1>shortest1)) {
        longest1 = ABS(h1) ;
        shortest1 = ABS(w1) ;
        if (h1<0) dy_longest1 = -1 ; else if (h1>0) dy_longest1 = 1 ;
        dx_longest1 = 0 ;            
    }

    numerator1 = longest1 / 2;

    w2 = s3[0] - s1[0];
    h2 = s3[1] - s1[1];
    x2 = s1[0];
    y2 = s1[1];

    int dx_shortest2 = 0, dy_shortest2 = 0, dx_longest2 = 0, dy_longest2 = 0 ;
    if (w2<0) dx_shortest2 = -1 ; else if (w2>0) dx_shortest2 = 1 ;
    if (h2<0) dy_shortest2 = -1 ; else if (h2>0) dy_shortest2 = 1 ;
    if (w2<0) dx_longest2 = -1 ; else if (w2>0) dx_longest2 = 1 ;
    int longest2 = ABS(w2) ;
    int shortest2 = ABS(h2) ;

    if (!(longest2>shortest2)) {
        longest2 = ABS(h2) ;
        shortest2 = ABS(w2) ;
        if (h2<0) dy_longest2 = -1 ; else if (h2>0) dy_longest2 = 1 ;
        dx_longest2 = 0 ;            
    }

    numerator2 = longest2 / 2;

    for(;;)
    {
        for(;;)
        {
            plot(vec3i(x1,y1,0), vector<float>());

            if (x1 == s2[0] && y1 == s2[1]) break;

            //add shortest to numerator each time we increment in the long direction.
            numerator1 += shortest1 ;

            if (numerator1 >= longest1) {
                //reset numerator when we've moved enough steps in the long direction. 
                numerator1 -= longest1 ;
                x1 += dx_shortest1 ;
                y1 += dy_shortest1 ;

                if (dy_shortest1 != 0) break;
            } 
            else 
            {
                x1 += dx_longest1 ;
                y1 += dy_longest1 ;

                if (dy_longest1 != 0) break;
            }
        }


        for(;;)
        {
            plot(vec3i(x2,y2,0), vector<float>());

            if (x2 == s3[0] && y2 == s3[1]) break;

            //add shortest to numerator each time we increment in the long direction.
            numerator2 += shortest2 ;

            if (numerator2 >= longest2) {
                //reset numerator when we've moved enough steps in the long direction. 
                numerator2 -= longest2 ;
                x2 += dx_shortest2 ;
                y2 += dy_shortest2 ;

                if (dy_shortest2 != 0) break;
            } 
            else 
            {
                x2 += dx_longest2 ;
                y2 += dy_longest2 ;

                if (dy_longest2 != 0) break;
            }
        }

        plot_horizontal(x1, x2, 0, 0, y2, vector<float>());

        if (x2 == s3[0] && y2 == s3[1] && x1 == s2[0] && y1 == s2[1]) break;

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
    if (polygon_mode == line)
    {
        plot_line(v1, v1_varying, v2, v2_varying);
        plot_line(v2, v2_varying, v3, v3_varying);
        plot_line(v3, v3_varying, v1, v1_varying);
    }
    else if (polygon_mode == fill)
    {
        //Sort the vertices in order of y coordinate
        if (v1[1] > v2[1])
            swap(v1,v2);
        if (v2[1] > v3[1])
            swap(v2,v3);
        if (v1[1] > v2[1])
            swap(v1,v2);

        //Find intersection of horizontal line from v2 to line between v1 and v3
        float x4 = v1[0] + (v2[1] - v1[1])/(v3[1] - v1[1]) * (v3[0] - v1[0]);
        vec3f v4(x4, v2[1], 0.0);

        //plot half triangle v1, v2, v4
        plot_half_triangle(to_pixel(v1), v1_varying, to_pixel(v2), v2_varying, to_pixel(v4), vector<float>(), vector<float>());
        //plot half triangle v3, v2, v4
        plot_half_triangle(to_pixel(v3), v3_varying, to_pixel(v2), v2_varying, to_pixel(v4), vector<float>(), vector<float>());
    }
    else
    {
        plot_point(v1, v1_varying);
        plot_point(v2, v2_varying);
        plot_point(v3, v3_varying);
    }

	// TODO Assignment 2: Calculate the average varying vector for flat shading and call plot_half_triangle as needed.
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

    for (std::vector<vec8f>::const_iterator it = geometry.begin(); it != geometry.end(); ++it)
    {
        std::vector<float> varying = std::vector<float>();
        vec8f p = shade_vertex(*it, varying);
        plot_point(p, varying);
    }

	// TODO Assignment 2: Implement frustum clipping and back-face culling
	// TODO Assignment 3: Update the normal matrix.
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

    for (std::vector<int>::const_iterator it = indices.begin(); it != indices.end(); ++it)
    {
        std::vector<float> v1 = std::vector<float>();
        vec3f p1 = shade_vertex(geometry[*it], v1);
        ++it;
        std::vector<float> v2 = std::vector<float>();
        vec3f p2 = shade_vertex(geometry[*it], v2);
        plot_line(p1, v1, p2, v2);
    }

	// TODO Assignment 2: Implement frustum clipping and back-face culling
	// TODO Assignment 3: Update the normal matrix.
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

    for (std::vector<int>::const_iterator it = indices.begin(); it != indices.end(); ++it)
    {
        std::vector<float> v1 = std::vector<float>();
        vec3f p1 = shade_vertex(geometry[*it], v1);
        ++it;
        std::vector<float> v2 = std::vector<float>();
        vec3f p2 = shade_vertex(geometry[*it], v2);
        ++it;
        std::vector<float> v3 = std::vector<float>();
        vec3f p3 = shade_vertex(geometry[*it], v3);
        plot_triangle(p1, v1, p2, v2, p3, v3);
    }

	// TODO Assignment 2: Implement frustum clipping and back-face culling
	// TODO Assignment 3: Update the normal matrix.
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
