/*
 * object.cpp
 *
 *  Created on: Jan 2, 2015
 *      Author: nbingham
 */

#include "object.h"
#include "canvas.h"

rigidhdl::rigidhdl()
{

}

rigidhdl::~rigidhdl()
{

}

/* draw
 *
 * Draw a rigid body.
 */
void rigidhdl::draw(canvashdl *canvas)
{
    canvas->draw_triangles(geometry, indices);
}

objecthdl::objecthdl()
{
	position = vec3f(0.0, 0.0, 0.0);
	orientation = vec3f(0.0, 0.0, 0.0);
	bound = vec6f(1.0e6, -1.0e6, 1.0e6, -1.0e6, 1.0e6, -1.0e6);
	scale = 1.0;
}

objecthdl::objecthdl(const objecthdl &o)
{
	position = o.position;
	orientation = o.orientation;
	bound = o.bound;
	scale = o.scale;
	rigid = o.rigid;
	for (map<string, materialhdl*>::const_iterator i = o.material.begin(); i != o.material.end(); i++)
		material.insert(pair<string, materialhdl*>(i->first, i->second->clone()));
}

objecthdl::~objecthdl()
{
	for (map<string, materialhdl*>::iterator i = material.begin(); i != material.end(); i++)
		if (i->second != NULL)
		{
			delete i->second;
			i->second = NULL;
		}

	material.clear();
}

/* draw
 *
 * Draw the model. Don't forget to apply the transformations necessary
 * for position, orientation, and scale.
 */
void objecthdl::draw(canvashdl *canvas)
{
    canvas->set_matrix(canvashdl::modelview_matrix);
    mat4f snapshot = canvas->matrices[canvas->active_matrix];
    canvas->translate(position);

    canvas->rotate(orientation[0], vec3f(1, 0, 0));
    canvas->rotate(orientation[1], vec3f(0, 1, 0));
    canvas->rotate(orientation[2], vec3f(0, 0, 1));

    canvas->scale(vec3f(scale,scale,scale));

    for (vector<rigidhdl>::iterator it = rigid.begin();
         it != rigid.end(); ++it)
    {
        it->draw(canvas);
    }
    
    canvas->matrices[canvas->active_matrix] = snapshot;

	// TODO Assignment 3: Pass the material as a uniform into the renderer
}

/* draw_bound
 *
 * Create a representation for the bounding box and
 * render it.
 */
void objecthdl::draw_bound(canvashdl *canvas)
{
    canvas->set_matrix(canvashdl::modelview_matrix);
    mat4f snapshot = canvas->matrices[canvas->active_matrix];
    canvas->translate(position);

    canvas->rotate(orientation[0], vec3f(1, 0, 0));
    canvas->rotate(orientation[1], vec3f(0, 1, 0));
    canvas->rotate(orientation[2], vec3f(0, 0, 1));

    canvas->scale(vec3f(scale, scale, scale));

    float l, r, b, t, n, f;
    l = bound[0];
    r = bound[1];
    b = bound[2];
    t = bound[3];
    n = bound[4];
    f = bound[5];

    std::vector<vec8f> v;
    v.push_back(vec8f(l, b, n, 0, 0, 0, 0, 0));
    v.push_back(vec8f(l, t, n, 0, 0, 0, 0, 0));
    v.push_back(vec8f(r, b, n, 0, 0, 0, 0, 0));
    v.push_back(vec8f(r, t, n, 0, 0, 0, 0, 0));

    v.push_back(vec8f(l, b, f, 0, 0, 0, 0, 0));
    v.push_back(vec8f(l, t, f, 0, 0, 0, 0, 0));
    v.push_back(vec8f(r, b, f, 0, 0, 0, 0, 0));
    v.push_back(vec8f(r, t, f, 0, 0, 0, 0, 0));

    std::vector<int> i;
    //front face
    i.push_back(0);
    i.push_back(1);

    i.push_back(2);
    i.push_back(3);

    i.push_back(1);
    i.push_back(3);
    
    i.push_back(0);
    i.push_back(2);
    
    //cross bars
    i.push_back(2);
    i.push_back(6);

    i.push_back(3);
    i.push_back(7);
    
    i.push_back(0);
    i.push_back(4);

    i.push_back(1);
    i.push_back(5);
    
    //back face
    i.push_back(4);
    i.push_back(5);

    i.push_back(6);
    i.push_back(7);

    i.push_back(5);
    i.push_back(7);
    
    i.push_back(4);
    i.push_back(6);
    
    canvas->draw_lines(v, i);
    canvas->matrices[canvas->active_matrix] = snapshot;

	// TODO Assignment 3: clear the material in the uniform list
}

/* draw_normals
 *
 * create a representation of the normals for this object.
 * If face is false, render the vertex normals. Otherwise,
 * calculate the normals for each face and render those.
 */
void objecthdl::draw_normals(canvashdl *canvas, bool face)
{
	/* TODO Assignment 1: Generate the geometry to display the normals and send the necessary
	 * transformations and geometry to the renderer
	 */

	// TODO Assignment 3: clear the material in the uniform list before rendering
}
