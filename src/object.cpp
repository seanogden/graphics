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
    canvas->set_matrix(canvashdl::modelview_matrix);
    mat4f snapshot = canvas->matrices[canvas->active_matrix];
    canvas->translate(position);

    canvas->rotate(orientation[0], vec3f(1, 0, 0));
    canvas->rotate(orientation[1], vec3f(0, 1, 0));
    canvas->rotate(orientation[2], vec3f(0, 0, 1));

    canvas->scale(vec3f(scale, scale, scale));

    std::vector<vec8f> v;
    std::vector<int> i;

    int j = 0;

    if (face)
    {
        //for face normals, iterate over indices so we can look at each
        //triangle as a whole.
        for (std::vector<rigidhdl>::iterator rit = rigid.begin();
                rit != rigid.end(); ++rit)
        {
            for (std::vector<int>::iterator it = rit->indices.begin();
                    it != rit->indices.end(); ++it)
            {
                vec3f p1 = rit->geometry[*(it++)];
                vec3f p2 = rit->geometry[*(it++)];
                vec3f p3 = rit->geometry[*it];
                vec3f norm = cross(p2 - p1, p3 - p1);

                //check if facing the right direction.
                float s = dot(norm, vec3f(p1[3],p1[4],p1[5]));
                if (s < 0)
                {
                    norm = -norm;
                }

                //Make the normal come from the centroid of the triangle.
                vec3f centroid = (p1 + p2 + p3)*float(1.0/3.0);
                v.push_back(centroid);
                v.push_back(centroid + norm);
                i.push_back(j++);
                i.push_back(j++);

            }
        }
    }
    else
    {
        //For vertex normals, iterate over all vertices.
        for (std::vector<rigidhdl>::iterator rit = rigid.begin();
                rit != rigid.end(); ++rit)
        {
            for (std::vector<vec8f>::iterator it = rit->geometry.begin();
                    it != rit->geometry.end(); ++it)
            {
                //make normal lines from vertex (p) in direction p[3:5]
                //  this requires adding p[0:2] to p[3:5] to start at the vertex and not
                //  the origin.
                vec8f p = *it;
                vec8f p2 = p + vec3f(p[3], p[4], p[5])/float(10.0);
                v.push_back(p);
                v.push_back(p2);
                i.push_back(j++);
                i.push_back(j++);
            }
        }
    }

    
    canvas->draw_lines(v, i);
    canvas->matrices[canvas->active_matrix] = snapshot;


	// TODO Assignment 3: clear the material in the uniform list before rendering
}
