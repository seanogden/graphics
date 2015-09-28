/*
 * scene.cpp
 *
 *  Created on: Dec 3, 2014
 *      Author: nbingham
 */

#include "scene.h"
#include "camera.h"
#include "object.h"
#include "light.h"

#include "primitive.h"
#include "model.h"

scenehdl::scenehdl()
{
	canvas = NULL;
	active_camera = -1;
	active_object = -1;
	render_normals = none;
	render_lights = false;
	render_cameras = false;
}

scenehdl::~scenehdl()
{

}


bool scenehdl::is_camera(objecthdl* obj)
{
    for (vector<camerahdl*>::iterator camit = cameras.begin();
            camit != cameras.end(); ++camit)
    {
        if (*camit != NULL && (*camit)->model == obj)
        {
            return true;
        }
    }

    return false;
}

bool scenehdl::is_active_camera(objecthdl* obj)
{
    return cameras[active_camera]->model == obj;
}


/* draw
 *
 * Update the locations of all of the lights, draw all of the objects, and
 * if enabled, draw the normals, the lights, the cameras, etc.
 */
void scenehdl::draw()
{
    cameras[active_camera]->project(canvas); 
    cameras[active_camera]->view(canvas); 

	for (vector<objecthdl*>::iterator it = objects.begin();
         it != objects.end(); ++it)
    {

        if (!render_cameras && is_camera(*it))
        {
            continue;
        }

        if (is_active_camera(*it))
        {
            continue;
        }

        (*it)->draw(canvas);

        if (render_normals)
        {
            (*it)->draw_normals(canvas);
        }
    }

    if (active_object_valid())
    {
        objects[active_object]->draw_bound(canvas);
    } 
   
	/* TODO Assignment 3: Clear the uniform variables and pass the vector of
	 * lights into the renderer as a uniform variable.
	 * TODO Assignment 3: Update the light positions and directions
	 * TODO Assignment 3: Render the lights
	 */
}

bool scenehdl::active_camera_valid()
{
	return (active_camera >= 0 && active_camera < cameras.size() && cameras[active_camera] != NULL);
}

bool scenehdl::active_object_valid()
{
	return (active_object >= 0 && active_object < objects.size() && objects[active_object] != NULL);
}
