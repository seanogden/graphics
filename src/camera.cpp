/*
 * camera.cpp
 *
 *  Created on: Dec 3, 2014
 *      Author: nbingham
 */

#include "camera.h"
#include "object.h"

camerahdl::camerahdl()
{
	position = vec3f(0.0, 0.0, 0.0);
	orientation = vec3f(0.0, 0.0, 0.0);
	model = NULL;
	type = "camera";
	focus = NULL;
	radius = 10.0f;
}

camerahdl::~camerahdl()
{

}

void camerahdl::view(canvashdl *canvas)
{

        canvas->set_matrix(canvashdl::modelview_matrix);
        canvas->load_identity();

//void canvashdl::look_at(vec3f eye, vec3f at, vec3f up)
    if (focus != NULL)
    {
        vec3f up = ror3(vec3f(0.0, 1.0, 0.0), orientation);
        vec3f pos = ror3(vec3f(0.0, 0.0, -radius), orientation);
        //given radius, orientation, and focus->position to calculate
        //camera position relative to focus->position, and feed to
        //look_at
        canvas->look_at(pos, focus->position, up);
    }
    else
    {
        canvas->rotate(-orientation[0], vec3f(1, 0, 0));
        canvas->rotate(-orientation[1], vec3f(0, 1, 0));
        canvas->rotate(-orientation[2], vec3f(0, 0, 1));
        canvas->translate(-position);
    }


}

orthohdl::orthohdl()
{
	left = -10.0;
	right = 10.0;
	bottom = -10.0;
	top = 10.0;
	front = 2.0;
	back = 101.0;
	type = "ortho";
}

orthohdl::~orthohdl()
{
}

void orthohdl::project(canvashdl *canvas)
{
    canvas->set_matrix(canvashdl::projection_matrix);
    canvas->load_identity();
    canvas->ortho(left, right, bottom, top, front, back);
}

frustumhdl::frustumhdl()
{
	left = -1.0;
	right = 1.0;
	bottom = -1.0;
	top = 1.0;
	front = 2.0;
	back = 101.0;
	type = "frustum";
}

frustumhdl::~frustumhdl()
{

}

void frustumhdl::project(canvashdl *canvas)
{
    canvas->set_matrix(canvashdl::projection_matrix);
    canvas->load_identity();
    canvas->frustum(left, right, bottom, top, front, back);
}

perspectivehdl::perspectivehdl()
{
	fovy = m_pi/4.0;
	aspect = 1.0;
	front = 2.0;
	back = 101.0;
	type = "perspective";
}

perspectivehdl::~perspectivehdl()
{

}

void perspectivehdl::project(canvashdl *canvas)
{
    canvas->set_matrix(canvashdl::projection_matrix);
    canvas->load_identity();
    canvas->perspective(fovy, aspect, front, back);
}
