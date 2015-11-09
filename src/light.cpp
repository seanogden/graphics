/*
 * light.cpp
 *
 *  Created on: Dec 17, 2014
 *      Author: nbingham
 */

#include "light.h"
#include "object.h"
#include "canvas.h"

#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))

lighthdl::lighthdl()
{
	model = NULL;
	type = "light";
}

lighthdl::lighthdl(const vec3f &ambient, const vec3f &diffuse, const vec3f &specular)
{
	this->ambient = ambient;
	this->diffuse = diffuse;
	this->specular = specular;
	model = NULL;
	type = "light";
}

lighthdl::~lighthdl()
{

}

directionalhdl::directionalhdl() : lighthdl(white*0.1f, white*0.5f, white)
{
	type = "directional";
}

directionalhdl::directionalhdl(const vec3f &direction, const vec3f &ambient, const vec3f &diffuse, const vec3f &specular) : lighthdl(ambient, diffuse, specular)
{
	type = "directional";
}

directionalhdl::~directionalhdl()
{

}

void directionalhdl::update(canvashdl *canvas)
{
    canvas->set_matrix(canvashdl::modelview_matrix);
    mat4f snapshot = canvas->matrices[canvas->active_matrix];
    canvas->translate(model->position);

    canvas->rotate(model->orientation[0], vec3f(1, 0, 0));
    canvas->rotate(model->orientation[1], vec3f(0, 1, 0));
    canvas->rotate(model->orientation[2], vec3f(0, 0, 1));

    canvas->update_normal_matrix();
    direction = canvas->matrices[canvashdl::normal_matrix] * vec4f(0.0, 0.0, -1.0, 0.0);
    canvas->matrices[canvas->active_matrix] = snapshot;
}

void directionalhdl::shade(vec3f &ambient, vec3f &diffuse, vec3f &specular, vec3f vertex, vec3f normal, float shininess) const
{
    /*  This code is mostly from the orange book. 
        reconciling terms from our code to the book code, 
        direction = location of the light source
        vertex = vertex being shaded. 
        */
    float ndotvp, ndothv, pf;
    vec3f vp = direction;
    vec3f hv = norm(direction - vertex/mag(vertex)); // direction from light position to vertex being shaded
    ndotvp = MAX(0.0, dot(normal, vp));
    ndothv = MAX(0.0, dot(normal, hv));

    if (ndotvp <= 0.0)
    {
        pf = 0.0;
    }
    else
    {
        pf = pow(ndothv, shininess);
    }

    ambient += this->ambient;
    diffuse += this->diffuse*ndotvp;
    specular += this->specular*pf;
}

pointhdl::pointhdl() : lighthdl(white*0.1f, white*0.5f, white)
{
	this->attenuation = vec3f(1.0, 0.14, 0.7);
	type = "point";
}

pointhdl::pointhdl(const vec3f &position, const vec3f &attenuation, const vec3f &ambient, const vec3f &diffuse, const vec3f &specular) : lighthdl(ambient, diffuse, specular)
{
	this->attenuation = attenuation;
	type = "point";
}

pointhdl::~pointhdl()
{

}

void pointhdl::update(canvashdl *canvas)
{
	/* TODO Assignment 3: Update the position of the light using the position of the attached model.
	 * The easiest thing is to do translations and rotations like you were going to render the object, and
	 * then just multiply the origin by the modelview matrix.
	 */
}

void pointhdl::shade(vec3f &ambient, vec3f &diffuse, vec3f &specular, vec3f vertex, vec3f normal, float shininess) const
{
	/* TODO Assignment 3: Implement a point light. See the OpenGL Orange Book in the references section
	 * of the course website. Its under the section about emulating the fixed function pipeline.
	 */
}

spothdl::spothdl() : lighthdl(white*0.1f, white*0.5f, white)
{
	this->attenuation = vec3f(1.0, 0.14, 0.7);
	this->cutoff = 0.5;
	this->exponent = 1.0;
	type = "spot";
}

spothdl::spothdl(const vec3f &attenuation, const float &cutoff, const float &exponent, const vec3f &ambient, const vec3f &diffuse, const vec3f &specular) : lighthdl(ambient, diffuse, specular)
{
	this->attenuation = attenuation;
	this->cutoff = cutoff;
	this->exponent = exponent;
	type = "spot";
}

spothdl::~spothdl()
{

}

void spothdl::update(canvashdl *canvas)
{
	/* TODO Assignment 3: Update both the direction and position of the light using the position and orientation
	 * of the attached model. See above.
	 */
}

void spothdl::shade(vec3f &ambient, vec3f &diffuse, vec3f &specular, vec3f vertex, vec3f normal, float shininess) const
{
	/* TODO Assignment 3: Implement a spot light. See the OpenGL Orange Book in the references section
	 * of the course website. Its under the section about emulating the fixed function pipeline.
	 */
}
