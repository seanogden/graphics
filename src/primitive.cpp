/*
 * primitive.cpp
 *
 *  Created on: Dec 21, 2014
 *      Author: nbingham
 */

#include "primitive.h"

/* boxhdl
 *
 * Generate the geometry and indices required to make a box.
 */
boxhdl::boxhdl(float width, float height, float depth)
{
	/* TODO Assignment 1: Generate the geometry and indices required to make a box.
	 * Calculate its bounding box.
	 */
	rigid.push_back(rigidhdl());
	rigid[0].geometry.reserve(8);
    rigid[0].geometry.push_back(vec8f(-width, -height,  depth, 0.0, 0.0, 0.0, 0.0, 0.0)); //0
    rigid[0].geometry.push_back(vec8f( width, -height,  depth, 0.0, 0.0, 0.0, 0.0, 0.0)); //1
    rigid[0].geometry.push_back(vec8f(-width,  height,  depth, 0.0, 0.0, 0.0, 0.0, 0.0)); //2
    rigid[0].geometry.push_back(vec8f( width,  height,  depth, 0.0, 0.0, 0.0, 0.0, 0.0)); //3
    rigid[0].geometry.push_back(vec8f(-width, -height, -depth, 0.0, 0.0, 0.0, 0.0, 0.0)); //4
    rigid[0].geometry.push_back(vec8f( width, -height, -depth, 0.0, 0.0, 0.0, 0.0, 0.0)); //5
    rigid[0].geometry.push_back(vec8f(-width,  height, -depth, 0.0, 0.0, 0.0, 0.0, 0.0)); //6
    rigid[0].geometry.push_back(vec8f( width,  height, -depth, 0.0, 0.0, 0.0, 0.0, 0.0)); //7

    //front
    rigid[0].indices.push_back(0);
    rigid[0].indices.push_back(1);
    rigid[0].indices.push_back(2);

    rigid[0].indices.push_back(1);
    rigid[0].indices.push_back(2);
    rigid[0].indices.push_back(3);

    //back
    rigid[0].indices.push_back(4);
    rigid[0].indices.push_back(5);
    rigid[0].indices.push_back(6);

    rigid[0].indices.push_back(5);
    rigid[0].indices.push_back(6);
    rigid[0].indices.push_back(7);

    //top
    rigid[0].indices.push_back(2);
    rigid[0].indices.push_back(3);
    rigid[0].indices.push_back(6);

    rigid[0].indices.push_back(6);
    rigid[0].indices.push_back(3);
    rigid[0].indices.push_back(7);

    //right side
    rigid[0].indices.push_back(1);
    rigid[0].indices.push_back(3);
    rigid[0].indices.push_back(5);

    rigid[0].indices.push_back(3);
    rigid[0].indices.push_back(5);
    rigid[0].indices.push_back(7);

    //left side
    rigid[0].indices.push_back(0);
    rigid[0].indices.push_back(4);
    rigid[0].indices.push_back(6);

    rigid[0].indices.push_back(0);
    rigid[0].indices.push_back(6);
    rigid[0].indices.push_back(2);
    
    //bottom
    rigid[0].indices.push_back(0);
    rigid[0].indices.push_back(1);
    rigid[0].indices.push_back(4);

    rigid[0].indices.push_back(1);
    rigid[0].indices.push_back(4);
    rigid[0].indices.push_back(5);



	bound = vec6f(-width, width, -height, height, -depth, depth);
	// TODO Assignment 3: Set up the material properties for this object
}

boxhdl::~boxhdl()
{

}

/* spherehdl
 *
 * Generate the geometry and indices required to make a sphere.
 */
spherehdl::spherehdl(float radius, int levels, int slices)
{
	rigid.push_back(rigidhdl());

	rigid[0].geometry.reserve(2 + (levels-1)*slices);
	rigid[0].geometry.push_back(vec8f(0.0, 0.0, radius, 0.0, 0.0, 1.0, 0.0, 0.0));
	for (int i = 1; i < levels; i++)
		for (int j = 0; j < slices; j++)
		{
			vec3f dir(sin(m_pi*(float)i/(float)levels)*cos(2.0*m_pi*(float)j/(float)slices),
					  sin(m_pi*(float)i/(float)levels)*sin(2.0*m_pi*(float)j/(float)slices),
					  cos(m_pi*(float)i/(float)levels));
			rigid[0].geometry.push_back(vec8f(radius*dir[0], radius*dir[1], radius*dir[2],
									 dir[0], dir[1], dir[2], 0.0, 0.0));
		}
	rigid[0].geometry.push_back(vec8f(0.0, 0.0, -radius, 0.0, 0.0, -1.0, 0.0, 0.0));

	for (int i = 0; i < slices; i++)
	{
		rigid[0].indices.push_back(1 + (i+1)%slices);
		rigid[0].indices.push_back(1 + i);
		rigid[0].indices.push_back(0);
	}

	for (int i = 0; i < levels-2; i++)
		for (int j = 0; j < slices; j++)
		{
			rigid[0].indices.push_back(1 + i*slices + j);
			rigid[0].indices.push_back(1 + i*slices + (j+1)%slices);
			rigid[0].indices.push_back(1 + (i+1)*slices + j);

			rigid[0].indices.push_back(1 + (i+1)*slices + j);
			rigid[0].indices.push_back(1 + i*slices + (j+1)%slices);
			rigid[0].indices.push_back(1 + (i+1)*slices + (j+1)%slices);
		}

	for (int i = 0; i < slices; i++)
	{
		rigid[0].indices.push_back(1 + (levels-1)*slices);
		rigid[0].indices.push_back(1 + (levels-2)*slices + i);
		rigid[0].indices.push_back(1 + (levels-2)*slices + (i+1)%slices);
	}

	bound = vec6f(-radius, radius, -radius, radius, -radius, radius);

	// TODO Assignment 3: Set up the material properties for this object
}

spherehdl::~spherehdl()
{

}

/* cylinderhdl
 *
 * Generate the geometry and indices required to make a cylinder.
 */
cylinderhdl::cylinderhdl(float radius, float height, int slices)
{
	rigid.push_back(rigidhdl());


	rigid[0].geometry.reserve(2 + 2*slices);
    rigid[0].geometry.push_back(vec8f(0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0));

    for (int j = 0; j < slices; ++j)
    {
        float y = 0.0;
        float x = sin(j*2*m_pi/slices);
        float z = cos(j*2*m_pi/slices);
        rigid[0].geometry.push_back(vec8f(x*radius, 0.0, z*radius,
                    x, 0.0, z, 0.0, 0.0));
    }

    for (int j = 0; j < slices; ++j)
    {
        float y = height/radius;
        float x = sin(j*2*m_pi/slices);
        float z = cos(j*2*m_pi/slices);
        rigid[0].geometry.push_back(vec8f(x*radius, y*radius, z*radius,
                    x, y, z, 0.0, 0.0));
    }


    rigid[0].geometry.push_back(vec8f(0.0, height, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0));

	for (int i = 0; i < slices; i++)
	{
		rigid[0].indices.push_back(1 + (i+1)%slices);
		rigid[0].indices.push_back(1 + i);
		rigid[0].indices.push_back(0);
	}

    for (int j = 0; j < slices; j++)
    {
        rigid[0].indices.push_back(1 + j);
        rigid[0].indices.push_back(1 + (j+1)%slices);
        rigid[0].indices.push_back(1 + slices + j);

        rigid[0].indices.push_back(1 + slices + j);
        rigid[0].indices.push_back(1 + (j+1)%slices);
        rigid[0].indices.push_back(1 + slices + (j+1)%slices);
    }

	for (int i = 0; i < slices; i++)
	{
		rigid[0].indices.push_back(1 + 2*slices);
		rigid[0].indices.push_back(1 + slices + i);
		rigid[0].indices.push_back(1 + slices + (i+1)%slices);
	}

	// (left, right, bottom, top, front, back)
	bound = vec6f(-radius, radius, 0.0, height, -radius, radius);


	// TODO Assignment 3: Set up the material properties for this object
}

cylinderhdl::~cylinderhdl()
{

}

/* pyramidhdl
 *
 * Generate the geometry and indices required to make a pyramid.
 */
pyramidhdl::pyramidhdl(float radius, float height, int slices)
{
	rigid.push_back(rigidhdl());


	rigid[0].geometry.reserve(1 + slices);
    rigid[0].geometry.push_back(vec8f(0.0, height, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0));

    for (int j = 0; j < slices; ++j)
    {
        float y = 0.0;
        float x = radius*sin(j*2*m_pi/slices);
        float z = radius*cos(j*2*m_pi/slices);
        rigid[0].geometry.push_back(vec8f(x*radius, 0.0, z*radius,
                    x, 0.0, z, 0.0, 0.0));
    }

    rigid[0].geometry.push_back(vec8f(0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0));

	for (int i = 0; i < slices; i++)
	{
		rigid[0].indices.push_back(1 + (i+1)%slices);
		rigid[0].indices.push_back(1 + i);
		rigid[0].indices.push_back(0);
	}

	for (int i = 0; i < slices; i++)
	{
		rigid[0].indices.push_back(1 + (i+1)%slices);
		rigid[0].indices.push_back(1 + i);
		rigid[0].indices.push_back(rigid[0].geometry.size() - 1);
    }
	
	// (left, right, bottom, top, front, back)
	bound = vec6f(-radius, radius, 0.0, height, -radius, radius);

	// TODO Assignment 3: Set up the material properties for this object
}

pyramidhdl::~pyramidhdl()
{

}
