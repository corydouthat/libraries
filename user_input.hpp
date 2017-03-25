// *****************************************************************************************************************************
// user_input.hpp
// User Input Functions using GLFW
// Author: Cory Douthat
// Copyright (c) 2015 Cory Douthat, All Rights Reserved.
// *****************************************************************************************************************************

#ifndef USER_INPUT_HPP_
#define USER_INPUT_HPP_

#define GLFW_DLL
#define GLFW_INCLUDE_GLU
#include "GLFW\glfw3.h"

#include "vec.hpp"

#define KEY_COUNT 349
#define MB_COUNT 8

// DATA
GLFWwindow *window_user_input = NULL;			// GLFW window pointer
unsigned int key_down_count[KEY_COUNT] = {0};	// Number of times each key has been pressed since last checked
unsigned int mouse_down_count[MB_COUNT];		// Number of times each mouse button has been pressed since last checked
double cursor_last_x = 0;						// Cursor x position when last checked
double cursor_last_y = 0;						// Cursor y position when last checked
double cursor_last_mb_x[MB_COUNT] = {0};		// Cursor x position when last checked or [] button pressed
double cursor_last_mb_y[MB_COUNT] = {0};		// Cursor x position when last checked or [] button pressed
double cursor_delta_mb_x[MB_COUNT] = {0};		// Change in cursor x position since last checked and while [] button pressed
double cursor_delta_mb_y[MB_COUNT] = {0};		// Change in cursor y position since last checked and while [] button pressed

// SETUP
void InitUserInput(GLFWwindow *window_in);

// CALLBACK FUNCTIONS
void KeyCallback(GLFWwindow *window_in,int key,int scancode,int action,int mods);
void MouseButtonCallback(GLFWwindow *window_in, int button, int action, int mods);
void CursorCallback(GLFWwindow *window_in, double xpos, double ypos);

// INFO
bool CheckKey(int key);
unsigned int KeyPressCount(int key);
bool CheckMousePress(int button);
unsigned int MousePressCount(int button);
Vec2f MouseMove();
Vec2f MouseMoveButtonDown(int button);

void KeyCallback(GLFWwindow *window,int key,int scancode,int action,int mods)
{
	// Increment key press counts
	if (key >= 0 && key < KEY_COUNT && action == GLFW_PRESS)
		key_down_count[key]++;
}

void MouseButtonCallback(GLFWwindow *window_in,int button,int action,int mods)
{
	// Increment button press counts
	if (button >= 0 && button < MB_COUNT && action == GLFW_PRESS)
		mouse_down_count[button]++;

	// Start button down motion tracking
	if (action == GLFW_PRESS)
		glfwGetCursorPos(window_in,&cursor_last_mb_x[button],&cursor_last_mb_y[button]);
}

void CursorCallback(GLFWwindow *window_in,double xpos,double ypos)
{
	// Button-down motion tracking
	for (int i = 0; i < MB_COUNT; i++)
	{
		if (glfwGetMouseButton(window_in,i) == GLFW_PRESS)
		{
			cursor_delta_mb_x[i] += xpos - cursor_last_mb_x[i];
			cursor_delta_mb_y[i] += ypos - cursor_last_mb_y[i];
			cursor_last_mb_x[i] = xpos;
			cursor_last_mb_y[i] = ypos;
		}
	}
}

// Initialize GLFW input callback functions
void InitUserInput(GLFWwindow *window_in)
{
	window_user_input = window_in;

	glfwSetKeyCallback(window_in, KeyCallback);
	glfwSetMouseButtonCallback(window_in, MouseButtonCallback);
	glfwSetCursorPosCallback(window_in, CursorCallback);

	// Reset Variables
	glfwGetCursorPos(window_in,&cursor_last_x,&cursor_last_y);
	for (int i = 0; i < KEY_COUNT; i++)
		key_down_count[i] = 0;
	for (int i = 0; i < MB_COUNT; i++)
	{
		mouse_down_count[i] = 0;
		cursor_last_mb_x[i] = cursor_last_x;
		cursor_last_mb_y[i] = cursor_last_y;
		cursor_delta_mb_x[i] = 0;
		cursor_delta_mb_y[i] = 0;
	}
}

// Check if key is currently pressed
bool CheckKey(int key)
{
	if (window_user_input == NULL)
		return false;

	return glfwGetKey(window_user_input,key) == GLFW_PRESS;
}

// Get key press count for a key and reset count to zero
unsigned int KeyPressCount(int key)
{
	if (key >= 0 && key <= KEY_COUNT)
	{
		unsigned int temp = key_down_count[key];
		key_down_count[key] = 0;
		return temp;
	}
	else
		return 0;
}

// Check if mouse button is currently pressed
bool CheckMousePress(int button)
{
	if (window_user_input == NULL)
		return false;

	return glfwGetMouseButton(window_user_input, button) == GLFW_PRESS;
}

// Get mouse button press count for a button and reset count to zero
unsigned int MousePressCount(int button)
{
	if (button >= 0 && button <= MB_COUNT)
	{
		unsigned int temp = mouse_down_count[button];
		mouse_down_count[button] = 0;
		return temp;
	}
	else
		return 0;
}

// Get change in mouse position and reset delta to zero
Vec2f MouseMove()
{
	if (window_user_input == NULL)
		return Vec2f(0.0f,0.0f);

	double xpos, ypos;
	Vec2f temp;
	glfwGetCursorPos(window_user_input, &xpos, &ypos);

	temp.x = (float)(xpos - cursor_last_x);
	temp.y = (float)(ypos - cursor_last_y);

	cursor_last_x = xpos;
	cursor_last_y = ypos;

	return temp;
}

// Get change in mouse position while button was held down and reset delta to zero
Vec2f MouseMoveButtonDown(int button)
{
	if (window_user_input == NULL)
		return Vec2f(0.0f,0.0f);

	if (button >= 0 && button < MB_COUNT)
	{
		// Save delta
		Vec2f temp((float)cursor_delta_mb_x[button],(float)cursor_delta_mb_y[button]);
		// Set delta to zero
		cursor_delta_mb_x[button] = 0;
		cursor_delta_mb_y[button] = 0;
		// Set last to current position
		glfwGetCursorPos(window_user_input,&cursor_last_mb_x[button],&cursor_last_mb_y[button]);

		return temp;
	}
	else
		return Vec2f();
}

#endif
