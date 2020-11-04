#pragma once
#include <glad.h>
#include <glfw3.h>

struct GLData
{
	// stores OpenGL data I will actully need
	// i.e vertices, their associated VAO
	// VBO, EBO, and the overall shader program

	float* vertices;
	int num_vertices;
	unsigned int* indices;
	int num_indices;
	unsigned int VAO;
	unsigned int VBO;
	unsigned int EBO;
	unsigned int shader;
};

void framebuffer_size_callback(GLFWwindow* window, int width, int height);

void processInput(GLFWwindow* window);

GLFWwindow* init_window();

int update_window(GLFWwindow* window, const GLData* data);

unsigned int gen_shader();

GLData* init_gl(float* vertices, int num_vertices, unsigned int* indices, int num_indices);

void update_gl_vertices(GLData* data);

void delete_gldata(GLData* data);