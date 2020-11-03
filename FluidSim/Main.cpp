#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <windows.h>
#include <synchapi.h> // sleep
#include "DataTypes.h"
#include "GLDisplay.h"



// grid spacing:
constexpr double dx = 1.;
// for speed mult is faster than div:
constexpr double inv_dx = 1.;
// step in time:
constexpr double dt = 1.;
// "kinematic viscosity"
constexpr double nu = 1.;
// for singularity at force...
constexpr double delta = 0.1;


Vec2D compute_force(const Vec2D& pt, const PtForce& f)
{
	// compute 1/r force i.e. f=q/r rhat direction
	double magnitude = f.magnitude / (delta + abs(sqrt(pow(f.position.x, 2) + pow(f.position.y, 2)) - sqrt(pow(pt.x, 2) + pow(pt.y, 2))));
	Vec2D direction = { pt.x - f.position.x, pt.y - f.position.y };
	return vecScale(direction, magnitude);
}

Vec2D vec_interpolate(VecVecField* u, Vec2D& pos)
{
	// bilinear interpolation
	int x1 = (int)pos.x;
	int x2 = (int)pos.x + 1;
	int y1 = (int)pos.y;
	int y2 = (int)pos.y + 1;
	Vec2D fQ11 = u->data[y1][x1];
	Vec2D fQ12 = u->data[y2][x1];
	Vec2D fQ21 = u->data[y1][x2];
	Vec2D fQ22 = u->data[y2][x2];
	return { (x2 - pos.x) * (fQ11.x * (y2 - pos.y) + fQ12.x * (pos.y - y1)) + (pos.x - x1) * (fQ21.x * (y2 - pos.y) + fQ22.x * (pos.y - y1)),
			 (x2 - pos.x) * (fQ11.y * (y2 - pos.y) + fQ12.y * (pos.y - y1)) + (pos.x - x1) * (fQ21.y * (y2 - pos.y) + fQ22.y * (pos.y - y1)) };
	// ^ may be faster since not creating so many structs... instead splitting manually
	/*return vecAdd(vecScale(vecAdd(vecScale(fQ11, (y2 - pos.y)), vecScale(fQ12, (pos.y - y1))), (x2 - pos.x)),
		vecScale(vecAdd(vecScale(fQ21, (y2 - pos.y)), vecScale(fQ22, (pos.y - y1))), (pos.x - x1)));*/
}

void advect(VecVecField* u, VecVecField* temp)
{
	for (int y = 1; y < GRID_SIZE - 1; y++)
	{
		for (int x = 1; x < GRID_SIZE - 1; x++)
		{
			Vec2D oldpos = { (double)x - u->data[y][x].x, (double)y - u->data[y][x].y };
			Vec2D oldval = vec_interpolate(u, oldpos);
			temp->data[y][x] = oldval;
		}
	}
	copy_temp_to_u(u, temp);
}

void iterative_poisson_diffusion(VecVecField* u, VecVecField* temp)
{
	// temp holds k+1, reads from kth u, copies back into u before next step
	double alpha = pow(dx, 2) / (nu * dt);
	double beta = 4 + alpha;
	double inv_beta = 1 / beta; // since mult faster than divide

	// we can more intelligently determine
	// when to stop by seeing if change overall is small
	// i.e. add up change for each cell, see if average change
	// within some limit
	// this however would be an extra expense
	int number_iterations = 40;
	for (int repeat = 0; repeat < number_iterations; repeat++)
	{
		// exclude edges?
		for (int y = 1; y < GRID_SIZE - 1; y++)
		{
			for (int x = 1; x < GRID_SIZE - 1; x++)
			{
				temp->data[y][x] = vecScale(vecAdd(vecAdd(vecAdd(vecAdd(u->data[y - 1][x], u->data[y + 1][x]), u->data[y][x - 1]), u->data[y][x + 1]), vecScale(u->data[y][x], alpha)), inv_beta);
			}
		}
		copy_temp_to_u(u, temp);
	}
}

inline double vec_div(VecVecField* v, int x, int y)
{
	return 0.5 * inv_dx * (v->data[y][x + 1].x - v->data[y][x - 1].x) + 0.5 * inv_dx * (v->data[y + 1][x].y - v->data[y - 1][x].y);
}

void iterative_poisson_pressure(ScalarVecField* p, ScalarVecField* temp, VecVecField* w)
{
	// now need calculate divergence of "w" i.e. u (u will be w and then we fix it to be divergence free)
	// looks like i j symmetric... but they are taking i->x j->y :(
	// div w = w.x(x+1,y) - w.x(x-1,y)/dx + w.y().../dy
	double alpha = -pow(dx, 2);
	double beta = 4;
	double inv_beta = 1 / beta;

	int number_iterations = 40;
	for (int repeat = 0; repeat < number_iterations; repeat++)
	{
		// exclude edges?
		for (int y = 1; y < GRID_SIZE - 1; y++)
		{
			for (int x = 1; x < GRID_SIZE - 1; x++)
			{
				temp->data[y][x] = inv_beta * (p->data[y - 1][x] + p->data[y + 1][x] + p->data[y][x + 1] + p->data[y][x - 1] + vec_div(w, x, y));
			}
		}
		copy_temp_to_p(p, temp);
	}
}

void subtract_pressure_gradient(VecVecField* w, ScalarVecField* p)
{
	for (int y = 1; y < GRID_SIZE - 1; y++)
	{
		for (int x = 1; x < GRID_SIZE - 1; x++)
		{
			double grad_p_x = 0.5 * inv_dx * (p->data[y][x + 1] - p->data[y][x - 1]);
			double grad_p_y = 0.5 * inv_dx * (p->data[y + 1][x] - p->data[y - 1][x]);
			w->data[y][x].y -= grad_p_y;
			w->data[y][x].x -= grad_p_x;
		}
	}
}

void add_force(VecVecField* u, const PtForce& f)
{
	// F*t = delta momentum = m * (delta v) i.e. take m=1 get delta v, that is u
	// just doing a simple F scalar not the vector case given.. (easier to input)
	double impulse_radius = 1.;
	for (int y = 1; y < GRID_SIZE - 1; y++)
	{
		for (int x = 1; x < GRID_SIZE - 1; x++)
		{
			Vec2D force = compute_force(u->data[y][x], f);
			u->data[y][x].x += force.x;
			u->data[y][x].y += force.y;
		}
	}
}

void enforce_boundary(VecVecField* u, ScalarVecField* p)
{
	// must have u on perimeter = 0, and normal p deriv = 0
	// update at each step
	for (int i = 1; i < GRID_SIZE - 1; i++)
	{
		// left:
		p->data[i][0] = p->data[i][1];
		u->data[i][0] = vecScale(u->data[i][1], -1.);
		// right:
		p->data[i][GRID_SIZE - 1] = p->data[i][GRID_SIZE - 2];
		u->data[i][GRID_SIZE - 1] = vecScale(u->data[i][GRID_SIZE - 2], 1.);
		// top:
		p->data[GRID_SIZE - 1][i] = p->data[GRID_SIZE - 2][i];
		u->data[GRID_SIZE - 1][i] = vecScale(u->data[GRID_SIZE - 2][i], -1.);
		// bot:
		p->data[0][i] = p->data[1][i];
		u->data[0][i] = vecScale(u->data[1][i], -1.);
	}
}

void step(VecVecField* u, VecVecField* temp_u, ScalarVecField* p, ScalarVecField* temp_p, const std::vector<PtForce>& forces)
{
	advect(u, temp_u);
	iterative_poisson_diffusion(u, temp_u);
	for (auto& f : forces)
	{
		// sum over all forces
		add_force(u, f);
	}
	iterative_poisson_pressure(p, temp_p, u);
	subtract_pressure_gradient(u, p);
	// done!
}

void display(ScalarVecField* p)
{
	std::cout.precision(2); // so shorter for output
	for (int y = 0; y < GRID_SIZE; y++)
	{
		for (int x = 0; x < GRID_SIZE; x++)
		{
			std::cout << p->data[y][x] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << std::endl;
}

void print_to_ppm(const ScalarVecField& p)
{
	// pgm not working well, try ppm instead
	int num_shades = 255; // what to set to?
	std::ofstream file;
	file.open("img.ppm", std::ios::trunc);
	file << "P3\n" << GRID_SIZE << " " << GRID_SIZE << "\n" << num_shades << "\n";
	// first determine range of values:
	double min_value = 0.;
	double max_value = 0.;
	for (int y = 0; y < GRID_SIZE; y++)
	{
		for (int x = 0; x < GRID_SIZE; x++)
		{
			if (p.data[y][x] > max_value)
				max_value = p.data[y][x];
			if (p.data[y][x] < min_value)
				min_value = p.data[y][x];
		}
	}
	//double range = max_value - min_value;
	for (int y = 0; y < GRID_SIZE; y++)
	{
		for (int x = 0; x < GRID_SIZE; x++)
		{
			// i.e. largest goes to max shade... not a perfect soln if a lot are much less
			// will handle negative as a different color
			if (p.data[y][x] > 0)
			{
				int shade = (num_shades * p.data[y][x]) / max_value;
				file << shade << " 0 0   ";
			}
			else
			{
				int shade = abs((num_shades * p.data[y][x]) / min_value);
				file << "0 " << shade << " 0   ";
			}

		}
		file << "\n";
	}
}

void print_vel_to_ppm(const VecVecField& u)
{
	// add velocity by multiplying each component! i.e. R=y, G=x
	// doesn't work well since takes max each time... maybe fix as really high?
	// now no negative! but doesn't contain direction...
	// pgm not working well, try ppm instead
	int num_shades = 255; // what to set to?
	std::ofstream file;
	file.open("velimg.ppm", std::ios::trunc);
	file << "P3\n" << GRID_SIZE << " " << GRID_SIZE << "\n" << num_shades << "\n";
	// first determine range of values:
	double max_valuex = 0.;
	double max_valuey = 0.;
	for (int y = 0; y < GRID_SIZE; y++)
	{
		for (int x = 0; x < GRID_SIZE; x++)
		{
			if (abs(u.data[y][x].x) > max_valuex)
				max_valuex = abs(u.data[y][x].x);
			if (abs(u.data[y][x].y) > max_valuey)
				max_valuey = abs(u.data[y][x].y);
		}
	}
	max_valuex = 2;
	max_valuey = 2;
	//double range = max_value - min_value;
	for (int y = 0; y < GRID_SIZE; y++)
	{
		for (int x = 0; x < GRID_SIZE; x++)
		{
			// i.e. largest goes to max shade... not a perfect soln if a lot are much less
			// will handle negative as a different color

			int shadex = (num_shades * abs(u.data[y][x].x)) / max_valuex;
			int shadey = (num_shades * abs(u.data[y][x].y)) / max_valuey;
			file << shadex << " 0 " << shadey << "   ";

		}
		file << "\n";
	}
}

inline float grid_to_screen(int gridpos)
{
	// convert grid pos in 0..SIZE
	// to screen coord in -1..1
	// was just int div error... std::cout << gridpos << " " << 2 * ((float)gridpos - 0.5 * SIZE) * INV_SIZE << std::endl;
	return 2 * ((float)gridpos - 0.5 * (GRID_SIZE-1.)) * INV_SIZE_SMALL; //(float)gridpos / (GRID_SIZE-1);//2 * ((float)gridpos - 0.5 * SIZE) * INV_SIZE;
}

inline int vertex_access(int y, int x)
{
	// access the vertex array
	return 3 * y * GRID_SIZE + 3 * x;
}

inline int index_access(int y, int x)
{
	// access the vertex array
	return 1 * y * GRID_SIZE + 1 * x;
}

float* generate_vertices()
{
	// need to normalize to -1..1
	// unfortunately doesn't match my setup, being a long 1d
	// array rather than 2d but it's ok just be careful
	float* vertices = new float[3 * GRID_SIZE * GRID_SIZE];
	for (int y = 0; y < GRID_SIZE; y++)
	{
		for (int x = 0; x < GRID_SIZE; x++)
		{
			float xpos = grid_to_screen(x);
			float ypos = grid_to_screen(y);
			std::cout << x << " " << y << " " << xpos << " " << ypos << std::endl; // above func wrong..
			// make sure x y z order is right
			// convert to vertex ordering: 3 (for now, # points, larger with color) 3*x + SIZE*y yup
			vertices[vertex_access(y, x) + 0] = xpos;
			vertices[vertex_access(y, x) + 1] = ypos;
			vertices[vertex_access(y, x) + 2] = 0;
		}
	}
	return vertices;
}

unsigned int* generate_indices()
{
	// need to normalize to -1..1
	// unfortunately doesn't match my setup, being a long 1d
	// array rather than 2d but it's ok just be careful
	// six lines per square from 2 triangles
	unsigned int* indices = new unsigned int[6*(GRID_SIZE-1)*(GRID_SIZE-1)];
	for (int y = 0; y < GRID_SIZE - 1; y++)
	{
		for (int x = 0; x < GRID_SIZE - 1; x++)
		{
			// need 6 per value of xy, following example given
			// but problem since don't correspond to 0..3, here SIZE
			// since many many more vertices
			// so do 0..3 (i.e a single square) + 3*x + 3*y*SIZE i.e. 
			// nope for a single square from bottom left xy, do 6 vertices, order not too important, but still have to group into 2 triangles:
			// sw = xy, 2*se = x+1y, 2*nw=xy+1, ne=x+1y+1
			indices[6*y * (GRID_SIZE - 1) + 6 * x + 0] = index_access(y, x); // sw
			indices[6 * y * (GRID_SIZE - 1) + 6 * x + 1] = index_access(y, x) + 1; // se
			indices[6 * y * (GRID_SIZE - 1) + 6 * x + 2] = index_access(y, x) + GRID_SIZE; // nw now not 3 scaled since looking at row not elem
			indices[6 * y * (GRID_SIZE - 1) + 6 * x + 3] = index_access(y, x) + 1; // se
			indices[6 * y * (GRID_SIZE - 1) + 6 * x + 4] = index_access(y, x) + GRID_SIZE + 1; // ne
			indices[6 * y * (GRID_SIZE - 1) + 6 * x + 5] = index_access(y, x) + GRID_SIZE; // nw
		}
	}
	return indices;
}

int main()
{
	// init to heap so doesn't overflow stack...
	VecVecField* u = new VecVecField();
	VecVecField* temp_u = new VecVecField();
	ScalarVecField* p = new ScalarVecField();
	ScalarVecField* temp_p = new ScalarVecField();
	// zero out:
	init_VecVecField(u);
	init_VecVecField(temp_u);
	init_ScalarVecField(p);
	init_ScalarVecField(temp_p);

	PtForce force = { {GRID_SIZE / 2. + 8, GRID_SIZE / 2}, 0.5 };
	PtForce force2 = { {GRID_SIZE / 2. - 8, GRID_SIZE / 2}, 0.5 };
	std::vector<PtForce> forces;
	forces.push_back(force);
	forces.push_back(force2);

	GLFWwindow* glwindow = init_window();
	
	int num_vertices = 3 * GRID_SIZE * GRID_SIZE;//12;
	float* vertices = generate_vertices();//new float[num_vertices];
	// box:
	/*vertices[0] = 0.5f;
	vertices[1] = 0.5f;
	vertices[2] = 0.0f;
	vertices[3] = 0.5f;
	vertices[4] = -0.5f;
	vertices[5] = 0.0f;
	vertices[6] = -0.5f;
	vertices[7] = -0.5f;
	vertices[8] = 0.0f;
	vertices[9] = -0.5f;
	vertices[10] = 0.5f;
	vertices[11] = 0.0f;*/
	// for square, need 6 per square, so will eventually have 6*SIZE*SIZE indices
	// and 12*SIZE**2 vertices? ignoring overlap since more difficult I guess
	int num_indices = 6 * (GRID_SIZE - 1) * (GRID_SIZE - 1);//6;
	unsigned int* indices = generate_indices();//new unsigned int[num_indices]; 
	//indices[0] = 0;
	//indices[1] = 1;
	//indices[2] = 3;
	//indices[3] = 1;
	//indices[4] = 2;
	//indices[5] = 3;
	GLData* gldata = init_gl(vertices, num_vertices, indices, num_indices);

	for (int i = 0; i < 10; i++)
	{
		// do 20 steps for now
		// do I want to display u or p? probably p since easier
		enforce_boundary(u, p); // should be in step but for display purposes it is here
		//display(p);
		print_to_ppm(*p);
		print_vel_to_ppm(*u);
		update_window(glwindow, gldata);
		//getchar();
		step(u, temp_u, p, temp_p, forces);
	}
	// when done just keep rendering openGL:
	while (!update_window(glwindow, gldata))
	{
		// need to cap framerate
		Sleep(500);
	}
	//getchar();
	delete u;
	delete temp_u;
	delete p;
	delete temp_p;
}