#include <iostream>
#include <fstream>
#include <cmath>

//# define SIZE 100
// grid size:
constexpr int SIZE = 50;
// grid spacing:
constexpr double dx = 1.;
constexpr double inv_dx = 1.; // for speed
// constexpr double dy = 1.; same as dx...
// step in time:
constexpr double dt = 1.;
// "kinematic viscosity"
constexpr double nu = 1.;
// for singularity at force...
constexpr double delta = 0.1;

typedef struct Vec2D
{
	double x;
	double y;
} Vec2D;

typedef struct VecVecField
{
	Vec2D data[SIZE][SIZE];
} VecVecField;

typedef struct ScalarVecField
{
	double data[SIZE][SIZE];
} ScalarVecField;

typedef struct PtForce
{
	Vec2D position;
	double magnitude;
} PtForce;

double vecMag(const Vec2D& v)
{
	return sqrt(pow(v.x, 2) + pow(v.y, 2));
}
Vec2D vecScale(const Vec2D& v, double s)
{
	return Vec2D{ v.x * s, v.y * s };
}

Vec2D vecAdd(const Vec2D& l, const Vec2D& r)
{
	return Vec2D{ l.x + r.x, l.y + r.y };
}

Vec2D compute_force(const Vec2D& pt, const PtForce& f)
{
	// compute 1/r potential i.e. f=q/r rhat direction
	double magnitude = f.magnitude / (delta + abs(sqrt(pow(f.position.x, 2) + pow(f.position.y, 2)) - sqrt(pow(pt.x, 2) + pow(pt.y, 2))));
	Vec2D direction = { pt.x - f.position.x, pt.y - f.position.y };
	return vecScale(direction, magnitude);
}

void init_ScalarVecField(ScalarVecField* p)
{
	for (int y = 0; y < SIZE; y++)
	{
		for (int x = 0; x < SIZE; x++)
		{
			p->data[y][x] = 0.;
		}
	}
}

void init_VecVecField(VecVecField* u)
{
	for (int y = 0; y < SIZE; y++)
	{
		for (int x = 0; x < SIZE; x++)
		{
			u->data[y][x] = { 0., 0. };
		}
	}
}

void copy_temp_to_u(VecVecField* u, VecVecField* temp)
{
	for (int y = 0; y < SIZE; y++)
	{
		for (int x = 0; x < SIZE; x++)
		{
			u->data[y][x] = temp->data[y][x];
		}
	}
}

void copy_temp_to_p(ScalarVecField* p, ScalarVecField* temp)
{
	for (int y = 0; y < SIZE; y++)
	{
		for (int x = 0; x < SIZE; x++)
		{
			p->data[y][x] = temp->data[y][x];
		}
	}
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
	for (int y = 1; y < SIZE - 1; y++)
	{
		for (int x = 1; x < SIZE - 1; x++)
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
	int number_iterations = 20;
	for (int repeat = 0; repeat < number_iterations; repeat++)
	{
		// exclude edges?
		for (int y = 1; y < SIZE - 1; y++)
		{
			for (int x = 1; x < SIZE - 1; x++)
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

	int number_iterations = 20;
	for (int repeat = 0; repeat < number_iterations; repeat++)
	{
		// exclude edges?
		for (int y = 1; y < SIZE - 1; y++)
		{
			for (int x = 1; x < SIZE - 1; x++)
			{
				temp->data[y][x] = inv_beta * (p->data[y - 1][x] + p->data[y + 1][x] + p->data[y][x + 1] + p->data[y][x - 1] + vec_div(w, x, y));
			}
		}
	}
	copy_temp_to_p(p, temp);
}

void subtract_pressure_gradient(VecVecField* w, ScalarVecField* p)
{
	for (int y = 1; y < SIZE - 1; y++)
	{
		for (int x = 1; x < SIZE - 1; x++)
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
	double impulse_radius = 1;
	for (int y = 1; y < SIZE - 1; y++)
	{
		for (int x = 1; x < SIZE - 1; x++)
		{
			Vec2D force = compute_force(u->data[y][x], f);
			u->data[y][x].x += force.x;
			u->data[y][x].y += force.y;
		}
	}
}

void enforce_boundary(VecVecField* u, ScalarVecField* p)
{
	// what to do about initial conditions vs just boundary?
	// must have u on perimeter = 0, and normal p deriv = 0
	// update at each step?
	for (int i = 1; i < SIZE - 1; i++)
	{
		// left:
		p->data[i][0] = p->data[i][1];
		u->data[i][0] = vecScale(u->data[i][1], -1.);
		// right:
		p->data[i][SIZE - 1] = p->data[i][SIZE - 2];
		u->data[i][SIZE - 1] = vecScale(u->data[i][SIZE - 2], 1.);
		// top:
		p->data[SIZE - 1][i] = p->data[SIZE - 2][i];
		u->data[SIZE - 1][i] = vecScale(u->data[SIZE - 2][i], -1.);
		// bot:
		p->data[0][i] = p->data[1][i];
		u->data[0][i] = vecScale(u->data[1][i], -1.);
	}
}

void step(VecVecField* u, VecVecField* temp_u, ScalarVecField* p, ScalarVecField* temp_p, const PtForce& f1, const PtForce& f2)
{
	advect(u, temp_u);
	iterative_poisson_diffusion(u, temp_u);
	add_force(u, f1);
	add_force(u, f2);
	iterative_poisson_pressure(p, temp_p, u);
	subtract_pressure_gradient(u, p);
	// done!
}


void display(ScalarVecField* p)
{
	std::cout.precision(2); // so shorter for output
	for (int y = 0; y < SIZE; y++)
	{
		for (int x = 0; x < SIZE; x++)
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
	file << "P3\n" << SIZE << " " << SIZE << "\n" << num_shades << "\n";
	// first determine range of values:
	double min_value = 0.;
	double max_value = 0.;
	for (int y = 0; y < SIZE; y++)
	{
		for (int x = 0; x < SIZE; x++)
		{
			if (p.data[y][x] > max_value)
				max_value = p.data[y][x];
			if (p.data[y][x] < min_value)
				min_value = p.data[y][x];
		}
	}
	//double range = max_value - min_value;
	for (int y = 0; y < SIZE; y++)
	{
		for (int x = 0; x < SIZE; x++)
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
	file << "P3\n" << SIZE << " " << SIZE << "\n" << num_shades << "\n";
	// first determine range of values:
	double max_valuex = 0.;
	double max_valuey = 0.;
	for (int y = 0; y < SIZE; y++)
	{
		for (int x = 0; x < SIZE; x++)
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
	for (int y = 0; y < SIZE; y++)
	{
		for (int x = 0; x < SIZE; x++)
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

	PtForce force = { {SIZE / 2. + 8, SIZE / 2}, 0.5 };
	PtForce force2 = { {SIZE / 2. - 8, SIZE / 2}, 0.5 };

	for (int i = 0; i < 10; i++)
	{
		// do 20 steps for now
		// do I want to display u or p? probably p since easier
		enforce_boundary(u, p); // should be in step but for display purposes it is here
		//display(p);
		print_to_ppm(*p);
		print_vel_to_ppm(*u);
		getchar();
		step(u, temp_u, p, temp_p, force, force2);
	}


	delete u;
	delete temp_u;
	delete p;
	delete temp_p;
}