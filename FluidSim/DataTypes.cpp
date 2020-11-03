#include "DataTypes.h"
#include <cmath>

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

void init_ScalarVecField(ScalarVecField* p)
{
	for (int y = 0; y < GRID_SIZE; y++)
	{
		for (int x = 0; x < GRID_SIZE; x++)
		{
			p->data[y][x] = 0.;
		}
	}
}

void init_VecVecField(VecVecField* u)
{
	for (int y = 0; y < GRID_SIZE; y++)
	{
		for (int x = 0; x < GRID_SIZE; x++)
		{
			u->data[y][x] = { 0., 0. };
		}
	}
}

void copy_temp_to_u(VecVecField* u, VecVecField* temp)
{
	for (int y = 0; y < GRID_SIZE; y++)
	{
		for (int x = 0; x < GRID_SIZE; x++)
		{
			u->data[y][x] = temp->data[y][x];
		}
	}
}

void copy_temp_to_p(ScalarVecField* p, ScalarVecField* temp)
{
	for (int y = 0; y < GRID_SIZE; y++)
	{
		for (int x = 0; x < GRID_SIZE; x++)
		{
			p->data[y][x] = temp->data[y][x];
		}
	}
}