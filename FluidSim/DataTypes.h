#pragma once

constexpr int GRID_SIZE = 150;
constexpr float INV_SIZE_SMALL = 1./(GRID_SIZE-1.);

struct Vec2D
{
	double x;
	double y;
};

struct VecVecField
{
	Vec2D data[GRID_SIZE][GRID_SIZE];
};

struct ScalarVecField
{
	double data[GRID_SIZE][GRID_SIZE];
};

struct PtForce
{
	Vec2D position;
	double magnitude;
};

double vecMag(const Vec2D& v);

Vec2D vecScale(const Vec2D& v, double s);

Vec2D vecAdd(const Vec2D& l, const Vec2D& r);

void init_ScalarVecField(ScalarVecField* p);

void init_VecVecField(VecVecField* u);

void copy_temp_to_u(VecVecField* u, VecVecField* temp);

void copy_temp_to_p(ScalarVecField* p, ScalarVecField* temp);