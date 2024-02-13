#pragma once

class RanPark {
public:
	RanPark(int);
	double uniform();
	double gaussian();
	void reset(int);
	void reset(int, double*);
	int state();

private:
	int seed, save;
	double second;
};