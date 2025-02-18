#pragma once

#include <vector>

class RanMars {
public:
	RanMars(int);
	~RanMars();
	double uniform();
	void shuffle(std::vector<int>& arr);
	double gaussian();
	double gaussian(double mu, double sigma);
	void get_state(double*);
	void set_state(double*);

private:
	int save;
	double second;
	double* u;
	int i97, j97;
	double c, cd, cm;
};