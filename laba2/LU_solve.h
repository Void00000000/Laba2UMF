#pragma once
#include <vector>
using namespace std;
// LU разложение
void calculation_lu(vector <double>& L, vector <double>& D, vector <double>& U,
	vector <double>& l, vector <double>& d, vector <double>& r, int n);

void calculation_y(vector <double>& L, vector <double>& b, vector <double>& q, int n);

void calucaltion_x(vector <double>& D, vector <double>& U, vector <double>& q, int n);