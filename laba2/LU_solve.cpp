#include <vector>
using namespace std;

// LU разложение с учётом ленточного формата хранения матрицы
void calculation_lu(vector <double>& L, vector <double>& D, vector <double>& U,
	vector <double>& l, vector <double>& d, vector <double>& r, int n) {
	D[0] = d[0];
	U[0] = r[0];
	for (int i = 1; i < n - 1; i++) {
		L[i - 1] = l[i - 1] / D[i - 1];
		D[i] = d[i] - L[i - 1] * U[i - 1];
		U[i] = r[i];
	}
	L[n - 2] = l[n - 2] / D[n - 2];
	D[n - 1] = d[n - 1] - L[n - 2] * U[n - 2];
}


// Вычисление вектора y прямым ходом, храниться он будет на месте вектора b
void calculation_y(vector <double>& L, vector <double>& b, vector <double>& q, int n) {
	q[0] = b[0];
	for (int i = 1; i < n; i++) 
		q[i] = b[i] - L[i - 1] * q[i - 1];
}


// Вычисление вектора x обратным ходом, храниться он будет на месте вектора b
void calucaltion_x(vector <double>& D, vector <double>& U, vector <double>& q, int n) {
	q[n - 1] /= D[n - 1];
	for (int i = n - 2; i >= 0; i--) {
		q[i] -= U[i] * q[i + 1];
		q[i] /= D[i];
	}
}