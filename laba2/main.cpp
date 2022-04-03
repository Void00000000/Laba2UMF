#include "FixedPointIteration.h"
int main() {
    double w = 1;
    Grid grid;
    FixedPointIteration solution(grid);
    double t1, t2;
    ofstream fout("output.txt");
    for (int i = 1; i < grid.getNT(); i++) {
        t1 = grid.getT(i - 1);
        t2 = grid.getT(i);
        solution.fixed_point_iteration(grid, w, t1, t2);
        solution.print_vector(fout, t2);
    }
    fout.close();
}