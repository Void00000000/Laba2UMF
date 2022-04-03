#pragma once
#include <fstream>
#include <vector>
using namespace std;

// ����������
struct sub_area {
    int nx1;  // ����� �������� � ������� X
    int nx2;  // ����� �������� � ������� X
    int ni;  // ����� ����������
};


class Area
{
protected:  // �������� �� ������ � ����� ������, �� � � ��� �����������
    int nw;  // ���������� �����������
    vector <double> Xw;  // �������� ������ ����������� 
    vector <sub_area> Mw;  // ������, �������� ��� ����������;
    double t0, t_end;
    void read_area() {  // ������ ��������� �������
        ifstream area_file;
        area_file.open("Area.txt");
        area_file >> nw;
        Xw.resize(nw + 1);
        Mw.resize(nw);
        for (int i = 0; i < nw + 1; i++)
            area_file >> Xw[i];
        for (int i = 0; i < nw; i++) {
            area_file >> Mw[i].ni >> Mw[i].nx1 >> Mw[i].nx2;
            Mw[i].nx1 -= 1;
            Mw[i].nx2 -= 1;
        }
        area_file.close();
    }
};
// ����� �����, ������� ��������� ����� ��������� �������
class Grid : public Area
{
private:
    vector <double> X;  // ������, �������� ���������� ���� �����
    vector <double> T;
    double nt;
    double qt;
    vector <int> IXw;  // ������, �������� ��������� ��������� ������ ����������� � X
    vector <int> nx;  // ���������� ���������
    vector <double> qx;  // ������������ ��������
    int n;  // ���������� �����

    void read_grid() {
        ifstream GridX;
        GridX.open("GridX.txt");
        nx.resize(nw);
        qx.resize(nw);
        for (int i = 0; i < nw; i++) {
            GridX >> nx[i] >> qx[i];
        }
        GridX.close();

        ifstream GridT;
        GridT.open("Time.txt");
        GridT >> t0 >> t_end >> nt >> qt;
        GridT.close();
    }


    void makeGrid(vector <double> &XT, double left, double right, int n, double qxt, int& i0, bool isTime) {
        double h0;
        static int j = 0;
        if (qxt - 1 < 1E-16)
            h0 = (right - left) / n;
        else
            h0 = (right - left) * (1 - qxt) / (1 - pow(qxt, n));

        XT[i0] = left;
        if (!isTime)
            IXw[j] = i0; j++;
        for (int i = i0 + 1; i < n + i0; i++) {
            XT[i] = XT[i - 1] + h0;
            h0 *= qxt;
        }
        i0 = n + i0;
    };


public:
    Grid() {   // ����������� ������
        read_area();
        read_grid();
        n = 0;
        for (int i = 0; i < nw; i++)
            n += nx[i];
        n++;
        X.resize(n);
        IXw.resize(nw + 1);
        int i0 = 0;
        for (int i = 0; i < nw; i++)
            makeGrid(X, Xw[i], Xw[i + 1], nx[i], qx[i], i0, false);
        X[n - 1] = Xw[nw];
        IXw[nw] = i0;

        T.resize(nt + 1);
        i0 = 0;
        makeGrid(T, t0, t_end, nt, qt, i0, true);
        T[nt] = t_end;
    }

    // ���������� ����� ����������, � ������� ��������� �������� ������� [x_p, x_p+1] 
    int inSubArea(int p) {
        int ixw1, ixw2;
        for (int i = 0; i < nw; i++) {
            ixw1 = IXw[Mw[i].nx1];
            ixw2 = IXw[Mw[i].nx2];
            if (p >= ixw1 && p <= ixw2 && p + 1 >= ixw1 && p + 1 <= ixw2)
                return Mw[i].ni;
        }
    }
    // ��������� �������
    double u0(double x) {
        return 10 * x + 9;
    }

    double sigma(int wi) {
        return 1;
    }

    double lambda(int wi, double qi, double t) {
        return 2 * qi + 3;
    }

    double f(int wi, int x, double t) {
        return -202;
    }

    double getX(int i) { return X[i]; }
    double getT(int i) { return T[i]; }
    int getN() { return n; }
    int getNT() { return nt + 1; }
};