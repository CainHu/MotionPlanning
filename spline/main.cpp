#include <iostream>
#include <map>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <set>
#include <cstdio>
#include <fstream>
#include "Eigen/Eigen"
#include "Eigen/LU"
#include "BSpline.h"
#include "BezierSpline.h"

using namespace std;

int main() {
    // m+1个端点的位置与时刻
    vector<double> pos = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};
    vector<double> t = {0, 4, 7, 9, 11, 13, 15, 17, 19, 22, 26};

    // 起始点与终止点的速度和加速度
    double vel_init = 0., acc_init = 0.;
    double vel_term = 0., acc_term = 0.;


    vector<double> diff_init = {vel_init, acc_init};
    vector<double> diff_term = {vel_term, acc_term};

    BSpline<5> bSpline;
    BezierSpline<5> bezierSpline;

    vector<double> tMap = {0, 0, 0, 0, 0, 0, 4, 7, 9, 11, 13, 15, 17, 19, 22, 26, 26, 26, 26, 26, 26};

    if (bSpline.fit(pos, tMap, diff_init, diff_term)){
        cout << "bSpline success!!!\n";
        cout << bSpline._controlPoints << endl;
    }

    if (bezierSpline.fit(pos, t, diff_init, diff_term)){
        cout << "bezierSpline success!!!\n";
        cout << bezierSpline._controlPoints << endl;
    }

    cout << " --------------------- BSpline --------------------- " << endl;
    vector<double> state(6);
    double step = 0., dStep = 26. / 1000.;
    for (int i = 0; i < 1000; ++i){
        bSpline.deBoor(step, tMap, state);
        step += dStep;
        for (const auto &_s : state){
            cout << _s << " ";
        }
        cout << endl;
    }

    cout << " --------------------- Bezier Spline --------------------- " << endl;
    step = 0., dStep = 26. / 1000.;
    for (int i = 0; i < 1000; ++i){
        bezierSpline(step, t, state);
        step += dStep;
        for (const auto &_s : state){
            cout << _s << " ";
        }
        cout << endl;
    }

    return 0;
}
