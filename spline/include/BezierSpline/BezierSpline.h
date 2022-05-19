//
// Created by Cain on 2022/4/24.
//

#ifndef DEMO5_BEZIERSPLINE_H
#define DEMO5_BEZIERSPLINE_H

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

using namespace std;

template<int nDegree>
class BezierSpline {
public:
    BezierSpline();
    bool fit(const vector<double> &wayPoints, const vector<double> &tMap,
             const vector<double> &diffInitial, const vector<double> &diffTerminal);
//    double bezier(double u, int n, const vector<double> &p);
    double bezier(double u, int n, const double *p);
    void operator()(double t, const vector<double> &tMap, vector<double> &state);
    static double combination(int m, int n);

    vector<vector<double>> _feature;
    vector<vector<double>> _featureAbs;
    Eigen::MatrixXd _matrixA;
    Eigen::VectorXd _vectorB;
    Eigen::VectorXd _controlPoints;
    Eigen::FullPivLU<Eigen::MatrixXd> _lu;
    int _nWayPoints, _nInitial, _nTerminal;
};

template<int nDegree>
BezierSpline<nDegree>::BezierSpline()
: _feature(nDegree + 1), _featureAbs(nDegree + 1), _nWayPoints(-1), _nInitial(-1), _nTerminal(-1) {
    double sign_n = 1., sign_m, cmn;
    for (int n = 0; n < _feature.size(); ++n){
        sign_m = sign_n;
        for (int m = 0; m <= n; ++m){
            cmn = combination(m, n);
            _feature[n].emplace_back(sign_m * cmn);
            _featureAbs[n].emplace_back(cmn);
            sign_m *= -1.;
        }
        sign_n *= -1.;
    }
}

template<int nDegree>
bool BezierSpline<nDegree>::fit(const vector<double> &wayPoints, const vector<double> &tMap,
                                const vector<double> &diffInitial, const vector<double> &diffTerminal) {
    int nWayPoints = wayPoints.size();
    int nInitial = diffInitial.size(), nTerminal = diffTerminal.size();
    if (nInitial + nTerminal != nDegree - 1){
        cerr << "The addition of the given number of differentials of the initial point and the end point is not equal to the polynomial degree - 1" << endl;
        return false;
    }

    if (tMap.size() != nWayPoints){
        cerr << "The length of tMap is not equal to the number of way points" << endl;
        return false;
    }

    int nControlPoints = (nWayPoints - 1) * (nDegree + 1);
    _matrixA.resize(nControlPoints, nControlPoints);
    _vectorB.resize(nControlPoints);
    _matrixA.setZero();
    _vectorB.setZero();

    // Matrix A
    // 0 ~ nDegree-1 阶连续, 右端点 - 左端点 = 0
    for (int i = 0; i < nWayPoints - 2; ++i){
        double rho = (tMap[i + 1] - tMap[i]) / (tMap[i + 2] - tMap[i + 1]), powRho = 1.;
        for (int j = 0; j < nDegree; ++j){
            for (int k = 0; k < _feature[j].size(); ++k){
                _matrixA(nDegree * i + j, (nDegree + 1) * (i + 1) - _feature[j].size() + k) = _feature[j][k];
                _matrixA(nDegree * i + j, (nDegree + 1) * (i + 1) + k) = -powRho * _feature[j][k];
            }
            powRho *= rho;
        }
    }

    // 第0 ~ nWayPoints-2 航点, 0 ~ nWayPoints-2 段的左端点
    int index = (nWayPoints - 2) * nDegree;
    for (int j = 0; j < nWayPoints - 1; ++j){
        for (int k = 0; k < _feature[0].size(); ++k){
            _matrixA(index + j, (nDegree + 1) * j + k) = _feature[0][k];
        }
    }

    // 第 nWayPoints - 1 个航点, 第 nWayPoints-2 段的右端点
    index += nWayPoints - 1;
    for (int k = 0; k < _feature[0].size(); ++k){
        _matrixA(index, (nDegree + 1) * (nWayPoints - 1) - _feature[0].size() + k) = _feature[0][k];
    }

    // 第 0 段的左端点 1 ~ nInitial 阶微分
    ++index;
    double rho = 1. / (tMap[1] - tMap[0]), powRho = 1., powN = 1.;
    int nDiff = nInitial;
    for (int j = 1; j <= nDiff; ++j){
        powRho *= rho;
        powN *= nDegree + 1 - j;
        for (int k = 0; k < _feature[j].size(); ++k){
            _matrixA(index + (j - 1), k) = powN * powRho * _feature[j][k];
        }
    }

    // 第 nWayPoints - 2 段的右端点 1 ~ nTerminal 阶微分
    index += nDiff;
    rho = 1. / (tMap[nWayPoints - 1] - tMap[nWayPoints - 2]), powRho = 1., powN = 1.;
    nDiff = nTerminal;
    for (int j = 1; j <= nDiff; ++j){
        powRho *= rho;
        powN *= nDegree + 1 - j;
        for (int k = 0; k < _feature[j].size(); ++k){
            _matrixA(index + (j - 1), (nDegree + 1) * (nWayPoints - 1) - _feature[j].size() + k) = powN * powRho * _feature[j][k];
        }
    }

    // Vector b
    index = (nWayPoints - 2) * nDegree;
    for (int i = 0; i < nWayPoints; ++i){
        _vectorB(index + i) = wayPoints[i];
    }

    index += nWayPoints;
    for (int i = 0; i < nInitial; ++i){
        _vectorB(index + i) = diffInitial[i];
    }

    index += nInitial;
    for (int i = 0; i < nTerminal; ++i){
        _vectorB(index + i) = diffTerminal[i];
    }

    _lu = _matrixA.fullPivLu();
    _controlPoints = _lu.solve(_vectorB);
//    _controlPoints = _matrixA.fullPivLu().solve(_vectorB);

    _nWayPoints = nWayPoints;
    _nInitial = nInitial;
    _nTerminal = nTerminal;

    return true;
}

template<int nDegree>
double BezierSpline<nDegree>::bezier(double u, int n, const double *p) {
    if (u <= 0.){
        return p[0];
    }
    if (u >= 1.){
        return p[n];
    }

    double v = 1. - u;
    vector<double> a(n + 1), b(n + 1);
    a[0] = 1.;
    b[0] = 1.;
    for (int i = 1; i <= n; ++i){
        a[i] = a[i - 1] * v;
        b[i] = b[i - 1] * u;
    }

    double value = 0.;
    for (int i = 0; i <= n; ++i){
        value += _featureAbs[n][i] * p[i] * a[n - i] * b[i];
    }

    return value;
}

template<int nDegree>
void BezierSpline<nDegree>::operator()(double t, const vector<double> &tMap, vector<double> &state) {
    if (tMap.size() != _nWayPoints){
        cerr << "The length of tMap is not equal to the number of way points" << endl;
        return;
    }

    t = min(max(tMap[0], t), tMap[tMap.size() - 1]);

    int l;
    for (l = 0; l < tMap.size() - 1; ++l){
        if (tMap[l + 1] >= t){
            break;
        }
    }

    int index = 0, indexLast;
    double p[(nDegree + 2) * (nDegree + 1) / 2];
    for (int j = 0; j <= nDegree; ++j){
        p[j] = _controlPoints(j + l * (nDegree + 1));
    }
    indexLast = index;
    index += nDegree + 1;
    for (int i = 1; i <= nDegree; ++i){
        for (int j = 0; j <= nDegree - i; ++j){
            // 其导数的控制点其控制点的差分, 省略了所乘的系数n_, 留到最右求值的时候在计算
            p[index + j] = p[indexLast + j + 1] - p[indexLast + j];
        }
        indexLast = index;
        index += nDegree + 1 - i;
    }

    index = 0;
    double rho = 1. / (tMap[l + 1] - tMap[l]), powRho = 1., powN = 1.;
    double u = (t - tMap[l]) * rho;
    state[0] = bezier(u, nDegree, p);
    index += nDegree + 1;
    for (int i = 1; i <= nDegree; ++i){
        // 导数需要考虑u与t的映射关系, 以及差分时省略的系数n_
        powRho *= rho;
        powN *= nDegree + 1 - i;
        state[i] = powN * powRho * bezier(u, nDegree - i, p + index);
        index += nDegree + 1 - i;
    }

}

template<int nDegree>
double BezierSpline<nDegree>::combination(int m, int n) {
    if (m > n || m < 0)
        return 0.;

    if (m == 0 || m == n){
        return 1.;
    }

    double value = 1.;
    for (int i = 0; i < m; ++i){
        value *= (double)(n - i) / (m - i);
    }
    return value;
}

#endif //DEMO5_BEZIERSPLINE_H
