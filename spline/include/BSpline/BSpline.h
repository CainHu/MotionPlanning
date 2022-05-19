//
// Created by Cain on 2022/4/20.
//

#ifndef DEMO5_BSPLINE_H
#define DEMO5_BSPLINE_H

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
class BSpline {
public:
    BSpline();
    bool fit(const vector<double> &wayPoints, const vector<double> &diffInitial, const vector<double> &diffTerminal);
    bool fit(const vector<double> &wayPoints, const vector<double> &tMap,
             const vector<double> &diffInitial, const vector<double> &diffTerminal);
    void deBoor(double t, const vector<double> &tMap, vector<double> &state);

    // M = feature / n!, P(t) = [t^n, t^(n-1), ..., 1] * M * [P0, P1, ..., Pn]^T
    Eigen::MatrixXd _feature;
    Eigen::MatrixXd _matrixA;
    Eigen::VectorXd _vectorB;
    Eigen::VectorXd _controlPoints;
    int _nWayPoints, _nInitial, _nTerminal;
    double operator()(double t, const vector<double> &tMap);
private:
    Eigen::FullPivLU<Eigen::MatrixXd> _lu;
    static double combination(int m, int n);
    static double arrangement(int m, int n);
};

template<int nDegree>
BSpline<nDegree>::BSpline() : _feature(nDegree + 1, nDegree + 1),
                              _nWayPoints(-1), _nInitial(-1), _nTerminal(-1) {
    vector<double> c(nDegree + 1, 1.);
    for (int i = 0; i <= nDegree; ++i){
        for (int j = 0; j <= nDegree ; ++j){
            double s = 0.;
            double sign = 1.;
            for (int k = j; k <= nDegree; ++k){
                s += c[k] * sign * combination(k - j, nDegree + 1);
                sign *= -1.;
            }
            _feature(i, j) = combination(i, nDegree) * s;
        }
        for (int k = 0; k <= nDegree; ++k){
            c[k] *= nDegree - k;
        }
    }
}

template<int k>
double BSpline<k>::combination(int m, int n) {
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

template<int nDegree>
double BSpline<nDegree>::arrangement(int m, int n) {
    if (m > n || m < 0)
        return 0.;

    double value = 1.;
    for (int i = 0; i < m; ++i){
        value *= (double)(n - i);
    }
    return value;
}

template<int nDegree>
bool BSpline<nDegree>::fit(const vector<double> &wayPoints, const vector<double> &diffInitial,
                           const vector<double> &diffTerminal) {
    int nWayPoints = wayPoints.size();
    if (nWayPoints < nDegree + 1){
        cerr << "The number of way points is less than polynomial degree + 1" << endl;
        return false;
    }

    int nInitial = diffInitial.size(), nTerminal = diffTerminal.size();
    if (nInitial + nTerminal != nDegree - 1){
        cerr << "The addition of the given number of differentials of the initial point and the end point is not equal to the polynomial degree - 1" << endl;
        return false;
    }

    if (nWayPoints != _nWayPoints){
        _matrixA.resize(nWayPoints + nDegree - 1, nWayPoints + nDegree - 1);
        _vectorB.resize(nWayPoints + nDegree - 1);
        _matrixA.setZero();
        _vectorB.setZero();

        // Matrix A
        for (int i = 0; i < nWayPoints; ++i){
            for (int j = 0; j < nDegree; ++j){
                _matrixA(i, i + j) = _feature(nDegree, j);
            }
        }

        for (int i = 0; i < nInitial; ++i){
            for (int j = 0; j < nDegree; ++j){
                _matrixA(i + nWayPoints, j) = _feature(nDegree - 1 - i, j);
            }
        }

        for (int i = 0; i < nTerminal; ++i){
            for (int j = 0; j < nDegree; ++j){
                _matrixA(i + nWayPoints + nInitial, j + nWayPoints - 1) = _feature(nDegree - 1 - i, j);
            }
        }

        _lu = _matrixA.fullPivLu();

    } else if (nInitial != _nInitial){
        // rewrite
        for (int i = 0; i < nInitial; ++i){
            for (int j = 0; j < nDegree; ++j){
                _matrixA(i + nWayPoints, j) = _feature(nDegree - 1 - i, j);
            }
        }

        for (int i = 0; i < nTerminal; ++i){
            for (int j = 0; j < nDegree; ++j){
                _matrixA(i + nWayPoints + nInitial, j + nWayPoints - 1) = _feature(nDegree - 1 - i, j);
            }
        }

        // clear
        for (int i = nInitial; i < _nInitial; ++i){
            for (int j = 0; j < nDegree; ++j){
                _matrixA(i + nWayPoints, j) = _feature(nDegree - 1 - i, j);
            }
        }

        for (int i = 0; i < _nTerminal - nTerminal; ++i){
            for (int j = 0; j < nDegree; ++j){
                _matrixA(i + nWayPoints + _nInitial, j + nWayPoints - 1) = _feature(nDegree - 1 - i, j);
            }
        }

        _lu = _matrixA.fullPivLu();
    }

    // Vector b
    double arrangeTemp, arrange = arrangement(nDegree, nDegree);
    for (int i = 0; i < nWayPoints; ++i){
        _vectorB(i) = wayPoints[i] * arrange;
    }

    arrangeTemp = arrange;
    for (int i = 0; i < nInitial; ++i){
        arrangeTemp = arrangeTemp / (i + 1);
        _vectorB(nWayPoints + i) = diffInitial[i] * arrangeTemp;
    }

    arrangeTemp = arrange;
    for (int i = 0; i < nTerminal; ++i){
        arrangeTemp = arrangeTemp / (i + 1);
        _vectorB(nWayPoints + nInitial + i) = diffTerminal[i] * arrangeTemp;
    }

    _controlPoints = _lu.solve(_vectorB);
    _nWayPoints = nWayPoints;
    _nInitial = nInitial;
    _nTerminal = nTerminal;

    return true;
}

template<int nDegree>
bool
BSpline<nDegree>::fit(const vector<double> &wayPoints, const vector<double> &tMap, const vector<double> &diffInitial,
                      const vector<double> &diffTerminal) {
    int nWayPoints = wayPoints.size();
    if (nWayPoints < nDegree + 1){
        cerr << "The number of way points is less than polynomial degree + 1" << endl;
        return false;
    }

    int nInitial = diffInitial.size(), nTerminal = diffTerminal.size();
    if (nInitial + nTerminal != nDegree - 1){
        cerr << "The addition of the given number of differentials of the initial point and the end point is not equal to the polynomial degree - 1" << endl;
        return false;
    }

    if (tMap.size() != nDegree + nWayPoints + nDegree){
        cerr << "The length of tMap is not equal to twice the polynomial degree + the number of way points" << endl;
        return false;
    }

    _matrixA.resize(nWayPoints + nDegree - 1, nWayPoints + nDegree - 1);
    _matrixA.setZero();

    _vectorB.resize(nWayPoints + nDegree - 1);
    _vectorB.setZero();

    vector<double> basisFull((nDegree + 2) * (nDegree + 1) / 2);
    vector<double> basisPartial(nDegree + 1, 0.);

    int basisFullIndex = 0;
    int l = nDegree;
    double t = tMap[l];
    basisPartial[nDegree] = 1.;
    basisFull[basisFullIndex++] = basisPartial[nDegree];
    for (int i = 1; i <= nDegree; ++i){
        double tor = (tMap[l + 1] - t) / (tMap[l + 1] - tMap[l - i + 1]);
        basisPartial[nDegree - i] = basisPartial[nDegree - i + 1] * tor;
        basisFull[basisFullIndex++] = basisPartial[nDegree - i];
        for (int j = 1; j < i; ++j){
            double tor_ = 1. - tor;
            tor = (tMap[l + j + 1] - t) / (tMap[l + j + 1] - tMap[l - i + j + 1]);
            basisPartial[nDegree - i + j] = basisPartial[nDegree - i + j] * tor_ + basisPartial[nDegree - i + j + 1] * tor;
            basisFull[basisFullIndex++] = basisPartial[nDegree - i + j];
        }
        basisPartial[nDegree] = basisPartial[nDegree] * (1. - tor);
        basisFull[basisFullIndex++] = basisPartial[nDegree];
    }

    int row = l - nDegree;
    for (int i = 0; i < nDegree; ++i){     // i = nDegree对应的coefficient为0 int i = 0; i <= nDegree; ++i
        _matrixA(row, row + i) = basisPartial[i];
    }

    auto tangentCoef = [&] (int nDiff, int tSetOff, int rowSetOff, int colSetOff){
        int rhoIndex = 0;
        vector<double> rho((nDegree + nDegree - nDiff + 1) * nDiff / 2);
        for (int k = 1; k <= nDiff; ++k){
            int order = nDegree - k + 1;
            for (int i = 0; i <= nDegree - k; ++i){
                rho[rhoIndex++] = order / (tMap[tSetOff + i + nDegree + 1] - tMap[tSetOff + i + k]);
            }
        }

        for (int k = 1; k <= nDiff; ++k){
            int row = rowSetOff + k - 1;
            for (int j = 0; j <= nDegree - k; ++j){
                _matrixA(row, colSetOff + j) = basisFull[(nDegree - k + 1) * (nDegree - k) / 2 + j];
            }
            for (int i = nDegree - k; i < nDegree; ++i){
                rhoIndex = (nDegree + i + 1) * (nDegree - i) / 2 - 1;
                double rhoLast, rhoCurr;
                int j = i + 1;
                rhoCurr = rho[rhoIndex];
                _matrixA(row, colSetOff + j) = rhoCurr * _matrixA(row, colSetOff + j - 1);
                rhoLast = rhoCurr;
                --j;
                for (; j > 0; --j){
                    rhoCurr = rho[rhoIndex + j - (i + 1)];
                    _matrixA(row, colSetOff + j) = rhoCurr * _matrixA(row, colSetOff + j - 1) - rhoLast * _matrixA(row, colSetOff + j);
                    rhoLast = rhoCurr;
                }
                _matrixA(row, colSetOff + j) = -rhoLast * _matrixA(row, colSetOff + j);
            }
        }
    };

    tangentCoef(diffInitial.size(), 0, nWayPoints, 0);

    for (l = nDegree + 1; l < nWayPoints + nDegree - 1; ++l){
        t = tMap[l];
        basisPartial[nDegree] = 1.;
        for (int i = 1; i <= nDegree; ++i){
            double tor = (tMap[l + 1] - t) / (tMap[l + 1] - tMap[l - i + 1]);
            basisPartial[nDegree - i] = basisPartial[nDegree - i + 1] * tor;
            for (int j = 1; j < i; ++j){
                double tor_ = 1. - tor;
                tor = (tMap[l + j + 1] - t) / (tMap[l + j + 1] - tMap[l - i + j + 1]);
                basisPartial[nDegree - i + j] = basisPartial[nDegree - i + j] * tor_ + basisPartial[nDegree - i + j + 1] * tor;
            }
            basisPartial[nDegree] = basisPartial[nDegree] * (1. - tor);
        }

        row = l - nDegree;
        for (int i = 0; i < nDegree; ++i){     // i = nDegree对应的coefficient为0 int i = 0; i <= nDegree; ++i
            _matrixA(row, row + i) = basisPartial[i];
        }
    }

    basisFullIndex = 0;
    l = nWayPoints + nDegree - 2;
    t = tMap[l + 1];
    basisPartial[nDegree] = 1.;
    basisFull[basisFullIndex++] = basisPartial[nDegree];
    for (int i = 1; i <= nDegree; ++i){
        double tor = (tMap[l + 1] - t) / (tMap[l + 1] - tMap[l - i + 1]);
        basisPartial[nDegree - i] = basisPartial[nDegree - i + 1] * tor;
        basisFull[basisFullIndex++] = basisPartial[nDegree - i];
        for (int j = 1; j < i; ++j){
            double tor_ = 1. - tor;
            tor = (tMap[l + j + 1] - t) / (tMap[l + j + 1] - tMap[l - i + j + 1]);
            basisPartial[nDegree - i + j] = basisPartial[nDegree - i + j] * tor_ + basisPartial[nDegree - i + j + 1] * tor;
            basisFull[basisFullIndex++] = basisPartial[nDegree - i + j];
        }
        basisPartial[nDegree] = basisPartial[nDegree] * (1. - tor);
        basisFull[basisFullIndex++] = basisPartial[nDegree];
    }

    row = l - nDegree;
    for (int i = 1; i <= nDegree; ++i){     // i = 0对应的coefficient为0 int i = 0; i <= nDegree; ++i
        _matrixA(row + 1, row + i) = basisPartial[i];
    }

    tangentCoef(diffTerminal.size(), l - nDegree, nWayPoints + diffInitial.size(), nWayPoints - 2);

    // Vector b
    for (int i = 0; i < nWayPoints; ++i){
        _vectorB(i) = wayPoints[i];
    }

    for (int i = 0; i < nInitial; ++i){
        _vectorB(nWayPoints + i) = diffInitial[i];
    }

    for (int i = 0; i < nTerminal; ++i){
        _vectorB(nWayPoints + nInitial + i) = diffTerminal[i];
    }

    _controlPoints = _matrixA.fullPivLu().solve(_vectorB);

    _nWayPoints = nWayPoints;
    _nInitial = nInitial;
    _nTerminal = nTerminal;

    return true;
}

template<int nDegree>
double BSpline<nDegree>::operator()(double t, const vector<double> &tMap) {
    t = min(max(tMap[nDegree], t), tMap[_controlPoints.size()]);

    int l;
    for (l = nDegree; l < _controlPoints.size() + 1; ++l){
        if (tMap[l + 1] >= t){
            break;
        }
    }

    vector<double> p(nDegree + 1);
    for (int i = 0; i < nDegree + 1; ++i){
        p[i] = _controlPoints(l - nDegree + i);
    }

    for (int j = 1; j < nDegree + 1; ++j){
        for (int i = l; i > l - (nDegree + 1) + j; --i){
            double tor = (t - tMap[i]) / (tMap[i + nDegree + 1 - j] - tMap[i]);
            p[i - l + nDegree] = (1. - tor) * p[i - l + nDegree - 1] + tor * p[i - l + nDegree];
        }
    }

    return p[nDegree];
}

template<int nDegree>
void BSpline<nDegree>::deBoor(double t, const vector<double> &tMap, vector<double> &state) {
    t = min(max(tMap[nDegree], t), tMap[_controlPoints.size()]);

    int l;
    for (l = nDegree; l < _controlPoints.size() + 1; ++l){
        if (tMap[l + 1] >= t){
            break;
        }
    }

    vector<double> p((nDegree + 2) * (nDegree + 1) / 2);
    for (int j = 0; j < nDegree + 1; ++j){
        p[j] = _controlPoints(l - nDegree + j);
    }
    int indexLast = 0;
    for (int k = 1; k < nDegree + 1; ++k){
        int index = (2 * nDegree + 3 - k) * k / 2;
        int order = nDegree - k + 1;
        for (int i = 0; i < nDegree + 1 - k; ++i){
            p[index + i] = (p[indexLast + i + 1] - p[indexLast + i]) * order / (tMap[l + i + 1] - tMap[l - nDegree + i + k]);
        }
        indexLast = index;
    }

    vector<double> basisPartial(nDegree + 1, 0.);
    basisPartial[nDegree] = 1.;
    state[nDegree] = basisPartial[nDegree] * p[p.size() - 1];
    for (int i = 1; i <= nDegree; ++i){
        double tor = (tMap[l + 1] - t) / (tMap[l + 1] - tMap[l - i + 1]);
        basisPartial[nDegree - i] = basisPartial[nDegree - i + 1] * tor;
        int index = (int)p.size() - (i + 2) * (i + 1) / 2;
        state[nDegree - i] = basisPartial[nDegree - i] * p[index];
        for (int j = 1; j < i; ++j){
            double tor_ = 1. - tor;
            tor = (tMap[l + j + 1] - t) / (tMap[l + j + 1] - tMap[l - i + j + 1]);
            basisPartial[nDegree - i + j] = basisPartial[nDegree - i + j] * tor_ + basisPartial[nDegree - i + j + 1] * tor;
            state[nDegree - i] += basisPartial[nDegree - i + j] * p[index + j];
        }
        basisPartial[nDegree] = basisPartial[nDegree] * (1. - tor);
        state[nDegree - i] += basisPartial[nDegree] * p[index + i];
    }

}


#endif //DEMO5_BSPLINE_H
