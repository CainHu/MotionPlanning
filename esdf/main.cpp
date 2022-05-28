#include <iostream>
#include "esdf.h"
#include "Eigen/Eigen"
#include <cstdio>
#include <ctime>
#include <cmath>

using namespace std;

int main() {
    clock_t t1, t2;

    int xBoundary = 10, yBoundary = 10;
    esdf::ESDF esdfDemo(xBoundary, yBoundary, 1);
    esdfDemo.addObstacle(0, 0, 0);
    esdfDemo.addObstacle(4, 4, 0);
    esdfDemo.addObstacle(4, 5, 0);
    esdfDemo.addObstacle(5, 4, 0);
    esdfDemo.addObstacle(5, 5, 0);
    t1 = clock();
    esdfDemo.updateESDF();
    t2 = clock();
    std::cout << "time: " << t2 - t1 << std::endl;

    esdf::Node *nodePtr;

    std::cout << "ESDF" << std::endl;
    for (int i = 0; i < xBoundary; ++i) {
        for (int j = 0; j < yBoundary; ++j) {
            int z = 0;
            nodePtr = esdfDemo.findFromIndex(i, j, z);
            printf("%.2f, ", nodePtr->_dis);
        }
        std::cout << std::endl;
    }

    esdfDemo.addObstacle(9, 8, 0);
    t1 = clock();
    esdfDemo.updateESDF();
    t2 = clock();
    std::cout << "time: " << t2 - t1 << std::endl;

    std::cout << "ESDF" << std::endl;
    for (int i = 0; i < xBoundary; ++i) {
        for (int j = 0; j < yBoundary; ++j) {
            int z = 0;
            nodePtr = esdfDemo.findFromIndex(i, j, z);
            printf("%.2f, ", nodePtr->_dis);
        }
        std::cout << std::endl;
    }

    esdfDemo.addObstacle(1, 8, 0);
    esdfDemo.delObstacle(9, 8, 0);
    esdfDemo.delObstacle(4, 4, 0);
    t1 = clock();
    esdfDemo.updateESDF();
    t2 = clock();
    std::cout << "time: " << t2 - t1 << std::endl;

    std::cout << "ESDF" << std::endl;
    for (int i = 0; i < xBoundary; ++i) {
        for (int j = 0; j < yBoundary; ++j) {
            int z = 0;
            nodePtr = esdfDemo.findFromIndex(i, j, z);
            printf("%.2f, ", nodePtr->_dis);
        }
        std::cout << std::endl;
    }

    return 0;
}
