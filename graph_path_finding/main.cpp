#include <iostream>
#include <cstdio>
#include <fstream>

//// D* Lite
//#include "DStarLite.h"
//int main() {
//    DStarLite<2> dStarLite;
//    for (int i = 0; i < 10; ++i){
//        for (int j = 0; j < 10; ++j){
//            int index[2];
//            index[0] = i, index[1] = j;
//            dStarLite._nodes.emplace_back(new Node<2>(index));
//        }
//    }
//
////    dStarLite._boundary.emplace_back(10);
////    dStarLite._boundary.emplace_back(10);
//    dStarLite._boundary[0] = 10;
//    dStarLite._boundary[1] = 10;
//
//    for (int i = 0; i < 9; ++i){
//        dStarLite._nodes[10 * i + 4]->_isObstacle = true;
//    }
//
//    dStarLite._start = dStarLite._nodes[0];
//    dStarLite._goal = dStarLite._nodes[99];
//
//    dStarLite.main();
//
//    cout << numeric_limits<double>::max() * 2 << ", " << numeric_limits<double>::max() + 1 << endl;
//
//    for (int i = 0; i < 10; ++i){
//        for (int j = 0; j < 10; ++j){
//            printf("%d, ", dStarLite._nodes[10 * i + j]->_rhs);
////            cout << dStarLite._nodes[10 * i + j]->_rhs << ", ";
//        }
//        cout << endl;
//    }
//
//    return 0;
//}

#include <iostream>
#include <functional>
#include <ctime>
#include "AStar.h"
#include "LazyThetaStar.h"

using namespace std;

int girdMap[20][20] = {
        {0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
        {1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0},
        {1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1},
        {1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1},
        {1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1},
        {1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1},
        {0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0},
        {1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1},
        {1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1},
        {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1},
        {1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
        {1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1},
        {0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1},
        {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1},
        {0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1},
        {0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1},
        {0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1}
};


int main(){
    // Lazy Theta*
    cout << "Lazy Theta*" << endl;

    lazy_theta_star::LazyThetaStar lazyThetaStar;

    lazyThetaStar.setBoundary(20, 20, 1);
    lazyThetaStar.setStart(19, 0, 0);
    lazyThetaStar.setGoal(19, 18, 0);

    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20; ++j) {
            if (girdMap[i][j] == 1) {
                lazyThetaStar.addObstacle(i, j, 0);
            }
        }
    }

//    lazyThetaStar.setBoundary(26, 26, 1);
//    lazyThetaStar.setStart(0, 3, 0);
//    lazyThetaStar.setGoal(25, 0, 0);
//    lazyThetaStar.addObstacleBox(3, 12, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(7, 8, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(8, 13, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(13, 11, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(20, 9, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(19, 14, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(0, 0, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(2, 5, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(5, 0, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(4, 17, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(1, 22, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(10, 2, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(9, 18, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(7, 23, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(16, 1, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(15, 6, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(14, 16, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(13, 22, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(21, 4, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(18, 21, 0, 3, 3, 1);
//    lazyThetaStar.addObstacleBox(23, 19, 0, 3, 3, 1);

    clock_t t1 = clock();
    lazyThetaStar.initialize();
    lazyThetaStar.findShortestPath();
    clock_t t2 = clock();

    cout << "Explored nodes: " << lazyThetaStar._nodeSet.size() << endl;
    cout << "Time: " << t2 - t1 << endl;

    auto lazyThetaStarNode = lazyThetaStar._nodeGoal;
    if (lazyThetaStarNode == nullptr){
        cout << "Do not find the path" << endl;
    } else {
        while (lazyThetaStarNode != lazyThetaStar._nodeStart){
            cout << "x = " << lazyThetaStarNode->_index[0] << ", y = " << lazyThetaStarNode->_index[1] << ", z = " << lazyThetaStarNode->_index[2] << endl;
            lazyThetaStarNode = lazyThetaStarNode->_parent;
        }
        lazyThetaStarNode = lazyThetaStar._nodeStart;
        cout << "x = " << lazyThetaStarNode->_index[0] << ", y = " << lazyThetaStarNode->_index[1] << ", z = " << lazyThetaStarNode->_index[2] << endl;
    }

    // A*
    cout << "A*" << endl;

    a_star::AStar aStar;

    aStar.setBoundary(20, 20, 1);
    aStar.setStart(19, 0, 0);
    aStar.setGoal(19, 18, 0);

    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20; ++j) {
            if (girdMap[i][j] == 1) {
                aStar.addObstacle(i, j, 0);
            }
        }
    }

//    aStar.setBoundary(26, 26, 1);
//    aStar.setStart(0, 3, 0);
//    aStar.setGoal(25, 25, 0);
//    aStar.addObstacleBox(3, 12, 0, 3, 3, 1);
//    aStar.addObstacleBox(7, 8, 0, 3, 3, 1);
//    aStar.addObstacleBox(8, 13, 0, 3, 3, 1);
//    aStar.addObstacleBox(13, 11, 0, 3, 3, 1);
//    aStar.addObstacleBox(20, 9, 0, 3, 3, 1);
//    aStar.addObstacleBox(19, 14, 0, 3, 3, 1);
//    aStar.addObstacleBox(0, 0, 0, 3, 3, 1);
//    aStar.addObstacleBox(2, 5, 0, 3, 3, 1);
//    aStar.addObstacleBox(5, 0, 0, 3, 3, 1);
//    aStar.addObstacleBox(4, 17, 0, 3, 3, 1);
//    aStar.addObstacleBox(1, 22, 0, 3, 3, 1);
//    aStar.addObstacleBox(10, 2, 0, 3, 3, 1);
//    aStar.addObstacleBox(9, 18, 0, 3, 3, 1);
//    aStar.addObstacleBox(7, 23, 0, 3, 3, 1);
//    aStar.addObstacleBox(16, 1, 0, 3, 3, 1);
//    aStar.addObstacleBox(15, 6, 0, 3, 3, 1);
//    aStar.addObstacleBox(14, 16, 0, 3, 3, 1);
//    aStar.addObstacleBox(13, 22, 0, 3, 3, 1);
//    aStar.addObstacleBox(21, 4, 0, 3, 3, 1);
//    aStar.addObstacleBox(18, 21, 0, 3, 3, 1);
//    aStar.addObstacleBox(23, 19, 0, 3, 3, 1);

    t1 = clock();
    aStar.initialize();
    aStar.findShortestPath();
    t2 = clock();

    cout << "Explored nodes: " << aStar._nodeSet.size() << endl;
    cout << "Time: " << t2 - t1 << endl;

    std::ofstream csvfile;
    csvfile.open("F:\\GithubProject\\MotionPlanning\\graph_path_finding\\data.csv", std::ios::out);

    auto aStarNode = aStar._nodeGoal;
    if (aStarNode == nullptr){
        cout << "Do not find the path" << endl;
    } else {
        while (aStarNode != aStar._nodeStart){
            csvfile << aStarNode->_index[0] << "," << aStarNode->_index[1] << "," << aStarNode->_index[2] << endl;
            cout << "x = " << aStarNode->_index[0] << ", y = " << aStarNode->_index[1] << ", z = " << aStarNode->_index[2] << endl;
            aStarNode = aStarNode->_parent;
        }
        aStarNode = aStar._nodeStart;
        csvfile << aStarNode->_index[0] << "," << aStarNode->_index[1] << "," << aStarNode->_index[2] << endl;
        cout << "x = " << aStarNode->_index[0] << ", y = " << aStarNode->_index[1] << ", z = " << aStarNode->_index[2] << endl;
    }

    return 0;
}
