#include <iostream>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include "HybridAStar.h"
#include <ctime>

int isPrime(int n){
    int i;
    for(i = 2; i <= (int)sqrt(n); i ++)
        if(n%i == 0) return 0;
    return 1;
}

int main(){
    hybrid_a_star::HybridAStar hybridAStar;

    // Case 1
    hybridAStar.setBoundary(26, 26, 1);
    hybridAStar.setStart(25, 0, 0, 25., 0., 0.);
    hybridAStar.setGoal(0, 25, 0, 0., 25., 0.);
    hybridAStar.addObstacleBox(3, 12, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(7, 8, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(8, 13, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(13, 11, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(20, 9, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(19, 14, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(0, 0, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(2, 5, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(5, 0, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(4, 17, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(1, 22, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(10, 2, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(9, 18, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(7, 23, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(16, 1, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(15, 6, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(14, 16, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(13, 22, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(21, 4, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(18, 21, 0, 3, 3, 1);
    hybridAStar.addObstacleBox(23, 19, 0, 3, 3, 1);

//    // Case 2
//    hybridAStar.setBoundary(20, 28, 1);
//    hybridAStar.setStart(1, 25, 0, 1., 25., 0.);
//    hybridAStar.setGoal(18, 2, 0, 18., 2., 0.);
//    hybridAStar.addObstacleBox(4, 6, 0, 4, 22, 1);
//    hybridAStar.addObstacleBox(12, 0, 0, 4, 22, 1);

    clock_t t1 = clock();
    hybridAStar.initialize();
    hybridAStar.findShortestPath();
    clock_t t2 = clock();

    auto node = hybridAStar._nodeGoal;
    if (node == nullptr){
        cout << "Do not find the path" << endl;
    } else {
        while (node != hybridAStar._nodeStart){
            cout << "x = " << node->_index[0] << ", y = " << node->_index[1] << ", z = " << node->_index[2] << endl;
            printf("pos: %.2f, %.2f, %.2f, vel: %.2f, %.2f, %.2f, acc: %.2f, %.2f, %.2f\n", node->_pos(0), node->_pos(1), node->_pos(2),
                   node->_vel(0), node->_vel(1), node->_vel(2), node->_acc(0), node->_acc(1), node->_acc(2));
            node = node->_parent;
        }
        node = hybridAStar._nodeStart;
        cout << "x = " << node->_index[0] << ", y = " << node->_index[1] << ", z = " << node->_index[2] << endl;
        printf("pos: %.2f, %.2f, %.2f, vel: %.2f, %.2f, %.2f, acc: %.2f, %.2f, %.2f\n", node->_pos(0), node->_pos(1), node->_pos(2),
               node->_vel(0), node->_vel(1), node->_vel(2), node->_acc(0), node->_acc(1), node->_acc(2));
    }

    cout << "Explored nodes: " << hybridAStar._nodeSet.size() << endl;
    cout << "Time: " << t2 - t1 << endl;

    return 0;
}