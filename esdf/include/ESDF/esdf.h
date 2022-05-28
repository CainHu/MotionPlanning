//
// Created by Cain on 2022/5/26.
//

#ifndef DEMO7_ESDF_H
#define DEMO7_ESDF_H

#include <iostream>
#include <map>
#include <cmath>
#include <algorithm>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <set>
#include <cstdio>
#include <queue>
#include "Eigen/Eigen"

namespace esdf {
    using namespace std;

    class Voxel {
    public:
        explicit Voxel() : _index(3) {}
        explicit Voxel(const vector<int> &index) : _index(index) {}
        explicit Voxel(const Eigen::Vector3i &indexPos) : _index(3) {
            _index[0] = indexPos[0];
            _index[1] = indexPos[1];
            _index[2] = indexPos[2];
        }
        explicit Voxel(int indexX, int indexY, int indexZ) : _index(3) {
            _index[0] = indexX;
            _index[1] = indexY;
            _index[2] = indexZ;
        }
        vector<int> _index;
    };

    class Node : public Voxel {
    public:
        explicit Node(double d=Node::inf, Node* coc=nullptr)
                : Voxel(), _dis(d), _coc(coc), _inQueue(false), _head(nullptr), _left(nullptr), _right(nullptr) {}

        explicit Node(const vector<int> &index, double d=Node::inf, Node* coc=nullptr)
                : Voxel(index), _dis(d), _coc(coc), _inQueue(false), _head(nullptr), _left(nullptr), _right(nullptr) {}

        explicit Node(const Eigen::Vector3i &index, double d=Node::inf, Node* coc=nullptr)
                : Voxel(index), _dis(d), _coc(coc), _inQueue(false), _head(nullptr), _left(nullptr), _right(nullptr) {}

        explicit Node(int indexX, int indexY, int indexZ, double d=Node::inf, Node* coc=nullptr)
                : Voxel(indexX, indexY, indexZ), _dis(d), _coc(coc), _inQueue(false), _head(nullptr), _left(nullptr), _right(nullptr) {}

        bool _inQueue;                  // 自身是否在已经在队列中
        Node *_coc;                     // 离自身最近的障碍物
        Node *_head, *_left, *_right;   // 双向链表所需要的指针
        double _dis;                    // 离自身最近的障碍物的距离
        constexpr static const double inf = 100000.;
    };

    class ESDF {
    public:
        ESDF(int xBoundary, int yBoundary, int zBoundary);
        Node* findFromIndex(int &indexX, int &indexY, int &indexZ);
        void neighbours(Node *&node, vector<Node*> &nbrs);
        void updateESDF();
        void addObstacle(int indexX, int indexY, int indexZ);
        void delObstacle(int indexX, int indexY, int indexZ);
        void deleteFromDLL(Node *&parent, Node *&son);
        void insertIntoDLL(Node *&parent, Node *&son);
        void pushIntoQueue(Node *&node);
        void popFromQueue(Node *&node);
        static double dist(Node *&lhs, Node *&rhs) ;

        int _xBoundary;
        int _yBoundary;
        int _zBoundary;
        vector<Node> _gridMap;
        queue<Node*> _queue;
    };
};

#endif //DEMO7_ESDF_H
