//
// Created by Cain on 2022/4/15.
//

#ifndef DEMO2_ASTAR_H
#define DEMO2_ASTAR_H

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
#include "Eigen/Eigen"

namespace a_star {
    class Voxel {
    public:
        explicit Voxel() : _index(3) {}
        explicit Voxel(const std::vector<int> &index) : _index(index) {}
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
        std::vector<int> _index;
    };

    class VoxelHash {
    public:
        std::size_t operator()(const Voxel *voxel) const {
            return hashVal(voxel);
        }
        static std::size_t hashVal(const Voxel *&voxel) {
            std::size_t seed = 0;
            for (const auto& elem : voxel->_index){
                seed ^= std::hash<int>()(elem) + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    struct VoxelPred : public std::binary_function<Voxel*, Voxel*, bool> {
    public:
        bool operator()(const Voxel *lhs, const Voxel *rhs) const {
            for (int i = 0; i < lhs->_index.size(); ++i){
                if (lhs->_index[i] != rhs->_index[i]){
                    return false;
                }
            }
            return true;
        }
    };

    class Obstacle : public Voxel {
        using Voxel::Voxel;
    };

    class ObstacleHash {
    public:
        std::size_t operator()(const Obstacle &obstacle) const {
            return hashVal(obstacle);
        }
        static std::size_t hashVal(const Obstacle &obstacle) {
            std::size_t seed = 0;
            for (const auto& elem : obstacle._index){
                seed ^= std::hash<int>()(elem) + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    struct ObstaclePred : public std::binary_function<Obstacle&, Obstacle&, bool> {
    public:
        bool operator()(const Obstacle &lhs, const Obstacle &rhs) const {
            for (int i = 0; i < lhs._index.size(); ++i){
                if (lhs._index[i] != rhs._index[i]){
                    return false;
                }
            }
            return true;
        }
    };

    class Node : public Voxel {
    public:
        explicit Node(const std::vector<int> &index, int id=0, double g=Node::inf, double h=0., Node* parent= nullptr)
                : Voxel(index), _id(id), _g(g), _h(h), _f(g + h), _parent(parent) {}

        explicit Node(int id=0, double g=Node::inf, double h=0., Node* parent= nullptr)
                : Voxel(), _id(id), _g(g), _h(h), _f(g + h), _parent(parent) {}

        explicit Node(const Eigen::Vector3i &index, int id=0, double g=Node::inf, double h=0., Node* parent= nullptr)
                : Voxel(index), _id(id), _g(g), _h(h), _f(g + h), _parent(parent) {}

        explicit Node(int indexX, int indexY, int indexZ, int id=0, double g=Node::inf, double h=0., Node* parent= nullptr)
                : Voxel(indexX, indexY, indexZ), _id(id), _g(g), _h(h), _f(g + h), _parent(parent) {}

        int _id;        // 1: in openSet, -1: in closedSet
        Node* _parent;
        double _g, _h, _f;
        constexpr static const double inf = 100000.;
    };

    class NodeHash {
    public:
        std::size_t operator()(const Node *node) const {
            return hashVal(node);
        }
        static std::size_t hashVal(const Node *&node) {
            std::size_t seed = 0;
            for (const auto& elem : node->_index){
                seed ^= std::hash<int>()(elem) + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    struct NodePred : public std::binary_function<Node*, Node*, bool> {
    public:
        bool operator()(const Node *lhs, const Node *rhs) const {
            for (int i = 0; i < lhs->_index.size(); ++i){
                if (lhs->_index[i] != rhs->_index[i]){
                    return false;
                }
            }
            return true;
        }
    };

    class AStar {
    public:
        explicit AStar(double weightTieBreaker=.01, double weightHeuristic=1.1, int nodeSize=19997, int obstacleSize=6661)
                : _weightTieBreaker(weightTieBreaker), _weightHeuristic(weightHeuristic), _nodeSize(nodeSize), _nodeIndex(0),
                  _nodeList(nodeSize + 1), _nodeSet(nodeSize), _obstacleSet(obstacleSize), _indexStart(3), _indexGoal(3),
                  _boundary(3), _nodeStart(nullptr), _nodeGoal(nullptr) {};

        void initialize();
        void findShortestPath();
        void neighbours(Node *node, std::vector<Node*> &neis);

        void setBoundary(const std::vector<int>& boundary);
        void setBoundary(int boundaryX, int boundaryY, int boundaryZ);

        void setStart(const std::vector<int>& posIndex);
        void setStart(int indexX, int indexY, int indexZ);
        void setGoal(const std::vector<int>& posIndex);
        void setGoal(int indexX, int indexY, int indexZ);

        bool collisionFree(const std::vector<int>& posIndex) const;
        bool collisionFree(const int &indexX, const int &indexY, const int &indexZ) const;
        bool isGoal(const std::vector<int>& posIndex) const;

        void addObstacle(const std::vector<int>& posIndex);
        void addObstacle(int indexX, int indexY, int indexZ);
        void addObstacleBox(int indexX, int indexY, int indexZ, int sizeX, int sizeY, int sizeZ);
        void clearObstacleSet();
        void clearNodeSet();
        void clearNodeList();

        static double cost(const std::vector<int> &lhs, const std::vector<int> &rhs) ;
        static double cost(const Node *&lhs, const Node *&rhs) ;
        static double cost(Node *&lhs, Node *&rhs) ;

        double tieBreaker(const std::vector<int>& posIndex) const;
        double tieBreaker(const Node *&node) const;
        double tieBreaker(Node *&node) const;

        double heuristic(const std::vector<int> &posIndex) const;
        double heuristic(const Node *&node) const;
        double heuristic(Node *&node) const;

        // Maybe faster
        static double fastCost(const std::vector<int> &lhs, const std::vector<int> &rhs);
        static double fastCost(Node *&lhs, Node *&rhs);
        double fastHeuristic(const std::vector<int> &posIndex) const;
        double fastHeuristic(Node *&node) const;

        void openSetInsert(Node *&node);
        void openSetRemove(Node *&node);
        void openSetUpdate(Node *&node);
        void openSetUpdate(Node *&node, Node *&parent, double g);
        void openSetPop(Node *&node);
        Node *openSetPop();

        double _weightTieBreaker;
        double _weightHeuristic;

        Node *_nodeStart, *_nodeGoal;
        std::vector<int> _indexStart, _indexGoal;
        std::vector<int> _boundary;
        Eigen::Vector3d _unitStartToGoal;

        const int _nodeSize;
        int _nodeIndex;
        std::vector<Node> _nodeList;

        std::unordered_set<Node*, NodeHash, NodePred> _nodeSet;
        std::unordered_set<Obstacle, ObstacleHash, ObstaclePred> _obstacleSet;
        std::multimap<double, Node*> _openSet;

    };

    inline void AStar::setBoundary(const std::vector<int> &boundary) {
        for (int i = 0; i < 3; ++i){
            _boundary[i] = boundary[i];
        }
    }

    inline void AStar::setBoundary(int boundaryX, int boundaryY, int boundaryZ) {
        _boundary[0] = boundaryX;
        _boundary[1] = boundaryY;
        _boundary[2] = boundaryZ;
    }

    inline void AStar::setStart(const std::vector<int> &posIndex) {
        for (int i = 0; i < 3; ++i){
            _indexStart[i] = posIndex[i];
        }
    }

    inline void AStar::setStart(int indexX, int indexY, int indexZ) {
        _indexStart[0] = indexX;
        _indexStart[1] = indexY;
        _indexStart[2] = indexZ;
    }

    inline void AStar::setGoal(const std::vector<int> &posIndex) {
        for (int i = 0; i < 3; ++i){
            _indexGoal[i] = posIndex[i];
        }
    }

    inline void AStar::setGoal(int indexX, int indexY, int indexZ) {
        _indexGoal[0] = indexX;
        _indexGoal[1] = indexY;
        _indexGoal[2] = indexZ;
    }

    inline bool AStar::collisionFree(const std::vector<int> &posIndex) const {
        Obstacle obstacle(posIndex);
        return _obstacleSet.find(obstacle) == _obstacleSet.end();
    }

    inline bool AStar::collisionFree(const int &indexX, const int &indexY, const int &indexZ) const {
        Obstacle obstacle(indexX, indexY, indexZ);
        return _obstacleSet.find(obstacle) == _obstacleSet.end();
    }

    inline bool AStar::isGoal(const std::vector<int> &posIndex) const {
        for (int i = 0; i < 3; ++i){
            if (posIndex[i] != _indexGoal[i]){
                return false;
            }
        }
        return true;
    }

    inline void AStar::addObstacle(const std::vector<int> &posIndex) {
        _obstacleSet.emplace(posIndex);
    }

    inline void AStar::addObstacle(int indexX, int indexY, int indexZ) {
        _obstacleSet.emplace(indexX, indexY, indexZ);
    }

    inline void AStar::addObstacleBox(int indexX, int indexY, int indexZ, int sizeX, int sizeY, int sizeZ) {
        for (int i = 0; i < sizeX; ++i){
            for (int j = 0; j < sizeY; ++j){
                for (int k = 0; k < sizeZ; ++k){
                    _obstacleSet.emplace(indexX + i, indexY + j, indexZ + k);
                }
            }
        }
    }

    inline void AStar::clearObstacleSet() {
        _obstacleSet.clear();
    }

    inline void AStar::clearNodeSet() {
        _nodeSet.clear();
    }

    inline void AStar::clearNodeList() {
        for (int i = 0; i < _nodeIndex; ++i){
            _nodeList[i]._id = 0;
            _nodeList[i]._g = Node::inf;
            _nodeList[i]._h = 0.;
            _nodeList[i]._parent = nullptr;
        }
        _nodeIndex = 0;
    }

    inline double AStar::cost(const std::vector<int> &lhs, const std::vector<int> &rhs) {
        double diff, dist2 = 0.;
        for (int i = 0; i < 3; ++i) {
            diff = lhs[i] - rhs[i];
            dist2 += diff * diff;
        }
        return sqrt(dist2);
    }

    inline double AStar::cost(const Node *&lhs, const Node *&rhs) {
        return cost(lhs->_index, rhs->_index);
    }

    inline double AStar::cost(Node *&lhs, Node *&rhs) {
        return cost(lhs->_index, rhs->_index);
    }

    inline double AStar::tieBreaker(const std::vector<int> &posIndex) const {
        Eigen::Vector3d start2pos(posIndex[0] - _indexStart[0],
                                  posIndex[1] - _indexStart[1],
                                  posIndex[2] - _indexStart[2]);
        return _weightTieBreaker * start2pos.cross(_unitStartToGoal).norm();
    }

    inline double AStar::tieBreaker(const Node *&node) const {
        return tieBreaker(node->_index);
    }

    inline double AStar::tieBreaker(Node *&node) const {
        return tieBreaker(node->_index);
    }

    inline double AStar::heuristic(const std::vector<int> &posIndex) const {
        return _weightHeuristic * cost(posIndex, _indexGoal);
    }

    inline double AStar::heuristic(const Node *&node) const {
        return heuristic(node->_index);
    }

    inline double AStar::heuristic(Node *&node) const {
        return heuristic(node->_index);
    }

    inline double AStar::fastCost(const std::vector<int> &lhs, const std::vector<int> &rhs) {
        static double c[4] = {0., 1., 1.4142135623730951, 1.7320508075688772};
        int index = 0;
        for (int i = 0; i < 3; ++i){
            index += abs(lhs[i] - rhs[i]);
        }
        return c[index];
    }

    inline double AStar::fastCost(Node *&lhs, Node *&rhs) {
        return fastCost(lhs->_index, rhs->_index);
    }

    inline double AStar::fastHeuristic(const std::vector<int> &posIndex) const {
        static double c[3] = {1., 0.41421356237309515, 0.31783724519578205};

        int diff[3];
        for (int i = 0; i < 3; ++i){
            diff[i] = abs(posIndex[i] - _indexGoal[i]);
        }
        sort(diff, diff + 3, std::greater<>());

        double h = 0.;
        for (int i = 0; i < 3; ++i){
            h += c[i] * diff[i];
        }
        return h;
    }

    inline double AStar::fastHeuristic(Node *&node) const {
        return fastHeuristic(node->_index);
    }

    inline void AStar::openSetInsert(Node *&node) {
        node->_id = 1;
        _openSet.emplace(node->_g + node->_h, node);
    }

    inline void AStar::openSetRemove(Node *&node) {
        node->_id = -1;
        auto iterRange = _openSet.equal_range(node->_g + node->_h);
        for (auto it = iterRange.first; it != iterRange.second; ++it){
            if (it->second == node){
                _openSet.erase(it);
                break;
            }
        }
    }

    inline void AStar::openSetUpdate(Node *&node) {
        auto iterRange = _openSet.equal_range(node->_f);
        for (auto it = iterRange.first; it != iterRange.second; ++it){
            if (it->second == node){
                _openSet.erase(it);
                break;
            }
        }
        node->_h = heuristic(node);
        node->_f = node->_g + node->_h;
        _openSet.emplace(node->_f, node);
    }

    inline void AStar::openSetUpdate(Node *&node, Node *&parent, double g) {
        auto iterRange = _openSet.equal_range(node->_f);
        for (auto it = iterRange.first; it != iterRange.second; ++it){
            if (it->second == node){
                _openSet.erase(it);
                break;
            }
        }
        node->_h = heuristic(node);
        node->_g = g;
        node->_f = node->_g + node->_h;
        node->_parent = parent;
        _openSet.emplace(node->_f, node);
    }

    inline void AStar::openSetPop(Node *&node) {
        auto iter = _openSet.begin();
        auto top = iter->second;
        top->_id = -1;
        _openSet.erase(iter);
    }

    inline Node *AStar::openSetPop() {
        auto iter = _openSet.begin();
        auto top = iter->second;
        top->_id = -1;
        _openSet.erase(iter);
        return top;
    }
}


#endif //DEMO2_ASTAR_H
