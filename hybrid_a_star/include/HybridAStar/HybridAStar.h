//
// Created by Cain on 2022/4/15.
//

#ifndef DEMO4_HYBRIDASTAR_H
#define DEMO4_HYBRIDASTAR_H

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

using namespace std;

namespace hybrid_a_star {
    class Obstacle : public vector<int> {
    public:
        using vector<int>::vector;
        Obstacle(int x, int y, int z){
            resize(3);
            (*this)[0] = x, (*this)[1] = y, (*this)[2] = z;
        }
        bool operator==(const Obstacle& rhs) const {
            if (this->size() != rhs.size()){
                return false;
            }
            for (int i = 0; i < this->size(); ++i){
                if ((*this)[i] != rhs[i]){
                    return false;
                }
            }
            return true;
        }
    };

    class ObstacleHash{
    public:
        std::size_t operator()(const Obstacle &obstacle) const {
            return hashVal(obstacle);
        }

        static std::size_t hashVal(const Obstacle &obstacle) {
            std::size_t seed = 0;
            for (const auto& elem : obstacle){
                seed ^= hash<int>()(elem) + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    class Voxel {
    public:
        explicit Voxel() : _index(4) {}
        explicit Voxel(const vector<int> &index) : _index(index) {}
        explicit Voxel(const vector<int> &indexPos, int indexT) : _index(4) {
            _index[0] = indexPos[0];
            _index[1] = indexPos[1];
            _index[2] = indexPos[2];
            _index[3] = indexT;
        }
        explicit Voxel(const Eigen::Vector3i &indexPos, int indexT) : _index(4) {
            _index[0] = indexPos[0];
            _index[1] = indexPos[1];
            _index[2] = indexPos[2];
            _index[3] = indexT;
        }
        explicit Voxel(int indexX, int indexY, int indexZ, int indexT) : _index(4) {
            _index[0] = indexX;
            _index[1] = indexY;
            _index[2] = indexZ;
            _index[3] = indexT;
        }
        vector<int> _index;
    };

    class VoxelHash {
    public:
        std::size_t operator()(const Voxel *voxel) const {
            return hashVal(voxel);
        }
        static std::size_t hashVal(const Voxel *&voxel) {
            std::size_t seed = 0;
            for (const auto& elem : voxel->_index){
                seed ^= hash<int>()(elem) + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    struct VoxelPred : public binary_function<Voxel*, Voxel*, bool> {
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

    class Node : public Voxel{
    public:
        explicit Node(const vector<int> &index, int id=0, double t=0., double g=Node::inf, double h=0., Node* parent= nullptr)
                : Voxel(index), _t(t), _id(id), _g(g), _h(h), _f(g + h), _parent(parent) {}

        explicit Node(int id=0, double t=0., double g=Node::inf, double h=0., Node* parent= nullptr)
                : Voxel(), _t(t), _id(id), _g(g), _h(h), _f(g + h), _parent(parent) {}

        explicit Node(const vector<int> &indexPos, int indexT, int id=0, double t=0., double g=Node::inf, double h=0., Node* parent= nullptr)
                : Voxel(indexPos, indexT), _t(t), _id(id), _g(g), _h(h), _f(g + h), _parent(parent) {}

        explicit Node(const Eigen::Vector3i &indexPos, int indexT, int id=0, double t=0., double g=Node::inf, double h=0., Node* parent= nullptr)
                : Voxel(indexPos, indexT), _t(t), _id(id), _g(g), _h(h), _f(g + h), _parent(parent) {}

        explicit Node(int indexX, int indexY, int indexZ, int indexT, int id=0, double t=0., double g=Node::inf, double h=0., Node* parent= nullptr)
                : Voxel(indexX, indexY, indexZ, indexT), _t(t), _id(id), _g(g), _h(h), _f(g + h), _parent(parent) {}

        int _id;        // 1: in openSet, -1: in closedSet
        Node* _parent;
        Eigen::Vector3d _pos;
        Eigen::Vector3d _vel;
        Eigen::Vector3d _acc;
        double _t;
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
                seed ^= hash<int>()(elem) + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    struct NodePred : public binary_function<Node*, Node*, bool> {
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

/*
 * g[n] = rho * n + Σ( (x[k] - r)'*Q*(x[k] - r) + u[k]'*R*u[k] )
 * c[n] = g[n] - g[n-1] = rho + (x[n] - r)'*Q*(x[n] - r) + u[n]'*R*u[n]
 *
 * min{u[.], N} g[N],
 * s.t x[0] = x0, x[N] = r, x[k+1] = A*x[k] + B*u[k], x[k] not in obstacle set for k = 1, ..., N-1
 * */
    class HybridAStar {
    public:
        explicit HybridAStar(double dt=1., double velMax=7., double weightPos=1., double weightVel=0., double weightTime=10.,
                             double weightTieBreaker=1., double weightHeuristic=5., int numSimSteps=50, int nodeSize=19997, int obstacleSize=6661);

        void initialize();
        void findShortestPath();
        void neighbours(Node *node, set<Node*> &neis);

        void setBoundary(const vector<int>& boundary);
        void setBoundary(int boundaryX, int boundaryY, int boundaryZ);

        void setStart(const vector<int>& posIndex, const vector<double>& posCoord);
        void setStart(int indexX, int indexY, int indexZ, double posX, double posY, double posZ);
        void setGoal(const vector<int>& posIndex, const vector<double>& posCoord);
        void setGoal(int indexX, int indexY, int indexZ, double posX, double posY, double posZ);

        bool collisionFree(const vector<int>& posIndex) const;
        bool collisionFree(const int &indexX, const int &indexY, const int &indexZ) const;
        bool collisionFree(const Eigen::Vector3d &pos) const;
        bool isGoal(const vector<int>& posIndex) const;

        void addObstacle(const Obstacle &obstacle);
        void addObstacle(int indexX, int indexY, int indexZ);
        void addObstacleBox(int indexX, int indexY, int indexZ, int sizeX, int sizeY, int sizeZ);
        void clearObstacleSet();
        void clearNodeSet();
        void clearNodeList();

        double cost(const Node *&node) const;
        double cost(const Eigen::Vector3d &acc) const;
        double cost(const Eigen::Vector3d &acc, const Eigen::Vector3d &pos) const;
        double cost(const Eigen::Vector3d &acc, const Eigen::Vector3d &pos, const Eigen::Vector3d &vel) const ;

        double tieBreaker(const Eigen::Vector3d& pos) const;

        double heuristic(const Node *node) const;
        double heuristic(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel) const;

        void openSetInsert(Node *&node);
        void openSetRemove(Node *&node);
        void openSetUpdate(Node *&node);
        void openSetUpdate(Node *&node, Node *&parent, double g);
        void openSetPop(Node *&node);
        Node *openSetPop();

        double _dt;
        double _velMax;
        double _weightPos, _weightVel, _weightTime, _weightTieBreaker, _weightHeuristic;
        int _numSimSteps;

        Node *_nodeStart, *_nodeGoal;
        Eigen::Vector3i _indexStart, _indexGoal;
        Eigen::Vector3d _posStart, _posGoal;
        Eigen::Vector3i _boundary;
        Eigen::Vector3d _unitStartToGoal;

        const int _nodeSize;
        int _nodeIndex;
        vector<Node> _nodeList;

        unordered_set<Node*, NodeHash, NodePred> _nodeSet;
        unordered_set<Obstacle, ObstacleHash> _obstacleSet;
        multimap<double, Node*> _openSet;

        static const vector<vector<double>> _actions;
    };

    inline HybridAStar::HybridAStar(double dt, double velMax, double weightPos, double weightVel, double weightTime,
                                    double weightTieBreaker, double weightHeuristic, int numSimSteps, int nodeSize, int obstacleSize)
            : _dt(dt), _velMax(velMax), _weightPos(weightPos), _weightVel(weightVel), _weightTime(weightTime),
              _weightTieBreaker(weightTieBreaker), _weightHeuristic(weightHeuristic), _numSimSteps(numSimSteps),
              _nodeSize(nodeSize), _nodeIndex(0), _nodeList(nodeSize + 1), _nodeSet(nodeSize), _obstacleSet(obstacleSize),
              _nodeStart(nullptr), _nodeGoal(nullptr) {}

    inline void HybridAStar::setBoundary(const vector<int>& boundary) {
        for (int i = 0; i < 3; ++i){
            _boundary[i] = boundary[i];
        }
    }

    inline void HybridAStar::setBoundary(int boundaryX, int boundaryY, int boundaryZ) {
        _boundary[0] = boundaryX;
        _boundary[1] = boundaryY;
        _boundary[2] = boundaryZ;
    }

    inline void HybridAStar::setStart(const vector<int>& posIndex, const vector<double>& posCoord) {
        for (int i = 0; i < 3; ++i){
            _indexStart[i] = posIndex[i];
            _posStart[i] = posCoord[i];
        }
    }

    inline void HybridAStar::setStart(int indexX, int indexY, int indexZ, double posX, double posY, double posZ) {
        _indexStart[0] = indexX;
        _indexStart[1] = indexY;
        _indexStart[2] = indexZ;
        _posStart[0] = posX;
        _posStart[1] = posY;
        _posStart[2] = posZ;
    }

    inline void HybridAStar::setGoal(const vector<int>& posIndex, const vector<double>& posCoord) {
        for (int i = 0; i < 3; ++i){
            _indexGoal[i] = posIndex[i];
            _posGoal[i] = posCoord[i];
        }
    }

    inline void HybridAStar::setGoal(int indexX, int indexY, int indexZ, double posX, double posY, double posZ) {
        _indexGoal[0] = indexX;
        _indexGoal[1] = indexY;
        _indexGoal[2] = indexZ;
        _posGoal[0] = posX;
        _posGoal[1] = posY;
        _posGoal[2] = posZ;
    }

    inline bool HybridAStar::isGoal(const vector<int> &posIndex) const {
        for (int i = 0; i < 3; ++i){
            if (posIndex[i] != _indexGoal[i]){
                return false;
            }
        }
        return true;
    }

    inline bool HybridAStar::collisionFree(const vector<int>& posIndex) const {
        return _obstacleSet.find((const Obstacle&)posIndex) == _obstacleSet.end();
    }

    inline bool HybridAStar::collisionFree(const int &indexX, const int &indexY, const int &indexZ) const {
        return _obstacleSet.find((const Obstacle) {indexX, indexY, indexZ}) == _obstacleSet.end();
    }

    inline bool HybridAStar::collisionFree(const Eigen::Vector3d &pos) const {
        int indexPosForCheck[3][2], numChecks[3];
        for (int i = 0; i < 3; ++i){
            indexPosForCheck[i][0] = lround(pos[i]), numChecks[i] = 1;
            indexPosForCheck[i][1] = lround(pos[i] + 1e-1);
            if (indexPosForCheck[i][1] != indexPosForCheck[i][0]){
                if (indexPosForCheck[i][1] >= _boundary[i]){
                    return false;
                }
                ++numChecks[i];
            } else {
                indexPosForCheck[i][1] = lround(pos[i] - 1e-1);
                if (indexPosForCheck[i][1] != indexPosForCheck[i][0]){
                    if (indexPosForCheck[i][1] < 0){
                        return false;
                    }
                    ++numChecks[i];
                }
            }
            if (indexPosForCheck[i][0] < 0 || indexPosForCheck[i][0] >= _boundary[i]){
                return false;
            }
        }

        for (int ix = 0; ix < numChecks[0]; ++ix) {
            for (int iy = 0; iy < numChecks[1]; ++iy) {
                for (int iz = 0; iz < numChecks[2]; ++iz) {
                    if (!collisionFree(indexPosForCheck[0][ix], indexPosForCheck[1][iy], indexPosForCheck[2][iz])){
                        return false;
                    }
                }
            }
        }

        return true;
    }

    inline void HybridAStar::addObstacle(const Obstacle &obstacle) {
        _obstacleSet.emplace(obstacle);
    }

    inline void HybridAStar::addObstacle(int indexX, int indexY, int indexZ) {
        _obstacleSet.emplace(indexX, indexY, indexZ);
    }

    inline void HybridAStar::addObstacleBox(int indexX, int indexY, int indexZ, int sizeX, int sizeY, int sizeZ) {
        for (int i = 0; i < sizeX; ++i){
            for (int j = 0; j < sizeY; ++j){
                for (int k = 0; k < sizeZ; ++k){
                    _obstacleSet.emplace(indexX + i, indexY + j, indexZ + k);
                }
            }
        }
    }

    inline void HybridAStar::clearObstacleSet() {
        _obstacleSet.clear();
    }

    inline void HybridAStar::clearNodeSet() {
        _nodeSet.clear();
    }

    inline void HybridAStar::clearNodeList() {
        for (int i = 0; i < _nodeIndex; ++i){
            _nodeList[i]._id = 0;
            _nodeList[i]._t = 0.;
            _nodeList[i]._g = Node::inf;
            _nodeList[i]._h = 0.;
            _nodeList[i]._parent = nullptr;
        }
        _nodeIndex = 0;
    }

    inline double HybridAStar::cost(const Node *&node) const {
        return _weightTime + node->_acc.squaredNorm() + _weightPos * (node->_pos - _posGoal).squaredNorm() + _weightVel * node->_vel.squaredNorm();
    }

    inline double HybridAStar::cost(const Eigen::Vector3d &acc) const {
        return _weightTime + acc.squaredNorm();
    }

    inline double HybridAStar::cost(const Eigen::Vector3d &acc, const Eigen::Vector3d &pos) const {
        return _weightTime + acc.squaredNorm() + _weightPos * (pos - _posGoal).squaredNorm();
    }

    inline double HybridAStar::cost(const Eigen::Vector3d &acc, const Eigen::Vector3d &pos, const Eigen::Vector3d &vel) const {
        return _weightTime + acc.squaredNorm() + _weightPos * (pos - _posGoal).squaredNorm() + _weightVel * vel.squaredNorm();
    }

    inline double HybridAStar::tieBreaker(const Eigen::Vector3d &pos) const {
        return _weightTieBreaker * (pos - _posStart).cross(_unitStartToGoal).squaredNorm();
    }

    inline void HybridAStar::openSetInsert(Node *&node) {
        node->_id = 1;
        _openSet.emplace(node->_g + node->_h, node);
    }

    inline void HybridAStar::openSetRemove(Node *&node) {
        node->_id = -1;
        auto iterRange = _openSet.equal_range(node->_g + node->_h);
        for (auto it = iterRange.first; it != iterRange.second; ++it){
            if (it->second == node){
                _openSet.erase(it);
                break;
            }
        }
    }

    inline void HybridAStar::openSetUpdate(Node *&node) {
        auto iterRange = _openSet.equal_range(node->_g + node->_h);
        for (auto it = iterRange.first; it != iterRange.second; ++it){
            if (it->second == node){
                _openSet.erase(it);
                break;
            }
        }
        node->_h = heuristic(node);
        node->_f = node->_g + node->_h;
        _openSet.emplace(node->_g + node->_h, node);
    }

    inline void HybridAStar::openSetUpdate(Node *&node, Node *&parent, double g) {
        auto iterRange = _openSet.equal_range(node->_g + node->_h);
        for (auto it = iterRange.first; it != iterRange.second; ++it){
            if (it->second == node){
                _openSet.erase(it);
                break;
            }
        }
        node->_h = heuristic(node);
        node->_g = g;
        node->_parent = parent;
        _openSet.emplace(node->_g + node->_h, node);
    }

    inline Node *HybridAStar::openSetPop() {
        auto iter = _openSet.begin();
        auto top = iter->second;
        top->_id = -1;
        _openSet.erase(iter);
        return top;
    }

    inline void HybridAStar::openSetPop(Node *&node) {
        auto iter = _openSet.begin();
        auto top = iter->second;
        top->_id = -1;
        _openSet.erase(iter);
    }
}

#endif //DEMO4_HYBRIDASTAR_H
