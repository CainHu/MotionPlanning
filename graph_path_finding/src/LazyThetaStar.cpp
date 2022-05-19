//
// Created by Cain on 2022/5/7.
//

#include "LazyThetaStar.h"

using namespace std;
using namespace lazy_theta_star;

void LazyThetaStar::initialize() {
    Eigen::Vector3d startToGoal(_indexGoal[0] - _indexStart[0],
                                _indexGoal[1] - _indexStart[1],
                                _indexGoal[2] - _indexStart[2]);
    _unitStartToGoal = startToGoal.normalized();

    _nodeStart = &_nodeList[_nodeIndex++];
    _nodeStart->_id = 1;
    _nodeStart->_index[0] = _indexStart[0];
    _nodeStart->_index[1] = _indexStart[1];
    _nodeStart->_index[2] = _indexStart[2];
    _nodeStart->_g = 0.;
    _nodeStart->_h = heuristic(_nodeStart);
    _nodeStart->_f = _nodeStart->_g + _nodeStart->_h;
    _nodeStart->_parent = _nodeStart;
    _nodeGoal = nullptr;

    _openSet.emplace(_nodeStart->_f , _nodeStart);
    _nodeSet.emplace(_nodeStart);
}

void LazyThetaStar::findShortestPath() {
    multimap<double, Node*>::iterator iter;
    Node *top;
    vector<Node*> neis(26);
    bool shouldFindNeis;
    double nei_g, nei_g_old;

    while (!_openSet.empty()){
        // expand f-cost最小的point，并加入close list
        iter = _openSet.begin();        // i的内容会因为insert而产生变化
        top = iter->second;   // 所以需要先把i的GridNodePtr取出来
        top->_id = -1;
        _openSet.erase(iter);

        if (!lineOfSight(top->_parent, top)) {
            shouldFindNeis = false;
            /* Path 1 */
            double top_g;
            top->_g = Node::inf;
            neighbours(top, neis);
            for (Node *nei : neis) {
                if (nei->_id == -1) {
                    top_g = nei->_g + cost(nei, top);
                    if (top_g < top->_g){
                        top->_g = top_g;
                        top->_parent = nei;
                    }
                }
            }
        } else {
            shouldFindNeis = true;
        }

        // 已经到了终点
        if (isGoal(top->_index)){
            _nodeGoal = top;
            cout << "Path found!" << endl;
            break;
        }

        if (shouldFindNeis) {
            neighbours(top, neis);
        }
        for (Node *nei : neis){
            if (nei->_id == -1) {  // 在closeSet, 不需要更新
                continue;
            }

            nei_g_old = nei->_g;

            /* Path 2 */
            nei_g = top->_parent->_g + cost(top->_parent, nei);
            if (nei_g < nei->_g) {
                nei->_parent = top->_parent;
                nei->_g = nei_g;
            }

//            if (nei->_g < nei_g_old) {  // 更新了g值
//                if (nei->_id == 0) {
//                    nei->_id = 1;
//                } else {
//                    auto iterRange = _openSet.equal_range(nei->_f);
//                    for (auto it = iterRange.first; it != iterRange.second; ++it){
//                        if (it->second == nei){
//                            _openSet.erase(it);
//                            break;
//                        }
//                    }
//                }
//                nei->_h = heuristic(nei);
//                nei->_f = nei->_g + nei->_h;
//                _openSet.emplace(nei->_f, nei);
//            }

            if (nei->_id == 0) {
                nei->_id = 1;
                nei->_h = heuristic(nei) + tieBreaker(nei);
                nei->_f = nei->_g + nei->_h;
                _openSet.emplace(nei->_f, nei);
            } else if (nei->_g < nei_g_old) {
                auto iterRange = _openSet.equal_range(nei->_f);
                for (auto it = iterRange.first; it != iterRange.second; ++it){
                    if (it->second == nei){
                        _openSet.erase(it);
                        break;
                    }
                }
                nei->_h = heuristic(nei) + tieBreaker(nei);
                nei->_f = nei->_g + nei->_h;
                _openSet.emplace(nei->_f, nei);
            }
        }
    }
}

void LazyThetaStar::neighbours(Node *node, std::vector<Node *> &neis) {
    neis.clear();

    vector<int> indexPos(3);
    Node *neighbor;
    for (int offsetX = -1; offsetX < 2; ++offsetX){
        indexPos[0] = node->_index[0] + offsetX;
        if (indexPos[0] < 0 || indexPos[0] >= _boundary[0]){   // 超过地图边界
            continue;
        }
        for (int offsetY = -1; offsetY < 2; ++offsetY){
            indexPos[1] = node->_index[1] + offsetY;
            if (indexPos[1] < 0 || indexPos[1] >= _boundary[1]){   // 超过地图边界
                continue;
            }
            for (int offsetZ = -1; offsetZ < 2; ++offsetZ){
                indexPos[2] = node->_index[2] + offsetZ;
                if (indexPos[2] < 0 || indexPos[2] >= _boundary[2]){   // 超过地图边界
                    continue;
                }

                if (offsetX == 0 && offsetY == 0 && offsetZ == 0){   // 同一个点
                    continue;
                }

                if (!collisionFree(indexPos)){  // 障碍物碰撞检测
                    continue;
                }

                _nodeList[_nodeIndex]._index[0] = indexPos[0];
                _nodeList[_nodeIndex]._index[1] = indexPos[1];
                _nodeList[_nodeIndex]._index[2] = indexPos[2];
                auto hashIter = _nodeSet.find(&_nodeList[_nodeIndex]);
                if (hashIter == _nodeSet.end()){     // 第一次发现这个点
                    if (_nodeIndex < _nodeSize){
                        neighbor = &_nodeList[_nodeIndex++];
                        _nodeSet.emplace(neighbor);
                    } else {
                        cout << "The nodeSet is full!" << endl;
                        continue;
                    }
                } else {            // 点已经在nodeSet中
                    neighbor = *hashIter;
                }
                neis.emplace_back(neighbor);
            }
        }
    }
}

bool LazyThetaStar::lineOfSight(const std::vector<int> &lhs, const std::vector<int> &rhs) const {
    // A Fast Voxel Traversal Algorithm for Ray Tracing
    // https://github.com/francisengelmann/fast_voxel_traversal
    vector<int> current_voxel = {lround(lhs[0]), lround(lhs[1]), lround(lhs[2])};
    if (!collisionFree(current_voxel)) {
        return false;
    }
    vector<int> last_voxel = {lround(rhs[0]), lround(rhs[1]), lround(rhs[2])};
    vector<int> ray = {rhs[0] - lhs[0], rhs[1] - lhs[1], rhs[2] - lhs[2]};

    double stepX = (ray[0] >= 0.) ? 1. : -1.;
    double stepY = (ray[1] >= 0.) ? 1. : -1.;
    double stepZ = (ray[2] >= 0.) ? 1. : -1.;

    double next_voxel_boundary_x = (current_voxel[0] + 0.5 * stepX);
    double next_voxel_boundary_y = (current_voxel[1] + 0.5 * stepY);
    double next_voxel_boundary_z = (current_voxel[2] + 0.5 * stepZ);

    double tMaxX = (ray[0] != 0.) ? (next_voxel_boundary_x - lhs[0]) / ray[0] : Node::inf;
    double tMaxY = (ray[1] != 0.) ? (next_voxel_boundary_y - lhs[1]) / ray[1] : Node::inf;
    double tMaxZ = (ray[2] != 0.) ? (next_voxel_boundary_z - lhs[2]) / ray[2] : Node::inf;

    double tDeltaX = (ray[0] != 0.) ? stepX / ray[0] : Node::inf;
    double tDeltaY = (ray[1] != 0.) ? stepY / ray[1] : Node::inf;
    double tDeltaZ = (ray[2] != 0.) ? stepZ / ray[2] : Node::inf;

    double tx = tMaxX, ty = tMaxY, tz = tMaxZ;
    while (last_voxel[0] != current_voxel[0] || last_voxel[1] != current_voxel[1] || last_voxel[2] != current_voxel[2]) {
        if (tx < ty) {
            if (tx < tz) {
                current_voxel[0] += stepX;
                tx += tDeltaX;
            } else {
                current_voxel[2] += stepZ;
                tz += tDeltaZ;
            }
        } else {
            if (ty < tz) {
                current_voxel[1] += stepY;
                ty += tDeltaY;
            } else {
                current_voxel[2] += stepZ;
                tz += tDeltaZ;
            }
        }
        if (!collisionFree(current_voxel)) {
            return false;
        }
    }
    return true;
}

bool LazyThetaStar::lineOfSight(Node *&lhs, Node *&rhs) const {
    return lineOfSight(lhs->_index, rhs->_index);
}