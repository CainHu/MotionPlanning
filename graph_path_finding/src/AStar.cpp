//
// Created by Cain on 2022/4/15.
//

#include "AStar.h"

using namespace std;
using namespace a_star;

void AStar::initialize() {
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

void AStar::findShortestPath() {
    multimap<double, Node*>::iterator iter;
    Node *top;
    vector<Node*> neis(26);
    double nei_g;

    while (!_openSet.empty()){
        // expand f-cost最小的point，并加入close list
        iter = _openSet.begin();        // i的内容会因为insert而产生变化
        top = iter->second;   // 所以需要先把i的GridNodePtr取出来
        top->_id = -1;
        _openSet.erase(iter);

        // 已经到了终点
        if (isGoal(top->_index)){
            _nodeGoal = top;
            cout << "Path found!" << endl;
            break;
        }

        neighbours(top, neis);
        for (Node *nei : neis){
            nei_g = top->_g + cost(top, nei);
            if (nei->_id == 0) {    // 新点, 直接更新
                nei->_id = 1;
                nei->_parent = top;
                nei->_g = nei_g;
                nei->_h = heuristic(nei) + tieBreaker(nei);
                nei->_f = nei->_g + nei->_h;
                _openSet.emplace(nei->_f, nei);
            } else if (nei_g < nei->_g){    // 代价更小才更新
                auto iterRange = _openSet.equal_range(nei->_f);
                for (auto it = iterRange.first; it != iterRange.second; ++it){
                    if (it->second == nei){
                        _openSet.erase(it);
                        break;
                    }
                }
                nei->_parent = top;
                nei->_g = nei_g;
                nei->_h = heuristic(nei) + tieBreaker(nei);
                nei->_f = nei->_g + nei->_h;
                _openSet.emplace(nei->_f, nei);
            }
        }
    }
}

void AStar::neighbours(Node *node, std::vector<Node *> &neis) {
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
                    if ((*hashIter)->_id == -1){    // 点在closeSet中
                        continue;
                    } else {
                        neighbor = *hashIter;
                    }
                }
                neis.emplace_back(neighbor);
            }
        }
    }
}