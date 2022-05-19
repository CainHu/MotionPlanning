//
// Created by Cain on 2022/4/15.
//

#include "HybridAStar.h"

using namespace hybrid_a_star;

const vector<vector<double>> HybridAStar::_actions = {
        {-2., -1., 0., 1., 2.},
        {-2., -1., 0., 1., 2.},
        {-2., -1., 0., 1., 2.}
};

void HybridAStar::initialize() {
    _unitStartToGoal = (_posGoal - _posStart).normalized();

    _nodeStart = &_nodeList[_nodeIndex++];
    _nodeStart->_id = 1;
    _nodeStart->_index[0] = _indexStart[0];
    _nodeStart->_index[1] = _indexStart[1];
    _nodeStart->_index[2] = _indexStart[2];
    _nodeStart->_index[3] = 0;
    _nodeStart->_pos[0] = _posStart[0];
    _nodeStart->_pos[1] = _posStart[1];
    _nodeStart->_pos[2] = _posStart[2];
    _nodeStart->_vel[0] = 0.;
    _nodeStart->_vel[1] = 0.;
    _nodeStart->_vel[2] = 0.;
    _nodeStart->_g = 0.;
    _nodeStart->_h = heuristic(_nodeStart);
    _nodeStart->_f = _nodeStart->_g + _nodeStart->_h;
    _nodeGoal = nullptr;

    _openSet.emplace(_nodeStart->_f , _nodeStart);
    _nodeSet.emplace(_nodeStart);
}

void HybridAStar::findShortestPath() {
    multimap<double, Node*>::iterator iter;
    Node *top;
    set<Node*> neis;

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
        for (auto &nei : neis){
            if (nei->_id == 0){
                nei->_id = 1;
            } else {
                auto iterRange = _openSet.equal_range(nei->_f);
                for (auto it = iterRange.first; it != iterRange.second; ++it){
                    if (it->second == nei){
                        _openSet.erase(it);
                        break;
                    }
                }
            }
            nei->_h = _weightHeuristic * heuristic(nei);
            nei->_f = nei->_g + nei->_h;
            _openSet.emplace(nei->_f, nei);
        }
    }
}

void HybridAStar::neighbours(Node *node, set<Node*> &neighbors) {
    neighbors.clear();

    double dt = _dt / _numSimSteps;
    Eigen::Vector3d pos, vel, acc;
    vector<int> indexPos(3);
    bool collision_free;
    Node *neighbor;
    double neighbor_g;
    for (auto &acc_x : _actions[0]){
        for (auto &acc_y : _actions[1]){
            for (auto &acc_z : _actions[2]){
                vel[0] = node->_vel[0] + acc_x * _dt;
                vel[1] = node->_vel[1] + acc_y * _dt;
                vel[2] = node->_vel[2] + acc_z * _dt;
                if (vel.squaredNorm() > _velMax * _velMax){
                    continue;
                }
                acc[0] = acc_x, acc[1] = acc_y, acc[2] = acc_z;
                pos[0] = node->_pos[0] + 0.5 * (vel[0] + node->_vel[0]) * _dt;
                pos[1] = node->_pos[1] + 0.5 * (vel[1] + node->_vel[1]) * _dt;
                pos[2] = node->_pos[2] + 0.5 * (vel[2] + node->_vel[2]) * _dt;

                // 超过地图边界
                indexPos[0] = lround(pos[0]);
                if (indexPos[0] < 0 || indexPos[0] >= _boundary[0]){
                    continue;
                }
                indexPos[1] = lround(pos[1]);
                if (indexPos[1] < 0 || indexPos[1] >= _boundary[1]){
                    continue;
                }
                indexPos[2] = lround(pos[2]);
                if (indexPos[2] < 0 || indexPos[2] >= _boundary[2]){
                    continue;
                }

                // 障碍物碰撞检测
                if (!collisionFree(indexPos)){
                    continue;
                }

                _nodeList[_nodeIndex]._index[0] = indexPos[0];
                _nodeList[_nodeIndex]._index[1] = indexPos[1];
                _nodeList[_nodeIndex]._index[2] = indexPos[2];
                _nodeList[_nodeIndex]._index[3] = node->_index[3] + 1;
                auto hashIter = _nodeSet.find(&_nodeList[_nodeIndex]);
                bool isNewNode = hashIter == _nodeSet.end();
                if (isNewNode){     // 第一次发现这个点
                    if (_nodeIndex < _nodeSize){
                        neighbor = &_nodeList[_nodeIndex];
                        neighbor_g = cost(acc, pos) + tieBreaker(pos) + node->_g;
                    } else {
                        cout << "nodeSet 已经满了" << endl;
                        continue;
                    }
                } else {            // 点已经在nodeSet中
                    neighbor = *hashIter;
                    if (neighbor->_id != -1){   // 在closeSet中
                        neighbor_g = cost(acc, pos) + tieBreaker(pos) + node->_g;
                        if (neighbor_g >= neighbor->_g){    // 新的代价没有比旧的小
                            continue;
                        }
                    } else {                    // 在closeSet中
                        continue;
                    }
                }

                // 延迟的碰撞与越界检测，只有点不在nodeSet中 或 点在nodeSet中但不在closeSet中且代价更小, 才能执行到此处
                collision_free = true;
                Eigen::Vector3d pos_(node->_pos);
                Eigen::Vector3d vel_(node->_vel);
                for (int step = 0; step < _numSimSteps - 1; ++step){    // 最后一个时刻已经检查过了
                    pos_[0] += (0.5 * acc[0] * dt + vel_[0]) * dt;
                    pos_[1] += (0.5 * acc[1] * dt + vel_[1]) * dt;
                    pos_[2] += (0.5 * acc[2] * dt + vel_[2]) * dt;
                    vel_[0] += acc[0] * dt;
                    vel_[1] += acc[1] * dt;
                    vel_[2] += acc[2] * dt;

                    // 超过地图边界
                    indexPos[0] = lround(pos_[0]);
                    if (indexPos[0] < 0 || indexPos[0] >= _boundary[0]){
                        collision_free = false;
                        break;
                    }
                    indexPos[1] = lround(pos_[1]);
                    if (indexPos[1] < 0 || indexPos[1] >= _boundary[1]){
                        collision_free = false;
                        break;
                    }
                    indexPos[2] = lround(pos_[2]);
                    if (indexPos[2] < 0 || indexPos[2] >= _boundary[2]){
                        collision_free = false;
                        break;
                    }

                    // 障碍物碰撞检测
                    collision_free = collisionFree(indexPos);
                    if (!collision_free){
                        break;
                    }
                }

                if (collision_free) {
                    if (isNewNode){
                        _nodeSet.emplace(neighbor);
                        ++_nodeIndex;
                    }
                    neighbor->_parent = node;
                    neighbor->_pos[0] = pos[0];
                    neighbor->_pos[1] = pos[1];
                    neighbor->_pos[2] = pos[2];
                    neighbor->_vel[0] = vel[0];
                    neighbor->_vel[1] = vel[1];
                    neighbor->_vel[2] = vel[2];
                    neighbor->_acc[0] = acc[0];
                    neighbor->_acc[1] = acc[1];
                    neighbor->_acc[2] = acc[2];
                    neighbor->_g = neighbor_g;
                    neighbors.emplace(neighbor);
                }
            }
        }
    }
}

/*
 * 每个维度均满足(省略下标):
 * Δp = pf - p0
 * Δv = vf - v0
 * u = a*t + b
 * T >= 0
 *
 * 上式等价于匀jerk运动， a为jerk, b为0时刻的acceleration
 * 所以
 * Δp = 1/6*a*T^3 + 1/2*b*T^2 + v0*T
 * Δv = 1/2*a*T^2 + b*T
 * 解出：
 * a = (6*(-2*dp + T*(Δv + 2*v0)))/T^3
 * b = (6*Δp - 2*T*(Δv + 3*v0))/T^2
 *
 * 代价函数:
 * j = ∫(1/3 + u^2)dt
 *   = T/3 + b^2*T + a*b*T^2 + (a^2*T^3)/3
 *
 * J = Σj
 *   = T + Σ(b^2*T + a*b*T^2 + (a^2*T^3)/3)
 *   = T + Σ(4*(3*Δp^2 - 3*Δp*T*(Δv + 2*v0) + T^2*(Δv^2 + 3*Δv*v0 + 3*v0^2)))/T^3
 *
 * J对T求导有:
 * dJ/dT = 1 - Σ((4*(9*Δp^2 - 6*Δp*T*(Δv + 2*v0) + T^2*(Δv^2 + 3*Δv*v0 + 3*v0^2)))/T^4)
 * 令:
 * dJ/dT = 0
 * 可得:
 * T^4 - Σ4*(9*Δp^2 - 6*Δp*T*(Δv + 2*v0) + T^2*(Δv^2 + 3*Δv*v0 + 3*v0^2)) = 0
 * => T^4 + 0*T^3 - Σ4*(Δv^2 + 3*Δv*v0 + 3*v0^2) * T^2 + Σ24*Δp*(Δv + 2*v0)*T - Σ36*Δp^2 = 0
 *    T^4 + c[3] * T^3 + c[2] * T^2 + c[1] * T + c[0] = 0
 *
 * 伴随矩阵:
 * Adj = [ 0 0 0 -c[0] ]
 *       [ 1 0 0 -c[1] ]
 *       [ 0 1 0 -c[2] ]
 *       [ 0 0 1 -c[3] ]
 * 特征方程:
 * |λI - Adj| = λ^4 + c[3] * λ^3 + c[2] * λ^2 + c[1] * λ + c[0]
 * 所以求T <=> 求Adj的实部大于0且虚部为0的根
 * */
double HybridAStar::heuristic(const Eigen::Vector3d &pos, const Eigen::Vector3d &vel) const {
    double optimal_cost = Node::inf;

    Eigen::Vector3d dp, dv;
    for (int it = 0; it < 3; ++it){
        dp(it) = _posGoal(it) - pos(it);
        dv(it) = 0. - vel(it);
    }

    vector<double> c(4);
    for (int it = 0; it < 3; ++it){
        c[0] += -36. * dp(it) * dp(it);
        c[1] += 24. * dp(it) * (dv(it) + 2. * vel(it));
        c[2] += -4. * (dv(it) * dv(it) + 3. * dv(it) * vel(it) + 3. * vel(it) * vel(it));
    }

    Eigen::Matrix<double, 4, 4> adj = Eigen::Matrix<double, 4, 4>::Zero();
    for (int it = 0; it < 3; ++it){
        adj(it + 1, it) = 1.;
    }
    for (int it = 0; it < 4; ++it){
        adj(it, 3) = -c[it];
    }

    Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> adj_eigenvalues = adj.eigenvalues();
    double T = -1.;
    for (int it = 0; it < 4; ++it){
        if (adj_eigenvalues(it).imag() == 0.){
            if (adj_eigenvalues(it).real() > 0.){
                T = adj_eigenvalues(it).real();
            }
        }
    }
    if (T > 0.){
        double T2 = T * T, T3 = T * T2;
        double J = T - c[0] / 3. / T3 - c[1] / 2. / T2 - c[2] / T;
        optimal_cost = J < optimal_cost ? J : optimal_cost;
    }

    return optimal_cost;
}

double HybridAStar::heuristic(const Node *node) const {
    return heuristic(node->_pos, node->_vel);
}
