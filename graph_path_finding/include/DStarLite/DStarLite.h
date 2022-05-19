//
// Created by Cain on 2022/4/12.
//

#ifndef DEMO2_DSTARLITE_H
#define DEMO2_DSTARLITE_H


#include <iostream>
#include <map>
#include <algorithm>
#include <cmath>
#include <vector>
#include <functional>

using namespace std;

template<int dim>
class Node {
public:
    Node(const int index[dim], bool obstacle=false, int g=inf, int rhs=inf, int h=0, Node *parent=nullptr)
            : _inOpenSet(false), _isObstacle(obstacle), _kmLast(-1), _g(g), _rhs(rhs), _h(h), _parent(parent), _index(dim), _coor(dim){
        for (int i = 0; i < dim; ++i){
            _index[i] = index[i];
            _coor[i] = index[i];
        }
    }
    Node(const int index[dim], const double coor[dim], bool obstacle=false, int g=inf, int rhs=inf, int h=0, Node *parent=nullptr)
            : _inOpenSet(false), _isObstacle(obstacle), _kmLast(-1), _g(g), _rhs(rhs), _h(h), _parent(parent), _index(dim), _coor(dim){
        for (int i = 0; i < dim; ++i){
            _index[i] = index[i];
            _coor[i] = coor[i];
        }
    }

    bool _inOpenSet;
    bool _isObstacle;
    vector<int> _index;
    vector<double> _coor;
    int _g;          // last rhs
    int _rhs;        // new rhs
    int _h;          // heuristic
    int _kmLast;
    pair<int, int> _key;
    Node *_parent;
    constexpr static const int inf = 1000000;

    static int heuristic(Node* node1, Node* node2);
};

template<int dim>
inline int Node<dim>::heuristic(Node *node1, Node *node2) {
    static int c[3] = {10, 4, 3};

    int diff[dim];
    for (int i = 0; i < dim; ++i){
        diff[i] = abs(node1->_index[i] - node2->_index[i]);
    }
    sort(diff, diff + dim, greater<>());

    int h = 0;
    for (int i = 0; i < dim; ++i){
        h += c[i] * diff[i];
    }
    return h;
}


template<int dim>
class DStarLite {
public:
    DStarLite() : _boundary(dim), _km(-1) {}
    pair<int, int> calculateKey(Node<dim>* s) const;
    void initialize();
    void updateVertex(Node<dim> *u);
    void computeShortestPath();
    void main();
    void predecessor(Node<dim> *node, vector<Node<dim>*> &pred);
    void successor(Node<dim> *node, vector<Node<dim>*> &succ);
    void neighbours(Node<dim> *node, vector<Node<dim>*> &nei);
    void insert(Node<dim> *node);
    void remove(Node<dim> *node);
    void update(Node<dim> *node, pair<int, int> key);
    static int cost(Node<dim> *s1, Node<dim> *s2);
    bool step(Node<dim>* start);
    bool step();
    void updateMap(vector<Node<dim>*> &changes);

    vector<int> _boundary;
    Node<dim> *_start, *_goal, *_last;
    int _km;
    vector<Node<dim>*> _nodes;
    multimap<pair<int, int>, Node<dim>*> _openSet;
};

template<int dim>
void DStarLite<dim>::initialize() {
    _km = 0;
    _last = _start;
    for (auto s : _nodes){
        s->_rhs = Node<dim>::inf;
        s->_g = s->_rhs;
    }
    _goal->_rhs = 0;
    _goal->_key = calculateKey(_goal);
    _openSet.emplace(pair<pair<int, int>, Node<dim>*>(_goal->_key, _goal));
}

template<int dim>
inline pair<int, int> DStarLite<dim>::calculateKey(Node<dim> *s) const {
    if (s->_kmLast != _km){
        s->_kmLast = _km;
        s->_h = Node<dim>::heuristic(_start, s);
    }
    int second = min(s->_g, s->_rhs);
    return {second + s->_h + _km, second};
}

template<int dim>
inline void DStarLite<dim>::updateVertex(Node<dim> *u) {
    if (u->_g != u->_rhs){
        if (u->_inOpenSet){
            update(u, calculateKey(u));
        } else {
            insert(u);
        }
    } else if (u->_inOpenSet){
        remove(u);
    }
}

template<int dim>
void DStarLite<dim>::computeShortestPath() {
    int vector_size = pow(3, dim) - 1;
    vector<Node<dim>*> preds(vector_size), succs(vector_size);
    while (_openSet.begin()->first < calculateKey(_start) || _start->_rhs > _start->_g){
        auto top = _openSet.begin()->second;
        auto kOld = _openSet.begin()->first;
        auto kNew = calculateKey(top);

        if (kOld < kNew){
            update(top, kNew);      // 若top的heuristic发生了变化时，则需要先update
        } else if (top->_g > top->_rhs){    // 过估计，检查能够通过让neighbours的parent为top而使代价更低
            top->_g = top->_rhs;
            remove(top);
            predecessor(top, preds);
            for (auto pred : preds){
                if (pred != _goal){
                    int rhs = cost(pred, top) + top->_g;
                    if (rhs < pred->_rhs){
                        pred->_rhs = rhs;
                        pred->_parent = top;
                    }
                }
                updateVertex(pred);
            }
        } else {    // 为top的neighbours重新选parent
            top->_g = Node<dim>::inf;
            predecessor(top, preds);
            for (auto pred : preds){
                if (pred->_parent == top){
                    if (pred != _goal){
                        pred->_rhs = Node<dim>::inf;
                        pred->_parent = nullptr;
                        successor(pred, succs);
                        for (auto succ : succs){
                            auto s_rhs = cost(pred, succ) + succ->_g;
                            if (s_rhs < pred->_rhs){
                                pred->_rhs = s_rhs;
                                pred->_parent = succ;
                            }
                        }
                    }
                }
                updateVertex(pred);
            }
        }
    }
}

template<int dim>
void DStarLite<dim>::main() {
    int step = 0;
    cout << _start->_index[0] << ", " << _start->_index[1] << endl;

    int vector_size = pow(3, dim) - 1;
    vector<Node<dim>*> preds(vector_size), succs(vector_size);
    vector<Node<dim>*> changes;
    Node<dim> *last = _start;
    initialize();
    computeShortestPath();
    while (_start != _goal){
        _start = _start->_parent;
        cout << _start->_index[0] << ", " << _start->_index[1] << endl;

        if (++step == 6){
            _nodes[2 * 10 + 4]->_isObstacle = false;
            _nodes[9 * 10 + 4]->_isObstacle = true;

            changes.emplace_back(_nodes[2 * 10 + 4]);
            changes.emplace_back(_nodes[9 * 10 + 4]);
        }

        if (!changes.empty()){
            _km += distance(last, _start);
            last = _start;

            for (auto change : changes){
                if (change->_isObstacle){   // 如果change的neighbours的parent为change, 需要为change的neighbours重新选parent
                    change->_rhs = Node<dim>::inf;
                    predecessor(change, preds);
                    for (auto pred : preds){
                        if (pred->_parent == change){
                            pred->_rhs = Node<dim>::inf;
                            pred->_parent = nullptr;
                            successor(pred, succs);
                            for (auto succ : succs){
                                auto u_rhs = cost(pred, succ) + succ->_g;
                                if (u_rhs < pred->_rhs){
                                    pred->_rhs = u_rhs;
                                    pred->_parent = succ;
                                }
                            }
                            updateVertex(pred);
                        }
                    }
                } else {    // 为change选代价最小的parent
                    change->_rhs = Node<dim>::inf;
                    change->_parent = nullptr;
                    successor(change, succs);
                    for (auto succ : succs){
                        auto changeNode_rhs = cost(change, succ) + succ->_g;
                        if (changeNode_rhs < change->_rhs){
                            change->_rhs = changeNode_rhs;
                            change->_parent = succ;
                        }
                    }
                }
                updateVertex(change);
            }
            computeShortestPath();
            changes.clear();
        }
    }
}

template<int dim>
inline int DStarLite<dim>::cost(Node<dim> *s1, Node<dim> *s2) {
    static int c[4] = {0, 10, 14, 17};
    if (s1->_isObstacle || s2->_isObstacle){
        return Node<dim>::inf;
    }
    int diff = 0;
    for (int i = 0; i < dim; ++i){
        diff += abs(s1->_index[i] - s2->_index[i]);
    }
    return c[diff];
}

template<int dim>
inline void DStarLite<dim>::predecessor(Node<dim> *node, vector<Node<dim>*> &pred) {
    neighbours(node, pred);
}

template<int dim>
inline void DStarLite<dim>::successor(Node<dim> *node, vector<Node<dim>*> &succ) {
    neighbours(node, succ);
}

template<int dim>
void DStarLite<dim>::neighbours(Node<dim> *node, vector<Node<dim>*> &nei) {
    nei.clear();

    function<void(int, int)> cycle = [&] (int index, int depth) {
        int index_next, depth_next = depth + 1;
        for (int i = -1; i < 2; ++i){
            int p = node->_index[depth] + i;
            index_next = index + p;
            if (p < 0 || p >= _boundary[depth]){
                continue;
            }
            if (depth_next == dim){
                auto neighbour = _nodes[index_next];
                if (neighbour->_isObstacle || neighbour == node){
                    continue;
                }
                nei.emplace_back(neighbour);
                continue;
            }
            cycle(_boundary[depth_next] * index_next, depth_next);
        }
    };

    cycle(0, 0);
}

template<int dim>
inline void DStarLite<dim>::insert(Node<dim> *node) {
    node->_inOpenSet = true;
    node->_key = calculateKey(node);
    _openSet.emplace(pair<pair<int, int>, Node<dim>*>(node->_key, node));
}

template<int dim>
void DStarLite<dim>::remove(Node<dim> *node) {
    auto iterRange = _openSet.equal_range(node->_key);
    for (auto iter = iterRange.first; iter != iterRange.second; ++iter){
        if (iter->second == node){
            _openSet.erase(iter);
            break;
        }
    }
    node->_inOpenSet = false;
}

template<int dim>
void DStarLite<dim>::update(Node<dim> *node, pair<int, int> key) {
    auto iterRange = _openSet.equal_range(node->_key);
    for (auto iter = iterRange.first; iter != iterRange.second; ++iter){
        if (iter->second == node){
            _openSet.erase(iter);
            break;
        }
    }
    node->_key = key;
    _openSet.emplace(pair<pair<int, int>, Node<dim>*>(key, node));
}

template<int dim>
bool DStarLite<dim>::step(Node<dim>* start) {
    if (_start == nullptr){
        return false;
    }
    if (start == nullptr){
        return false;
    }
    _start = start;
    return true;
}

template<int dim>
bool DStarLite<dim>::step() {
    if (_start == nullptr){
        return false;
    }
    Node<dim> *start = _start->_parent;
    if (start == nullptr){
        return false;
    }
    _start = start;
    return true;
}

template<int dim>
void DStarLite<dim>::updateMap(vector<Node<dim> *> &changes) {
    _km += distance(_last, _start);
    _last = _start;

    int vector_size = pow(3, dim) - 1;
    vector<Node<dim>*> preds(vector_size), succs(vector_size);
    for (auto change : changes){
        if (change->_isObstacle){   // 如果change的neighbours的parent为change, 需要为change的neighbours重新选parent
            change->_rhs = Node<dim>::inf;
            predecessor(change, preds);
            for (auto pred : preds){
                if (pred->_parent == change){
                    pred->_rhs = Node<dim>::inf;
                    pred->_parent = nullptr;
                    successor(pred, succs);
                    for (auto succ : succs){
                        auto u_rhs = cost(pred, succ) + succ->_g;
                        if (u_rhs < pred->_rhs){
                            pred->_rhs = u_rhs;
                            pred->_parent = succ;
                        }
                    }
                    updateVertex(pred);
                }
            }
        } else {    // 为change选代价最小的parent
            change->_rhs = Node<dim>::inf;
            change->_parent = nullptr;
            successor(change, succs);
            for (auto succ : succs){
                auto changeNode_rhs = cost(change, succ) + succ->_g;
                if (changeNode_rhs < change->_rhs){
                    change->_rhs = changeNode_rhs;
                    change->_parent = succ;
                }
            }
        }
        updateVertex(change);
    }
    computeShortestPath();
}


#endif //DEMO2_DSTARLITE_H
