//
// Created by Cain on 2022/5/26.
//

#include "esdf.h"

using namespace esdf;

ESDF::ESDF(int xBoundary, int yBoundary, int zBoundary)
: _xBoundary(xBoundary), _yBoundary(yBoundary), _zBoundary(zBoundary), _gridMap(xBoundary * yBoundary * zBoundary) {
    for (int i = 0; i < _xBoundary; ++i) {
        for (int j = 0; j < _yBoundary; ++j) {
            for (int k = 0; k < _zBoundary; ++k) {
                auto node = findFromIndex(i, j, k);
                node->_index[0] = i;
                node->_index[1] = j;
                node->_index[2] = k;
            }
        }
    }
}

Node* ESDF::findFromIndex(int &indexX, int &indexY, int &indexZ) {
    return &_gridMap[indexX + _xBoundary * indexY + _xBoundary * _yBoundary * indexZ];
}

void ESDF::neighbours(Node *&node, vector<Node*> &nbrs) {
    static const vector<vector<int>> offSets = {
            {-1, -1, -1},
            {-1, -1, 0},
            {-1, -1, 1},
            {-1, 0, -1},
            {-1, 0, 0},
            {-1, 0, 1},
            {-1, 1, -1},
            {-1, 1, 0},
            {-1, 1, 1},

            {0, -1, -1},
            {0, -1, 0},
            {0, -1, 1},
            {0, 0, -1},
            {0, 0, 1},
            {0, 1, -1},
            {0, 1, 0},
            {0, 1, 1},

            {1, -1, -1},
            {1, -1, 0},
            {1, -1, 1},
            {1, 0, -1},
            {1, 0, 0},
            {1, 0, 1},
            {1, 1, -1},
            {1, 1, 0},
            {1, 1, 1},
    };

    nbrs.clear();
    int indexX, indexY, indexZ;
    for (const auto &offSet : offSets) {
        indexX = node->_index[0] + offSet[0];
        if (indexX < 0 || indexX >= _xBoundary) {
            continue;
        }
        indexY = node->_index[1] + offSet[1];
        if (indexY < 0 || indexY >= _yBoundary) {
            continue;
        }
        indexZ = node->_index[2] + offSet[2];
        if (indexZ < 0 || indexZ >= _zBoundary) {
            continue;
        }
        nbrs.emplace_back(findFromIndex(indexX, indexY, indexZ));
    }
}

void ESDF::updateESDF() {
    vector<Node*> curNbrs(26);

    while (!_queue.empty()) {
        Node *cur;
        popFromQueue(cur);

        neighbours(cur, curNbrs);
        bool flag = false;
        for (auto & nbr : curNbrs) {
            if (nbr->_coc == nullptr) {     // nbr???coc???????????????
                continue;
            }
            double distTemp = dist(nbr->_coc, cur);
            if (distTemp < cur->_dis) {
                cur->_dis = distTemp;
                deleteFromDLL(cur->_coc, cur);  // cur->_coc????????????
                cur->_coc = nbr->_coc;
                insertIntoDLL(cur->_coc, cur);  // cur->_coc????????????
                flag = true;
            }
        }
        if (flag) {
            pushIntoQueue(cur);
            continue;
        }

        for (auto & nbr : curNbrs) {
            double distTemp = dist(cur->_coc, nbr);     // queue????????????coc????????????
            if (distTemp < nbr->_dis) {
                nbr->_dis = distTemp;
                if (nbr->_coc != nullptr) { // nbr->_coc???????????????
                    deleteFromDLL(nbr->_coc, nbr);
                }
                nbr->_coc = cur->_coc;
                insertIntoDLL(nbr->_coc, nbr);  // nbr->_coc????????????
                pushIntoQueue(nbr);
            }
        }
    }
}

void ESDF::addObstacle(int indexX, int indexY, int indexZ) {
    Node *cur = findFromIndex(indexX, indexY, indexZ);
    if (cur->_coc != nullptr) { // cur->_coc???????????????
        deleteFromDLL(cur->_coc, cur);
    }
    cur->_coc = cur;
    cur->_dis = 0.;
    insertIntoDLL(cur->_coc, cur);  // cur->_coc????????????
    pushIntoQueue(cur);
}

void ESDF::delObstacle(int indexX, int indexY, int indexZ) {
    vector<Node*> voxNbrs(26);

    Node *cur = findFromIndex(indexX, indexY, indexZ);
    Node *vox, *voxNext = cur->_head;
    while (voxNext != nullptr) {    // ??????????????????????????????, ???while?????????????????????, ??????cur->_head = nullptr
        /*
         * ??????insertIntoDLL?????????vox->_right, ????????????????????????vox->_right
         * */
        vox = voxNext;
        voxNext = vox->_right;

        deleteFromDLL(vox->_coc, vox);  // vox->_coc???cur, ??????vox->_coc????????????
        vox->_coc = nullptr;
        vox->_dis = Node::inf;
        neighbours(vox, voxNbrs);
        for (auto & voxNbr : voxNbrs) {
            if (voxNbr->_coc == nullptr || voxNbr->_coc == cur) {  // nbr???coc???????????????, ??????????????????voxNbr->_coc??????cur?????????
                continue;
            }
            double distTemp = dist(voxNbr->_coc, vox);
            if (distTemp < vox->_dis) {
                vox->_dis = distTemp;
                vox->_coc = voxNbr->_coc;
            }
        }
        if (vox->_coc != nullptr) {
            /*
             * insertIntoDLL?????????vox->_right
             * */
            insertIntoDLL(vox->_coc, vox);  // vox->_coc????????????
            pushIntoQueue(vox);
        }
    }
}

void ESDF::deleteFromDLL(Node *&parent, Node *&son) {
    /*
     * ?????????????????????????????????parent==nullptr?????????, ?????????????????????????????????
     * */
//    if (parent == nullptr) {
//        return;
//    }

    if (son->_left != nullptr) {
        son->_left->_right = son->_right;
    } else {
        parent->_head = son->_right;    // ???DLL????????????????????????, ?????????parent->_head???????????????nullptr, ????????????son->_right = son->_left = nullptr
    }
    if (son->_right != nullptr) {
        son->_right->_left = son->_left;
    }
    son->_left = nullptr;
    son->_right = nullptr;
}

void ESDF::insertIntoDLL(Node *&parent, Node *&son) {
    /*
     * ?????????????????????????????????insertIntoDLL???, parent????????????nullptr, ?????????????????????????????????
     * */
//    if (parent->_coc == nullptr) {
//        return;
//    }

    if (parent->_head == nullptr) {
        parent->_head = son;
    } else {
        parent->_head->_left = son;
        son->_right = parent->_head;
        /*
         * ??????insertIntoDLL??????????????????deleteFromDLL, ???????????????son->_left = nullptr, ??????????????????????????????
         * */
//        if (son->_left != nullptr) {
//            son->_left = nullptr;
//        }
        parent->_head = son;
    }
}

void ESDF::pushIntoQueue(Node *&node) {
    if (!node->_inQueue) {
        node->_inQueue = true;
        _queue.emplace(node);
    }
}

void ESDF::popFromQueue(Node *&node) {
    node = _queue.front();
    node->_inQueue = false;
    _queue.pop();
}

double ESDF::dist(Node *&lhs, Node *&rhs) {
    double diff, value = 0.;
    for (int i = 0; i < 3; ++i) {
        diff = lhs->_index[i] - rhs->_index[i];
        value += diff * diff;
    }
    return sqrt(value);
}