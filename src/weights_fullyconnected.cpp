#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <tuple>
#include <unordered_set>
#include <vector>
#include <queue>
#include <map>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <random>
#include <chrono>
using namespace std;

map<pair<int, int>, double> calc_weights(int nodecount2, const std::vector<int> &nodes2, int n, int N, double remcap, int Delta)
{
    set<int> orderednodes2;
    for (int ele : nodes2)
    {
        orderednodes2.insert(ele);
    }
    map<pair<int, int>, double> addcost;
    set<int> checkednodes;
    int bestcount = 0;
    int current;
    int neighbour;
    vector<int> visited(N, 0);
    int bestsizednode = 0;

    for (int i = 0; i < nodecount2; i++)
    {
        if (visited[nodes2[i]] == 0)
        {
            int checkedcount = 1;
            queue<int> myqueue;
            myqueue.push(nodes2[i]);
            visited[nodes2[i]] = i + 1;
            while (!myqueue.empty())
            {
                current = myqueue.front();
                myqueue.pop();
                neighbour = current + 1;
                if (orderednodes2.count(neighbour) > 0 and neighbour % n != 0 and visited[neighbour] == 0)
                {
                    myqueue.push(neighbour);
                    visited[neighbour] = i + 1;
                    checkedcount++;
                }
                neighbour = current - 1;
                if (orderednodes2.count(neighbour) > 0 and current % n != 0 and visited[neighbour] == 0)
                {
                    myqueue.push(neighbour);
                    visited[neighbour] = i + 1;
                    checkedcount++;
                }
                neighbour = current + n;
                if (orderednodes2.count(neighbour) > 0 and neighbour < N and visited[neighbour] == 0)
                {
                    myqueue.push(neighbour);
                    visited[neighbour] = i + 1;
                    checkedcount++;
                }
                neighbour = current - n;
                if (orderednodes2.count(neighbour) > 0 and neighbour >= 0 and visited[neighbour] == 0)
                {
                    myqueue.push(neighbour);
                    visited[neighbour] = i + 1;
                    checkedcount++;
                }
            }
            if (checkedcount > bestcount)
            {
                bestcount = checkedcount;
                bestsizednode = i;
            }
        }
    }
    int nodecount = bestcount;
    vector<int> nodes;
    for (int i = 0; i < N; i++)
    {
        if (visited[i] == bestsizednode + 1)
        {
            nodes.push_back(i);
        }
    }
    if (0.1 * nodecount > remcap or 0.9 * nodecount < remcap)
    {
        addcost[make_pair(-1, -1)] = -1;
    }
    else
    {
        addcost[make_pair(-1, -1)] = nodecount;
    }

    set<int> orderednodes;
    for (int ele : nodes)
    {
        orderednodes.insert(ele);
    }

    for (int i = 0; i < nodecount; i++)
    {
        if (orderednodes.count(nodes[i] + 1) > 0 and (nodes[i] + 1) % n > 0)
        {
            addcost[make_pair(nodes[i], nodes[i] + 1)] = 2;
        }
        if (orderednodes.count(nodes[i] + n) > 0 and (nodes[i] + n) < N)
        {
            addcost[make_pair(nodes[i], nodes[i] + n)] = 2;
        }
    }
    auto rng = std::default_random_engine{};
    queue<int> evenqueue;
    queue<int> oddqueue;
    for (int i = 0; i < nodecount; i++)
    {
        int mode = 0;
        vector<int> prevvec(N, -1);
        priority_queue<pair<int, int>> prio;
        vector<int> children(N, 1);

        evenqueue.push(nodes[i]);
        prevvec[nodes[i]] = nodes[i];
        int count = 1;
        int level = 1;

        while (!evenqueue.empty() or !oddqueue.empty())
        {

            vector<int> tempvec;
            if (mode == 0)
            {
                while (!evenqueue.empty())
                {
                    current = evenqueue.front();
                    evenqueue.pop();

                    neighbour = current + 1;
                    if (neighbour % n > 0 and orderednodes.count(neighbour) > 0 and prevvec[neighbour] == -1)
                    {
                        prevvec[neighbour] = current;
                        tempvec.push_back(neighbour);
                        prio.push(make_pair(level, neighbour));
                        count++;
                    }
                    neighbour = current - 1;
                    if (current % n > 0 and orderednodes.count(neighbour) > 0 and prevvec[neighbour] == -1)
                    {
                        prevvec[neighbour] = current;
                        tempvec.push_back(neighbour);
                        prio.push(make_pair(level, neighbour));
                        count++;
                    }
                    neighbour = current + n;
                    if (neighbour < N and orderednodes.count(neighbour) > 0 and prevvec[neighbour] == -1)
                    {
                        prevvec[neighbour] = current;
                        tempvec.push_back(neighbour);
                        prio.push(make_pair(level, neighbour));
                        count++;
                    }
                    neighbour = current - n;
                    if (neighbour >= 0 and orderednodes.count(neighbour) > 0 and prevvec[neighbour] == -1)
                    {
                        prevvec[neighbour] = current;
                        tempvec.push_back(neighbour);
                        prio.push(make_pair(level, neighbour));
                        count++;
                    }
                }
                std::shuffle(std::begin(tempvec), std::end(tempvec), rng);
                for (const auto &e : tempvec)
                {
                    oddqueue.push(e);
                }
                mode = 1;
                level++;
            }
            else
            {
                while (!oddqueue.empty())
                {
                    current = oddqueue.front();
                    oddqueue.pop();

                    neighbour = current + 1;
                    if (neighbour % n > 0 and orderednodes.count(neighbour) > 0 and prevvec[neighbour] == -1)
                    {
                        prevvec[neighbour] = current;
                        tempvec.push_back(neighbour);
                        prio.push(make_pair(level, neighbour));
                        count++;
                    }
                    neighbour = current - 1;
                    if (current % n > 0 and orderednodes.count(neighbour) > 0 and prevvec[neighbour] == -1)
                    {
                        prevvec[neighbour] = current;
                        tempvec.push_back(neighbour);
                        prio.push(make_pair(level, neighbour));
                        count++;
                    }
                    neighbour = current + n;
                    if (neighbour < N and orderednodes.count(neighbour) > 0 and prevvec[neighbour] == -1)
                    {
                        prevvec[neighbour] = current;
                        tempvec.push_back(neighbour);
                        prio.push(make_pair(level, neighbour));
                        count++;
                    }
                    neighbour = current - n;
                    if (neighbour >= 0 and orderednodes.count(neighbour) > 0 and prevvec[neighbour] == -1)
                    {
                        prevvec[neighbour] = current;
                        tempvec.push_back(neighbour);
                        prio.push(make_pair(level, neighbour));
                        count++;
                    }
                }
                std::shuffle(std::begin(tempvec), std::end(tempvec), rng);
                for (const auto &e : tempvec)
                {
                    evenqueue.push(e);
                }
                mode = 0;
                level++;
            }
        }
        if (count < nodecount)
        {
            cout << count << "n:" << i << "," << nodes[i] << endl;
            addcost[make_pair(-1, -1)] = -1;
        }

        while (!prio.empty())
        {
            current = prio.top().second;
            int currlvl = prio.top().first;
            prio.pop();
            int mymin = min(current, prevvec[current]);
            int mymax = max(current, prevvec[current]);
            if (Delta + currlvl > currlvl)
            {
                if (currlvl > 1)
                {

                    addcost[make_pair(mymin, mymax)] = addcost[make_pair(mymin, mymax)] + children[current];
                }
                else
                {
                    addcost[make_pair(mymin, mymax)] = addcost[make_pair(mymin, mymax)] + children[current] - 1;
                }
            }
            children[prevvec[current]] = children[prevvec[current]] + children[current];
        }
    }
    return addcost;
}

PYBIND11_MODULE(calc_weights, m)
{
    m.doc() = "pybind11 example plugin";
    m.def("calc_weights", &calc_weights, "calc_weights");
}