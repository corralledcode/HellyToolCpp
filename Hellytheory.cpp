//
// Created by peterglenn on 3/18/24.
//

#ifndef HELLYTOOLCPP_HELLYTHEORY_H
#define HELLYTOOLCPP_HELLYTHEORY_H

#include <tuple>
#include "Graph.cpp"
#include "Theory.cpp"

class EdgesforHelly : public Edges {
public:
    std::vector<std::tuple<vertextype,vertextype,vertextype>> triangles;
    EdgesforHelly( std::vector<Edge> edgein) : Edges(edgein) {
        //
    }
    EdgesforHelly( int s ) : Edges(s) {
        //
    }
    EdgesforHelly() : Edges() {
        //
    }
    ~EdgesforHelly() {
        //
    }
    void process() override;

};

inline void EdgesforHelly::process() {
    Edges::process();
    auto b = paused();
    int sz = size();
    pause();
    triangles.clear();
    for (vertextype v1 = vertextype(0); v1 <= maxvertex; ++v1) {
        for (vertextype v2 = v1+1; v2 <= maxvertex; ++v2) {
            for (vertextype v3 = v2+1; v3 <= maxvertex; ++v3 ) {
                bool b1,b2,b3 {false};
                for (int i = 0; i < sz; ++i) {
                    Edge e = getdata(i);
                    b1 = b1 || (v1 == e.first && v2 == e.second);
                    b2 = b2 || (v1 == e.first && v3 == e.second);
                    b3 = b3 || (v2 == e.first && v3 == e.second);
                }
                if (b1 && b2 && b3) {
                    triangles.push_back({v1, v2, v3});
                    //std::cout << v1 << v2 << v3 << "\n";
                }
            }
        }
    }
}

class Hellytheory : public Theory {
public:
    Vertices* V;
    EdgesforHelly* E;
    Cover* C;
    bool checkrs();

    Hellytheory() : Theory() {};
};

inline bool Hellytheory::checkrs() {
    V->process();
    E->process();
    C->process();
    int tsz = E->triangles.size();
    bool rscover = true;
    for (int n = 0; (n<tsz && rscover);++n) {
       std::vector<Edges> Emeet {};
       const auto [tfirst,tsecond,tthird] = E->triangles[n];
       //std::cout << "tfirst,tsecond,tthird: "<< tfirst << tsecond << tthird << "\n";
       int csz = C->size();
       for (int i = 0; i < csz; ++i) {
           Edges es = C->getdata(i);
           bool twoofthree = false;
           int esz = es.size();
           int j = 0;
           while (j < esz && !twoofthree) {
               Edge e = es.getdata(j);
               int cnt = 0;
               if (e.first == tfirst || e.second == tfirst)
                   ++cnt;
               if (e.first == tsecond || e.second == tsecond)
                   ++cnt;
               if (e.first == tthird || e.second == tthird)
                   ++cnt;
               twoofthree = (cnt >= 2);
               if (twoofthree) {
                   std::cout << "Edge: " << e.first << " " << e.second << "\n";
                   Emeet.push_back(es);
               }
               ++j;
           }
//           if (twoofthree)
//               Emeet.push_back(es);
       }
       bool sharedvertex = false;
       int vsz = V->size();
       int k = 0;
       int Emeets = Emeet.size();
       while (k < vsz && !sharedvertex) {
           vertextype v = V->getdata(k);
           std::cout << "vertex " << v << "\n";
           int l = 0;
           bool shared = true;
           while (shared && (l < Emeets)) { // check every cover sharing two of three vertices
               Edges es = Emeet[l];
               int m = 0;
               shared = false;
               while (!shared && (m < es.size())) { // for each such cover see if it contains v
                   shared = shared || ((es.getdata(m).first == v) || (es.getdata(m).second == v));
                   //std::cout << "shared " << shared << " " << es.getdata(m).first << es.getdata(m).second << "\n";
                   ++m;
               }
               ++l;
           }
           sharedvertex = sharedvertex || shared;
           ++k;
       }
       rscover = rscover && sharedvertex;
       std::cout << "Triangle <" << tfirst << ", " << tsecond << ", " << tthird << ">";
       if (rscover)
           std::cout << " all meet\n";
       else
           std::cout << " fail to meet at this triangle\n";
    }
    return rscover;
}

#endif //HELLYTOOLCPP_HELLYTHEORY_H

