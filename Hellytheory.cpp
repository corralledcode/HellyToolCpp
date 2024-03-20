//
// Created by peterglenn on 3/18/24.
//

#ifndef HELLYTOOLCPP_HELLYTHEORY_H
#define HELLYTOOLCPP_HELLYTHEORY_H

#include <tuple>
#include "Graph.cpp"
#include "Theory.cpp"
#include "Formatgraph.cpp"

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
    virtual void process() override;

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
                bool b1 {false};
                bool b2 {false};
                bool b3 {false};
                for (int i = 0; i < sz; ++i) {
                    Edge e = getdata(i);
                    b1 = b1 || (v1 == e.first && v2 == e.second);
                    b2 = b2 || (v1 == e.first && v3 == e.second);
                    b3 = b3 || (v2 == e.first && v3 == e.second);
                }
                if (b1 && b2 && b3) {
                    triangles.push_back({v1, v2, v3});
                    //std::cout << "classic EdgesforHelly found triangle: " << v1 << v2 << v3 << "\n";
                }
            }
        }
    }
}

class EdgesforHellyfast : public EdgesforHelly {
public:
    std::vector<std::tuple<vertextype,vertextype,vertextype>> triangles;
    EdgesforHellyfast( std::vector<Edge> edgein) : EdgesforHelly(edgein) {
        //
    }
    EdgesforHellyfast( int s ) : EdgesforHelly(s) {
        //
    }
    EdgesforHellyfast() : EdgesforHelly() {
        //
    }
    ~EdgesforHellyfast() {
        //
    }
    void process() override;

};

inline void EdgesforHellyfast::process() {
    Edges::process();
    auto b = paused();
    pause();
    triangles.clear();
    int esz = size();
    for (int n = 0; n < esz; ++n) {
        Edge e = getdata(n);
        //std::cout << e.first << " " << e.second << "\n";
        bool b1,b2 {false};
        for (int m = n+1; m < esz; ++m) {
            Edge e2 = getdata(m);
            //std::cout << e2.first << " " << e2.second << "\n";
            vertextype soughtv1, soughtv2;
            bool b1 = true;
            if (e.first == e2.first) {
                soughtv1 = e.second;
                soughtv2 = e2.second;
            } else {
                if (e.first == e2.second) {
                    soughtv1 = e.second;
                    soughtv2 = e2.first;
                } else {
                    if (e.second == e2.first) {
                        soughtv1 = e.first;
                        soughtv2 = e2.second;
                    } else {
                        if (e.second == e2.second) {
                            soughtv1 = e.first;
                            soughtv2 = e2.first;
                        } else {
                            b1 = false;
                        }
                    }
                }
            }

            if (b1) {
                //std::cout << "soughtv2, soughtv2 " << soughtv1 << " " << soughtv2 << "\n";
                for (int j = m+1; j < esz; ++j) {
                    Edge e3 = getdata(j);
                    if ((e3.first == soughtv1 && e3.second == soughtv2)
                        || (e3.first == soughtv2 && e3.second == soughtv1)) {
                        //std::cout << "found " << e3.first << " " << e3.second << "\n";
                        //std::cout << "pushing back " << e.first << e.second
                        //          << (e.first == e2.first ? e2.second : e2.first) << "\n";
                        triangles.push_back({e.first, e.second, (e.first == e2.first ? e2.second : e2.first)});

                    }
                }
            }
        }
    }
    if (!b)
        resume();
}

class Hellytheory : public Theory {
public:
    Vertices* V;
    EdgesforHellyfast* E;
    Cover* C;
    Formatvertices* FV;
    Formatedges* FE;
    Formatcover* FC;

    bool checkcover();
    bool checkcoverlegal();
    bool checkrs();

    Hellytheory() : Theory() {};
};

inline bool Hellytheory::checkcover() {
    int esz = E->size();
    bool allcovered = true;
    for (int n = 0; n < esz; ++n) {
        bool covered = false;
        Edge e = E->getdata(n);
        int csz = C->size();
        for (int j = 0; j < csz && !covered; ++j) {
            Edges es = C->getdata(j);
            int cesz = es.size();
            for (int k = 0; k < cesz && !covered; ++k) {
                Edge ese = es.getdata(k);
                covered = (covered || (ese == e));
            }
        }
#ifdef VERBOSE
        if (!covered)
            std::cout << "Edge " << FV->lookup(e.first) << "," << FV->lookup(e.second) << " not covered.\n";
#endif
        allcovered = allcovered && covered;
    }
    return allcovered;
}

inline bool Hellytheory::checkcoverlegal() {
    int csz = C->size();
    int esz = E->size();
    bool alllegal = true;
    for (int n = 0; n < csz; ++n) {
        Edges ce = C->getdata(n);
        bool legal = false;
        int cesz = ce.size();
        for (int m = 0; m < cesz && !legal; ++m) {
            for (int j = 0; j < esz; ++j) {
                legal = legal || ce.getdata(m) == E->getdata(j);
            }
#ifdef VERBOSE
            if (!legal)
                std::cout << "Edge " << FV->lookup(ce.getdata(m).first)
                            << " " << FV->lookup(ce.getdata(m).second) << " in cover is not legal.\n";
#endif

        }
        alllegal = alllegal && legal;

    }
    return alllegal;
}

inline bool Hellytheory::checkrs() {
    V->resume();
    E->resume();
    C->resume();
    std::cout << "Triangles size " << E->triangles.size() << "\n";
#ifdef CHECKS
    std::cout << (checkcover() ? "All edges covered.\n" : "All edges not covered.\n");
    std::cout << (checkcoverlegal() ? "" : "Illegal edges found in the cover.\n");
#endif
    int tsz = E->triangles.size();
    bool rscover = true;
    for (int n = 0; (n<tsz && rscover);++n) {
       std::vector<Edges> Emeet {};
       const auto [tfirst,tsecond,tthird] = E->triangles[n];
#ifdef VERBOSE
       std::cout << "Triangle: <"<< FV->lookup(tfirst) << ", "<< FV->lookup(tsecond)
                << ", " << FV->lookup(tthird) << ">\n";
#endif
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
#ifdef VERBOSE
                   std::cout << "   meets Cover: ";
                   for (int l = 0; l < esz; ++l) {
                       Edge e2 = es.getdata(l);
                       std::cout << FV->lookup(e2.first) << " " << FV->lookup(e2.second) << ", ";
                   }
                   std::cout << "\b\b  \n";
#endif
                   Emeet.push_back(es);
               }
               ++j;
           }
       }
       bool sharedvertex = false;
       int vsz = V->size();
       int k = 0;
       int Emeets = Emeet.size();
       while (k < vsz && !sharedvertex) {
           vertextype v = V->getdata(k);
#ifdef VERBOSE
//           std::cout << "vertex " << FV->lookup(v) << "\n";
#endif
           int l = 0;
           bool shared = true;
           while (shared && (l < Emeets)) { // check every cover sharing two of three vertices
               Edges es = Emeet[l];
               int m = 0;
               shared = false;
               int esz = es.size();
               while (!shared && (m < esz)) { // for each such cover see if it contains v
                   Edge e = es.getdata(m);
                   shared = shared || ((e.first == v) || (e.second == v));
                   ++m;
               }
               ++l;
           }
           sharedvertex = sharedvertex || shared;
           ++k;
       }
       rscover = rscover && sharedvertex;
#ifdef VERBOSE
       if (rscover)
           std::cout << "   All meet.\n";
       else
           std::cout << "   Fail to meet at this triangle\n";
#endif
    }
    return rscover;
}

#endif //HELLYTOOLCPP_HELLYTHEORY_H

