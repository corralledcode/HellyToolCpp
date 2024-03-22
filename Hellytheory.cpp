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
    EdgesforHelly* E;
    Cover* C;
    Formatvertices* FV;
    Formatedges* FE;
    Formatcover* FC;

    bool checkcover();
    bool checkcoverlegal();
    bool checkrs();
    void findrscovers(Cover hintCover);
    std::vector<Cover> recursefindcovers( std::vector<Edges>* completeedgesets, Cover* hintCover, int progress);
    void enumeratecompletesets(std::vector<Edges>* Evar, Cover hintCover);
    void findncovers(std::vector<Cover>* Cvar, Cover hintCover, int n);
    void findncoversbare(std::vector<Cover>* Cvar, Cover hintCover, std::vector<Edges>* completeEdgesets, int n);

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
    bool allalllegal = true;
    for (int n = 0; n < csz; ++n) {
        Edges ce = C->getdata(n);
        bool alllegal = true;
        int cesz = ce.size();
        for (int m = 0; m < cesz && alllegal; ++m) {
            bool legal = false;
            for (int j = 0; j < esz && !legal; ++j) {
                legal = legal || ce.getdata(m) == E->getdata(j);
            }
            alllegal = alllegal && legal;
#ifdef VERBOSE
            if (!legal)
                std::cout << "Edge " << FV->lookup(ce.getdata(m).first)
                            << " " << FV->lookup(ce.getdata(m).second) << " in cover is not legal.\n";
#endif

            alllegal = alllegal && legal;
        }
        allalllegal = allalllegal && alllegal;

    }
    return allalllegal;
}

inline bool Hellytheory::checkrs() {
    V->resume();
    E->resume();
    C->resume();
#ifdef VERBOSE
    std::cout << "Triangles size " << E->triangles.size() << "\n";
#endif
#ifdef CHECKS
#ifdef VERBOSE
    std::cout << (checkcover() ? "All edges covered.\n" : "All edges not covered.\n");
    std::cout << (checkcoverlegal() ? "" : "Illegal edges found in the cover.\n");
#endif
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

inline void Hellytheory::enumeratecompletesets(std::vector<Edges>* Evar, Cover hintCover) {
    int Vsz = V->size();
    int Esz = E->size();
    int s = 0;
    int n = V->size();
    int bin = 0;
    if (n > 0)
        bin = 1;
    bin = bin << n;
    //if ((*Evar).size() < bin) {
    //    std::cout << "Evar passed is not large enough.\n";
    //    return;
    //}
    std::vector<Edges> Es{};
    Es.clear();
    //std::vector<vertextype> vertexpool {};
    //for (int i = 0; i < hintCover.size(); ++i) {
    //    Es.push_back(hintCover.getdata(i));
     /*
        for (int j = 0; j < hintCover.getdata(i).size(); ++j) {
            vertexpool.push_back(hintCover.getdata(i).getdata(j).first);
            vertexpool.push_back(hintCover.getdata(i).getdata(j).second);
        }*/
    //}

    for (int c = 0; c <= bin; ++c) {
        //Cover Cvr = (*Cvar)[c];
        std::vector<Edge> es{};
        es.clear();
        std::vector<vertextype> vertices{};
        vertices.clear();
        for (int i = 0; i < n; ++i) {
            if ((c & (1 << i)) > 0)
                vertices.push_back(V->getdata(i));
        }
        bool allfound = true;
        for (int j = 0; (j < vertices.size()) && allfound; ++j) {
            for (int l = j+1; (l < vertices.size()) && allfound; ++l) {
                bool found = false;
                for (int k = 0; (k < Esz) && !found; ++k) {
                    Edge e = E->getdata(k);
                    if ((e.first == vertices[j]) && (e.second == vertices[l])
                        || ((e.second == vertices[j]) && (e.first == vertices[l]))) {
                        es.push_back(e);
                        found = true;
                    }
                }
                allfound = allfound && found;
            }
        }

        if (allfound) {
//            for (int l = 0; l < es.size(); ++l) {
//                std::cout << es[l].first << " " << es[l].second << ", ";
//            }
//            std::cout << "\n";
            bool allall = true;
            int esz = es.size();
            for (int i = 0; i < esz && allall; ++i) {
                int hCsz = hintCover.size();
                bool all = false;
                for (int k = 0; k < hCsz && !all; ++k) {
                    bool found = false;
                    Edges he = hintCover.getdata(k);
                    int hesz = he.size();
                    for (int j = 0; j < hesz && !found; ++j) {
                        found = found || (es[i] == he.getdata(j));
                    }
                    all = all || found;
                }
                allall = allall && all;
            }
            if (!allall)
                Es.push_back(es);
            else { /*
                for (int j = 0; j < es.size(); ++j)
                    std::cout << es[j].first << " " << es[j].second << ",";
                std::cout << "\n"; */
            }
        }
        else { /*
            for (int j = 0; j < es.size(); ++j)
                std::cout << es[j].first << " " << es[j].second << ",";
            std::cout << "not allfound \n"; */
        }
    }
    int Essz = Es.size();
    std::cout << Essz << "Essz\n";
    //std::vector<Edges> CompleteEs{};
    for (int k = 0; k < Essz; ++k) { // check that the cover is a complete set cover
        Edges es2 = Es[k];
        int es2z = es2.size();
        std::cout << es2z << "es2z\n";
        for (int l = 0; l < es2z; ++l) {
            std::cout << es2.getdata(l) << "\n";
        }
        if (es2z > 2) {
            Evar->push_back(es2);

            /*
            std::vector<vertextype> v{};
            v.clear();
            for (int l = 0; l < es2z; ++l) {
                Edge e2 = es2.getdata(l);

                bool foundfirst = false;
                bool foundsecond = false;
                for (int j = 0; j < v.size(); ++j) {
                    foundfirst = foundfirst || v[j] == e2.first;
                    foundsecond = foundsecond || v[j] == e2.second;
                }
                if (!foundfirst)
                    v.push_back(e2.first);
                if (!foundsecond)
                    v.push_back(e2.second);
                //std::cout << "adding vertex " << FV->lookup (e2.first) << " " << FV->lookup(e2.second)<<"\n";

                v.push_back(e2.first);
                v.push_back(e2.second);
            } */
         /*
            bool allcomplete = true;
            int vsz = v.size();
            for (int m = 0; m < vsz && allcomplete; ++m) { // for every vertex pair...
                for (n = m + 1; n < vsz && allcomplete; ++n) {
                    if (v[m] != v[n]) { //redundant
                        bool complete = false;
                        for (int l = 0; l < es2z && !complete; ++l) {
                            Edge e3 = es2.getdata(l);
                            complete = (complete || (e3.first == v[m] && e3.second == v[n])
                                       || (e3.second == v[m] && e3.first == v[n]));
                        }
                        allcomplete = allcomplete && complete;
                    }
                }
            }
            if (allcomplete) {
                /*int hsz = hintCover.size();
                // if any set of edges es in hintCover contains  es2
                bool contains = false;
                for (int i = 0; i < hsz && !contains; ++i) {
                    Edges he = hintCover.getdata(i);
                    Evar->push_back(he);
                    int esz = es2.size();
                    bool allmatch = true;
                    for (int j = 0; j < esz; ++j) {
                        bool match = false;
                        for (int k = 0; k < he.size(); ++k) {
                            match = match || (es2.getdata(j) == hintCover.getdata(i)[k]);
                        }
                        allmatch = allmatch && match;
                    }
                    contains = contains || allmatch;

                }
                if (!contains)
                    Evar->push_back(es2);
                //std::cout << "allcomplete at ";
                //for (int n = 0; n < es2.size(); ++n)
                //    std::cout << es2.getdata(n) << " ";
                //std::cout << "\n";
            } */
        }
    }
    std::cout << Evar->size() << "evarsz\n";
    // data is returned in Evar
}

inline std::vector<std::vector<bool>> recursebool(int n) {
    std::vector<std::vector<bool>> B {};
    if (n <= 0) {
        std::vector<bool> b {};
        B.push_back(b);
        return B;
    }

    std::vector<bool> workingbool {};

    --n;
    B = recursebool(n);
    int Bsz = B.size();
    std::vector<std::vector<bool>> B2 {};
    for (int m = 0; m < Bsz; ++m) {
        workingbool = B[m];
        workingbool.push_back(false);
        B2.push_back(workingbool);
        workingbool = B[m];
        workingbool.push_back(true);
        B2.push_back(workingbool);
    }
    return B2;
}


inline std::vector<Cover> Hellytheory::recursefindcovers(std::vector<Edges>* completeedgesets, Cover* hintCover, int progress) {
    std::vector<Cover> Cvrs {};
    if (progress <= 0) {
        Cover c {*hintCover};
        Cvrs.push_back(c);
        //std::cout << "Ground zero recursefindcovers, Cvrs.size()"<<Cvrs.size()<<"\n";
        return Cvrs;
    }

    --progress;
    Cvrs = recursefindcovers(completeedgesets,hintCover,progress);
    std::vector<Cover> CvrsReturn {};
    CvrsReturn.clear();
    int Csz = Cvrs.size();
    Edges es;
    es = (*completeedgesets)[progress];
    // for each cover
    // branch into one invocation with es added and one invocation without es added

    for (int m = 0; m < Csz; ++m) {
        //std::cout << "m, Csz = " << m << ","<<Csz<<"\n";
        //C = &(Cvrs[m]);
        //FC->simplifycover();

        Cover c = Cvrs[m];
        std::vector<Edges> ces {};
        ces.clear();
        int csz = c.size();
        for (int k = 0; k < csz; ++k) {
            ces.push_back(c.getdata(k));
        }
        ces.push_back(es);

        Cover c2 {};
        c2.readvector(ces);

        CvrsReturn.push_back(c);
        CvrsReturn.push_back(c2);
    }
    return CvrsReturn;
}

inline void Hellytheory::findncovers(std::vector<Cover>* Cvar, Cover hintCover, int n) {
    std::vector<Edges> Es;
    int vsz = V->size();
    if (vsz > 31) {
        std::cout << "Unable to process more than 31 vertices.\n";
        return;
    }
    //Es.size() = (1 << vsz);

    std::cout << "Using maxcliquesize = " << n << "\n";

    Es.clear();
    //std::cout << "resize " << (1<<n) << "\n";
    //Es.resize(1 << n);
    enumeratecompletesets(&Es, hintCover );
    for (int m = 0; m < Es.size(); ++m) {
        if (Es[m].size() > 0)
            std::cout << "Edges [";
        else
            std::cout << "Blank edge\n";
        for (int j = 0; j < Es[m].size(); ++j) {
            std::cout << "[" << FV->lookup(Es[m][j].first) << ", " << FV->lookup(Es[m][j].second) << "], ";
        }
        if (Es[m].size() > 0)
            std::cout << "\b\b] \n";
    }
    std::vector<Cover> Cvrs = recursefindcovers(&Es,&hintCover, Es.size());
    for (int n = 0; n < Cvrs.size();++n) {
        //Cvrs[n].simplifycover();
        Cvar->push_back(Cvrs[n]);
    }
}


inline void Hellytheory::findncoversbare(std::vector<Cover>* Cvar, Cover hintCover, std::vector<Edges>* completeEdgesets, int n) {

    // if n == 0 return hintCover

    // n = n-1

    // for every edgeset Es in completeedgesets
    //    if Es size == n
    //        split
    // for every complete set Es of vertices/edges
    //    add Es to a local copy of hintCover
    //    create a vector Cstemp of Cover type
    //    call findncoversbare( Cstemp, new hintCover, n)
    // add each Cstemp to Cvar





    // for each vertex v1
    //    for each adjacent vertex v2
    //        Etemp.push_back (v1,v2)
    // for (i = 0; i < 2^^size(Etemp); ++i)
    //
    //        tempCover = hintCover;
    //          tempCover.
    //        findncoversbare( Cvar, hintCover
    //        hintCover.add edge v1,v2
    //        findncoversbare( Ctemp, hintCover, n )
    //        combine Cvar with Ctemp
    // "return" Cvar
/*
    if (n == 0) {
        (*Cvar).clear();
        (*Cvar).push_back(hintCover);
        return;
    }

    n = n-1;

    Cover Ctemp(hintCover);
    Edges Etemp(0);
    int vsz = V->size();
    int esz = E->size();
    for (int m = 0; m < vsz; ++m) {
        vertextype v1 = V->getdata(m);
        for (int n = 0; n < vsz; ++n) {
            if (m != n) {
                vertextype v2 = V->getdata(n);
                for (int k=0;k<esz;++k) {
                    Edge e = E->getdata(k);
                    bool adjacent = ((e.first == v1 && e.second == v2) || (e.first == v2 && e.second == v2));
                    if (adjacent)

                    findncoversbare(*Ctemp, hintCover, n);
                }
            }
        }
    }
*/

    //std::vector<Cover> allc {};
    /*if (n <= 0) {
        Cvar.clear();
        Cvar.push_back(hintCover);
        std::cout << "Zero " << Cvar.size() << "\n";
        return;
    }
    n = n-1;*/

/*
    s=0;
    int k = 0;
    int Cesz = CompleteEs.size();
    for (int i = 0; i < Cesz; ++i) {
        std::cout << "Round 1: CompleteEs[i].size() = " << CompleteEs[i].size() << "\n";
    }
    std::cout << "Cesz " << Cesz << "\n";


                Es.push_back(CompleteEs[i]);
//                std::cout << "i, CompleteEs[i].size()=" << i << ", " << CompleteEs[i].size() << "\n";
//                std::cout << "\n";
            }
            //for (int n = 0; n < Es.size(); ++n) {
            // simplifycover

        }
        (*Cvar)[k].readvector(Es);
        ++k;
    }*/
}

inline void Hellytheory::findrscovers( Cover hintCover ) {
    V->resume();
    E->resume();
    C->resume();

    //std::vector<Edges> Evar {};
    std::vector<Cover> Cvar {};
    //Evar.resize(1 << V->size());
    //enumeratecompletesets(&Evar, hintCover );
    findncovers(&Cvar, hintCover,V->size());

    /*
    for (int n = 0; n < Evar.size(); ++n) {
        if (Evar[n].size() > 0)
            std::cout << "Edges ";
        for (int j = 0; j < Evar[n].size(); ++j) {
            std::cout << "[" << Evar[n][j].first << ", " << Evar[n][j].second << "], ";
        }
        if (Evar[n].size() > 0)
            std::cout << "\b\b  \n";
    }
*/
    int Cvsz = Cvar.size();
    bool Brs[Cvsz];
    for (int m = 0; m < Cvsz; ++m) {
        C = &(Cvar[m]);
        //FC->matchiedata();
        std::cout << "Cover " << m << "\n";
        for (int n = 0; n < C->size(); ++n) {
            Edges es = C->getdata(n);
            int essz = es.size();
            if (essz > 0)
                std::cout << "...Edges [";
            for (int j=0; j < essz; ++j) {
                std::cout << "["<< FV->lookup(es.getdata(j).first) << "," << FV->lookup(es.getdata(j).second) << "], ";
            }
            if (essz > 0)
                std::cout << "\b\b]\n";
        }

        if (checkcover()) {
            Brs[m] = checkrs();
            std::cout << (Brs[m] ? "checkrs() returns true." : "checkrs() returns false.");
            std::cout << "\n";
        } else {
            std::cout << "Does not cover all edges\n";
            Brs[m] = false;
        }
    }

    std::cout << "Recap: ";
    int cnt = 0;
    for (int m = 0;m<Cvsz;++m) {
        if (Brs[m]) {
            cnt++;
        }
    }
    std::cout << cnt << " passed.\n";













/*
    auto emptyCover = new Cover(1);
    auto emptyEdges = new Edges(0);
    emptyCover->setdata(*emptyEdges,0);
    //start by enumerating the complete sets
    int sz = 1;
    int vsz = V->size();
    for (int c = 0; c < vsz; ++c)
        sz = sz * 2;
    std::vector<Cover> vsubsets;
    vsubsets.resize(sz);
    findncovers( &vsubsets, *emptyCover, vsz ); // replace with E->maxclique() once coded

    for (int n = 0; n < vsubsets.size(); ++n) {
        std::cout << "Cover " << n;
        for (int j = 0; j < vsubsets[n].size(); ++j)
            for (int k = 0; k < vsubsets[n].getdata(j).size(); ++k) {
                Edge e = vsubsets[n].getdata(j).getdata(k);
                std::cout << FV->lookup(e.first) << " " << FV->lookup(e.second) << " ";
            }
        std::cout << "\n";
    }

    bool brs[sz];
    for (int n = 0; n < sz; ++n) {
        C = &vsubsets[n];
        std::cout << "Cover " << n << ": [";
        for (int j = 0; j < C->size(); ++j) {
            for (int k = 0; k < C->getdata(j).size(); ++k) {
                Edge e = C->getdata(j).getdata(k);
                std::cout << "[" << FV->lookup(e.first) << "," << FV->lookup(e.second) << "], ";
            }
            std::cout << "\b\b], [";
        }
        std::cout << "\b\b\n";
        if (checkcover()) {
            brs[n] = checkrs();
            std::cout << "is RS cover? ";
            if (brs[n])
                std::cout << "True";
            else
                std::cout << "False";
            std::cout  << "\n";
        } else {
            brs[n] = false;
            //std::cout << "does not cover all edges.\n";
        }
    }*/
}

#endif //HELLYTOOLCPP_HELLYTHEORY_H
