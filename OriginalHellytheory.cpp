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
    int maxcliquesize = 0;

    bool checkcover();
    bool checkcoverlegal();
    bool checkrs();
    void findrscovers(Cover hintCover);
    std::vector<Cover> recursefindcovers( std::vector<Edges>* completeedgesets, Cover* hintCover, int progress);
    std::vector<Edges> enumeratevertexsubsets(  std::vector<Edge>* EE, Cover hintCover );
    void enumeratecompletesets(std::vector<Edges>* Evar, Cover hintCover, int n);
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
    for (int n = 0; (n < tsz && rscover); ++n) {
        std::vector<Edges> Emeet{};
        const auto [tfirst, tsecond, tthird] = E->triangles[n];
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

inline bool adjacent( Edge e1, Edge e2 ) {
    return ((e1.first == e2.first) || (e1.first == e2.second) || (e1.second == e2.second));
}

inline bool completeedgeset( std::vector<Edge>* EE) {
    const int EEsz = (*EE).size();
    std::vector<vertextype> v{};
    v.resize(2 * EEsz);
    int vcnt = 0;
    for (; vcnt < EEsz; ++vcnt) {
        v[2 * vcnt] = (*EE)[vcnt].first;
        v[2 * vcnt + 1] = (*EE)[vcnt].second;
    }
    bool allfound = true;
    for (int m = 0; (m < (2*EEsz)) && allfound; ++m)
        for (int n = m + 1; (n < (2*EEsz)) && allfound; ++n) {
            if (v[m] != v[n]) {
                bool found = false;
                for (int k = 0; (k < EEsz) && !found; ++k) {
                    Edge e = (*EE)[k];
                    Edge e2{};
                    e2.first = v[m];
                    e2.second = v[n];
                    found = (found || ((e.first == e2.first)&&(e.second == e2.second))
                                   || ((e.first == e2.second)&&(e.second == e2.first)));
                }
                allfound = (allfound && found);
            }
        }
    return allfound;
}


inline std::vector<Edges> Hellytheory::enumeratevertexsubsets(  std::vector<Edge>* EE, Cover hintCover) {
    Edges es{};
    int EEsz = EE->size();
    int Esz = E->size();
    int Vsz = V->size();
    std::vector<Edges> EsReturn{};
    if (EEsz <= 0) {
        EsReturn.push_back(es);
        return EsReturn;
    }

    // for each edge <n,k> in EE
    //    remove <n,k> from EE and call recursively
    //    let j be the number of edges e in E that meet <n,k>
    //    for all powersets of j (use int up to say twenty decimal digits)
    //       if p is complete
    //          EsReturn = TmpEs with each Es in TmpEs branched by count of p
    //
    //


    Edge e = (*EE)[0];
    //std::cout << "Edge e = " << e << "\n";
    // for every vertex v that forms a triangle with e
    //    add the edges v,e.first and v,e.second


    std::vector<Edge> Etmp{};
    Etmp.resize(Esz);
    Etmp[0] = e;
    int cnt = 1;
    int vcnt = 0;
    std::vector<vertextype> v{};
    v.resize(Vsz - 2);
    for (int n = 0; n < Vsz; ++n) {  // prime opportunity to use adjacency matrix
        vertextype newv = V->getdata(n);
        if ((newv != e.first) && (newv != e.second)) {
            Edge e1{};
            e1.first = ((newv < e.first) ? newv : e.first);
            e1.second = ((newv > e.first) ? newv : e.first);
            Edge e2{};
            e2.first = ((newv < e.second) ? newv : e.second);
            e2.second = ((newv > e.second) ? newv : e.second);
            bool found1 = false;
            bool found2 = false;
            for (int m = 0; (m < Esz) && !(found1 && found2); ++m) {
                vertextype newv1 = E->getdata(m).first;
                vertextype newv2 = E->getdata(m).second;
                found1 = (found1 || (E->getdata(m) == e1));
                found2 = (found2 || (E->getdata(m) == e2));
            }
            if (found1 && found2) {
                //std::cout << "newv passes " << newv << "\n";
                Etmp[cnt] = e1;
                Etmp[cnt + 1] = e2;
                cnt += 2;
                v[vcnt] = newv;
                vcnt++;
            }
        }
    }


    //std::vector<Edge> Etmp2 {};
    Etmp.resize(cnt + vcnt); //sloppy
    //Etmp2.resize(cnt+1);
    //std::cout << "Etmp2.resize first round = " << cnt+1 << "\n";
    int cnt2 = cnt;
    for (int k = 0; (k < vcnt); ++k)
        for (int l = k + 1; (l < vcnt); ++l) {
            bool found = false;
            for (int m = 0; (m < Esz) && !found; ++m) {
                Edge etemp;
                etemp.first = v[k];
                etemp.second = v[l];
                if ((*E)[m] == etemp) {
                    Etmp[cnt2] = etemp;
                    ++cnt2;
                    found = true;
                }
            }
        }
    Etmp.resize(cnt2);

    //std::cout << "Etmp2.resize = " << cnt2 << "\n";
    std::vector<Edge> Etemp;
    int cnt3 = 0;
    Etemp.resize(EEsz);
    for (int s = 0; s < EEsz; ++s) {
        bool found = false;
        for (int m = 0; (m < cnt2) && !found; ++m) {
            found = (found || (Etmp[m] == (*EE)[s]));
        }
        if (!found) {
            Etemp[cnt3] = (*EE)[s];
            ++cnt3;
        } else {
            //std::cout << "removing edge " << (*EE)[s] << "\n"; // to aid in the recursive step, remove vertices here considered
        }
    }
    Etemp.resize(cnt3);

    int bin = ((1 << cnt2) - 1);
    //std::cout << "using bin = " << bin << "\n";
    std::vector<std::vector<Edge>> Ps{};
    Ps.clear();
    Ps.resize(bin + 1);
    int i = 0; // index for Ps
    for (int c = 0; c <= bin; ++c) {
        std::vector<Edge> P{};
        P.clear();
        P.resize(cnt2);
        int j = 0;
        for (int k = 0; k < cnt2; ++k) {
            if ((c & (1 << k)) > 0) {
                P[j] = Etmp[k];
                ++j;
            }
        }
        P.resize(j);



        //check that P is complete


        //std::cout << "Checking P for completeness: " << P.size() << "\n";

        if (completeedgeset(&P)) {
            Ps[i] = P;
            ++i;
        }
    }


    Ps.resize(i);

    // for each edgeset EPs in Ps
    //    for each edgeset Es in hintcover
    //        if EPs is contained in Es
    //           remove EPs from Ps

    std::vector<std::vector<Edge>> PsTemp {};
    int idx = 0; // index for PsTemp
    PsTemp.resize(i);
    int hCsz = hintCover.size();
    //std::cout << "HintCover size " << hCsz << "\n";
    for (int t = 0; t < i; ++t) {                // for each edgeset EPs in Ps
        std::vector<Edge> EPs = Ps[t];

        bool allallfound = false;
        for (int s = 0; !allallfound && (s < hCsz); ++s) { // for each edgeset Es in hintCover
            Edges Es = hintCover.getdata(s);
            bool allfound = true;                            // is every edge in EPs found in Es?
            for (int u = 0; allfound && (u < EPs.size()); ++u) {
                // is this edge in Es?
                bool found = false;
                for (int v = 0; !found && (v < Es.size()); ++v) {
                    found = found || (EPs[u] == Es[v]);
                }
                allfound = allfound && found;
            }
            allallfound = allallfound || allfound;
        }
        if (!allallfound || (EPs.size() == 0)) {
            PsTemp[idx] = EPs;
            ++idx;
            //std::cout << " not allfound ";
            //for (int a = 0; a < EPs.size(); ++a)
            //    std::cout << FV->lookup(EPs[a].first) << FV->lookup(EPs[a].second) << "\n";
        } else {
            //std::cout << "all found ";
            //for (int a = 0; a < EPs.size(); ++a )
            //    std::cout << EPs[a].first << EPs[a].second << "\n";
        }
    }

    PsTemp.resize(idx);
    for (int u = 0; u < idx; ++u) {
        Ps[u] = PsTemp[u];  // lazy about normalizing variable names
        //for (int v = 0; v < Ps[u].size(); ++v)
        //    std::cout << Ps[u][v].first << Ps[u][v].second << " PS \n";
        //std::cout << "\n";
    }
    //std::cout << "\n";
    Ps.resize(idx);

    std::vector<Edges> TmpEs{};
    TmpEs = enumeratevertexsubsets(&Etemp, hintCover);

    // for each Es in TmpEs
    //    for each c
    //       EsReturn.copy out Es followed by c

    int pos = 0;
    EsReturn.resize(TmpEs.size() * Ps.size());
    //std::cout << "EsReturn.resize " << TmpEs.size() << " " << Ps.size() << "\n";
    for (int n = 0; n < TmpEs.size(); ++n) {
        for (int k = 0; k < Ps.size(); ++k) {
            // copy out Es followed by c
            std::vector<Edge> EsTmp{};
            EsTmp.clear();
            EsTmp.resize(TmpEs[n].size() + Ps[k].size());
            int l;
            for (l = 0; l < TmpEs[n].size(); ++l) {
                EsTmp[l] = TmpEs[n][l];
            }
            int m2 = 0;
            for (int m = 0; m < Ps[k].size(); ++m) {
                bool dupe = false;
                for (int l2 = 0; l2 < TmpEs[n].size(); ++l2)
                    dupe = dupe || (TmpEs[n][l2] == Ps[k][m]);
                if (!dupe) {
                    EsTmp[l + m2] = Ps[k][m];
                    ++m2;
                }
            }
            if (completeedgeset(&EsTmp)) {
                EsTmp.resize(l + m2);
                Edges E2{};
                E2.readvector(EsTmp);
                //E2.removeduplicates();
                EsReturn[pos] = E2;
                ++pos;
            }
        }
    }
    EsReturn.resize(pos);
    std::vector<Edges> EsReturn2 {};
    EsReturn2.resize(EsReturn.size());
    int pos4 = 0; // index for EsReturn2
    //std::cout << "EsReturn.size() = " << EsReturn.size() << "\n";
    for (int r = 0; r< EsReturn.size(); ++r) {
        bool dupe = false;
        for (int s = 0; (s< r) && !dupe; ++s){
            dupe = dupe || (EsReturn[r] == EsReturn[s]);
        }
        if (!dupe) {
            EsReturn2[pos4].setsize(EsReturn[r].size());
            //std::cout << "EsReturn.size() " << EsReturn[r].size() << "\n";
            for (int i = 0; i < EsReturn[r].size(); ++i)
                EsReturn2[pos4].setdata(EsReturn[r].getdata(i),i);
            ++pos4;
        }
    }
    EsReturn2.resize(pos4);
    //std::cout << "EsReturn2 size " << EsReturn2.size() << "\n";
    return EsReturn2;

}

inline void Hellytheory::enumeratecompletesets(std::vector<Edges>* Evar, Cover hintCover, int n) {
    int Vsz = V->size();
    int Esz = E->size();
    int s = 0;
    std::vector<Edges> Es{};
    Es.clear();

    // first enumerate all edges not covered by the hint
    // next consider every n-or-less-sized subset of vertices

    std::vector<Edge> es {};
    for (int k = 0; k < Esz; ++k) {  // for each edge
        Edge etmp = E->getdata(k);
        bool covered = false;
        for (int m = 0; !covered && (m < hintCover.size()); ++m) {
            for (int j = 0; (j < hintCover.getdata(m).size()) && !covered; ++j) {
                covered = (covered || (etmp == hintCover.getdata(m).getdata(j)));
            }
        }
        if (!covered) {
            es.push_back(E->getdata(k));
        }
    }

    Es = enumeratevertexsubsets(&es, hintCover);
    for (int m=0;m<Es.size();++m)
        if (Es[m].size() >= 3)
            Evar->push_back(Es[m]);
    return;
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
    if (progress > 19)
        std::cout << "Recursion depth high = " << progress << "\n";
    Cvrs = recursefindcovers(completeedgesets,hintCover,progress);
    if (progress > 19)
        std::cout << "Note deep recursion depth ...returns..." << progress << "\n";
    std::vector<Cover> CvrsReturn {};
    CvrsReturn.clear();
    int Csz = Cvrs.size();
    Edges es;
    es = (*completeedgesets)[progress];
    // for each cover
    // branch into one invocation with es added and one invocation without es added

    CvrsReturn.resize(Csz*2);  // looks like about 5% speed-up exlicitly managing this vector's size as opposed to push_back calls
    int cri = 0;
    for (int m = 0; m < Csz; ++m) {

        //std::cout << "m, Csz = " << m << ","<<Csz<<"\n";
        //C = &(Cvrs[m]);
        //FC->simplifycover();

        Cover c = Cvrs[m];
        C = &c;
        if (!checkrs()) // 25 seconds becomes 1 second (96% speed-up)
            continue; // envision cutting off an entire branch...

        //c.simplifycover();       // massive slowdown (15%) in some cases, and speedup (10%) in others albeit with many duplicates
        //Cvrs[m].simplifycover();


        std::vector<Edges> ces {};
        ces.clear();
        int csz = c.size();
        std::vector<Edges> ces2 {};
        ces2.resize(csz + 1);
        ces.resize(csz +1);
        int l2 = 0;
        for (int k = 0; k < csz; ++k) {
            Edges es2 = c.getdata(k);
            if (es2.size() >= 3) {
                bool dupe = false;
                for (int l = 0; l < k; ++l)
                    dupe = (dupe || (ces[l] == c[k]));
                if (!dupe) {
                    ces[l2] = es2;
                    ces2[l2] = es2;
                    ++l2;
                }
            }
        }
        ces2.resize(l2);
        bool dupe = false;
        //for (int l = 0; l < l2; ++l)
        //    dupe = (dupe || (ces[l] == es));
        //bool special = false;
        //if (!dupe) {
            if (es.size() >= 3) {
                ces[l2] = es;
                ++l2;
        //        special = true;
            }
        //}
        ces.resize(l2);

        Cover c2 {};
        //if (special) {
            c2.readvector(ces);
        //    c2.simplifycover();    // massive slowdown (33% if both call it)
            CvrsReturn[cri]=c2;;
            ++cri;
        //}
        Cover c3 {};
        c3.readvector(ces2);
        //c3.simplifycover();
        CvrsReturn[cri] = c3;
        ++cri;
    }
    CvrsReturn.resize(cri);
    return CvrsReturn;
}


inline void Hellytheory::findncovers(std::vector<Cover>* Cvar, Cover hintCover, int n) {
    std::vector<Edges> Es;
    int vsz = V->size();
    //if (vsz > 24) {
    //    std::cout << "Unable to process more than 24 vertices.\n";
    //    return;
    //}
    //Es.size() = (1 << vsz);

    //std::cout << "Using maxcliquesize = " << n << "\n";


    // auto-add to hintCover all pure triangles

    std::vector<Edges> hintEs {};
    Cover hintCover2 {};
    for (int k = 0; k < vsz; ++k) {
        std::vector<Edge> N = E->vertexneighborsasedges(V->getdata(k));
        int cnt = 0;
        if (N.size() == 2) {
            Edge e{};
            vertextype v1 = N[0].first;
            vertextype v2 = N[0].second;
            vertextype v3 = ((N[1].first == N[0].first) || (N[1].first == N[0].second)) ? N[1].second : N[1].first;
            if (N[0].first == N[1].first) {
                e.first = N[0].second;
                e.second = N[1].second;
            } else {
                if (N[0].first == N[1].second) {
                    e.first = N[0].second;
                    e.second = N[1].first;
                } else {
                    if (N[0].second == N[1].second) {
                        e.first = N[0].first;
                        e.second = N[1].first;
                    } else {
                        if (N[0].second == N[1].first) {
                            e.first = N[0].first;
                            e.second = N[0].second;
                        }
                    }
                }
            }
            // if N is already in hintCover then ignore
            for (int h = 0; h < hintCover.size(); ++h) {
                if (hintCover.getdata(h).size() == 3) {
                    Edge e2 = hintCover.getdata(h)[0];
                    Edge e3 = hintCover.getdata(h)[1];
                    cnt += e2.first == v1;
                    cnt += e2.first == v2;
                    cnt += e2.first == v3;
                    cnt += e2.second == v1;
                    cnt += e2.second == v2;
                    cnt += e2.second == v3;
                    cnt += e3.first == v1;
                    cnt += e3.first == v2;
                    cnt += e3.first == v3;
                    cnt += e3.second == v1;
                    cnt += e3.second == v2;
                    cnt += e3.second == v3;
                }
            }
            if (cnt < 4) {
                N.push_back(e);
                Edges Etmp{};
                Etmp.readvector(N);
                hintEs.push_back(Etmp);
            }
        }
    }

    for (int m = 0; m < hintCover.size(); ++m) {
        hintEs.push_back(hintCover.getdata(m));
    }
    hintCover2.readvector(hintEs);

    Es.clear();
    //std::cout << "resize " << (1<<n) << "\n";
    //Es.resize(1 << n);


    enumeratecompletesets(&Es, hintCover2, n);

    //TO DO:
    // for each Edges EE in Es
    //    if its size is >= 3
    //       make a cover C and give it hintcover and EE
    //       call simplifycover
    // remove duplicates from this new Es.


    Edges EE{};
    Cover c{};
    std::vector<Edges> Es2{};
//    for (int i = 0; i < Es.size(); ++i) {
//        if (Es[i].size() >= 3) {
//            Es2.push_back(Es[i]);
//            for (int j = 0; j < hintCover.size(); ++j) {
//                Es2.push_back(hintCover.getdata(j));
//            }
//            c.readvector(Es2);
//            //c.simplifycover();

//        }
//    }
    //Es.clear();
    //for (int k = 0; k < c.size(); ++k)
    //    Es.push_back(c.getdata(k));

    for (int m = 0; m < Es.size(); ++m) {
        if (Es[m].size() > 0)
            std::cout << "Edges [";
        else
            std::cout << "Blank edge at m= " << m << "\n";
        for (int j = 0; j < Es[m].size(); ++j) {
            std::cout << "[" << FV->lookup(Es[m][j].first) << ", " << FV->lookup(Es[m][j].second) << "], ";
        }
        if (Es[m].size() > 0)
            std::cout << "\b\b] \n";
    }
    std::vector<Cover> Cvrs = recursefindcovers(&Es, &hintCover2, Es.size());
    for (int n = 0; n < Cvrs.size(); ++n) {
        //Cvrs[n].simplifycover();
        Cvar->push_back((Cvrs[n]));
    }
}

inline void Hellytheory::findrscovers( Cover hintCover ) {
    V->resume();
    E->resume();
    C->resume();
    maxcliquesize = 5;

    //std::vector<Edges> Evar {};
    std::vector<Cover> Cvar {};
    //Evar.resize(1 << V->size());
    //enumeratecompletesets(&Evar, hintCover );
    findncovers(&Cvar, hintCover,maxcliquesize);


    int Cvsz = Cvar.size();
    bool Brs[Cvsz];
    for (int m = 0; m < Cvsz; ++m) {
        C = &(Cvar[m]);
        //FC->matchiedata();
        if (checkcover()) {
            Brs[m] = checkrs();
        } else
            Brs[m] = false;

        if (Brs[m]) {
            std::cout << "Found in edge-complete cover " << m << "\n";
            for (int n = 0; n < C->size(); ++n) {
                Edges es = C->getdata(n);
                int essz = es.size();
                if (essz > 0)
                    std::cout << "...Edges [";
                for (int j = 0; j < essz; ++j) {
                    std::cout << "[" << FV->lookup(es.getdata(j).first) << "," << FV->lookup(es.getdata(j).second)
                              << "], ";
                }
                if (essz > 0)
                    std::cout << "\b\b]\n";
            }
            std::cout << "checkrs() returns true.\n";
        } else {
            //std::cout << "Does not cover all edges\n";
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
    std::cout << cnt << " of " << Cvsz << " passed.\n";

    std::cout << "Checking for duplicates...\n"; // so far provided the recursion isn't using simplifycover,
                                                 // the removal of duplicates turns up empty
                                                 // Moreover, much of the time as recorded is due to this loop;
                                                 // feel free to comment out, therefore, all this duplicate check

    std::vector<Cover> Cvar2 {};
    Cvar2.resize(Cvsz);
    int idx2 = 0;
    for (int m = 0; m < Cvsz; ++m) {
        bool found = false;
        if (!Brs[m])
            continue;
        for (int l = 0; !found && (l < m); ++l) {
            found = (found || (Cvar[m] == Cvar[l]));
        }
        if (!found) {
            Cvar2[idx2] = Cvar[m];
            ++idx2;
        }
    }
    Cvar2.resize(idx2);
    std::cout << "Among those that are RS-covers, found: " << cnt-idx2 << " duplicates...\n";
    if (cnt > idx2) {
        std::cout << "Recap; no dupes:\n";
        for (int m = 0; m < idx2; ++m) {
            std::cout << "Found in edge-complete cover: "<< m << "\n";
            Cover c = Cvar2[m];
            for (int n = 0; n < c.size(); ++n) {
                Edges es = c.getdata(n);
                int essz = es.size();
                if (essz > 0)
                    std::cout << "...Edges [";
                for (int j = 0; j < essz; ++j) {
                    std::cout << "[" << FV->lookup(es.getdata(j).first) << "," << FV->lookup(es.getdata(j).second)
                              << "], ";
                }
                if (essz > 0)
                    std::cout << "\b\b]\n";
            }

        }
    }





}

#endif //HELLYTOOLCPP_HELLYTHEORY_H
