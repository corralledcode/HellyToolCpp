//
// Created by peterglenn on 3/13/24.
//

#ifndef GRAPH_H
#define GRAPH_H



#include <utility>
#include <vector>
#include <sys/types.h>
#include "Batchprocesseddata.cpp"
#include <string>

using vertextype = uint;
using strtype = std::string;


class Vertices : public Batchprocesseddata<vertextype> {
public:
    vertextype maxvertex = 0;
    Vertices(std::vector<vertextype>& verticesin) : Batchprocesseddata<vertextype>(verticesin) {
        pause();
        maxvertex = verticesin.size();
        //for (int n=0;n<maxvertex;++n){
        //    setdata(verticesin[n],n);
        //}
        resume();
    }
    Vertices(int s) : Batchprocesseddata<vertextype>(s) {
        pause();
        setsize(s);
        resume();
    }
    Vertices() : Batchprocesseddata<vertextype>() {}
    ~Vertices() {}
    virtual void process() override;
    virtual void setsize(const int s) override {
        bool p = paused();
        pause();
        Batchprocesseddata<vertextype>::setsize(s);
        auto sz = size();
        for (int i = 0; i < sz; ++i)
            setdata( i, i);
        if (!p)
            resume();
    }
};

class Edge : public std::pair<vertextype,vertextype> {
public:
    bool operator<(const Edge& other) const
    {
        return (first < other.first)
               || ((first == other.first) && (second < other.second));
    }
    bool operator>(const Edge& other) const
    {
        return (first > other.first)
               || ((first == other.first) && (second > other.second));
    }
    bool operator==(const Edge& other) const
    {
        return ((first == other.first && second == other.second)
                || (first == other.second && second == other.first));
    }
};


class Edges : public Batchprocesseddata<Edge> {
public:
    vertextype maxvertex=0;
    vertextype* edgematrix; // adjacency matrix
    vertextype* edgeadjacency; // two dimensional array
    int maxdegree=0;  // the size of each edge's adjacency list
    Edges(std::vector<Edge> edgesin) : Batchprocesseddata<Edge>(edgesin) {
        maxvertex = computemaxvertex();
        maxdegree = computemaxdegree();
        auto sz = size();
        edgematrix = new vertextype[maxdegree*sz];
        edgeadjacency = new vertextype[sz*sz];
        computeadjacencymatrix();
        computeedgematrix();
    }
    Edges(int s) : Batchprocesseddata<Edge>(s) {
        edgematrix = new vertextype[maxdegree*s];
        edgeadjacency = new vertextype[s*s];
    };
    Edges() : Batchprocesseddata<Edge>() {
        edgematrix = new vertextype[0];
        edgeadjacency = new vertextype[0];
    };

    ~Edges() {
        pause();
        //delete[] edgematrix;
        //delete[] edgeadjacency;
    }
    virtual void process() override;
    //virtual void sortdata() override;
    bool operator<(const Edges& other) const
    {
        int n = 0;
        int sz = (size() <= other.size()) ? size() : other.size();
        while ((n<sz) && (getdata(n) == other.getdata(n)))
            n++;
        if (n < sz)
            return getdata(n) < other.getdata(n);
        else
            return false;
    }
    bool operator>(const Edges& other) const
    {
        int n = 0;
        int sz = size() <= other.size() ? size() : other.size();
        while ((n<sz) && (getdata(n) == other.getdata(n)))
            n++;
        if (n < sz)
            return getdata(n) > other.getdata(n);
        else
            return false;
    }
    bool operator==(const Edges& other) const {
        int sz = size();
        if (sz != other.size())
            return false;
        int n = 0;
        bool allfound = true;
        while ((n < sz) && allfound) {
            bool found = false;
            for (int m = 0; (m < sz) && !found; ++m) {
                found = (found || (getdata(n) == other.getdata(m))); // or call sortdata
            }
            allfound = allfound && found;
            ++n;
        }
        return allfound;
    }

    std::vector<vertextype> vertexneighbors(vertextype v) {
        const int sz = size();
        std::vector<vertextype> adjacent {};
        for (int n = 0; n < sz; ++n) {
            Edge e = getdata(n);
            if (e.first != e.second) {
                if (e.first == v)
                    adjacent.push_back(e.second);
                if (e.second == v)
                    adjacent.push_back(e.first);
            }
        }
        return adjacent;
    }

    std::vector<Edge> vertexneighborsasedges(vertextype v) {
        const int sz = size();
        std::vector<Edge> adjacentedges {};
        for (int n = 0; n < sz; ++n) {
            Edge e = getdata(n);
            if (e.first != e.second) {
                if (e.first == v)
                    adjacentedges.push_back(e);
                if (e.second == v)
                    adjacentedges.push_back(e);
            }
        }
        return adjacentedges;
    }

    int vertexdegree(vertextype v) {
        const int sz = size();
        std::vector<vertextype> adjacent {};
        for (int n = 0; n < sz; ++n) {
            Edge e = getdata(n);
            if (e.first != e.second) {
                if (e.first == v)
                    adjacent.push_back(e.second);
                if (e.second == v)
                    adjacent.push_back(e.first);
            }
        }
        return adjacent.size();
    }


private:
    void computeadjacencymatrix();
    void computeedgematrix();
    int computemaxdegree();
    int computemaxvertex() {
        int tempmax = vertextype(0);
        int sz = size();
        for (int n = 0; n < sz; ++n) {
            tempmax = getdata(n).first > tempmax ? getdata(n).first : tempmax;
            tempmax = getdata(n).second > tempmax ? getdata(n).second : tempmax;
        }
        maxvertex = tempmax;
        return maxvertex;
    }
};


class Cover : public Batchprocesseddata<Edges> {
public:
    Cover(std::vector<Edges> edgesin) : Batchprocesseddata<Edges>(edgesin) {
        //
    }
    Cover(int s) : Batchprocesseddata<Edges>(s) {
        //
    }
    Cover() : Batchprocesseddata<Edges>() {
        //
    }
    ~Cover() {
        //
    }

    void process() override {
        bool p = paused();
        pause();
        Batchprocesseddata<Edges>::process();
        int sz = size();
        for (int n = 0;n < sz;++n) {
            getdata(n).process();
        }
        if (!p)
            resume();
    };
    bool coversedges( Edges E ) const {
        int sz = E.size();
        int sz2 = size();
        auto covered = true;
        int n = 0;
        while (covered && n < sz) {
            covered = false;
            int i = 0;
            Edge d = E.getdata(n);
            while (!covered && i < sz2) {
                Edges CE = getdata(i);
                int j = 0;
                int sz3 = E.size();
                while (!covered && (j < sz3)) {
                    covered = (d == CE.getdata(j));
                    ++j;
                }
                ++i;
            }
            ++n;
        }
        return covered;
    }
    void simplifycover();

    bool operator==(const Cover& other) const {
        int sz = size();
        if (sz != other.size())
            return false;
        int n = 0;
        bool allfound = true;
        while ((n < sz) && allfound) {
            bool found = false;
            for (int m = 0; (m < sz) && !found; ++m) {
                found = (found || (getdata(n) == other.getdata(m))); // or call sortdata
            }
            allfound = allfound && found;
            ++n;
        }
        return allfound;
    }


};

class Graph : public Batchprocessed {
public:
    Vertices* vertices;
    Edges* edges;

    Graph(Vertices* verticesin, int vs, Edges* edgesin, int es) : Batchprocessed() {pause(); vertices = verticesin; edges = edgesin; resume();}
    Graph() : Batchprocessed() {}

    void process() override;

};


inline void Vertices::process() {
    //
    Batchprocesseddata<vertextype>::process();
}

inline void Edges::process() {
    //delete[] edgematrix; // getting "double free" message when not commented out
    //delete[] edgeadjacency;
    auto sz = size();
    for (int n = 0; n < sz; ++n) {
        Edge e = getdata(n);
        maxvertex = ((e.first > maxvertex) ? e.first : maxvertex);
        maxvertex = ((e.second > maxvertex) ? e.second : maxvertex);
    }
    maxdegree = computemaxdegree();
    //std::cout << "maxvertex: " << maxvertex << " maxdegree: " << maxdegree << "size()" << size() << "\n";
    edgeadjacency = new vertextype[sz * maxdegree];
    edgematrix = new vertextype[sz * sz];
    Batchprocesseddata<Edge>::process();
}

/* now handled by operator overloading of < and > (see code above)
inline void Edges::sortdata() {
    bool ch = true;
    bool p = paused();
    int sz = size();
    pause();
    while (ch) {
        ch = false;
        for (int i = 0; i < sz - 1; ++i)
            if (getdata(i).first > getdata(i).second) {
                Edge tmp = getdata(i);
                Edge tmp2;
                tmp2.first = tmp.second;
                tmp2.second = tmp.first;
                setdata(tmp2, i);
            }
    }
    if (!p)
        resume(); // implicitly calls inherited sortdata

}
*/

inline int Edges::computemaxdegree() {
    auto szE = size();
    int m = 0;
    auto szV = maxvertex+1;
    auto tally = new int[szV];
    for (int n = 0; n < szV; ++n) {
        tally[n] = 0;
    }
    for (int n = 0; n < szE; ++n) {
        ++tally[getdata(n).first];
        ++tally[getdata(n).second];
    }
    for (int i = 0; i < szV; ++i)
        m = (tally[i] > m) ? tally[i] : m;
    delete[] tally;
    maxdegree = m;
    return maxdegree;
}



inline void Edges::computeadjacencymatrix() {
    //
}

inline void Edges::computeedgematrix() {
    //
}

inline void Cover::simplifycover() { // this should be made instead into a Hellytheory method
    int sz = size();
    std::vector<Edges> Es {};
    Es.clear();
    //sortdata();
    for (int n = 0; n < sz; ++n) {
        Edges es = getdata(n);
        int essz = es.size();
        // is every edge e in es covered by some one es2?
        bool allallcovered = false;
        for (int k = 0; k < sz && !allallcovered; ++k) {
            if (k != n) {
                Edges es2 = getdata(k);
                // does es2 contain every edge e in es?
                int es2sz = es2.size();
                bool allcovered = true;
                for (int m = 0; m < essz && allcovered; ++m) {
                    Edge e = es.getdata(m);
                    bool covered = false;
                    for (int j = 0; j < es2sz && !covered; ++j) {
                        Edge e2 = es2.getdata(j);
                        covered = covered || (e == e2);
                    }
                    allcovered = allcovered && covered;
                }
                allallcovered = allallcovered || allcovered;
            }
        }
        if (!allallcovered) {
            Es.push_back(es);
        }
#ifdef VERBOSE
        else {
            std::cout << "Removing covered edges ";
            for (int m = 0; m < es.size(); ++m)
                std::cout << "[" << es.getdata(m).first << ", " << es.getdata(m).second << "], ";
            std::cout << "\b\b  \n";
        }
#endif

    }
    readvector(Es);
}

#endif //GRAPH_H
