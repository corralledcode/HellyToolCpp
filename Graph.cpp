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
    Vertices(std::vector<vertextype>& verticesin) : Batchprocesseddata<vertextype>(verticesin.size()) {
        pause();
        maxvertex = verticesin.size();
        for (int n=0;n<maxvertex;++n){
            setdata(verticesin[n],n);
        }
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
        bool b = paused();
        pause();
        Batchprocesseddata<vertextype>::setsize(s);
        auto sz = size();
        for (int i = 0; i < sz; ++i)
            setdata(vertextype(i),i);
        if (!b)
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
    vertextype maxvertex;
    vertextype* edgematrix; // adjacency matrix
    vertextype* edgeadjacency; // two dimensional array
    int maxdegree;  // the size of each edge's adjacency list
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
        delete[] edgematrix;
        delete[] edgeadjacency;
    }
    virtual void process() override;
    virtual void sortdata() override;
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
        auto b = true;
        int n = 0;
        while (b && n < sz) {
            b = b && (getdata(n) == other.getdata(n)); // rewrite to be immune to unordered situation... or call sortdata
            ++n;
        }
        return b;
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
        Batchprocesseddata<Edges>::process();
        int sz = size();
        for (int n = 0;n < sz;++n) {
            getdata(n).process();
        }
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
    delete[] edgematrix;
    delete[] edgeadjacency;
    auto sz = size();
    int maxvertex = 0;
    for (int n = 0; n < sz; ++n) {
        maxvertex = getdata(n).first > maxvertex ? getdata(n).first : maxvertex;
        maxvertex = getdata(n).second > maxvertex ? getdata(n).second : maxvertex;
    }
    maxdegree = computemaxdegree();
    edgeadjacency = new vertextype[size() * maxdegree];
    edgematrix = new vertextype[size() * size()];
    Batchprocesseddata<Edge>::process();
}

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

inline int Edges::computemaxdegree() {
    auto szE = size();
    int m = 0;
    auto szV = maxvertex;
    auto tally = new int[szV];
    for (int n = 0; n < szE; ++n) {
        tally[getdata(n).first]++;
        tally[getdata(n).second]++;
    }
    for (int i = 0; i < szV; ++i)
        m = tally[i] > m ? tally[i] : m;
    delete[] tally;
    maxdegree = m;
    return maxdegree;
}



inline void Edges::computeadjacencymatrix() {
    //
}

#endif //GRAPH_H
