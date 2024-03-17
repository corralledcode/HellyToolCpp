//
// Created by peterglenn on 3/13/24.
//

#ifndef HELLYTOOLCPP_FORMATGRAPH_CPP
#define HELLYTOOLCPP_FORMATGRAPH_CPP

#include <string>
#include <iostream>
#include <ostream>
#include <vector>
#include <c++/11/regex>
#include "Graph.cpp"
#include "Formatdata.cpp"

//{"([a-zA-Z]{1}[\\d_]*)"}

class Formatvertices : public Formatdata<vertextype,strtype,Vertices,Batchprocesseddata<strtype>> {
public:
    Formatvertices( Vertices* iin, Batchprocesseddata<strtype>* ein)
        : Formatdata<vertextype,strtype,Vertices,Batchprocesseddata<strtype>>(*iin,*ein) {
        //
    };
    Formatvertices( Vertices* iin ) : Formatdata<vertextype, strtype, Vertices, Batchprocesseddata<strtype>>(*iin) {
        //
    }
    explicit Formatvertices(std::istream &is) : Formatdata<vertextype,strtype,Vertices,Batchprocesseddata<strtype>>(is) {
        //
    }
    Formatvertices(int s) : Formatdata<vertextype,strtype,Vertices,Batchprocesseddata<strtype>>(s) {
        //
    }
    Formatvertices() : Formatdata<vertextype,strtype,Vertices,Batchprocesseddata<strtype>>() {
        //
    }
    ~Formatvertices() {
        //Formatdata<vertextype,strtype,Vertices,Batchprocesseddata<strtype>>::~Formatdata();
    }

protected:
    void parseexternal(const std::string str, std::vector<strtype>* out) override {
        out->push_back(str);
    }

    void parseinternal(const std::string str, std::vector<vertextype>* out) override {
        out->push_back(Formatdata<vertextype,strtype,Vertices,Batchprocesseddata<strtype>>::lookup(str));
    }
};

class Edgestr : public std::pair<strtype,strtype> {
public:
    bool operator<(const Edgestr& other) const
    {
        return (first < other.first)
               || ((first == other.first) && (second < other.second));
    }
    bool operator>(const Edgestr& other) const
    {
        return (first > other.first)
               || ((first == other.first) && (second > other.second));
    }
    bool operator==(const Edgestr& other) const
    {
        return ((first == other.first && second == other.second)
                || (first == other.second && second == other.first));
    }
};

inline std::ostream& operator<<(std::ostream& os, const Edgestr& e) {
    return os << "[" << e.first << ", " << e.second << "]";
}

inline std::ostream& operator<<(std::ostream& os, const Edge& e) {
    return os << "[" << e.first << ", " << e.second << "]";
}

inline std::ostream& operator<<(std::ostream& os, const Edges& e) {
    int sz = e.size();
    os << "[";
    for (int n = 0; n < sz; ++n) {
        os << "[" << e.getdata(n).first << "," << e.getdata(n).second << "]";
        if (n < sz-1)
            os << ",";
    }
    return os;
}


/*
inline std::istream& operator>>(std::istream& is, Edgestr& e) {

}
*/


class Formatedges : public Formatdata<Edge,Edgestr,Edges,Batchprocesseddata<Edgestr>> {
public:
    Formatedges( Edges* iin, Batchprocesseddata<Edgestr>* ein )
        : Formatdata<Edge,Edgestr,Edges,Batchprocesseddata<Edgestr>>( *iin, *ein ) {
        //
    }
    Formatedges( Edges* iin ) : Formatdata<Edge,Edgestr,Edges,Batchprocesseddata<Edgestr>>(*iin) {
        //
    }
    Formatedges(std::istream &is) : Formatdata<Edge,Edgestr,Edges, Batchprocesseddata<Edgestr>>(is) {
        //
    }
    Formatedges(int s) : Formatdata<Edge,Edgestr,Edges, Batchprocesseddata<Edgestr>>(s) {
        //
    }
    Formatedges() : Formatdata<Edge,Edgestr,Edges,Batchprocesseddata<Edgestr>>() {
        //
    }
    ~Formatedges() {
    }
    Formatvertices* FV;  // Formatedges is not responsible for memory management of this pointer
    void setvertices( Formatvertices* FVin ) {
        bool b = paused();
        pause();
        int sz = size();
        FV = FVin;
        for (int n = 0; n < sz; ++n) {
            Edge e;
            Edgestr es = edata->getdata(n);
            e.first = FV->lookup(es.first);
            e.second = FV->lookup(es.second);
            idata->setdata(e,n);
        }
    }

protected:
    void parseexternal(const std::string str, std::vector<Edgestr>* edge) override {
        Edgestr e;
        std::regex pat {"([a-zA-Z]{1}[\\d_]*)"};
        int n = 0;
        std::vector<strtype> v;
        for (std::sregex_iterator p(str.begin(),str.end(),pat); p != std::sregex_iterator{};++p)
            v.push_back((*p)[1]);
        std::sort(v.begin(), v.end());
        for (int m = 0; m < v.size(); ++m) {
            for (int n = m+1; n < v.size(); ++n) {
                e.first = v[m];
                e.second = v[n];
                edge->push_back(e);
                //std::cout << "v[m]: " << v[m] << " v[n]: " << v[n] << "\n";
            }
        }
    }

    void parseinternal(const std::string str, std::vector<Edge>* edge) override {
        Edgestr es;
        std::regex pat {"([a-zA-Z]{1}[\\d_]*)"};
        int n = 0;
        std::vector<strtype> v;
        for (std::sregex_iterator p(str.begin(),str.end(),pat); p != std::sregex_iterator{};++p)
            v.push_back((*p)[1]);
        std::sort(v.begin(), v.end());
        for (int m = 0; m < v.size(); ++m) {
            for (int n = m+1; n < v.size(); ++n) {
                es.first = v[m];
                es.second = v[n];
                //std::cout << ":: v[m]: " << v[m] << " v[n]: " << v[n] << "\n";
                edge->push_back(lookup(es));
            }
        }
    }
};

class Formatcover : public Formatdata<Edges,Edgestr*,Cover,Batchprocesseddata<Edgestr*>> {
public:
    Formatcover( Cover* iin, Batchprocesseddata<Edgestr*>* ein )
            : Formatdata<Edges,Edgestr*,Cover,Batchprocesseddata<Edgestr*>>( *iin, *ein ) {
        //
    }
    Formatcover( Cover* iin ) : Formatdata<Edges,Edgestr*,Cover,Batchprocesseddata<Edgestr*>>( *iin ) {
        //
    }
    Formatcover(std::istream &is) : Formatdata<Edges,Edgestr*,Cover,Batchprocesseddata<Edgestr*>>(is) {
        //
    }
    Formatcover(int s) : Formatdata<Edges,Edgestr*,Cover,Batchprocesseddata<Edgestr*>>(s) {
        //
    }
    Formatcover() : Formatdata<Edges,Edgestr*,Cover,Batchprocesseddata<Edgestr*>>() {
        //
    }
    ~Formatcover() {
    }
    Formatvertices* FV;  // Formatedges is not responsible for memory management of this pointer
    void setvertices( Formatvertices* FVin ) {
        bool b = paused();
        pause();
        int sz = size();
        FV = FVin;
        for (int n = 0; n < sz; ++n) {
            Edges E = idata->getdata(n);
            Edgestr* ES = edata->getdata(n);
            int sz2 = E.size();
            for (int i = 0; i < sz2; ++i) {
                Edge e = E.getdata(i);
                Edgestr es = ES[i];
                e.first = FV->lookup(es.first);
                e.second = FV->lookup(es.second);
                E.setdata(e, n);
            }
        }
    }

protected:
    void parseexternal(const std::string str, std::vector<Edgestr*>* edges ) override {
        std::regex pat {"([a-zA-Z]{1}[\\d_]*)"};
        int n = 0;
        std::vector<strtype> v;
        for (std::sregex_iterator p(str.begin(),str.end(),pat); p != std::sregex_iterator{};++p)
            v.push_back((*p)[1]);
        std::sort(v.begin(), v.end());
        int j = 0;
        int sz = v.size();
        Edgestr* e[sz*(sz-1)];
        for (int k = 0; k < sz*(sz-1); ++k) {
            e[k] = new Edgestr;
        }
        for (int m = 0; m < sz; ++m) {
            for (int n = m+1; n < sz; ++n) {
                e[j]->first = v[m];
                e[j]->second = v[n];
                ++j;
                //std::cout << "v[m]: " << v[m] << " v[n]: " << v[n] << "\n";
            }
        }
        edges->push_back(*e);
    }

    void parseinternal(const std::string str, std::vector<Edges>* edges) override {
        std::regex pat {"([a-zA-Z]{1}[\\d_]*)"};
        int n = 0;
        std::vector<strtype> v;
        for (std::sregex_iterator p(str.begin(),str.end(),pat); p != std::sregex_iterator{};++p)
            v.push_back((*p)[1]);
        std::sort(v.begin(), v.end());
        int j = 0;
        int sz = v.size();
        Edgestr* e[sz*(sz-1)];
        for (int k = 0; k < sz*(sz-1); ++k) {
            e[k] = new Edgestr;
        }
        for (int m = 0; m < v.size(); ++m) {
            for (int n = m+1; n < v.size(); ++n) {
                e[j]->first = v[m];
                e[j]->second = v[n];
                //std::cout << ":: v[m]: " << v[m] << " v[n]: " << v[n] << "\n";
                ++j;
            }
        }
        // lookup which edge matches e (until refactoring the lookup stuff)

        for (int m = 0; m < idata->size(); ++m) {
            Edgestr* etemp = edata->getdata(m);
            if (j == edata->size()) {
                bool match = false;
                bool allmatch = true;
                for (int k = 0; k < j; ++k) {
                    for (int l = 0; l < j; ++l) {
                        match = match ||
                                (etemp[k].first == e[l]->first && etemp[k].second == e[l]->second) ||
                                (etemp[k].second == e[l]->first && etemp[k].first == e[l]->second);
                    }
                    allmatch = allmatch && match;
                }
                if (allmatch)
                    edges->push_back(idata->getdata(m));
            }
        }
    }

    void cleanup() {

    }

};





#endif //HELLYTOOLCPP_FORMATGRAPH_CPP
