//
// Created by peterglenn on 3/13/24.
//

#ifndef HELLYTOOLCPP_FORMATDATA_CPP
#define HELLYTOOLCPP_FORMATDATA_CPP

#include <string>
#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <vector>
#include <regex>
#include "Batchprocesseddata.cpp"

template<typename I, typename E, typename IBPD, typename EBPD>
class Formatdata : public Batchprocessed {
private:
    int size_ = 0;
protected:
    virtual void parseexternal(const std::string str, std::vector<E>*) = 0;
    virtual void parseinternal(const std::string str, std::vector<I>*) = 0;

public:
    IBPD* idata = nullptr;
    EBPD* edata = nullptr;

    E lookup(I di);
    I lookup(E de);
    int size() {
        if (paused())
            return size_;
        else
            return idata->size()>=edata->size() ? idata->size() : edata->size();
    }

    virtual void readdata(std::istream& is) {
        auto p = paused();
        pause();
        if (!idata || !edata)
            throw std::exception();
        int s=0;
        std::vector<E> eres;
        std::string tmp = "";
        std::string tmp2 = "";
        while ((is >> tmp) && (tmp != "END")) {
            parseexternal(tmp,&eres);
            tmp2 += tmp + " ";
            tmp = "";
            s++;
        }
        edata->readvector(eres);
        size_ = edata->size(); // may be larger than s due to condensed parsing
        idata->setsize(size_);
        std::vector<I> ires;
        std::regex pat {"([\\w]+)"};
        int n = 0;
        for (std::sregex_iterator p(tmp2.begin(),tmp2.end(),pat); p != std::sregex_iterator{};++p) {
            parseinternal((*p)[1],&ires);
            ++n;
        }
        idata->readvector(ires);
        if (!p) {
            resume();
            edata->resume();
            idata->resume();
        }
    }

    virtual void writedata(std::ostream& os) {
        auto sz = size();
        for (int n = 0; n < sz; ++n) {
            os << "{" << idata->getdata(n) << ":" << edata->getdata(n) << "}, ";
        }
        os << "\n";
    }

    Formatdata(IBPD& iin, EBPD& ein) : Batchprocessed() {
        size_ = iin.size();
        idata = &iin;
        edata = &ein;
    }

    Formatdata(IBPD& iin) : Batchprocessed() {
        idata = &iin;
    }
    Formatdata(int s) : Batchprocessed() {
        size_ = s;
    }
    Formatdata() : Batchprocessed() {
        idata = nullptr;
        edata = nullptr;
        size_ = 0;
    }
    Formatdata(std::istream& is) : Batchprocessed() {
        if (idata && edata)
            readdata(is);
        else
            throw std::exception();
    }
    ~Formatdata() {
        size_ = 0;
        //Batchprocessed::~Batchprocessed();
    }
    void process() override {
        if (idata && edata) {
            idata->process();
            edata->process();
        }
        Batchprocessed::process();
    };

};

template<typename I, typename E, typename IBPD, typename EBPD>
I Formatdata<I,E,IBPD,EBPD>::lookup(E ve) {
    bool found = false;
    int n = 0;
    int sz = size();
    if (sz <= 0)
        throw std::out_of_range("I lookup: Cannot find if size is zero");
    E e = edata->getdata(0);
    //std::cout << "ve: " << ve << "e: " << e << "\n";
    found = (ve == e);
    while (!found && n < sz-1) {
        n++;
        e = edata->getdata(n);
        //std::cout << "ve: " << ve << "e: " << e << " sz: " << sz << "\n";
        found = (ve == e);
    }
    if (found)
        return idata->getdata(n);
    else
        throw std::out_of_range("I lookup: Unknown data element");
}

template<typename I, typename E, typename IBPD, typename EBPD>
E Formatdata<I,E,IBPD,EBPD>::lookup(I vi) {
    bool found = false;
    int n = 0;
    int sz = size();
    if (sz <= 0)
        throw std::out_of_range("E lookup: Cannot find if size is zero");
    I i = idata->getdata(0);
    found = (vi == i);
    while (!found && n < sz-1) {
        n++;
        i = idata->getdata(n);
        found = (vi == i);
    }
    if (found)
        return edata->getdata(n);
    else
        throw std::out_of_range("E lookup: Unknown data element");
}

#endif //HELLYTOOLCPP_FORMATDATA_CPP
