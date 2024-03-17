#include <iostream>
#include <fstream>
#include <regex>
#include "Graph.cpp"
#include "Formatgraph.cpp"


int main(int argc, char *argv[]) {
    std::cout << "Hello, World!" << std::endl;
/*
    auto bpd = new Batchprocesseddata<int>(10);
    bpd->process();
    bpd->setdata(13,0);
    for (int i = 0; i < bpd->size(); ++i) {
        std::cout << bpd->getdata(i) << '\n';
    }

    uint ar[] {0,1,2,3,4,5,6,7,8,9};

    auto V = new Vertices(ar,10);

    V->pause();
    for (int i = 0; i < 10; ++i) {
        V->setdata(10-i,i);
    }
    std::cout << '\n';
    V->resume();
    for (int i = 0; i < 10; ++i) {
        std::cout << V->getdata(i) << ", ";
    }
    std::cout << "\n";

    auto FV = new Formatvertices(10);
    FV->pause();
    for (int n = 0; n < FV->size();++n) {
        auto lp = new lookuppair<vertextype>;
        std::cin >> lp->sr;
        lp->ir = n;
        FV->setdata({lp->ir,lp->sr},n);
    }
    FV->resume();

    for (int n = 0; n < FV->size();++n)
        std::cout << '{' << FV->getdata(n).sr << ',' << FV->getdata(n).ir << '}' << ", ";
    std::cout << '\n';
    for (int n = 0; n < FV->size();++n)
        std::cout << FV->lookup(FV->getdata(n).sr ) << ", ";
*/

    std::ifstream ifs;
    if (argc > 1) {
        std::string filename = argv[1];
        std::cout << "Opening file " << filename << "\n";
        ifs.open(filename);
        if (!ifs)
            std::cout << "Couldn't open file for reading \n";
    }

    //std::cout << "How many vertices? ";
    //std::cin >> sV;
    auto V = new Vertices();
    auto EV = new Batchprocesseddata<strtype>();
    auto FV = new Formatvertices(V,EV);
    FV->pause();
    FV->readdata(ifs);
    FV->resume();
    int sV = FV->size();
    for( int n = 0; n<sV; ++n)
        std::cout << n << "{" << FV->idata->getdata(n) << ":" << FV->edata->getdata(n) << "}, ";
    std::cout << "\n";
    //std::cout << "How many edges? ";
    //std::cin >> sE;
/*
    auto E = new Edges(10);
    E->pause();
    Edge e;
    e.first = 0;
    e.second = 1;
    E->setdata(e,0);
    e.first = 3;
    e.second = 2;
    E->setdata(e,1);
    for (int n = 0;n<E->size();++n) {
        std::cout << "n: " << n << " " << E->getdata(n) << "\n";
    }
    E->resume();
    for (int n = 0;n<E->size();++n) {
        std::cout << "n: " << n << " " << E->getdata(n) << "\n";
    }*/


    auto E = new Edges();
    auto EE = new Batchprocesseddata<Edgestr>();
    auto FE = new Formatedges(E,EE);
    FE->pause();
    FE->readdata(ifs);
    FE->setvertices(FV);
    FE->resume();
    int sE = FE->size();

    for (int m = 0; m < sE; ++m) {
        std::cout << "{" << FE->idata->getdata(m).first << ", " << FE->idata->getdata(m).second;
        std::cout << ":" << FE->edata->getdata(m).first << ", " << FE->edata->getdata(m).second << "}, ";
    }
    std::cout << "\n";
    //std::cout << "computemaxdegree: " << E->computemaxdegree() << "\n";

    auto C = new Cover();
    auto EC = new Batchprocesseddata<Edgestr*>();
    auto FC = new Formatcover(C,EC);
    FC->pause();
    FC->readdata(ifs); //replace with std::cin for a console input
    FC->setvertices(FV);
    FC->resume();
    int sC = FC->size();
    for (int m = 0; m < sE; ++m) {
        int sc = FC->idata->getdata(m).size();
        std::cout << "C[" << m << "]: ";
        for (int n = 0; n < sc; ++n) {
            std::cout << FC->edata->getdata(m)[n];
        }
    }
    std::cout << "\n";



    return 0;
}
