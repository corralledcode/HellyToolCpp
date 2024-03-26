//define VERBOSE   // compile with verbose output enabled
//define VERBOSE2
#define CHECKS    // compile with (slow) checks in checkrs routine and elsewhere
#define MONITORRECURSIONDEPTH
//define VERBOSETRIANGLE
#define VERBOSESIMPLIFYCOVER
#include <iostream>
#include <fstream>
#include <regex>
#include <chrono>
#include "Graph.cpp"
#include "Formatgraph.cpp"
#include "Hellytheory.cpp"



int main(int argc, char *argv[]) {
    std::cout << "Hello, World!" << std::endl;

    std::ifstream ifs;
    std::istream* s = &std::cin;
    if (argc > 1) {
        std::string filename = argv[1];
        std::cout << "Opening file " << filename << "\n";
        ifs.open(filename);
        if (!ifs)
            std::cout << "Couldn't open file for reading \n";
        s = &ifs;
    } else {
        std::cout << "Enter a filename or enter T for terminal mode: ";
        std::string filename;
        std::cin >> filename;
        if (filename != "T") {
            ifs.open(filename);
            if (!ifs)
                std::cout << "Couldn't open file for reading \n";
            s = &ifs;
        }
    }

    auto V = new Vertices();
    auto EV = new Batchprocesseddata<strtype>();
    auto FV = new Formatvertices(V,EV);
    FV->pause();
    FV->readdata(*s);
    FV->resume();
    int sV = FV->size();
    for( int n = 0; n<sV; ++n)
        std::cout << n << "{" << FV->idata->getdata(n) << ":" << FV->edata->getdata(n) << "}, ";
    std::cout << "\b\b  \n";

    auto E = new EdgesforHelly();
    auto EE = new Batchprocesseddata<Edgestr>();
    auto FE = new Formatedges(E,EE);
    FE->pause();
    FE->readdata(*s);
    FE->setvertices(FV);
    FE->resume();
    int sE = FE->size();

    for (int m = 0; m < sE; ++m) {
        std::cout << "{" << FE->idata->getdata(m).first << ", " << FE->idata->getdata(m).second;
        std::cout << ":" << FE->edata->getdata(m).first << ", " << FE->edata->getdata(m).second << "}, ";
    }
    std::cout << "\b\b  \n";
    //std::cout << "computemaxdegree: " << E->computemaxdegree() << "\n";

    auto C = new Cover();
    auto EC = new Batchprocesseddata<std::vector<Edgestr>>();
    auto FC = new Formatcover(C,EC);
    FC->pause();
    FC->readdata(*s);
    FC->setvertices(FV);
    //FC->simplifycover();
    FC->resume();
    int sC = FC->size();
    for (int m = 0; m < sC; ++m) {
        int sc = FC->idata->getdata(m).size();
        std::cout << "C[" << m << "]: ";
        for (int n = 0; n < sc; ++n) {
            std::cout << FC->idata->getdata(m)[n] << ": " << FC->edata->getdata(m)[n] << ", ";
        }
        std::cout << "\b\b  \n";
    }
    std::cout << "\n";

    auto T = new Hellytheory();
    V->pause();
    E->pause();
    C->pause();
    FV->pause();
    FE->pause();
    FC->pause();
    T->V = V;
    T->E = E;
    T->C = C;
    T->FV = FV;
    T->FE = FE;
    T->FC = FC;

    auto starttime = std::chrono::high_resolution_clock::now();
    bool rscover = T->checkrs();
    std::cout << "checkrs returns " << rscover << "\n";
    auto stoptime = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::microseconds>(stoptime - starttime);
    std::cout << "Time elapsed: " << float(duration.count())/1000000 << "\n";

    Cover hintCover {};
    std::vector<Edges> es {};
    for (int i = 0; i < (C->size()); ++i) {
        es.push_back(C->getdata(i));
    }
    hintCover.readvector(es);
    hintCover.simplifycover();
    starttime = std::chrono::high_resolution_clock::now();

    T->findrscovers(hintCover);
    stoptime = std::chrono::high_resolution_clock::now();
    duration = duration_cast<std::chrono::microseconds>(stoptime - starttime);
    std::cout << "Time elapsed: " << float(duration.count())/1000000 << "\n";

    return 0;
}
