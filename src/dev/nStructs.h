#ifndef NS_H
#define NS_H
#include <vector>
#include <random>
#include <fstream>
#include <array>
#include "NeuroLibrary.h"
#include "typedefs.h"

using std::vector;
using std::cout;
using std::endl;

struct ijr {
    // func(x) = func(i) + r*(func(j)-func(i));
    size i;
    size j;
    double r;
    ijr(){};
    ijr(double i0, double j0, double r0) : i(i0), j(j0), r(r0) { };
    ijr operator+(struct ijr ijr_input);
    ijr operator+(double r0);
    ijr operator+(size i0);
};
typedef struct ijr IJR;

struct jumpyNeuronData {
    std::vector<double> t,v;
    jumpyNeuronData(size rSize);
    void initialize(size rSize);
};
typedef struct jumpyNeuronData jND;

struct CrossData {
    // index 0 for vinit
    // counting from 1.
    size nCross;
    std::vector<size> iCross; // next Cross index
    std::vector<IJR> vCross;
    std::vector<double> tCross; // t of crossing back
    std::vector<double> t,v,gE,gI,m,n,h,hE,hI;
    CrossData(size nt, double vinit);
    void initialize(size nt, double vinit);
};
typedef struct CrossData Cross;

struct BilinearRelationships {
    // dt to the second input;
    std::vector<IJR> dTijr; 
    std::vector<size> idt; 
    std::vector<unsigned int> ID;
    BilinearRelationships(size corrSize);
};
typedef struct BilinearRelationships biR;

struct Inputs {
    // 0 for sample tstep;
    // 1 for demanded tstep;
    std::vector<IJR> Vijr;
    // dT for cross and K
    std::vector<IJR> dTijr,Tmpijr;
    // current index of effecting cross
    std::vector<size> cCross;
    std::vector<double> t,dt,tMax;
    std::vector<biR> bir;
    std::vector<size> ID;
    Inputs(size rSize);
    void initialize(size rSize);
    void junk_this();
    void assert_size();
    void print_this(int i);
};
typedef struct Inputs Input;
