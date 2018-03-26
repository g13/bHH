#include "nStructs.h"

IJR ijr::operator+(IJR ijr_input) {
    IJR ijr_output;
    double r0 = r + ijr_input.r;
    size ir = static_cast<size>(r0);
    ijr_output.i = i + ijr_input.i + ir;
    ijr_output.r = r0-ir; 
    ijr_output.j = ijr_output.i+1;
    return ijr_output;
}
IJR ijr::operator+(double r0) {
    IJR ijr_output;
    double r1 = r + r0;
    size ir = static_cast<size>(r1);
    ijr_output.i = i + ir;
    ijr_output.r = r1-ir; 
    ijr_output.j = ijr_output.i+1;
    return ijr_output;
}
IJR ijr::operator+(size i0) {
    IJR ijr_output;
    ijr_output.i = i + i0;
    ijr_output.r = r; 
    ijr_output.j = ijr_output.i+1;
    return ijr_output;
}

jumpyNeuronData::jumpyNeuronData(size rSize) {
    t.reserve(2*rSize);
    v.reserve(2*rSize);
}
void jumpyNeuronData::initialize(size rSize) {
    t.clear();
    v.clear();
    t.reserve(2*rSize);
    v.reserve(2*rSize);
}

CrossData::CrossData(size nt, double vinit) {
    nCross = 0;
    iCross.reserve(nt/2);
    tCross.reserve(nt/2);
    vCross.reserve(nt/2);
    tCross.push_back(0);
    v.reserve(nt);
    t.reserve(nt);
    iCross.push_back(1);
    v.push_back(vinit);
    t.push_back(0);
}
void CrossData::initialize(size nt, double vinit) {
    nCross = 0;
    iCross.clear();
    iCross.reserve(nt/2);
    iCross.push_back(1);
    tCross.clear();
    tCross.reserve(nt/2);
    tCross.push_back(0);
    vCross.clear();
    vCross.reserve(nt/2);
    v.clear();
    v.reserve(nt);
    v.push_back(vinit);
    t.clear();
    t.reserve(nt);
    t.push_back(0);
}

BilinearRelationships::BilinearRelationships(size corrSize) {
    dTijr.reserve(corrSize);
    idt.reserve(corrSize);
    ID.reserve(corrSize);
}

Inputs::Inputs(size rSize){
    t.reserve(rSize); 
    dt.reserve(rSize); 
    tMax.reserve(rSize); 
    cCross.reserve(rSize);
    dTijr.reserve(rSize); 
    Vijr.reserve(rSize);
    ID.reserve(rSize);
    Tmpijr.reserve(rSize); 
    bir.reserve(rSize);
}
void Inputs::initialize(size rSize) {
    t.clear(); 
    dt.clear(); 
    tMax.clear(); 
    cCross.clear();
    dTijr.clear(); 
    Vijr.clear();
    ID.clear();
    Tmpijr.clear(); 
    bir.clear();
    t.reserve(rSize); 
    dt.reserve(rSize); 
    tMax.reserve(rSize); 
    cCross.reserve(rSize);
    dTijr.reserve(rSize); 
    Vijr.reserve(rSize);
    ID.reserve(rSize);
    Tmpijr.reserve(rSize); 
    bir.reserve(rSize);
}
void Inputs::junk_this() {
    t.push_back(0); 
    ID.push_back(-1);
    dt.push_back(0);
    tMax.push_back(0); 
    cCross.push_back(0);
    dTijr.push_back(IJR(0,1,0)); 
    Vijr.push_back(IJR(0,1,0));
    Tmpijr.push_back(IJR(0,1,0)); 
    bir.push_back(BilinearRelationships(0));
}   
void Inputs::assert_size() {
    size Size = t.size(); 
    assert(dt.size() == Size); 
    assert(tMax.size() == Size); 
    assert(cCross.size() == Size);
    assert(dTijr.size() == Size); 
    if (Vijr.size() != Size) {
        cout << Vijr.size() << "!=" << Size << endl;
        assert(Vijr.size() == Size);
    }
    assert(Tmpijr.size() == Size); 
    assert(bir.size() == Size);
}   
void Inputs::print_this(int i) {
    cout << "t: " << t[i] << endl;
    cout << "ID: " << ID[i] << endl;
    cout << "V index: " << Vijr[i].i << ", " << Vijr[i].j << ", " << Vijr[i].r << endl;
    cout << "dT index: " << dTijr[i].i << ", " << dTijr[i].j << ", " << dTijr[i].r << endl;
}
