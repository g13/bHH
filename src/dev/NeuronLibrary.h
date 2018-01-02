#ifndef NEURONLIBRARY_H
#define NEURONLIBRARY_H
#include "matFunc.h"
#include <cmath>

typedef struct NeuronLibrary 
{
    size ndt,nt,nt0,nv,nE,nI,ith,l0;

    double  ****EPSP,   ***EPSP0,  ****IPSP,   ***IPSP0,
           *EPSP_ptr, *EPSP0_ptr, *IPSP_ptr, *IPSP0_ptr,
             ***kVEE,    ***kVEI,   ***kVIE,    ***kVII,
           *kVEE_ptr,  *kVEI_ptr, *kVIE_ptr,  *kVII_ptr,
             *vRange,   *dtRange, *fE, *fI, dur, tstep,
             **vLeak, *vLeak_ptr;
    //double vReset, vThres;
    size *idtRange;
    const char *file;
    bool ei;

    void readLib(const char *filename) {
        MATFile *pmat;
        const char **var;
        mxArray *tmp;
        size arraySize, dimSize[4];
        file = filename;
        int i,j,k,l;
        int n;
    
        openMat(pmat, file);
        
        var = (const char **)matGetDir(pmat, &n);
        if (var == NULL) {
            std::cout << "Error reading content of file: " << file;
            mxFree(var);
            abort();
        } 
        std::cout << "Reading file " << file << "... " << std::endl;
        std::cout << n << " variable(s):" << std::endl;
        for (i=0; i<n; i++) {
            std::cout << var[i];
            if (i<n-1) std::cout << ", ";
        }
        std::cout << std::endl;
        mxFree(var);
    
        readArray(tmp,"sEPSP", dimSize, arraySize, pmat, file);
        EPSP_ptr = new double[arraySize];
        memcpy((void *) EPSP_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer4d(EPSP,EPSP_ptr,dimSize);
        mxDestroyArray(tmp);
    
        readArray(tmp,"sIPSP", dimSize, arraySize, pmat, file);
        IPSP_ptr = new double[arraySize];
        memcpy((void *) IPSP_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer4d(IPSP,IPSP_ptr,dimSize);
        mxDestroyArray(tmp);
    
        nt = dimSize[3];
    
        readArray(tmp,"sEPSP0", dimSize, arraySize, pmat, file);
        EPSP0_ptr = new double[arraySize];
        memcpy((void *) EPSP0_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer3d(EPSP0,EPSP0_ptr,dimSize);
        mxDestroyArray(tmp);
    
        readArray(tmp,"sIPSP0", dimSize, arraySize, pmat, file);
        IPSP0_ptr = new double[arraySize];
        memcpy((void *) IPSP0_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer3d(IPSP0,IPSP0_ptr,dimSize);
        mxDestroyArray(tmp);
    
        nt0 = dimSize[2];
    
        readArray(tmp,"kVEE", dimSize, arraySize, pmat, file);
        kVEE_ptr = new double[arraySize];
        memcpy((void *) kVEE_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer3d(kVEE,kVEE_ptr,dimSize);
        mxDestroyArray(tmp);
    
        readArray(tmp,"kVEI", dimSize, arraySize, pmat, file);
        kVEI_ptr = new double[arraySize];
        memcpy((void *) kVEI_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer3d(kVEI,kVEI_ptr,dimSize);
        mxDestroyArray(tmp);
    
        readArray(tmp,"kVIE", dimSize, arraySize, pmat, file);
        kVIE_ptr = new double[arraySize];
        memcpy((void *) kVIE_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer3d(kVIE,kVIE_ptr,dimSize);
        mxDestroyArray(tmp);
    
        readArray(tmp,"kVII", dimSize, arraySize, pmat, file);
        kVII_ptr = new double[arraySize];
        memcpy((void *) kVII_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer3d(kVII,kVII_ptr,dimSize);
        mxDestroyArray(tmp);
    
        readArray(tmp,"vleakage0", dimSize, arraySize, pmat, file);
        vLeak_ptr = new double[arraySize];
        memcpy((void *) vLeak_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer2d(vLeak,vLeak_ptr,dimSize);
        mxDestroyArray(tmp);

        readArray(tmp,"vRange", dimSize, arraySize, pmat, file);
        vRange = new double[arraySize];
        memcpy((void *) vRange, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        //std::cout << "vRange:" << std::endl; disp1d(vRange,dimSize[0]);
        mxDestroyArray(tmp);
        
        nv = arraySize;
    
        readArray(tmp,"dtRange", dimSize, arraySize, pmat, file);
        dtRange = new double[arraySize];
        memcpy((void *) dtRange, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        //std::cout << "dtRange:" << std::endl; disp1d(dtRange,dimSize[0]);
        mxDestroyArray(tmp);
        
        ndt = arraySize;
    
        readArray(tmp,"fE", dimSize, arraySize, pmat, file);
        fE = new double[arraySize];
        memcpy((void *) fE, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        //std::cout << "fE:" << std::endl; disp1d(fE,dimSize[0]);
        mxDestroyArray(tmp);
    
        nE = arraySize;
    
        readArray(tmp,"fI", dimSize, arraySize, pmat, file);
        fI = new double[arraySize];
        memcpy((void *) fI, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        //std::cout << "fI:" << std::endl; disp1d(fI,dimSize[0]);
        mxDestroyArray(tmp);

        nI = arraySize;
    
        readVar(ith,"i", pmat, file);
        readVar(dur,"dur", pmat, file);
        readVar(l0,"l0", pmat, file);
        //readVar(vThres,"vThres", pmat, file);
        //readVar(vReset,"vReset", pmat, file);
        readVar(ei,"ei", pmat, file);
    
        closeMat(pmat, file);
        std::cout << "Done" << std::endl;
        std::cout << std::endl;
        tstep = dur/(nt-1);
        std::cout << "nt = " << nt << std::endl;
        std::cout << "ndt = " << ndt << std::endl;
        std::cout << "nv = " << nv << std::endl;
        std::cout << "nE = " << nE << std::endl;
        std::cout << "nI = " << nI << std::endl;
        std::cout << "nt0 = " << nt0 << std::endl;
        std::cout << "tstep = " << tstep << std::endl;

        std::cout << "vRange:" << std::endl; disp1d(vRange,nv);
        std::cout << "dtRange:" << std::endl; disp1d(dtRange,ndt);
        std::cout << "fE:" << std::endl; disp1d(fE,nI);
        std::cout << "fI:" << std::endl; disp1d(fI,nI);
        idtRange = new size[ndt];
        for (i=0; i<ndt; i++){
            idtRange[i] = static_cast<size>(round(dtRange[i]/tstep));
            assert(idtRange[i]*tstep-dtRange[i] < 1e-12);
        }
    }
    // clear mem
    void clearLib() {
        size dimSize[4];
        dimSize[3] = nt;
        dimSize[2] = nE;
        dimSize[1] = ndt;
        dimSize[0] = nv;
        del4d(EPSP,dimSize);
        delete []EPSP_ptr;

        dimSize[2] = nI;
        del4d(IPSP,dimSize);
        delete []IPSP_ptr;

        dimSize[2] = nt0;
        dimSize[1] = nE;
        dimSize[0] = nv;
        del3d(EPSP0,dimSize);
        delete []EPSP0_ptr;

        dimSize[1] = nI;
        del3d(IPSP0,dimSize);
        delete []IPSP0_ptr;

        dimSize[2] = nt;
        dimSize[1] = ndt;
        del3d(kVEE,dimSize);
        delete []kVEE_ptr;
        del3d(kVEI,dimSize);
        delete []kVEI_ptr;
        del3d(kVIE,dimSize);
        delete []kVIE_ptr;
        del3d(kVII,dimSize);
        delete []kVII_ptr;
        del2d(vLeak);
        delete []vLeak_ptr;

        delete []vRange;
        delete []dtRange;
        delete []idtRange;
        delete []fE;
        delete []fI;
    }
} NeuroLib;
#endif
