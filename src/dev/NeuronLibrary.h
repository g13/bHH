#ifndef NEURONLIBRARY_H
#define NEURONLIBRARY_H
#include "matFunc.h"
#include <cmath>

typedef struct NeuronLibrary 
{
    size ndt,nt,nt0,nv,nE,nI,ith;

    double  ****EPSP,   ***EPSP0,  ****IPSP,   ***IPSP0,
           *EPSP_ptr, *EPSP0_ptr, *IPSP_ptr, *IPSP0_ptr,
             ****kEE,    ****kEI,   ****kIE,    ****kII,
           ******kV,
           *kEE_ptr,  *kEI_ptr, *kIE_ptr,  *kII_ptr,
           *kV_ptr,
             *vRange,   *dtRange, *fE, *fI, dur, tstep,
             **vLeak, *vLeak_ptr,
            ***E_tmax, ***I_tmax, ***tMax,
            *E_tmax_ptr, *I_tmax_ptr, *tMax_ptr;
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

        readArray(tmp,"E_tmax", dimSize, arraySize, pmat, file);
        E_tmax_ptr = new double[arraySize];
        memcpy((void *) E_tmax_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer3d(E_tmax,E_tmax_ptr,dimSize);
        mxDestroyArray(tmp);
    
        readArray(tmp,"I_tmax", dimSize, arraySize, pmat, file);
        I_tmax_ptr = new double[arraySize];
        memcpy((void *) I_tmax_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer3d(I_tmax,I_tmax_ptr,dimSize);
        mxDestroyArray(tmp);

        tMax_ptr = new double[(nE+nI)*ndt*nv];
        for (k=0; k<nv; k++) {
            for (j=0; j<ndt; j++) {
                int ii = j*(nE+nI)+k*ndt*(nE+nI);
                for (i=0; i<nE; i++) {
                    tMax_ptr[i+ii] = E_tmax[i+j*nE+k*ndt*nE];
                }
                for (i=0; i<nI; i++) {
                    tMax_ptr[i+nE+ii] = E_tmax[i+j*nI+k*ndt*nI];
                }
            }
        }
        pointer3d(tMax,tMax_ptr,dimSize);
    
        dimSize[2] = nE;
        dimSize[1] = ndt;
        dimSize[0] = nv;
        del3d(E_tmax,dimSize);
        delete []E_tmax_ptr;
        dimSize[2] = nI;
        del3d(I_tmax,dimSize);
        delete []I_tmax_ptr;

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
    
        readArray(tmp,"kV", dimSize, arraySize, pmat, file);
        kV_ptr = new double[arraySize];
        memcpy((void *) kV_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer6d(kV,kV_ptr,dimSize);
        mxDestroyArray(tmp);
    
        readArray(tmp,"kEE", dimSize, arraySize, pmat, file);
        kEE_ptr = new double[arraySize];
        memcpy((void *) kEE_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer4d(kEE,kEE_ptr,dimSize);
        mxDestroyArray(tmp);
    
        readArray(tmp,"kEI", dimSize, arraySize, pmat, file);
        kEI_ptr = new double[arraySize];
        memcpy((void *) kEI_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer4d(kEI,kEI_ptr,dimSize);
        mxDestroyArray(tmp);
    
        readArray(tmp,"kIE", dimSize, arraySize, pmat, file);
        kIE_ptr = new double[arraySize];
        memcpy((void *) kIE_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer4d(kIE,kIE_ptr,dimSize);
        mxDestroyArray(tmp);
    
        readArray(tmp,"kII", dimSize, arraySize, pmat, file);
        kII_ptr = new double[arraySize];
        memcpy((void *) kII_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer4d(kII,kII_ptr,dimSize);
        mxDestroyArray(tmp);
    
        readArray(tmp,"vleakage0", dimSize, arraySize, pmat, file);
        vLeak_ptr = new double[arraySize];
        memcpy((void *) vLeak_ptr, (void *)(mxGetPr(tmp)),arraySize*sizeof(double));
        pointer2d(vLeak,vLeak_ptr,dimSize);
        mxDestroyArray(tmp);

        readVar(ith,"i", pmat, file);
        readVar(dur,"dur", pmat, file);
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
        for (i=1; i<ndt; i++) {
            if (nt - idtRange[i] < idtRange[ndt-1] - idtRange[i-1]) {
                std::cout << "check the " << i << "th entry of dtRange" << std::endl;
                assert(nt - idtRange[i] >= idtRange[ndt-1] - idtRange[i-1]);
            }
        }
    }
    // clear mem
    void clearLib() {
        size dimSize[6];
        dimSize[5] = nt;
        dimSize[4] = ndt;
        dimSize[3] = nE+nI;
        dimSize[2] = nE+nI;
        dimSize[1] = ndt;
        dimSize[0] = nv;
        del6d(kV,dimSize);
        delete []kV_ptr;

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

        dimSize[3] = nt;
        dimSize[2] = ndt;
        dimSize[1] = ndt;
        del4d(kEE,dimSize);
        delete []kEE_ptr;
        del4d(kEI,dimSize);
        delete []kEI_ptr;
        del4d(kIE,dimSize);
        delete []kIE_ptr;
        del4d(kII,dimSize);
        delete []kII_ptr;
        del2d(vLeak);
        delete []vLeak_ptr;

        dimSize[2] = nE+nI;
        dimSize[1] = ndt;
        dimSize[0] = nv;
        del3d(tMax,dimSize);
        delete []tMax_ptr;

        delete []vRange;
        delete []dtRange;
        delete []idtRange;
        delete []fE;
        delete []fI;
    }
} NeuroLib;
#endif
