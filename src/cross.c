#include "cross.h"

#ifdef _OPENMP
extern fftw_real *recVol, *recR12, *recR6, *recEle;
#endif

void xVol(PRO *rec, PRO *lig, PARM sys, const int l,fftw_real ligGrd[l][l][2*(l/2+1)]){
#ifdef _OPENMP
        static bool first=true;
        fftw_real (*recGrd)[l][2*(l/2+1)]=(void *)recVol;
        if (first){
#else
	fftw_real recGrd[l][l][2*(l/2+1)];
#endif
	SetZero(l,recGrd);
	RGrdVol(l, &(recGrd[0][0][0]), rec, sys,RdfVol);
	fft3d_r2c(l,l,l,recGrd);
#ifdef _OPENMP
		first=false;
        }
#endif
        SetZero(l,ligGrd);
	RGrdVol(l, &(ligGrd[0][0][0]), lig, sys,RdfVol);
	fft3d_r2c(l,l,l,ligGrd);
        fft3d_add_inplace(l,l,l,recGrd,ligGrd);
        fft3d_c2r(l,l,l,ligGrd);
}

void xR12(PRO *rec, PRO *lig, PARM sys, const int l,fftw_real ligGrd[l][l][2*(l/2+1)]){
#ifdef _OPENMP
        static bool first=true;
        fftw_real (*recGrd)[l][2*(l/2+1)]=(void *)recR12;
        if (first){
#else
	fftw_real recGrd[l][l][2*(l/2+1)];
#endif
	SetZero(l,recGrd);
        RGrdR12(l, &(recGrd[0][0][0]), rec, sys,RdfR12m);
        fft3d_r2c(l,l,l,recGrd);
#ifdef _OPENMP
		first=false;
        }
#endif
        SetZero(l,ligGrd);
        LGrdR12(l, &(ligGrd[0][0][0]), lig);
        fft3d_r2c(l,l,l,ligGrd);
        fft3d_add_inplace(l,l,l,recGrd,ligGrd);
        fft3d_c2r(l,l,l,ligGrd);
}

void xR6(PRO *rec, PRO *lig, PARM sys, const int l,fftw_real ligGrd[l][l][2*(l/2+1)]){
#ifdef _OPENMP
	static bool first=true;
	fftw_real (*recGrd)[l][2*(l/2+1)]=(void *)recR6;
	if (first){
#else
	fftw_real recGrd[l][l][2*(l/2+1)];
#endif
        SetZero(l,recGrd);
        RGrdR6(l, &(recGrd[0][0][0]), rec, sys,RdfR6m);
        fft3d_r2c(l,l,l,recGrd);
#ifdef _OPENMP
		first=false;
	}
#endif
        SetZero(l,ligGrd);
        LGrdR6(l, &(ligGrd[0][0][0]), lig);
        fft3d_r2c(l,l,l,ligGrd);
        fft3d_add_inplace(l,l,l,recGrd,ligGrd);
        fft3d_c2r(l,l,l,ligGrd);
}

void xEle(PRO *rec, PRO *lig, PARM sys, const int l,fftw_real ligGrd[l][l][2*(l/2+1)]){
#ifdef _OPENMP
        static bool first=true;
        fftw_real (*recGrd)[l][2*(l/2+1)]=(void *)recEle;
        if (first){
#else
        fftw_real recGrd[l][l][2*(l/2+1)];
#endif
        SetZero(l,recGrd);
        RGrdEle(l, &(recGrd[0][0][0]), rec, sys,RdfDebye);
        fft3d_r2c(l,l,l,recGrd);
#ifdef _OPENMP
                first=false;
        }
#endif
        SetZero(l,ligGrd);
        LGrd2ndEle(l, &(ligGrd[0][0][0]), lig);
        fft3d_r2c(l,l,l,ligGrd);
        fft3d_add_inplace(l,l,l,recGrd,ligGrd);
        fft3d_c2r(l,l,l,ligGrd);
}

void Cross(PRO *rec, PRO *lig, PARM sys, const int l,fftw_real sav[l][l][2*(l/2+1)], const double ang[3]){
	bool vol[l][l][l];
	static double erncut=10.0;
	fftw_real local[l][l][2*(l/2+1)];
	fftw_real ligGrd[l][l][2*(l/2+1)];
	SetZero(l,local);
	xVol(rec,lig,sys,l,ligGrd);
	vRep(l,vol,0.001,ligGrd,1.0,"vol",sys.kBT);
	//exit(0);
	xEle(rec,lig,sys,l,ligGrd);
	AddTo1(l,local,1.0,ligGrd,sys.escl);
	vRep(l,vol,0.001,ligGrd,sys.escl,"ele",sys.kBT);
	
	xR6(rec,lig,sys,l,ligGrd);
        AddTo1(l,local,1.0,ligGrd,sys.vscl);
        vRep(l,vol,0.001,ligGrd,sys.vscl,"vdw",sys.kBT);

	xR12(rec,lig,sys,l,ligGrd);
	AddTo1(l,local,1.0,ligGrd,sys.vscl);
	vRep(l,vol,0.001,ligGrd,sys.vscl,"v12",sys.kBT);
	vRep(l,vol,0.001,local,1.0,"v+e",sys.kBT);
	if (erncut>0){
		if (sys.erncut>0){
			erncut=GetSpotCut(l,vol,local,100,100,sys.kBT);
		}else{
			erncut=sys.erncut;
		}
		fprintf(stderr,"erncut:%f\n",erncut);
	}
	SelErn(l,vol,local,ang,erncut);
	BltzSumFilt(l,vol,sav,local,sys.kBT);
}
