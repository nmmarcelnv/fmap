#include "cross.h"

#ifdef _OPENMP
extern fftw_real *recVol, *recR12, *recR6, *recEle;
#endif

void xVol(PRO *rec, PRO *lig, PARM sys, const int l,fftw_real ligGrd[l][l][2*(l/2+1)]){
#ifdef _OPENMP
        static bool first=true;
        fftw_real (*recGrd)[l][2*(l/2+1)]=(void *)recVol;
        if (first){
		Times(true,1,1);
#else
	fftw_real recGrd[l][l][2*(l/2+1)];
#endif
	SetZero(l,recGrd);
	Times(true,1,2);
	RGrdVol(l, &(recGrd[0][0][0]), rec, sys,RdfVol);
	Times(false,1,2);
	fft3d_r2c(l,l,l,recGrd);
#ifdef _OPENMP
		Times(false,1,1);
		first=false;
        }
#endif
        SetZero(l,ligGrd);
	Times(true,1,3);
	RGrdVol(l, &(ligGrd[0][0][0]), lig, sys,RdfVol);
	Times(false,1,3);
	fft3d_r2c(l,l,l,ligGrd);
        fft3d_add_inplace(l,l,l,recGrd,ligGrd);
        fft3d_c2r(l,l,l,ligGrd);
}

void xR12(PRO *rec, PRO *lig, PARM sys, const int l,fftw_real ligGrd[l][l][2*(l/2+1)]){
#ifdef _OPENMP
        static bool first=true;
        fftw_real (*recGrd)[l][2*(l/2+1)]=(void *)recR12;
        if (first){
		Times(true,1,1);
#else
	fftw_real recGrd[l][l][2*(l/2+1)];
#endif
	SetZero(l,recGrd);
	Times(true,1,2);
        RGrdR12(l, &(recGrd[0][0][0]), rec, sys,RdfR12m);
	Times(false,1,2);
        fft3d_r2c(l,l,l,recGrd);
#ifdef _OPENMP
		Times(false,1,1);
		first=false;
        }
#endif
        SetZero(l,ligGrd);
	Times(true,1,3);
        LGrdR12(l, &(ligGrd[0][0][0]), lig);
	Times(false,1,3);
        fft3d_r2c(l,l,l,ligGrd);
        fft3d_add_inplace(l,l,l,recGrd,ligGrd);
        fft3d_c2r(l,l,l,ligGrd);
}

void xR6(PRO *rec, PRO *lig, PARM sys, const int l,fftw_real ligGrd[l][l][2*(l/2+1)]){
#ifdef _OPENMP
	static bool first=true;
	fftw_real (*recGrd)[l][2*(l/2+1)]=(void *)recR6;
	if (first){
		Times(true,1,1);
#else
	fftw_real recGrd[l][l][2*(l/2+1)];
#endif
        SetZero(l,recGrd);
	Times(true,1,2);
        RGrdR6(l, &(recGrd[0][0][0]), rec, sys,RdfR6m);
	Times(false,1,2);
        fft3d_r2c(l,l,l,recGrd);
#ifdef _OPENMP
		Times(false,1,1);
		first=false;
	}
#endif
        SetZero(l,ligGrd);
	Times(true,1,3);
        LGrdR6(l, &(ligGrd[0][0][0]), lig);
	Times(false,1,3);
        fft3d_r2c(l,l,l,ligGrd);
        fft3d_add_inplace(l,l,l,recGrd,ligGrd);
        fft3d_c2r(l,l,l,ligGrd);
}

void xEle(PRO *rec, PRO *lig, PARM sys, const int l,fftw_real ligGrd[l][l][2*(l/2+1)]){
#ifdef _OPENMP
        static bool first=true;
        fftw_real (*recGrd)[l][2*(l/2+1)]=(void *)recEle;
        if (first){
		Times(true,1,1);
#else
        fftw_real recGrd[l][l][2*(l/2+1)];
#endif
        SetZero(l,recGrd);
	Times(true,1,2);
        RGrdEle(l, &(recGrd[0][0][0]), rec, sys,RdfDebye);
	Times(false,1,2);
        fft3d_r2c(l,l,l,recGrd);
#ifdef _OPENMP
		Times(false,1,1);
                first=false;
        }
#endif
        SetZero(l,ligGrd);
	Times(true,1,3);
        LGrd2ndEle(l, &(ligGrd[0][0][0]), lig);
	Times(false,1,3);
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
			fprintf(stderr,"erncut:%f\n",erncut);
		}else{
			erncut=sys.erncut;
		}
	}
	SelErn(l,vol,local,ang,erncut);
	BltzSumFilt(l,vol,sav,local,sys.kBT);
}
