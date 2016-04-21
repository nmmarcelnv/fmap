ICC=icc
MPICC=mpicc
#CFLAGS=-I/panfs/storage.local/imb/home/sqin/fftw/fftw3.3/include -c -D DEBUG -Wall -pg
#CFLAGS=-I/panfs/storage.local/imb/home/sqin/fftw/fftw3.3/include -c -Wall -O3
#CFLAGS=-I/panfs/storage.local/imb/home/sqin/fftw/fftw3.3/include -c -Wall -pg -O3 -ftree-vectorizer-verbose=2
CFLAGS= -I${TACC_FFTW3_INC} -c -Wall -O3 -std=c99 -openmp -DUSE_MPI #-pg -openmp_report=2
CFLAGSMIC= -I${TACC_FFTW3_INC} -c -Wall -O3 -std=c99 -openmp -mmic #-pg -openmp_report=2
LDFLAGS= -I${TACC_FFTW3_INC} -Wl,-rpath,${TACC_FFTW3_LIB}  -L${TACC_FFTW3_LIB} -lfftw3 -lfftw3_omp -lm -Wall -std=c99 -openmp -O3 -DUSE_MPI #-pg
LDFLAGSMIC= -I${TACC_FFTW3_INC} -lm -Wall -std=c99 -openmp -O3 -mmic #-pg
SOURCES= cross.c  fftw.c  grdL.c  grdR.c  parm.c  pdb.c  pro.c  rdf.c  rot.c  score.c cputime.c
SOURCESdd= ele.drt.c  vdw.drt.c  vol.drt.c linked.c drt.c drtlnk.c angijk.c zdread.c pdb.c rot.c score.c parm.c cputime.c 
SOURCESddatm= ele.drt.c  vdw.drt.c  vol.drt.c linked.c drt.atm.c drtlnk.atm.c angijk.atm.c zdread.c pdb.c rot.c score.c parm.c cputime.c
OBJECTS=$(SOURCES:.c=.o)
OBJECTSdd=$(SOURCESdd:.c=.o)
OBJECTSddMIC=$(SOURCESdd:.c=.om)
OBJECTSddatm=$(SOURCESddatm:.c=.o)
OBJECTSddatmMIC=$(SOURCESddatm:.c=.om)
EXECUTABLE=fmap fmap.cen fmapdd fmapdd.atm fmapdd.showpqr matevscl fmappdbs 
EXECUTABLEMIC=fmapdd.mic fmapdd.atm.mic fmapdd.showpqr.mic

all: $(SOURCES) $(SOURCESdd) $(EXECUTABLE) $(EXECUTABLEMIC)
	
#$(EXECUTABLE): $(OBJECTS) $(EXECUTABLE).o
#	$(CC) $(OBJECTS) $(EXECUTABLE).o $(LDFLAGS) -o $@

fmap: $(OBJECTS) fmap.o
	$(MPICC) $(OBJECTS) fmap.o $(LDFLAGS) -o $@
fmap.cen: $(OBJECTS) fmap.c
	grep -v -e ToCtd -e CalCtd fmap.c >fmap.cen.c
	$(MPICC) $(OBJECTS) fmap.cen.c $(LDFLAGS) -o $@

fmappdbs:  $(OBJECTS) fmappdbs.c
	$(MPICC) $(OBJECTS) fmappdbs.c $(LDFLAGS) -o $@

fmapdd: $(OBJECTSdd)  fmapdd.c
	$(MPICC) $(OBJECTSdd) fmapdd.c $(LDFLAGS) -o $@
fmapdd.mic: $(OBJECTSddMIC)  fmapdd.c
	$(ICC) $(OBJECTSddMIC) fmapdd.c $(LDFLAGSMIC) -o $@

fmapdd.showpqr: $(OBJECTSdd)  fmapdd.showpqr.c
	$(MPICC) $(OBJECTSdd) fmapdd.showpqr.c $(LDFLAGS) -o $@
fmapdd.showpqr.mic: $(OBJECTSddMIC)  fmapdd.showpqr.c
	$(ICC) $(OBJECTSddMIC) fmapdd.showpqr.c $(LDFLAGSMIC) -o $@

fmapdd.atm: $(OBJECTSddatm)  fmapdd.c
	$(MPICC) $(OBJECTSddatm) fmapdd.c $(LDFLAGS) -o $@
fmapdd.atm.mic: $(OBJECTSddatmMIC)  fmapdd.c
	$(ICC) $(OBJECTSddatmMIC) fmapdd.c $(LDFLAGSMIC) -o $@
matevscl: matevscl.c $(OBJECTS)
	$(MPICC) $(CFLAGS) $< -o matevscl.o
	$(MPICC) pdb.o score.o parm.o matevscl.o $(LDFLAGS) -o $@

.c.o:
	$(MPICC) $(CFLAGS) $< -o $@
	$(ICC) $(CFLAGSMIC) $< -o $(@)m

clean:
	rm -f *.o *.om
