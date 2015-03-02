CC=gcc
ICC=icc
#CFLAGS=-I/panfs/storage.local/imb/home/sqin/fftw/fftw3.3/include -c -D DEBUG -Wall -pg
#CFLAGS=-I/panfs/storage.local/imb/home/sqin/fftw/fftw3.3/include -c -Wall -O3
#CFLAGS=-I/panfs/storage.local/imb/home/sqin/fftw/fftw3.3/include -c -Wall -pg -O3 -ftree-vectorizer-verbose=2
CFLAGS= -I${TACC_FFTW3_INC} -c -Wall -O3 -std=c99 -fopenmp #-pg -openmp_report=2
LDFLAGS= -I${TACC_FFTW3_INC} -Wl,-rpath,${TACC_FFTW3_LIB}  -L${TACC_FFTW3_LIB} -lfftw3 -lfftw3_omp -lm -Wall -fopenmp -O3 #-pg
SOURCES= cross.c  fftw.c  grdL.c  grdR.c  parm.c  pdb.c  pro.c  rdf.c  rot.c  score.c
SOURCESdd= ele.drt.c  vdw.drt.c  vol.drt.c linked.c drt.c drtlnk.c angijk.c zdread.c pdb.c rot.c score.c parm.c 
SOURCESddatm= ele.drt.c  vdw.drt.c  vol.drt.c linked.c drt.atm.c drtlnk.atm.c angijk.atm.c zdread.c pdb.c rot.c score.c parm.c
OBJECTS=$(SOURCES:.c=.o)
OBJECTSdd=$(SOURCESdd:.c=.o)
OBJECTSddatm=$(SOURCESddatm:.c=.o)
EXECUTABLE=fmap fmapdd fmapdd.atm fmapdd.showpqr

all: $(SOURCES) $(SOURCESdd) $(EXECUTABLE)
	
#$(EXECUTABLE): $(OBJECTS) $(EXECUTABLE).o
#	$(CC) $(OBJECTS) $(EXECUTABLE).o $(LDFLAGS) -o $@

fmap: $(OBJECTS) fmap.o
	$(CC) $(OBJECTS) fmap.o $(LDFLAGS) -o $@

fmapdd: $(OBJECTSdd)  fmapdd.c
	$(CC) $(OBJECTSdd) fmapdd.c $(LDFLAGS) -o $@

fmapdd.showpqr: $(OBJECTSdd)  fmapdd.showpqr.c
	$(CC) $(OBJECTSdd) fmapdd.showpqr.c $(LDFLAGS) -o $@

fmapdd.atm: $(OBJECTSddatm)  fmapdd.c
	$(CC) $(OBJECTSddatm) fmapdd.c $(LDFLAGS) -o $@
.c.o:
	$(ICC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o
