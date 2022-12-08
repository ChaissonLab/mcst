BLASR=blasr/common
LINKS=-lz -lpthread
BINARIES=printGaps

all: $(BINARIES)

screenInversions: ScreenInversions.cpp InversionAlign.h
	g++ -O3 $< -o $@ -I $(BLASR) -lpthread

printGaps: PrintGaps.cpp 
	g++ -O2  -I $(CONDA_PREFIX)/include $< -o $@ -L$(CONDA_PREFIX)/lib -lhts  -Wl,-rpath,$(CONDA_PREFIX)/lib
clean:
	rm -f $(BINARIES)
