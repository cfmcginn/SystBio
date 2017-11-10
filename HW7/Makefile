#std-c++11
umeCXX = g++
CXXFLAGS = -Wall -O2 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

ROOT = `root-config --cflags --glibs`
BOOST='/Users/cfmcginn/Packages/Boost/boost_1_65_1/'

MKDIR_BIN = mkdir -p $(PWD)/bin
MKDIR_PDFDIR = mkdir -p $(PWD)/pdfDir
MKDIR_GIFDIR1 = mkdir -p $(PWD)/gifDir1
MKDIR_GIFDIR2 = mkdir -p $(PWD)/gifDir2
MKDIR_OUTPUT = mkdir -p $(PWD)/output

all: mkdirBin mkdirOutput mkdirPdfdir mkdirGifdir1 mkdirGifdir2 odeInt nullClines campSim campSim2

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdfdir:
	$(MKDIR_PDFDIR)

mkdirGifdir1:
	$(MKDIR_GIFDIR1)

mkdirGifdir2:
	$(MKDIR_GIFDIR2)

odeInt: src/odeInt.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -I $(BOOST) -o bin/odeInt.exe src/odeInt.C

nullClines: src/nullClines.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -I $(BOOST) -o bin/nullClines.exe src/nullClines.C

campSim: src/campSim.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -I $(BOOST) -o bin/campSim.exe src/campSim.C

campSim2: src/campSim2.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -I $(BOOST) -o bin/campSim2.exe src/campSim2.C

clean:
	rm *~ || true
	rm *# || true
	rm include/*~ || true
	rm include/#*# || true
	rm src/*~ || true
	rm src/#*# || true
	rm bin/*.exe || true
	rmdir bin || true