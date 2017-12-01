#std-c++11
CXX = g++
CXXFLAGS = -Wall -O2 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

ROOT = `root-config --cflags --glibs`
BOOST='/Users/cfmcginn/Packages/Boost/boost_1_65_1/'

MKDIR_BIN = mkdir -p $(PWD)/bin
MKDIR_PDFDIR = mkdir -p $(PWD)/pdfDir
MKDIR_OUTPUT = mkdir -p $(PWD)/output

all: mkdirBin mkdirOutput mkdirPdfdir hawkDoveSim simRNG popTakeOver stochIBM quickAmplitude

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdfdir:
	$(MKDIR_PDFDIR)

hawkDoveSim: src/hawkDoveSim.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -I $(BOOST) -o bin/hawkDoveSim.exe src/hawkDoveSim.C

simRNG: src/simRNG.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -I $(BOOST) -o bin/simRNG.exe src/simRNG.C

popTakeOver: src/popTakeOver.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -I $(BOOST) -o bin/popTakeOver.exe src/popTakeOver.C

stochIBM: src/stochIBM.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -I $(BOOST) -o bin/stochIBM.exe src/stochIBM.C

quickAmplitude: src/quickAmplitude.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -I $(BOOST) -o bin/quickAmplitude.exe src/quickAmplitude.C

clean:
	rm *~ || true
	rm *# || true
	rm include/*~ || true
	rm include/#*# || true
	rm src/*~ || true
	rm src/#*# || true
	rm bin/*.exe || true
	rmdir bin || true