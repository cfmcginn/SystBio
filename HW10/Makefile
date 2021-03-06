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
MKDIR_GIFDIR = mkdir -p $(PWD)/gifDir
MKDIR_OUTPUT = mkdir -p $(PWD)/output


all: mkdirBin mkdirOutput mkdirPdfdir mkdirGifdir moranModel  betHedging plotMoranModel moranModel_Path_NOROOT moranModel_Path_CSVPlot

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdfdir:
	$(MKDIR_PDFDIR)

mkdirGifdir:
	$(MKDIR_GIFDIR)

moranModel: src/moranModel.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -I $(BOOST) -o bin/moranModel.exe src/moranModel.C

moranModel_Path: src/moranModel_Path.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -I $(BOOST) -o bin/moranModel_Path.exe src/moranModel_Path.C

moranModel_Path_NOROOT: src/moranModel_Path_NOROOT.C
	$(CXX) $(CXXFLAGS) -I $(PWD) -I $(BOOST) -o bin/moranModel_Path_NOROOT.exe src/moranModel_Path_NOROOT.C

moranModel_Path_CSVPlot: src/moranModel_Path_CSVPlot.C
	$(CXX) $(CXXFLAGS) -I $(PWD) -I $(BOOST)  $(ROOT) -o bin/moranModel_Path_CSVPlot.exe src/moranModel_Path_CSVPlot.C

betHedging: src/betHedging.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -I $(BOOST) -o bin/betHedging.exe src/betHedging.C

plotMoranModel: src/plotMoranModel.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -I $(BOOST) -o bin/plotMoranModel.exe src/plotMoranModel.C


clean:
	rm *~ || true
	rm *# || true
	rm include/*~ || true
	rm include/#*# || true
	rm src/*~ || true
	rm src/#*# || true
	rm bin/*.exe || true
	rmdir bin || true