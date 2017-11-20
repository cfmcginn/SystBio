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
MKDIR_OUTPUT = mkdir -p $(PWD)/output

all: mkdirBin mkdirOutput mkdirPdfdir logMap logMap2 allee

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdfdir:
	$(MKDIR_PDFDIR)

logMap: src/logMap.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -I $(BOOST) -o bin/logMap.exe src/logMap.C

logMap2: src/logMap2.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -I $(BOOST) -o bin/logMap2.exe src/logMap2.C

allee: src/allee.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -I $(BOOST) -o bin/allee.exe src/allee.C

clean:
	rm *~ || true
	rm *# || true
	rm include/*~ || true
	rm include/#*# || true
	rm src/*~ || true
	rm src/#*# || true
	rm bin/*.exe || true
	rmdir bin || true