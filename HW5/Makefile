#std-c++11
umeCXX = g++
CXXFLAGS = -Wall -O2 -Wextra -Wno-unused-local-typedefs  -Werror -Wno-deprecated-declarations -std=c++11
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

ROOT = `root-config --cflags --glibs`

MKDIR_BIN = mkdir -p $(PWD)/bin
MKDIR_PDFDIR = mkdir -p $(PWD)/pdfDir
MKDIR_OUTPUT = mkdir -p $(PWD)/output

all: mkdirBin mkdirOutput mkdirPdfdir stochProteins stochBistable stochPotential stochBistableGil stochPotentialProb

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdfdir:
	$(MKDIR_PDFDIR)

stochProteins: src/stochProteins.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/stochProteins.exe src/stochProteins.C

stochBistable: src/stochBistable.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/stochBistable.exe src/stochBistable.C

stochBistableGil: src/stochBistableGil.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/stochBistableGil.exe src/stochBistableGil.C

stochPotential: src/stochPotential.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/stochPotential.exe src/stochPotential.C

stochPotentialProb: src/stochPotentialProb.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/stochPotentialProb.exe src/stochPotentialProb.C

clean:
	rm *~ || true
	rm *# || true
	rm include/*~ || true
	rm include/#*# || true
	rm src/*~ || true
	rm src/#*# || true
	rm bin/*.exe || true
	rmdir bin || true