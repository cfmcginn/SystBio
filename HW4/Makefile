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

all: mkdirBin mkdirOutput mkdirPdfdir geneDuplication

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdfdir:
	$(MKDIR_PDFDIR)

geneDuplication: src/geneDuplication.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/geneDuplication.exe src/geneDuplication.C

clean:
	rm *~ || true
	rm *# || true
	rm include/*~ || true
	rm include/#*# || true
	rm src/*~ || true
	rm src/#*# || true
	rm bin/*.exe || true
	rmdir bin || true