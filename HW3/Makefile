#std-c++11
CXX = g++
CXXFLAGS = -Wall -O2 -Wextra -Wno-unused-local-typedefs  -Werror -Wno-deprecated-declarations -std=c++11
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

ROOT = `root-config --cflags --glibs`

MKDIR_BIN = mkdir -p $(PWD)/bin

all: mkdirBin emulateCircadian nodeEdgeModel

mkdirBin:
	$(MKDIR_BIN)

emulateCircadian: src/emulateCircadian.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/emulateCircadian.exe src/emulateCircadian.C

nodeEdgeModel: src/nodeEdgeModel.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(PWD) -o bin/nodeEdgeModel.exe src/nodeEdgeModel.C

clean:
	rm *~ || true
	rm *# || true
	rm include/*~ || true
	rm include/#*# || true
	rm src/*~ || true
	rm src/#*# || true
	rm bin/*.exe || true
	rmdir bin || true