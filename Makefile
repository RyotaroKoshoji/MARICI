CXX = mpiicpc
CXXFLAGS = -std=c++20 -qmkl=sequential -parallel -O3 -lstdc++fs -ip -ipo -no-prec-div -fp-model=fast=2 -march=core-avx2 -diag-disable=10441
INCLUDEDIR = ./include $(HOME)/.Library/spglib/include
SRCDIR := ./src
BUILDDIR := ./build
OBJDIR := $(BUILDDIR)/obj
TARGET := $(BUILDDIR)/Marici.exe


INCLUDE := $(addprefix -iquote, $(INCLUDEDIR))
SRC := $(shell find $(SRCDIR) -name *.cpp)
OBJ := $(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

MKDIR = mkdir -p
RM = rm -rf


$(TARGET) : $(OBJ)
	$(CXX) -o $@ $^  -std=c++20 -qmkl=sequential -parallel -O3 -lstdc++fs -ip -ipo -no-prec-div -fp-model=fast=2 -march=core-avx2 -diag-disable=10441 libsymspg.so

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(MKDIR) $(OBJDIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(INCLUDE)

.PHONY : clean
clean :
	$(RM) $(OBJDIR)
	$(RM) $(TARGET)

.PHONY : rebuild
rebuild :
	make clean
	make

