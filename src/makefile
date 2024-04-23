CXX      ?= g++
CXXFLAGS ?= -std=c++20
CPPFLAGS ?= -O3 -Wall -I"../include"
LINK.o := $(LINK.cc) # implicit flag to enable the linking

# Set the default Eigen include path
# EIGEN_INC ?= ${mkEigenInc}
# 
# # Override the default Eigen include path based on a flag
# ifdef EIGEN_MAC
# 	EIGEN_INC := /opt/homebrew/include/eigen3
# endif
# 
# CPPFLAGS += -I$(EIGEN_INC)


SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
HEADERS = $(wildcard *.hpp)

exe_sources = $(filter main%.cpp, $(SRCS))
EXEC = $(exe_sources:.cpp=)

.PHONY: all debug clean distclean

all: $(EXEC)

debug: CPPFLAGS += -DDEBUG
debug: all

$(EXEC): $(OBJS)

$(OBJS): $(SRCS) $(HEADERS)

clean:
	-\rm -f $(EXEC) $(OBJS)

distclean: clean
	$(RM) -f ./doc $(DEPEND)
	$(RM) *.out *.bak *~
