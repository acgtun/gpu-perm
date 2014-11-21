MPIPATH = /usr/usc/openmpi/default/bin/
MPILIB = /usr/usc/openmpi/default/lib/

#MPIPATH = /usr/lib64/openmpi/bin/
#MPILIB = /usr/lib64/openmpi/lib/

CC = $(MPIPATH)mpicc
CXX = $(MPIPATH)mpicxx

CXXFLAGS := -O3 -Wall -fmessage-length=0

SRCS := hash.cpp iofile.cpp option.cpp match.cpp

OBJS := $(SRCS:.cpp=.o)

default: buildindex saligner paligner 

.cpp.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<

buildindex: index.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^
	
paligner: pmatch.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^
	
saligner: smatch.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^
                

all: buildindex saligner paligner 

clean:
	rm -rf aligner saligner paligner buildindex *.o *.exe
                                        
