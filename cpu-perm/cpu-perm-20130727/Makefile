####################################################
#                                                  #
#  GPU-PerM                                        #
#  Haifeng Chen (haifengc at usc dot edu)          #
#  University of Southern California               #
#  July 1, 2013                                    #
#                                                  #
####################################################

CXXFLAGS := -O3 -Wall -fmessage-length=0

TARGET := cpu-perm

SRCS := hash.cpp iofile.cpp option.cpp refin.cpp allocate.cpp bitscode.cpp match_kernel.cpp seed.cpp cpu-perm.cpp

OBJS := $(SRCS:.cpp=.o)

all: $(TARGET)

$(OBJS): $(SRCS)
	$(CXX) $(CXXFLAGS) $(PROFILE) -c $(SRCS)
                
$(TARGET): $(OBJS)
	$(CXX) $(PROFILE) -o $(TARGET) $(OBJS)

clean:
	rm -rf $(TARGET) *.o *.exe
                                        
