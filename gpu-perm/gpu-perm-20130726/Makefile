####################################################
#                                                  #
#  GPU-PerM                                        #
#  Haifeng Chen (haifengc at usc dot edu)          #
#  University of Southern California               #
#  July 1, 2013                                    #
#                                                  #
####################################################

CUDA_INSTALL_PATH := /usr/usc/cuda/5.0
NVCC := $(CUDA_INSTALL_PATH)/bin/nvcc

CXXFLAGS = -O3 --compile  -x c++
NVCCFLAGS := --device-c -O3 -gencode arch=compute_35,code=sm_35  -x cu -rdc=true -g -G
LINK := --relocatable-device-code=true -gencode arch=compute_35,code=sm_35 -link

TARGET := gpu-perm

CPP_SRC := hash.cpp iofile.cpp option.cpp refin.cpp
CU_SRC := allocate.cu bitscode.cu match_kernel.cu seed.cu gpu-perm.cu
CPP_OBJS := $(patsubst %.cpp, %.o, $(CPP_SRC))
CU_OBJS := $(patsubst %.cu, %.o, $(CU_SRC))

%.o: %.cpp
	$(NVCC) -c $(CXXFLAGS) -o $@ $<
%.o: %.cu
	$(NVCC) -c $(NVCCFLAGS) -o $@ $<

$(TARGET): $(CPP_OBJS) $(CU_OBJS)
	$(NVCC) $(LINK) -o $(TARGET) $(CU_OBJS) $(CPP_OBJS)
	
clean:
	rm -rf $(TARGET) *.o
