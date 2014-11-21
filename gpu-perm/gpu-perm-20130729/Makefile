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
NVCCFLAGS := --device-c -O3 -gencode arch=compute_30,code=sm_30  -x cu -rdc=true -lineinfo
LINK := --relocatable-device-code=true -gencode arch=compute_30,code=sm_30 -link

TARGET := gpu-perm

CPP_SRC := hash.cpp iofile.cpp option.cpp refin.cpp
CU_SRC := alloc_kernel.cu bitscode.cu match.cu seed.cu gpu-perm.cu
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
