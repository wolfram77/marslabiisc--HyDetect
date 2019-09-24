#!/bin/bash
#GCC Compilers:
CC  = nvcc 
CPP = nvcc 
#CFLAGS   = -Ofast -fopenmp -std=c99 
CFLAGS=-Xcompiler "-Wall -Wconversion -Wsign-conversion -Wextra -Wshadow -fopenmp" $(NVCC_ARCH) -Xptxas="-v"
#CPPFLAGS = -Ofast -fopenmp -std=c++0x
CPPFLAGS=-Xcompiler "-Wall -Wconversion -Wsign-conversion -Wextra -Wshadow -fopenmp" $(NVCC_ARCH) -Xptxas="-v"

##############

CUDAVERSION=8.0
CFLAGS1= -I/usr/local/cuda-$(CUDAVERSION)/include -O3

CC1=/usr/local/cuda-$(CUDAVERSION)/bin/nvcc
DFLAGS= -D RUNONGPU
CUDAFLAGS= -arch sm_35

DEPS = communityGPU.h  graphGPU.h  graphHOST.h openaddressing.h

OBJ = binWiseGaussSeidel.o communityGPU.o preprocessing.o  aggregateCommunity.o coreutility.o independentKernels.o gatherInformation.o graphHOST.o graphGPU.o main.o assignGraph.o computeModularity.o computeTime.o


LIBS1= -L/usr/local/cuda-$(CUDAVERSION)/lib64 -lcudart
EXEC=run_CU_community

################

METIS_HOME = /afs/msrc.pnl.gov/files/home/hala533/metis-5.0.2
METIS_INCLUDE = -I$(METIS_HOME)/include
METIS_LIB = -L$(METIS_HOME)/libmetis -lmetis -lm

LDFLAGS  = $(CPPFLAGS)
INCLUDES = . $(METIS_INCLUDE)
LIBS     = -lm


TARGET_1 = driverForGraphClustering
TARGET_2 = convertFileToBinary
TARGET_3 = convertFileToEdgeList
TARGET_4 = driverForColoringExperiments

TARGET_5 = driverForRmat
TARGET_6 = driverForRGG
TARGET_7 = driverForPartitioningWithMetis

TARGET_8 = convertSnapFileToBinary

TARGET = $(TARGET_1) $(EXEC)

OBJECTS  =graphReaderFromGr.o\
 binWiseGaussSeidel.o communityGPU.o preprocessing.o  aggregateCommunity.o\
gpuonly.o\
coreutility.o independentKernels.o gatherInformation.o graphHOST.o graphGPU.o  assignGraph.o computeModularity.o\
computeTime.o RngStream.o utilityFunctions.o  \
writeGraphDimacsFormat.o buildNextPhase.o \
coloringDistanceOne.o utilityClusteringFunctions.o \
parallelLouvainMethod.o parallelLouvainWithColoring.o \
louvainMultiPhaseRun.o parseInputParameters.o vertexFollowing.o cpuonly.o\
verticesToMoveToGPU.o verticesToMoveToCPU.o modifyCPUstructure.o modifyGPUstructure.o \
	

all: $(TARGET_1) message

$(TARGET_1): $(OBJECTS) $(TARGET_1).o
	 $(CPP) $(LDFLAGS)   -o  $(TARGET_1) $(TARGET_1).o $(OBJECTS) $(LIBS1)  
.cu.o: $(DEPS)
	$(CC1) -o $@ -c $< $(CFLAGS1) $(DFLAGS) $(CUDAFLAGS)

#%.o: %.cpp $(DEPS)
#	$(CC1) -o $@ -c $< $(CFLAGS1)

#############
#$(EXEC): $(OBJ)
#	$(CC1) -o $@ $^ $(LIBS1) 
###########
.c.o:
	$(CC) $(CFLAGS)  -c $< -I$(INCLUDES) -o $@

.cpp.o:
	$(DEPS)$(CPP) $(CPPFLAGS) -c $< -I$(INCLUDES) -o $@

clean:
	rm -f $(TARGET1).o $(OBJECTS)

wipe:
	rm -f $(TARGET1).o $(OBJECTS) $(TARGET1) *~ *.bak

message:
	echo "Executables: " $(TARGET1) " have been created"

#CUDAVERSION=8.0
#CFLAGS1= -I/usr/local/cuda-$(CUDAVERSION)/include -O3 

#CC1=/usr/local/cuda-$(CUDAVERSION)/bin/nvcc 
#DFLAGS= -D RUNONGPU
#CUDAFLAGS= -arch sm_35 

#DEPS = communityGPU.h  graphGPU.h  graphHOST.h openaddressing.h

#OBJ = binWiseGaussSeidel.o communityGPU.o preprocessing.o  aggregateCommunity.o coreutility.o independentKernels.o gatherInformation.o graphHOST.o graphGPU.o main.o assignGraph.o computeModularity.o computeTime.o


#LIBS1= -L/usr/local/cuda-$(CUDAVERSION)/lib64 -lcudart 


#EXEC=run_CU_community
all:$(EXEC)

#$(EXEC): $(OBJ)
#	$(CC1) -o $@ $^ $(LIBS1) 

%.o: %.cu $(DEPS)
	$(CC1) -o $@ -c $< $(CFLAGS1) $(DFLAGS) $(CUDAFLAGS)

%.o: %.cpp $(DEPS)
	$(CC1) -o $@ -c $< $(CFLAGS1) 


clean:
	rm -f *.o *~ $(EXEC)


