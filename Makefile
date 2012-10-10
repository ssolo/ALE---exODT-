
CC=g++
mCC=mpic++

FLAGS = -O3  -fmerge-all-constants -funroll-loops -DNDEBUG -Wall 
DEV_FLAGS =  -g -Wall
boost_DIR=/usr/local/
bpp_DIR= /usr/local/
boost_libs= -lboost_mpi -lboost_serialization 


MPI_INCLUDE=-I/usr/lib/openmpi/include/

ifndef OSTYPE
  OSTYPE = $(shell uname -s|awk '{print tolower($$0)}')
  export OSTYPE
endif

ifeq ($(OSTYPE),darwin)
	bpp_libs = $(bpp_DIR)lib/libbpp-core.a $(bpp_DIR)lib/libbpp-phyl.a $(bpp_DIR)lib/libbpp-seq.a 
	boost_libs = $(boost_DIR)lib/libboost_serialization.a $(boost_DIR)lib/libboost_wserialization.a $(boost_DIR)lib/libboost_mpi.a
else
	CC=g++-4.4
	bpp_libs = -lbpp-core -lbpp-phyl -lbpp-seq
	boost_DIR= $(HOME)/
	STATIC= libexODT.a	
endif

DYNAMIC =  -L. -L$(bpp_DIR)lib $(bpp_libs) $(boost_libs) 
LINK = $(DYNAMIC) -lexODT
INCLUDE = -I$(boost_DIR)include -I$(bpp_DIR)include  


ALE.o: ALE.h ALE.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE) -c -o ALE.o ALE.cpp
ALE_util.o: ALE.h ALE_util.h ALE_util.cpp Makefile
	 $(CC) $(FLAGS) $(INCLUDE)  -c -o ALE_util.o ALE_util.cpp
ML_util.o: ALE.h ML_util.h ML_util.cpp Makefile
	 $(CC) $(FLAGS) $(INCLUDE)  -c -o ML_util.o ML_util.cpp
exODT.o: ALE.h exODT.h exODT.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o exODT.o exODT.cpp
model.o: ALE.h exODT.h model.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o model.o model.cpp
traceback.o: ALE.h exODT.h traceback.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o traceback.o traceback.cpp

sample.o: ALE.h exODT.h sample.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o sample.o sample.cpp

traceback_lowmem.o: ALE.h exODT.h traceback_lowmem.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o traceback_lowmem.o traceback_lowmem.cpp


mpi_tree.o: ALE.h exODT.h mpi_tree.h mpi_tree.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o mpi_tree.o mpi_tree.cpp
rate_estimate.o: ALE.h exODT.h mpi_tree.h rate_estimate.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o rate_estimate.o rate_estimate.cpp

libexODT.a: ALE.o ALE_util.o ML_util.o exODT.o model.o traceback.o mpi_tree.o rate_estimate.o sample.o Makefile
	ar rcs libexODT.a ALE.o exODT.o model.o traceback.o sample.o mpi_tree.o rate_estimate.o ALE_util.o ML_util.o

test:	libexODT.a test.cpp Makefile
	$(CC) $(FLAGS) $(INCLUDE) $(LINK) -o test test.cpp $(STATIC)
	strip test

metat:	libexODT.a metat.cpp Makefile
	$(CC) $(FLAGS) $(INCLUDE) $(LINK) -o metat metat.cpp $(STATIC)
	strip metat

transform_tree:	libexODT.a transform_tree.cpp Makefile
	$(CC) $(FLAGS) $(INCLUDE) $(LINK) -o transform_tree transform_tree.cpp $(STATIC)
	strip transform_tree
transform_transfer:	libexODT.a transform_transfer.cpp Makefile
	$(CC) $(FLAGS) $(INCLUDE) $(LINK) -o transform_transfer transform_transfer.cpp $(STATIC)
	strip transform_transfer


test_mpi:	libexODT.a test_mpi.cpp Makefile
	$(mCC) $(FLAGS) $(INCLUDE) $(MPI_INCLUDE) $(LINK) -o test_mpi test_mpi.cpp $(STATIC)
	strip test_mpi

test_mpi3:	libexODT.a test_mpi3.cpp Makefile
	$(mCC) $(FLAGS) $(INCLUDE) $(MPI_INCLUDE) $(LINK) -o test_mpi3 test_mpi3.cpp $(STATIC)
	strip test_mpi3


decompose: decompose.cpp ALE.o Makefile
	$(CC) $(FLAGS) $(INCLUDE) $(LINK) ALE.o -o decompose decompose.cpp $(STATIC)
	strip decompose	

process_trees: process_trees.cpp ALE.o Makefile
	$(mCC) $(FLAGS) $(INCLUDE) $(MPI_INCLUDE) $(LINK) -o process_trees process_trees.cpp $(STATIC)
	strip process_trees	

