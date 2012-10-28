
CC=g++ -pipe
mCC=mpic++

FLAGS = -O3  -fmerge-all-constants -funroll-loops -DNDEBUG -Wall -std=gnu++11
DEV_FLAGS =  -g -Wall

boost_DIR=/usr/local/
#bpp_DIR= $(HOME)/newest_bpp/
bpp_DIR= /usr/local/
boost_libs= -lboost_mpi -lboost_serialization 


MPI_INCLUDE=-I/usr/lib/openmpi/include/ 

ifndef OSTYPE
  OSTYPE = $(shell uname -s|awk '{print tolower($$0)}')
  export OSTYPE
endif

ifeq ($(OSTYPE),darwin)
	bpp_libs = $(bpp_DIR)lib/libbpp-core.a $(bpp_DIR)lib/libbpp-phyl.a $(bpp_DIR)lib/libbpp-seq.a 
	boost_libs = 
	STATIC = libexODT.a $(boost_DIR)lib/libboost_serialization.a  $(boost_DIR)lib/libboost_mpi.a	
else
	bpp_libs = -lbpp-core -lbpp-seq -lbpp-phyl
	boost_libs = -lboost_mpi -lboost_serialization 
	#boost_DIR= $(HOME)/
	STATIC= libexODT.a	
endif

DYNAMIC =  -L. -L$(bpp_DIR)lib $(bpp_libs) $(boost_libs) 
LINK = $(DYNAMIC) #-lexODT
INCLUDE = -I$(boost_DIR)include -I$(bpp_DIR)include  $(MPI_INCLUDE)


ALE.o: ALE.h ALE.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o ALE.o ALE.cpp

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

mpi_tree.o: ALE.h exODT.h mpi_tree.h mpi_tree.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o mpi_tree.o mpi_tree.cpp

rate_estimate.o: ALE.h exODT.h mpi_tree.h rate_estimate.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o rate_estimate.o rate_estimate.cpp

sample.o: ALE.h exODT.h sample.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o sample.o sample.cpp	


libexODT.a: ALE.o ALE_util.o ML_util.o exODT.o model.o traceback.o mpi_tree.o rate_estimate.o sample.o Makefile
	ar rcs libexODT.a ALE.o exODT.o model.o traceback.o mpi_tree.o rate_estimate.o ALE_util.o ML_util.o sample.o

test:	libexODT.a test.cpp Makefile
	$(CC) test.cpp -o test $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)
	strip test

nni:	libexODT.a nni.cpp Makefile
	$(CC) nni.cpp -o nni $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)
	strip nni

metat:    libexODT.a metat.cpp Makefile
	$(CC) metat.cpp -o metat  $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)
	strip metat

metat_sample:    libexODT.a metat_sample.cpp Makefile
	$(CC) metat_sample.cpp -o metat_sample $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)
	strip metat_sample

test_mpi:	libexODT.a test_mpi.cpp Makefile
	$(mCC) test_mpi.cpp -o test_mpi $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)
	strip test_mpi

test_mpi2:	libexODT.a test_mpi2.cpp Makefile
	$(mCC) test_mpi2.cpp -o test_mpi2 $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)
	strip test_mpi2

test_mpi3:	libexODT.a test_mpi3.cpp Makefile
	$(mCC) test_mpi3.cpp -o test_mpi3 $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)
	strip test_mpi3

decompose: decompose.cpp ALE.o Makefile
	$(CC) decompose.cpp -o decompose $(FLAGS) $(INCLUDE) ALE.o $(LINK)   
	strip decompose	

process_trees: process_trees.cpp ALE.o Makefile
	$(CC) process_trees.cpp -o process_trees $(FLAGS) $(INCLUDE) ALE.o $(LINK)  
	strip process_trees	

small_ales: small_ales.cpp ALE_util.o ALE.o Makefile
	$(CC) small_ales.cpp -o small_ales $(FLAGS) $(INCLUDE) ALE_util.o ALE.o $(LINK)   
	strip small_ales 

strip_strap: strip_strap.cpp Makefile
	$(CC) strip_strap.cpp -o strip_strap $(FLAGS) $(INCLUDE) $(LINK) 
	strip strip_strap

malate:	libexODT.a malate.cpp Makefile
	$(CC) malate.cpp -o malate $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)
	strip malate

test_sample: 	libexODT.a test_sample.cpp Makefile
	$(CC) test_sample.cpp -o test_sample  $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)
	strip test_sample

rate_sample: 	libexODT.a rate_sample.cpp Makefile
	$(CC) rate_sample.cpp -o rate_sample  $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)
	strip rate_sample
