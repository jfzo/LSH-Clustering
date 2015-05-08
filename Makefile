#COMPILER = g++-mp-4.7
COMPILER = g++
#FLAGS = -Wno-sign-compare -std=c++11 -c -g -Wall 
FLAGS = -fpermissive -Wno-sign-compare -g -Wall -std=c++11 -O3 -pg
OPTIONS = -Wall -std=c++11 
JACCSIM = -D SIMTYPE=1
COSSIM = -D SIMTYPE=2
SIMMOD1 = -D SIMESTMOD=1
SIMMOD2 = -D SIMESTMOD=2
# Use -pg for debugging.
#FLAGS = -Wno-sign-compare -c -g -Wall

# Useit for Linux
#LIBS = -lstxxl_debug  -lboost_system -ll -lboost_filesystem -lgomp

# Useit for MacOS
LIBS = -lstxxl_debug  -lboost_system -ll -lboost_filesystem

%.o : %.cpp
	$(COMPILER) $(FLAGS) -c $<

all: jaccardest_mod1 cosest_mod1 jaccardest_mod2 cosest_mod2


clean: 
	rm -Rf *.o bin/*


jaccardest_mod1: src/lshsimilarity.cpp
	$(COMPILER) $(OPTIONS) $(SIMMOD1) $(JACCSIM) -o bin/jaccLshsimMode1 src/lshsimilarity.cpp

cosest_mod1: src/lshsimilarity.cpp
	$(COMPILER) $(OPTIONS) $(SIMMOD1) $(COSSIM)  -o bin/cosLshsimMode1 src/lshsimilarity.cpp

jaccardest_mod2: src/lshsimilarity.cpp
	$(COMPILER) $(OPTIONS) $(SIMMOD2) $(JACCSIM) -o bin/jaccLshsimMode2 src/lshsimilarity.cpp

cosest_mod2: src/lshsimilarity.cpp
	$(COMPILER) $(OPTIONS) $(SIMMOD2) $(COSSIM)  -o bin/cosLshsimMode2 src/lshsimilarity.cpp


