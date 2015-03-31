#COMPILER = g++-mp-4.7
COMPILER = g++
#FLAGS = -Wno-sign-compare -std=c++11 -c -g -Wall 
FLAGS = -fpermissive -Wno-sign-compare -g -Wall -std=c++11 -O3 -pg
# Use -pg for debugging.
#FLAGS = -Wno-sign-compare -c -g -Wall

# Useit for Linux
#LIBS = -lstxxl_debug  -lboost_system -ll -lboost_filesystem -lgomp

# Useit for MacOS
LIBS = -lstxxl_debug  -lboost_system -ll -lboost_filesystem

%.o : %.cpp
	$(COMPILER) $(FLAGS) -c $<

clean: 
	rm -Rf *.o

anglelsh: src/perfectlshangle_custom.cpp
	$(COMPILER) -Wall -std=c++11 -o bin/perfectlshangle_custom src/perfectlshangle_custom.cpp

minwiselsh: src/perfectlshmwise_custom.cpp
	$(COMPILER) -Wall -std=c++11 -o bin/perfectlshmwise_custom src/perfectlshmwise_custom.cpp


