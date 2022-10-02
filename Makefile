#CXXFLAGS = -std=c++0x -MMD -O0 -Wall -g
CXXFLAGS = -std=c++0x -MMD -O3 -Wall

Type : main.o split.o type.o alignment.o
	$(CXX) -o $@ $^

-include *.d

