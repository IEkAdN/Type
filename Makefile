#CXXFLAGS = -std=c++0x -MMD -O0 -Wall -g
CXXFLAGS = -std=c++0x -MMD -O3 -Wall

Type : main.o hoge.o fuga.o lib.o
	$(CXX) -o $@ $^

-include *.d

