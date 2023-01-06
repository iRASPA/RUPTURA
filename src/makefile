CXXFLAGS=-g -O3 -std=c++17 -march=native -ffast-math -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion -Wnull-dereference -Wdouble-promotion -Wformat=2 -Werror -fomit-frame-pointer -ftree-vectorize -fno-stack-check -funroll-loops

default: ruptura;

random_numbers.o: random_numbers.cpp random_numbers.h
	$(CXX) $(CXXFLAGS) -c random_numbers.cpp

special_functions.o: special_functions.cpp special_functions.h
	$(CXX) $(CXXFLAGS) -c special_functions.cpp

isotherm.o: isotherm.cpp isotherm.h
	$(CXX) $(CXXFLAGS) -c isotherm.cpp

multi_site_isotherm.o: multi_site_isotherm.cpp multi_site_isotherm.h isotherm.h
	$(CXX) $(CXXFLAGS) -c multi_site_isotherm.cpp

component.o: component.cpp component.h isotherm.h
	$(CXX) $(CXXFLAGS) -c component.cpp

mixture_prediction.o: mixture_prediction.cpp mixture_prediction.h
	$(CXX) $(CXXFLAGS) -c mixture_prediction.cpp

inputreader.o: inputreader.cpp inputreader.h
	$(CXX) $(CXXFLAGS) -c inputreader.cpp

breakthrough.o: breakthrough.cpp breakthrough.h
	$(CXX) $(CXXFLAGS) -c breakthrough.cpp

fitting.o: fitting.cpp fitting.h
	$(CXX) $(CXXFLAGS) -c fitting.cpp

main.o : main.cpp 
	$(CXX) $(CXXFLAGS) -c main.cpp

ruptura: random_numbers.o special_functions.o isotherm.o multi_site_isotherm.o component.o mixture_prediction.o inputreader.o breakthrough.o fitting.o main.o
	$(CXX) main.o fitting.o breakthrough.o inputreader.o mixture_prediction.o component.o multi_site_isotherm.o isotherm.o special_functions.o random_numbers.o -o ruptura

clean:
	rm -f *.pcm *.o *.a ruptura
