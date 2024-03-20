CXX = g++
CXXFLAGS = -O2 -std=c++17

all: Modelo_SIR.gif

%.x: %.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

data.txt: Modelo_SIR.x
	./$^ > $@

Modelo_SIR.gif: data.txt
	gnuplot $^

.PHONY: clean

clean:
	rm -f *.x *.out *.gp *.pdf *.txt *.gif 
