CXX = g++
CXXFLAGS = -O2 -std=c++17

all: Modelo_SIR.pdf

%.x: %.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

data.txt: Modelo_SIR.x
	./$^ > $@

Modelo_SIR.pdf: data.txt plot.gp
	gnuplot plot.gp

.PHONY: clean

clean:
	rm -f *.x *.out *.pdf *.txt *.png *.gif
