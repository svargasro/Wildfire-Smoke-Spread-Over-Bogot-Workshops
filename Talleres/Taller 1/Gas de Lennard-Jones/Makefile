CC = g++
CXXFLAGS =-O2 -std=c++17
TARGET = 4.1_GasDeLennard-Jones

all: $(TARGET).gif

$(TARGET).gif.x: $(TARGET).cpp
	$(CC) $(CXXFLAGS) $^ -o $@

$(TARGET).txt: $(TARGET).gif.x
	./$^ > $@

$(TARGET).gif: $(TARGET).txt
	gnuplot $<

.PHONY: clean

clean:
	rm -f *.x *.gp *.pdf *.txt *.gif
