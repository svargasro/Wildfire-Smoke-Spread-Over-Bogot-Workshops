CC = g++
CXXFLAGS = -O2 -std=c++17
TARGET = Pendulos

all: $(TARGET).gif

$(TARGET).x: $(TARGET).cpp
	$(CC) $(CXXFLAGS) $^ -o $@

$(TARGET).txt: $(TARGET).x
	./$^ > $@

$(TARGET).gif: $(TARGET).txt
	gnuplot $^

.PHONY: clean

clean:
	rm -f *.x *.gp *.pdf *.txt *.gif
