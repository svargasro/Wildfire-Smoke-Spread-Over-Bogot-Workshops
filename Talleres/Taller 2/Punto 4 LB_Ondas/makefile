GCC = g++
CFLAGS = -O3 --std=c++17 -Wall
FSANITIZE = -fsanitize=address -fsanitize=leak -fsanitize=undefined

TARGET = $(patsubst %.cpp,%,$(wildcard *.cpp))

all: $(TARGET).jpg

$(TARGET).out: $(TARGET).cpp
	$(GCC) $(CFLAGS) $(TARGET).cpp -o $(TARGET).out

$(TARGET).dat: $(TARGET).out
	./$(TARGET).out 

$(TARGET).jpg: $(TARGET).dat 
	gnuplot $(TARGET).gp

.PHONY: clean

clean: 
	rm *.out *.txt *.x *.dat 