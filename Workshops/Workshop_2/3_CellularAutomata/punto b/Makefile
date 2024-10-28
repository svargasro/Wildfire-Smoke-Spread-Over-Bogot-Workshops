GCC = g++
GPYTHON = python3
CFLAGS = -O3 --std=c++17 -Wall
FSANITIZE = -fsanitize=address -fsanitize=leak -fsanitize=undefined

TARGET = $(patsubst %.cpp,%,$(wildcard *.cpp))

all: $(TARGET).png

$(TARGET).out: $(TARGET).cpp
	$(GCC) $(CFLAGS) $(TARGET).cpp -o $(TARGET).out

$(TARGET).txt: $(TARGET).out
	./$(TARGET).out 1> $(TARGET).txt 2> parametros.txt

$(TARGET).png: $(TARGET).txt parametros.txt
	$(GPYTHON) $(TARGET).py

.PHONY: clean

clean: 
	rm -f *.out *.txt *.x 